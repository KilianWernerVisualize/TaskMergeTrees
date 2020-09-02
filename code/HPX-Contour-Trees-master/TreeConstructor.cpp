#include "TreeConstructor.h"
#include "DataManager.h"
#include "Log.h"

#include "Arc.h"
#include "Augmentation.h"
#include "Boundary.h"
#include "SweepQueue.h"
#include "profiler.h"

#include <boost/range/irange.hpp>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/execution_policy.hpp>
#include <hpx/include/parallel_sort.hpp>
#include <unistd.h>

#ifdef FILEOUT
#include <fstream>
#include <iostream>
#include <string>
#endif

HPX_REGISTER_COMPONENT_MODULE();

typedef hpx::components::simple_component<TreeConstructor> TreeConstructorComponent;

HPX_REGISTER_COMPONENT(TreeConstructorComponent, TreeConstructor);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::init_action, treeConstructor_init_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::construct_action, treeConstructor_construct_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::scanMinima_action, treeConstructor_scanMinima_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::destroy_action, treeConstructor_destroy_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::startSweep_action, treeConstructor_startSweep_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::continueLocalSweep_action, treeConstructor_continueLocalSweep_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::terminateSweep_action, treeConstructor_terminateSweep_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::continueSweep_action, treeConstructor_continueSweep_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::continueSweepPeer_action, treeConstructor_continueSweepPeer_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::sendSweepCount_action, treeConstructor_sendSweepCount_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::countSweeps_action, treeConstructor_countSweeps_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::finish_action, treeConstructor_finish_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::terminate_action, treeConstructor_terminate_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::registerArc_action, treeConstructor_registerArc_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::postProc_action, treeConstructor_postProc_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::stitch_action, treeConstructor_stitch_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::startTrunk_action, treeConstructor_startTrunk_action);
#ifdef FILEOUT
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::convertToGlobal_action, treeConstructor_convertToGlobal_action);

uint64_t TreeConstructor::convertToGlobal(uint64_t v){
    return this->dataManager->convertToGlobal(v);
}

#endif

/**
 * @brief Initialize data.
 * @param treeConstructors
 * @param input
 */
void TreeConstructor::init(const std::vector<hpx::id_type>& treeConstructors, const std::string& input, const Options& options)
{
    // Store list of tree constructor components and determine index of this component
    this->treeConstructors = treeConstructors;
    for (uint32_t i = 0; i < this->treeConstructors.size(); ++i) {
        if (this->treeConstructors[i] == this->get_id()) {
            this->index = i;
            break;
        }
    }

    // Load data
    this->dataManager = DataManager::load(input, this->index, this->treeConstructors.size());
    Profiler profiler(__func__);
    this->options = options;

    // Init structures
    this->sweeps.store(0l);
    this->numMinima = 0l;

    numVertices = this->dataManager->getNumVerticesLocal(true); // TODO

    this->swept.init(numVertices, INVALID_VERTEX, this->dataManager);
    this->UF.init(numVertices, INVALID_VERTEX, this->dataManager);
    this->arcMap.init(numVertices, nullptr, this->dataManager);

    this->sweepsLocal.store(0l);
    sweepsLocal_last = 0l;
    scanCount = 0ul;

    this->done.reset();
    this->leavesDone.reset();
    this->termination.reset();
    this->postProcessing.reset();
    this->trunkStarted.reset();

    graph.init("graph.csv");


#ifdef FILEOUT
    this->clock.restart();
#endif

    hpx::apply(this->executor_low, TreeConstructor::sendSweepCount_action(), this->get_id());
}

/**
 * @brief Shutdown.
 */
void TreeConstructor::destroy()
{
    if (this->dataManager){
#ifdef FILEOUT
        std::streampos size;

        std::ofstream file ("output_"+ std::to_string(this->index) +".bin", std::ios::out|std::ios::binary);

        if (file.is_open()){

            if (this->index == 0){
                uint64_t x = (uint64_t)this->dataManager->getDimensions().x;
                uint64_t y = (uint64_t)this->dataManager->getDimensions().y;
                LogInfo().tag(std::to_string(this->index)) << "max_x : " << x << " max_y: " << y;
                file.write((char*)&x, sizeof(uint64_t));
                file.write((char*)&y, sizeof(uint64_t));
            }

            for (fileEntry entry : this->arcFileEntries){
                file.write((char*)&entry.type, sizeof(char));
                file.write((char*)&entry.time, sizeof(timestamp));
                file.write((char*)&entry.id1, sizeof(std::uint64_t));
                file.write((char*)&entry.id2, sizeof(std::uint64_t));
            }
            file.close();
        }
#endif
#ifdef VTIOUT
        vtkSmartPointer<vtkImageData> imageData =
                vtkSmartPointer<vtkImageData>::New();
        imageData->SetDimensions(this->dataManager->getLocalDimensions().x, this->dataManager->getLocalDimensions().y, this->dataManager->getLocalDimensions().z);
#if VTK_MAJOR_VERSION <= 5
        imageData->SetNumberOfScalarComponents(1);
        imageData->SetScalarTypeToDouble();
#else
        imageData->AllocateScalars(VTK_UNSIGNED_LONG, 1);
#endif
        int* dims = imageData->GetDimensions();

        uint64_t* data = static_cast<uint64_t*>(imageData->GetScalarPointer());

        int j = 0;

        for (int i = 0; (i < this->dataManager->getNumVerticesLocal(true)); i++){

            uint64_t v = this->dataManager->getLocalVertex(i);
            if (this->dataManager->isGhost(v))
                continue;
            /*
            uint64_t host = swept[localVertices[i]] >> BLOCK_INDEX_SHIFT;
            if (host < this->treeConstructors.size())
                data[i] = hpx::async<TreeConstructor::convertToGlobal_action>(this->treeConstructors[host], swept[localVertices[i]]).get() & VERTEX_INDEX_MASK;
                */
            if (this->dataManager->isLocal(swept[v])){
                data[j] = this->dataManager->getValue(swept[v]).value;
            }
            else {
                data[j] = 0ul;
            }
            j++;
        }

        vtkSmartPointer<vtkXMLImageDataWriter> writer =
                vtkSmartPointer<vtkXMLImageDataWriter>::New();
        std::string filename = "output_" + std::to_string(this->index) + ".vti";
        writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(imageData->GetProducerPort());
#else
        writer->SetInputData(imageData);
#endif
        writer->Write();
#endif
        delete this->dataManager;
    }
}

void TreeConstructor::finish(){
    this->done.set();
}

void TreeConstructor::terminate(){
    this->done.set();
    this->termination.set();
}

void TreeConstructor::postProc(){
    if (std::atomic_fetch_add(&this->remoteLeavesDone, -1) == 1){
        LogInfo().tag(std::to_string(this->index)) << "Post Processing starts!";
        this->postProcessing.set();
    }
}

void TreeConstructor::stitch(std::vector<Value> danglings){

    Profiler profiler(__func__);
    if (danglings.size() == 0)
        return;

    Arc* trunkStart;
    fetchCreateArc(trunkStart, danglings[0].vertex);

    //trunkStart->body->augmentation.inherit(trunkStart->inheritedAugmentations);

    if (danglings.size() > 1){
        for (int i = 0; (i < danglings.size() - 2); i++){
            if (this->dataManager->isLocal(danglings[i].vertex)){
                Arc* newArc;
                fetchCreateArc(newArc, danglings[i].vertex);
                newArc->saddle = danglings[i+1].vertex;
#ifdef FILEOUT
                fileEntry entry;
                entry.id1 = (hpx::async<TreeConstructor::convertToGlobal_action>(this->treeConstructors[newArc->extremum >> BLOCK_INDEX_SHIFT], newArc->extremum).get()) & VERTEX_INDEX_MASK;
                entry.id2 = (hpx::async<TreeConstructor::convertToGlobal_action>(this->treeConstructors[newArc->saddle >> BLOCK_INDEX_SHIFT], newArc->saddle).get()) & VERTEX_INDEX_MASK;
                entry.time = this->clock.elapsed_microseconds();
                entry.type = 2;
                this->arcFileLock.lock();
                this->arcFileEntries.push_back(entry);
                this->arcFileLock.unlock();
#endif
            }
            //if (this->arcMap.contains(danglings[i].vertex)){
            //    Arc* arc = arcMap[danglings[i].vertex];
#ifndef FLATAUGMENTATION
                for (auto iter = arc->body->augmentation.vertices.begin(); (iter != arc->body->augmentation.vertices.end()); iter.operator ++()){
                    swept[(*iter)->key.vertex] = INVALID_VERTEX;
#else
//                for (Value c : arc->body->augmentation.vertices){
//                    swept[c.vertex] = INVALID_VERTEX;
#endif
                    /*
                    LogDebug().tag(std::to_string(this->index)) << "Dangling vertex "
                                                                << ( (*iter)->key.vertex  >> BLOCK_INDEX_SHIFT) << "::" << ((*iter)->key.vertex  & VERTEX_INDEX_MASK) << " from dangling arc "
                                                                << ( arc->extremum  >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum  & VERTEX_INDEX_MASK);*/
                //}
  //              arc->body->augmentation.clear();
            //}
        }
    }

    hpx::util::high_resolution_timer timer;


    hpx::parallel::for_each(
        hpx::parallel::execution::parallel_unsequenced_policy(),
        boost::irange(0ul, numVertices).begin(), boost::irange(0ul, numVertices).end(),
        [&](std::uint64_t i)
    {
        std::uint64_t v = this->dataManager->getLocalVertex(i);

        if (swept[v] == INVALID_VERTEX){
            if (this->dataManager->getValue(v) < danglings[0]){
                swept[v] = INVALID_VERTEX; //To-Do Trunkstart!

            } else if (this->dataManager->getValue(v) > danglings[danglings.size()-1]){
                swept[v] = danglings[danglings.size()-1].vertex;
                Arc* arc;
                fetchCreateArc(arc, danglings[danglings.size()-1].vertex);
                arc->body->lock.lock();
                arc->body->augmentation.sweep(this->dataManager->getValue(v));
                arc->body->lock.unlock();
            } else {
                int64_t min = 0;
                int64_t max = danglings.size()-1;
                while (true){
                    int64_t j = (max+min)/2;
                    if (this->dataManager->getValue(v) > danglings[j]){
                        if (this->dataManager->getValue(v) < danglings[j+1]){
                            swept[v] = danglings[j].vertex;
                            Arc* arc;
                            fetchCreateArc(arc, danglings[j].vertex);
                            arc->body->lock.lock();
                            arc->body->augmentation.sweep(this->dataManager->getValue(v));
                            arc->body->lock.unlock();
                            break;
                        } else {
                            min = j+1;
                            if (min > max){
                                swept[v] = danglings[min].vertex;
                                Arc* arc;
                                fetchCreateArc(arc, danglings[min].vertex);
                                arc->body->lock.lock();
                                arc->body->augmentation.sweep(this->dataManager->getValue(v));
                                arc->body->lock.unlock();
                                break;
                            }
                        }
                    } else {
                        max = j-1;
                        if (min > max){
                            swept[v] = danglings[min].vertex;
                            Arc* arc;
                            fetchCreateArc(arc, danglings[min].vertex);
                            arc->body->lock.lock();
                            arc->body->augmentation.sweep(this->dataManager->getValue(v));
                            arc->body->lock.unlock();
                            break;
                        }
                    }
                }
            }
        }
    });
            LogInfo().tag(std::to_string(this->index)) << "Stitch: " << timer.elapsed() << " s";
#ifdef VTIOUT
#ifndef FLATAUGMENTATION
    for (Arc* arc : this->arcMap.local){
        if (arc != nullptr){
            for (auto iter = arc->body->augmentation.vertices.begin(); (iter != arc->body->augmentation.vertices.end()); iter.operator ++()){
                swept[(*iter)->key.vertex] = arc->extremum;
            }
        }
    }

    for (auto entry : this->arcMap.remote){
        if (entry.second != nullptr){
            Arc* arc = entry.second;
            for (auto iter = arc->body->augmentation.vertices.begin(); (iter != arc->body->augmentation.vertices.end()); iter.operator ++()){
                swept[(*iter)->key.vertex] = arc->extremum;
            }
        }
    }
#else
    for (Arc* arc : this->arcMap.local){
        if (arc != nullptr){
            for (Value c : arc->body->augmentation.vertices){
                swept[c.vertex] = arc->extremum;
            }
        }
    }

    for (auto entry : this->arcMap.remote){
        if (entry.second != nullptr){
            Arc* arc = entry.second;
            for (Value c : arc->body->augmentation.vertices){
                swept[c.vertex] = arc->extremum;
            }
        }
    }
#endif
#endif
    postProcessing.set();
}

void TreeConstructor::sendSweepCount(){
    if (this->termination.occurred())
        return;

    int64_t totalCount = this->sweepsLocal.load() - this->sweepsLocal_last;
    this->sweepsLocal_last += totalCount;


    if (totalCount != 0l){
        hpx::apply<TreeConstructor::countSweeps_action>(this->treeConstructors[0], totalCount);
    }

    hpx::apply(this->executor_low, TreeConstructor::sendSweepCount_action(), this->get_id());
}

void TreeConstructor::countSweepsLocal(int64_t diff){
    std::int64_t currentcount = std::atomic_fetch_add(&this->sweepsLocal, diff) + diff;
    if (currentcount % 100000 == 0 || (currentcount < 1000 && currentcount > -1000) || (currentcount < -19700000)){
        LogDebug().tag(std::to_string(this->index))<< "Current Count: " << currentcount;
    }

    //LogDebug().tag(std::to_string(this->index)) << "Counting sweeps locally " << sweepsLocal[hpx::get_worker_thread_num()] << "@" << hpx::get_worker_thread_num();
}

void TreeConstructor::countSweeps(std::int64_t diff){
    //this->countLock.lock();
    if (diff == std::numeric_limits<std::int64_t>::max()){
        if (std::atomic_fetch_add(&this->remoteLeavesDone, 1) == (this->treeConstructors.size()-1)){
            LogDebug().tag(std::to_string(this->index)) << "LeavesDone!";
            this->leavesDone.set();
        }
        diff = 0l;
    }
    std::int64_t old = std::atomic_fetch_add(&this->sweeps, diff);
    std::int64_t currentcount = old + diff;
    if (currentcount % 100000 == 0 || (currentcount < 1000 && currentcount > -1000)){
        LogDebug().tag(std::to_string(this->index))<< "Current CountSweeps: " << currentcount;
    }

    LogDebug().tag(std::to_string(this->index)) << "Increasing Counter " << old << " by " << diff;
    if (diff > 0l){
        if (std::atomic_fetch_add(&this->remoteLeavesDone, 1) == (this->treeConstructors.size()-1)){
            LogDebug().tag(std::to_string(this->index)) << "LeavesDone!";
            this->leavesDone.set();
        }
    }

    if ((old + diff <= (options.trunkskip ? 0l : -1l)) && leavesDone.occurred()) {
        LogDebug().tag(std::to_string(this->index)) << "Terminating!";
        for (hpx::id_type peer : this->treeConstructors){
            hpx::apply<TreeConstructor::terminate_action>(peer);
        }
    } else if (old + diff <= 1l && options.trunkskip && leavesDone.occurred()){
        LogDebug().tag(std::to_string(this->index)) << "Finishing!";
        for (hpx::id_type peer : this->treeConstructors){
            hpx::apply(this->executor_high, TreeConstructor::finish_action(), peer);
        }
    }

    //this->countLock.unlock();
}

void TreeConstructor::registerArc(std::vector<Value> extrema){
    Profiler profiler(__func__);
    this->mapLock.lock();
    danglingArcs.insert(danglingArcs.end(), extrema.begin(), extrema.end());
    this->mapLock.unlock();
}

void TreeConstructor::startTrunk(Value v){
    this->trunkStart = v;
    this->trunkStarted.set();
}

void TreeConstructor::scanMinima(){
    for (uint64_t i = 0; (i < 1000 && scanCount < this->dataManager->getNumVerticesLocal(true)); ++i) {
        uint64_t v = this->dataManager->getLocalVertex(scanCount);
        if (this->dataManager->isMinimum(v)) {

            uint64_t neighbors[6];
            uint32_t numNeighbors = this->dataManager->getNeighbors(v, neighbors);

            assert(v != INVALID_VERTEX);
            if (!this->dataManager->isGhost(v)){
                hpx::apply(TreeConstructor::startSweep_action(), this->get_id(), v, true);
                ++numMinima;
                if (numMinima % 5000 == 0){
                    LogInfo().tag(std::to_string(this->index))<< "Scheduled Sweeps: " << numMinima;
                }
            }
        }
        scanCount++;
    }

    if (scanCount < this->dataManager->getNumVerticesLocal(true)){
        hpx::apply(this->executor_low, TreeConstructor::scanMinima_action(), this->treeConstructors[0]);
        return;
    }

    LogInfo().tag(std::to_string(this->index))<< "ALL SWEEPS SCHEDULED!: " << numMinima;

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "Minima (local): " << numMinima;
#endif
    if (numMinima == 0l)
        numMinima = std::numeric_limits<std::int64_t>::max();

    if (!options.trunkskip){
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "Sending " << numMinima << " leaves!";
#endif
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], std::numeric_limits<std::int64_t>::max());
    } else {
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "Sending " << numMinima << " leaves!";
#endif
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], numMinima);
    }
}


/**
 * @brief Construction main entry point.
 */
uint64_t TreeConstructor::construct()
{
    graph.taskStart("construct", "construct", "construct", "construct");
    Profiler profiler(__func__);
    hpx::util::high_resolution_timer timer;

    // Find minima in local data
    std::vector<uint64_t> localVertices = this->dataManager->getLocalVertices();

    for (uint64_t v : localVertices) {
        if (this->dataManager->isMinimum(v)) {

            uint64_t neighbors[6];
            uint32_t numNeighbors = this->dataManager->getNeighbors(v, neighbors);

            assert(v != INVALID_VERTEX);
            if (!this->dataManager->isGhost(v)){
                graph.taskStart("construct", "construct", "start"+std::to_string(v), "startsweep");
                hpx::apply(this->executor_start_sweeps, TreeConstructor::startSweep_action(), this->get_id(), v, true);
                ++numMinima;
            }
        }
    }

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "Minima (local): " << numMinima;
#endif
    if (numMinima == 0l)
        numMinima = std::numeric_limits<std::int64_t>::max();

    if (!options.trunkskip){
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "Sending " << numMinima << " leaves!";
#endif
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], std::numeric_limits<std::int64_t>::max());
    } else {
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "Sending " << numMinima << " leaves!";
#endif
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], numMinima);
    }

    this->termination.wait();

    graph.taskStop("construct", "construct");
    if (this->index == 0){
        this->postProc();
        this->postProcessing.wait();

        int numArcs = 0;
        int numArcss[hpx::get_num_worker_threads()];
        std::vector<std::vector<Value>> danglingArcss;
        danglingArcss.resize(hpx::get_num_worker_threads());
        for (int i = 0; (i < hpx::get_num_worker_threads()); i++){
            numArcss[i] = 0;
        }

        hpx::parallel::for_each(
            hpx::parallel::execution::parallel_unsequenced_policy(),
            this->arcMap.local.begin(), this->arcMap.local.end(),
            [&](Arc* arc)
        {
            if (arc){
                if (arc->saddle != INVALID_VERTEX){
                    numArcss[hpx::get_worker_thread_num()]++;
                } else {
                    danglingArcss[hpx::get_worker_thread_num()].push_back(this->dataManager->getValue(arc->extremum));
                }
            }
        });

        this->trunkStarted.wait();
        danglingArcs.push_back(this->trunkStart);

        for (int i = 0; (i < hpx::get_num_worker_threads()); i++){
            numArcs += numArcss[i];
            danglingArcs.insert(danglingArcs.end(), danglingArcss[i].begin(), danglingArcss[i].end());
        }
        numArcs += danglingArcs.size();
        hpx::parallel::sort(hpx::parallel::execution::parallel_unsequenced_policy(), this->danglingArcs.begin(), this->danglingArcs.end(),
                            [&](const Value a1, const Value a2) {
                                return a1 < a2;
                            });

        if (this->danglingArcs[0] != trunkStart)
            LogError().tag(std::to_string(this->index)) << "WTF";

        for (int i = 1; (i < this->treeConstructors.size()); i++){
            hpx::apply<TreeConstructor::stitch_action>(this->treeConstructors[i], danglingArcs);
        }

        stitch(danglingArcs);

        LogInfo().tag(std::to_string(this->index)) << "Tree Constructed: " << timer.elapsed() << " s";

        LogInfo().tag(std::to_string(this->index))<< "Arcs: " << numArcs;

        return numArcs;

    } else {
        std::vector<std::vector<Value>> danglingArcs;
        danglingArcs.resize(hpx::get_num_worker_threads());
        int numArcs = 0;
        int numArcss[hpx::get_num_worker_threads()];
        for (int i = 0; (i < hpx::get_num_worker_threads()); i++){
            numArcss[i] = 0;
        }

        hpx::parallel::for_each(
            hpx::parallel::execution::parallel_unsequenced_policy(),
            this->arcMap.local.begin(), this->arcMap.local.end(),
            [&](Arc* arc)
        {
            if (arc){
                if (arc->saddle == INVALID_VERTEX){
                    danglingArcs[hpx::get_worker_thread_num()].push_back(this->dataManager->getValue(arc->extremum));
                }  else {
                    numArcss[hpx::get_worker_thread_num()]++;
                }
            }
        });

        std::vector<hpx::future<void>> registered;

        for (int i = 0; (i < hpx::get_num_worker_threads()); i++){
            numArcs += numArcss[i];
            registered.push_back(hpx::async(this->executor_high, TreeConstructor::registerArc_action(), this->treeConstructors[0], danglingArcs[i]));
        }

        hpx::wait_all(registered);

        LogInfo().tag(std::to_string(this->index))<< "Arcs: " << numArcs;
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "Sending " << registered.size() << " danglings!";
#endif
        hpx::async<TreeConstructor::postProc_action>(this->treeConstructors[0]);

        this->postProcessing.wait();

        return numArcs;
    }
    /*
     * TODO:
     * Send all Arc Information to Node0 (or do this on Arc completion)
     */

    /*
    std::vector<Arc*> arcList;
        for (Arc* arc : arcMap)
            if (arc != nullptr && arc->saddle == INVALID_VERTEX)
                arcList.push_back(arc);

        std::sort(arcList.begin(), arcList.end(),
            [&](const Arc* a1, const Arc* a2) {
                return this->dataManager->getValue(a1->extremum) < this->dataManager->getValue(a2->extremum);
            });

    #pragma omp parallel for
        for (uint64_t i = 0; i < arcList.size() - 1; ++i)
            arcList[i]->saddle = arcList[i + 1]->extremum;

        for (uint64_t v : localVertices){
            swept[v & VERTEX_INDEX_MASK] = INVALID_VERTEX;
        }

    #pragma omp parallel for
        for (uint64_t i = 0; (i < localVertices.size()); i++){
            if (arcMap[i]){
                Arc* arc = arcMap[i];
                for (auto iter = arc->body->augmentation.vertices.begin(); (iter != arc->body->augmentation.vertices.end()); iter.operator ++()){
                    if (swept[(*iter)->key.vertex  & VERTEX_INDEX_MASK] != INVALID_VERTEX){
                        if (arcMap[swept[(*iter)->key.vertex  & VERTEX_INDEX_MASK] & VERTEX_INDEX_MASK]->extremum < arc->extremum)
                            swept[(*iter)->key.vertex  & VERTEX_INDEX_MASK] = arc->extremum;
                    } else {
                        swept[(*iter)->key.vertex  & VERTEX_INDEX_MASK] = arc->extremum;
                    }
                }
            }
        }

        hpx::parallel::for_each(
            hpx::parallel::execution::parallel_unsequenced_policy(),
            boost::irange(0, (int)localVertices.size()).begin(), boost::irange(0, (int)localVertices.size()).end(),
            [&](std::size_t i)
            {
                if (swept[i] == INVALID_VERTEX){
                    if (this->dataManager->getValue(i) < this->dataManager->getValue(arcList[0]->extremum)){
                       swept[i] = trunkStart->extremum;
                    } else if (this->dataManager->getValue(i) > this->dataManager->getValue(arcList[arcList.size()-1]->extremum)){
                        swept[i] = arcList[arcList.size()-1]->extremum;
                    } else {
                        int64_t min = 0;
                        int64_t max = arcList.size()-1;
                        while (true){
                            int64_t j = (max+min)/2;
                            if (this->dataManager->getValue(i) > this->dataManager->getValue(arcList[j]->extremum)){
                                if (this->dataManager->getValue(i) < this->dataManager->getValue(arcList[j+1]->extremum)){
                                    swept[i] = arcList[j]->extremum;
                                    break;
                                } else {
                                    min = j+1;
                                    if (min > max){
                                        swept[i] = arcList[min]->extremum;
                                        break;
                                    }
                                }
                            } else {
                                max = j-1;
                                if (min > max){
                                    swept[i] = arcList[min]->extremum;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        );
*/

    //Log().tag(std::to_string(this->index)) << "Local Arcs Constructed: " << timer.elapsed() << " s";
    /*
    for (int i = 0; (i < localVertices.size()); i++){
        Arc* arc = arcMap[i];
        if (arc)
            Log()
    }

    for (int i = 0; (i < localVertices.size()); i++){
        Log().tag(std::to_string(this->index)) << "Index " << i << " valued " << this->dataManager->getValue(i).value << " swept by " << (swept[i] & VERTEX_INDEX_MASK);
    }*/
}

void TreeConstructor::sweepLoop(Arc* arc){
    //Profiler profiler(__func__);
    SweepQueue& queue = arc->body->queue;
    Boundary& boundary = arc->body->boundary;
    uint64_t& v = arc->extremum;

    while (!queue.empty()) {
    //for (int j = 0; (j < 100); j++){
        uint64_t c = queue.pop();

        //Log().tag(std::to_string(this->index)) << " cont touching " << (c & VERTEX_INDEX_MASK);

        if (c == INVALID_VERTEX) // empty (queue skips elements which have been swept in the mean time)  -> could move skip check to empty()
            break;

        const Value value = this->dataManager->getValue(c);

        if (this->dataManager->isGhost(c)){
//#ifdef ENABLE_DEBUG_LOGGING
//            LogDebug().tag(std::to_string(this->index)) << "ContinueSweep removing from Boundary of: " << (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
//                                                        << " the vertex: " << (c >> BLOCK_INDEX_SHIFT) << "::" << (c& VERTEX_INDEX_MASK);
//#endif
            boundary.remove(value); // remove from boundary

            //mark as swept
            swept[c] = v;
//#ifdef ENABLE_DEBUG_LOGGING
//            LogDebug().tag(std::to_string(this->index)) << "Augmenting: " << (c >> BLOCK_INDEX_SHIFT) << "::" << (c & VERTEX_INDEX_MASK)
//                                                        << " to: " << (v>> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK);
//#endif

            uint64_t neighbors[6];
            uint32_t numNeighbors = this->dataManager->getNeighbors(c, neighbors);

            for (int i = 0; (i < numNeighbors); i++){
                if ((!this->dataManager->isGhost(neighbors[i])) && (this->swept[neighbors[i]] == INVALID_VERTEX)){
                    queue.push(neighbors[i]);
                }
            }
//#ifdef ENABLE_DEBUG_LOGGING
//            LogDebug().tag(std::to_string(this->index)) << " sweeping ghost " << (c >> BLOCK_INDEX_SHIFT) << "::" << (c & VERTEX_INDEX_MASK) << " from min: " << (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK);
//#endif
            continue;
        }
//#ifdef ENABLE_DEBUG_LOGGING
//        LogDebug().tag(std::to_string(this->index)) << "Sweep from " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
//                                                   << " touches " <<  (c >> BLOCK_INDEX_SHIFT) << "::" << (c & VERTEX_INDEX_MASK);
//#endif
        // check if can be swept
        if (this->touch(c, v)) {
//#ifdef ENABLE_DEBUG_LOGGING
//            LogDebug().tag(std::to_string(this->index)) << "Sweep from " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
//                                                       << " sweeps " <<  (c >> BLOCK_INDEX_SHIFT) << "::" << (c & VERTEX_INDEX_MASK);
//#endif
            //Log().tag(std::to_string(this->index)) << " cont sweeping " << (c & VERTEX_INDEX_MASK);
//#ifdef ENABLE_DEBUG_LOGGING
//            LogDebug().tag(std::to_string(this->index)) << "SweepLoop removing from Boundary of: " << (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
//                                                                    << " the vertex: " << (value.vertex >> BLOCK_INDEX_SHIFT) << "::" << (value.vertex& VERTEX_INDEX_MASK);
//#endif
            boundary.remove(value); // remove from boundary

            // put into our augmentation and mark as swept
            this->swept[c] = v;
            arc->body->augmentation.sweep(value);
//#ifdef ENABLE_DEBUG_LOGGING
//            LogDebug().tag(std::to_string(this->index)) << "Augmenting: " << (c >> BLOCK_INDEX_SHIFT) << "::" << (c & VERTEX_INDEX_MASK)
//                                                                        << " to: " << (v>> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK);
//#endif

            uint64_t neighbors[6];
            uint32_t numNeighbors = this->dataManager->getNeighbors(c, neighbors);

            for (int i = 0; (i < numNeighbors); i++){
                if (this->dataManager->isGhost(neighbors[i]) && (this->swept[neighbors[i]] == INVALID_VERTEX)){
                    arc->body->remoteCallLock.lock();
#ifdef LAZYREMOTECALLS
                    arc->body->remoteCalls.push_back(RemoteCall(neighbors[i], c));
#else
                    arc->body->remoteCalls.push_back(hpx::async<TreeConstructor::continueSweep_action>(this->treeConstructors[this->dataManager->getBlockIndex(neighbors[i])], v, this->dataManager->convertToGlobal(neighbors[i]), this->dataManager->convertToGlobal(c), false));
#endif
                    arc->body->remoteCallLock.unlock();
                } else {
                    queue.push(neighbors[i]);
                }
            }
        } else {
            // not allowed to sweep yet -> add as potential boundary candidate
//#ifdef ENABLE_DEBUG_LOGGING
//            LogDebug().tag(std::to_string(this->index)) << "Adding to Boundary of: " << (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
//                                                        << " the vertex: " << (value.vertex >> BLOCK_INDEX_SHIFT) << "::" << (value.vertex& VERTEX_INDEX_MASK);
//#endif
            boundary.add(value);
        }
    }
}

std::vector<RemoteAnswer> TreeConstructor::continueSweep(uint64_t v, uint64_t c, uint64_t from, bool leaf)
{
    uint64_t sent = dataManager->convertToLocal(c);
    uint64_t swept_from = dataManager->convertToLocal(from);
#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "receiving " << (c & VERTEX_INDEX_MASK) << " @global: " << this->dataManager->getGlobalCoords(c) << " taking as: " << (sent & VERTEX_INDEX_MASK) << " @local: " << this->dataManager->getLocalCoords(sent);
#endif
    bool leaveItToContinuePeer = false;

    //Fetch Arc
    Arc* arc;
    this->mapLock.lock();
    if (fetchCreateArc(arc, v)){
        this->mapLock.unlock();
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << " cont anew " << (v >> BLOCK_INDEX_SHIFT) << (v & VERTEX_INDEX_MASK);
#endif
    } else {
        this->mapLock.unlock();
        arc->body->lock.lock();
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << " cont not anew " << (v >> BLOCK_INDEX_SHIFT) << (v & VERTEX_INDEX_MASK) << " state " << arc->body->state;
#endif
        if (arc->body->state == 0) {
            leaveItToContinuePeer = true;
        } else {
            arc->body->lock.unlock();
        }
    }

    Boundary& boundary = arc->body->boundary;
    SweepQueue& queue = arc->body->queue;

    //Sweep Startpoint
    queue.push(swept_from);  //pushing ghost to queue instantly sweeps it and issues sweep from there.

    if (leaveItToContinuePeer){
        arc->body->queue.push(sent);
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << " injecting for future continueSweepPeer " << (sent & VERTEX_INDEX_MASK);
#endif
        arc->body->lock.unlock();
        RemoteAnswer answer(this->get_id());
        std::vector<RemoteAnswer> answerMap;
        answerMap.push_back(answer);
        return answerMap;
    }

    //Initial Fill of Queue
    arc->body->lock.lock();
    if (swept[sent] != INVALID_VERTEX){
        assert(swept[sent] == v);
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << " ignoring already known " << (sent & VERTEX_INDEX_MASK);
#endif
        arc->body->lock.unlock();
        RemoteAnswer answer(this->get_id());
        std::vector<RemoteAnswer> answerMap;
        answerMap.push_back(answer);
        return answerMap;
    }

    /*
    if (this->index == 1u){
        arc->body->lock.unlock();
        RemoteAnswer answer(this->get_id());
        return answer;
    }
    */
    arc->body->queue.push(sent);
    //No boundary or augmentation initialization needed

    //Inject in possibly running sweep
    if (arc->body->state == 1) {
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << " injecting " << (sent & VERTEX_INDEX_MASK);
#endif
        arc->body->lock.unlock();
        RemoteAnswer answer(this->get_id());
        std::vector<RemoteAnswer> answerMap;
        answerMap.push_back(answer);
        return answerMap;
    }

    arc->body->state = 1;
    arc->body->lock.unlock();

    uint64_t min;

    while (true){
        this->sweepLoop(arc);
        arc->body->lock.lock();
        if (arc->body->queue.empty()){
            arc->body->state = 2;

            if (arc->body->boundary.empty()){
                min = INVALID_VERTEX;
            } else {
                min = boundary.min();
            }
            break;
        }
        arc->body->lock.unlock();
    }

    //Couldn't we do this only on terminate?
    if ((arc->saddle == INVALID_VERTEX || this->dataManager->getValue(arc->saddle) > this->dataManager->getValue(min)) || swept[arc->saddle] == v)
        arc->saddle = min;

    min = arc->saddle;

    RemoteAnswer myWinner(this->get_id());
    myWinner.counter = ++arc->body->counter;
    myWinner.minimum = min;
    myWinner.value = this->dataManager->getValue(min);

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "ContSweeploop empty for " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK) << " with min: " << (min >> BLOCK_INDEX_SHIFT) << "::" << (min & VERTEX_INDEX_MASK) << " awaiting remotes.";
#endif

    arc->body->lock.unlock();

    AnswerMap answerMap = awaitRemoteAnswers(arc);

    if (myWinner.counter > answerMap[this->get_id()].counter){
        answerMap[this->get_id()] = myWinner;
    }

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "cont " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)  << " returns answerMap with size: " << answerMap.size();
#endif

    std::vector<RemoteAnswer> answers;
    answers.reserve(answerMap.size());

    for (auto entry: answerMap){
        answers.push_back(entry.second);
    }

    return answers;
}

std::vector<RemoteAnswer> TreeConstructor::continueSweepPeer(uint64_t v, const std::vector<std::uint64_t> children)
{
    //Fetch Arc
    Arc* arc;
    this->mapLock.lock();
    fetchCreateArc(arc, v);
    this->mapLock.unlock();

    Boundary& boundary = arc->body->boundary;
    SweepQueue& queue = arc->body->queue;

    //Should be identical?
    arc->body->children = children;

    //Sweep Startpoint
    this->swept[v] = v;

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "Arc from " <<  (arc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum & VERTEX_INDEX_MASK)
                                                << " continuePeer finds childcount " << children.size();
#endif

    //No Initial Fill except boundary intersections
    for (uint64_t child : children){
        this->swept[child] = v;
        this->UF[child] = v;
    }

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "Arc from " << (arc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum & VERTEX_INDEX_MASK)
                                                << " continuePeer starts with min: "  <<  (boundary.min() >> BLOCK_INDEX_SHIFT) << "::" << (boundary.min() & VERTEX_INDEX_MASK);
#endif

    // Mutual intersections of child boundaries -> this is where sweeps must continue
    mergeBoundaries(arc);

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "Arc from " << (arc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum & VERTEX_INDEX_MASK)
                                                << " continuePeer extracted min from childboundaries:"  <<  (boundary.min() >> BLOCK_INDEX_SHIFT) << "::" << (boundary.min() & VERTEX_INDEX_MASK);
#endif

    arc->body->augmentation.inherit(arc->body->inheritedAugmentations);

    //Inject in possibly running sweep
    arc->body->lock.lock();
    if (arc->body->state == 1){
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << " injecting old boundary intersection.";
#endif

        arc->body->lock.unlock();
        RemoteAnswer answer(this->get_id());
        std::vector<RemoteAnswer> answerMap;
        answerMap.push_back(answer);
        return answerMap;
    }

    arc->body->state = 1;
    arc->body->lock.unlock();

    uint64_t min;

    while (true){
        this->sweepLoop(arc);
        arc->body->lock.lock();
        if (arc->body->queue.empty()){
            arc->body->state = 2;

            if (arc->body->boundary.empty()){
                min = INVALID_VERTEX;
            } else {
                min = boundary.min();
            }
            break;
        }
        arc->body->lock.unlock();
    }

    if ((arc->saddle == INVALID_VERTEX || (this->dataManager->getValue(arc->saddle) > this->dataManager->getValue(min))) || swept[arc->saddle] == v)
        arc->saddle = min;

    min = arc->saddle;

    RemoteAnswer myWinner(this->get_id());
    myWinner.counter = ++arc->body->counter;
    myWinner.minimum = min;
    myWinner.value = this->dataManager->getValue(min);

    arc->body->lock.unlock();

    AnswerMap answerMap = awaitRemoteAnswers(arc);

    if (myWinner.counter > answerMap[this->get_id()].counter){
        answerMap[this->get_id()] = myWinner;
    }

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "cont returns answerMap with size: " << answerMap.size();
#endif
    std::vector<RemoteAnswer> answers;
    answers.reserve(answerMap.size());

    for (auto entry: answerMap){
        answers.push_back(entry.second);
    }

    return answers;
}

bool TreeConstructor::terminateSweep(uint64_t v, Value s, std::vector<hpx::id_type> peers){

    this->mapLock.lock();
    Arc*& arc = this->arcMap[v];
    assert (arc != nullptr);
    this->mapLock.unlock();

    if ((static_cast<uint32_t>(s.vertex >> BLOCK_INDEX_SHIFT)) == this->index){
        if (s.vertex != arc->saddle){
            LogError().tag(std::to_string((this->index))) << "woopsie";
            arc->saddle = s.vertex;
        }
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "TerminateSweep removing from Boundary of: " << (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
                                                    << " the vertex: " << (arc->saddle >> BLOCK_INDEX_SHIFT) << "::" << (arc->saddle & VERTEX_INDEX_MASK);
#endif
        arc->body->boundary.remove(s);
        this->assignSaddle(arc, s, peers);

//#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "Arc from " <<  (arc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum & VERTEX_INDEX_MASK) << " to "  << (arc->saddle >> BLOCK_INDEX_SHIFT) << "::" << (arc->saddle & VERTEX_INDEX_MASK);
//#endif

        graph.taskStop("terminate"+std::to_string(v), "terminate");
        // check if we can continue sweep at this saddle
        if (this->checkSaddle(arc->saddle)) {
            if (!done.occurred()){
                assert(true);
                assert(arc->saddle != INVALID_VERTEX);
#ifdef ENABLE_DEBUG_LOGGING
                LogDebug().tag(std::to_string(this->index)) << "cont starting next arc  " << (arc->saddle & VERTEX_INDEX_MASK);
#endif
                return true;
            } else {
                return false;
            }
        } else {
#ifdef ENABLE_DEBUG_LOGGING
            LogDebug().tag(std::to_string(this->index)) << "cont Checking saddle false" << (arc->saddle & VERTEX_INDEX_MASK);
#endif
            return false;
        }
    } else {
        std::vector<hpx::id_type> empty;
        arc->saddle = s.vertex;
        this->assignSaddle(arc, s, empty);
#ifdef ENABLE_DEBUG_LOGGING
        LogDebug().tag(std::to_string(this->index)) << "Informed about Arc from " <<  (arc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum & VERTEX_INDEX_MASK) << " to "  << (arc->saddle >> BLOCK_INDEX_SHIFT) << "::" << (arc->saddle & VERTEX_INDEX_MASK);
#endif
        graph.taskStop("terminate"+std::to_string(v), "terminate");
        return false;
    }
}


/**
 * @brief Fetches Arc from Map NOT! threadsafe
 * @param arc
 * @param v
 * @return True if Arc had to be created, False if Arc could be fetched
 */
bool TreeConstructor::fetchCreateArc(Arc*& arc, uint64_t v){
    Profiler profiler(__func__);
    Arc*& tmparc = this->arcMap[v];
    if (tmparc == nullptr){
        tmparc = new Arc(v, this->dataManager, &swept);
        arc = tmparc;
        return true;
    } else {
        arc = tmparc;
        return false;
    }
}

void TreeConstructor::mergeBoundaries(Arc*& arc){
    Profiler profiler(__func__);
    for (uint32_t i = 0; i < arc->body->children.size(); ++i){
        mapLock.lock();
        Arc* childptr = arcMap[arc->body->children[i]];
        mapLock.unlock();
        if (childptr == nullptr)
            continue;
        for (uint32_t j = i+1; j < arc->body->children.size(); ++j){
            mapLock.lock();
            Arc* child2ptr = arcMap[arc->body->children[j]];
            mapLock.unlock();
            if (child2ptr == nullptr)
                continue;
            arc->body->queue.push(childptr->body->boundary.intersect(child2ptr->body->boundary));
        }

        //Is this Threadsafe?
        arc->body->boundary.unite(childptr->body->boundary);
    }
}

AnswerMap TreeConstructor::awaitRemoteAnswers(Arc *&arc){
    AnswerMap answerMap;
    hpx::lcos::local::mutex remoteAnswerLock;

    std::vector<hpx::future<std::vector<RemoteAnswer>>> myAnswers;

#ifndef LAZYREMOTECALLS
    arc->body->remoteCallLock.lock();
    arc->body->remoteCalls.swap(myAnswers);
    arc->body->remoteCallLock.unlock();
#else
    std::vector<RemoteCall> myCalls;
    arc->body->remoteCallLock.lock();
    arc->body->remoteCalls.swap(myCalls);
    arc->body->remoteCallLock.unlock();

    myAnswers.reserve(myCalls.size());
    for (RemoteCall c : myCalls){
        if (c.peer != hpx::naming::id_type(0ul, (hpx::naming::id_type::management_type)0))
            myAnswers.push_back(hpx::async<TreeConstructor::continueSweepPeer_action>(c.peer, arc->extremum, arc->body->children));
        else
            myAnswers.push_back(hpx::async<TreeConstructor::continueSweep_action>(this->treeConstructors[this->dataManager->getBlockIndex(c.neighbor)], arc->extremum, this->dataManager->convertToGlobal(c.neighbor), this->dataManager->convertToGlobal(c.origin), false));
    }
#endif

    hpx::wait_each([&](hpx::future<std::vector<RemoteAnswer>> af) -> void
    {
        std::vector<RemoteAnswer> a = af.get();
        std::lock_guard<hpx::lcos::local::mutex> tmplock(remoteAnswerLock);

        for (RemoteAnswer entry : a){
            if (entry.counter >= answerMap[entry.location].counter){
                answerMap[entry.location] = entry;
            }
        }
    }, myAnswers);

    return std::move(answerMap);
}

/**
 * @brief Starts a new sweep at the given local vertex.
 * @param v
 * @param leaf
 */
void TreeConstructor::startSweep(uint64_t v, bool leaf)
{

#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "Starting Sweep from " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK);
#endif
    //Fetch Arc
    Arc* arc;
    this->mapLock.lock();
    bool created = fetchCreateArc(arc, v);
    if (created) {
        arc->activate(this->dataManager, &swept);
        arc->body->leaf = true;
    } else {

        //arc->activate(this->dataManager, &swept);
    }
    this->mapLock.unlock();

    Boundary& boundary = arc->body->boundary;
    SweepQueue& queue = arc->body->queue;

    //Sweep Startpoint
    this->swept[v] = v;

    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
        graph.taskStop("start"+std::to_string(v), "startsweep");
        return;
    }

//    // Mutual intersections of child boundaries -> this is where sweeps must continue
    mergeBoundaries(arc);

//    for (int i = 0; (i < arc->body->children.size()); i++){
//        std::uint64_t v = arc->body->children[i];
//        arcMap[v]->deactivate();
//    }

    arc->body->augmentation.inherit(arc->body->inheritedAugmentations); // gathered and merged here because no lock required here
    arc->body->augmentation.sweep(this->dataManager->getValue(v)); // also add saddle/local minimum to augmentation
//#ifdef ENABLE_DEBUG_LOGGING
//    LogDebug().tag(std::to_string(this->index)) << "StartSweep Augmenting: " << (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
//                                                 << " to: " << (arc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum & VERTEX_INDEX_MASK);
//#endif

    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
        graph.taskStop("start"+std::to_string(v), "startsweep");
        return;
    }

    //No injection possible, only authority for v is here
    arc->body->state = 1;

    for (int i = 0; (i < 6); i++){
        if (this->dataManager->isGhost(this->dataManager->getNeighbor(v, i)) && (this->swept[this->dataManager->getNeighbor(v, i)] == INVALID_VERTEX)){ 
            arc->body->remoteCallLock.lock();

    //Now we can authorize remote sweeps:
    //Initial Fill of Queue

#ifdef LAZYREMOTECALLS
            arc->body->remoteCalls.push_back(RemoteCall(neighbors[i], v));
#else
            arc->body->remoteCalls.push_back(hpx::async<TreeConstructor::continueSweep_action>(this->treeConstructors[this->dataManager->getBlockIndex(this->dataManager->getNeighbor(v, i))], v, this->dataManager->convertToGlobal(this->dataManager->getNeighbor(v,i)), this->dataManager->convertToGlobal(v), false));
#endif
            arc->body->remoteCallLock.unlock();
        } else {
            queue.push(this->dataManager->getNeighbor(v, i));
        }
    }


#ifdef ENABLE_DEBUG_LOGGING
    LogDebug().tag(std::to_string(this->index)) << "Reminding" << arc->body->peers.size() <<  "peers of " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK);
#endif
    //Call ContinueSweepPeer on all remote portions of v's cc
    arc->body->remoteCallLock.lock();
    for (hpx::id_type peer : arc->body->peers){
        if (peer != this->get_id())
#ifndef LAZYREMOTECALLS
            arc->body->remoteCalls.push_back(hpx::async<TreeConstructor::continueSweepPeer_action>(peer, v, arc->body->children));
#else
            arc->body->remoteCalls.push_back(RemoteCall(0ul, 0ul, peer));
#endif
    }
    arc->body->remoteCallLock.unlock();

    continueLocalSweep(v);
}


void TreeConstructor::continueLocalSweep(uint64_t v){

    Arc* arc;
    this->mapLock.lock();
    if (fetchCreateArc(arc, v)) {
        arc->body->leaf = true;
    }
    this->mapLock.unlock();

    Boundary& boundary = arc->body->boundary;
    SweepQueue& queue = arc->body->queue;

    uint64_t min;


        this->sweepLoop(arc);
        arc->body->lock.lock();
        if (arc->body->queue.empty()){
            arc->body->state = 2;
            if (arc->body->boundary.empty()){
                min = INVALID_VERTEX;
            } else {
                min = boundary.min();
            }
        } else {
            if (done.occurred()){
                arc->body->lock.unlock();
                awaitRemoteAnswers(arc);
                hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
                hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
                graph.taskStop("start"+std::to_string(v), "startsweep");
                return;
            } else {
                arc->body->lock.unlock();
                hpx::apply(TreeConstructor::continueLocalSweep_action(), this->get_id(), v);
                return;
            }
        }


    RemoteAnswer myWinner(this->get_id());
    myWinner.counter = ++arc->body->counter;
    myWinner.minimum = min;
    myWinner.value = this->dataManager->getValue(min);

    arc->body->lock.unlock();

//#ifdef ENABLE_DEBUG_LOGGING
//    LogDebug().tag(std::to_string(this->index)) << "Sweeploop empty for " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK) << " with min: " << (min >> BLOCK_INDEX_SHIFT) << "::" << (min & VERTEX_INDEX_MASK) << " awaiting remotes.";
//#endif
    AnswerMap answerMap = awaitRemoteAnswers(arc);
    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
        graph.taskStop("start"+std::to_string(v), "startsweep");
        return;
    }

    if (myWinner.counter > answerMap[this->get_id()].counter){
        answerMap[this->get_id()] = myWinner;
    }
//#ifdef ENABLE_DEBUG_LOGGING
//    LogDebug().tag(std::to_string(this->index)) << "Remotes of " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK) << " returned!";
//#endif
    RemoteAnswer winner(this->get_id());

    for (auto& entry : answerMap){
        RemoteAnswer a = entry.second;
        if ((a.counter > 0) && (a.value < winner.value)){
            winner = a;
        }
    }

//#ifdef ENABLE_DEBUG_LOGGING
//    LogDebug().tag(std::to_string(this->index)) << "Remotes of " <<  (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK) << " resulted in winner: " << (winner.minimum >> BLOCK_INDEX_SHIFT) << "::" << (winner.minimum & VERTEX_INDEX_MASK);
//#endif
    if (winner.minimum == INVALID_VERTEX){
        //Last sweep
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        this->countSweepsLocal(-1l);
        return;
    } else {
        arc->saddle = winner.minimum;
//#ifdef ENABLE_DEBUG_LOGGING
//        LogDebug().tag(std::to_string(this->index)) << "StartSweep removing from Boundary of: " << (v >> BLOCK_INDEX_SHIFT) << "::" << (v & VERTEX_INDEX_MASK)
//                                                    << " the vertex: " << (arc->saddle >> BLOCK_INDEX_SHIFT) << "::" << (arc->saddle & VERTEX_INDEX_MASK);
//#endif
        boundary.remove(winner.value);
    }

    std::vector<hpx::id_type> peers;
    peers.reserve(answerMap.size());
    for (auto& entry : answerMap){
        peers.push_back(entry.first);
    }
    //peers.push_back(this->get_id());
    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
        return;
    }

    this->assignSaddle(arc, winner.value, peers);

#ifdef FILEOUT
    fileEntry entry;
    entry.id1 = (hpx::async<TreeConstructor::convertToGlobal_action>(this->treeConstructors[arc->extremum >> BLOCK_INDEX_SHIFT], arc->extremum).get()) & VERTEX_INDEX_MASK;
    entry.id2 = (hpx::async<TreeConstructor::convertToGlobal_action>(this->treeConstructors[arc->saddle  >> BLOCK_INDEX_SHIFT], arc->saddle).get()) & VERTEX_INDEX_MASK;
    entry.time = this->clock.elapsed_microseconds();
    entry.type = 1;
    this->arcFileLock.lock();
    this->arcFileEntries.push_back(entry);
    this->arcFileLock.unlock();
#endif

    std::vector<hpx::id_type> empty;
    std::vector<hpx::future<bool>> terminations;
    int winnerLoc;
    for (hpx::id_type peer : peers){
        if (peer == this->get_id())
            continue;
        if (peer == winner.location){
            graph.taskStart("start"+std::to_string(v), "startsweep", "terminate"+std::to_string(v), "terminate");
            terminations.push_back(hpx::async<TreeConstructor::terminateSweep_action>(winner.location, v, winner.value, peers));
            winnerLoc = terminations.size()-1;
        } else {
            graph.taskStart("start"+std::to_string(v), "startsweep", "terminate"+std::to_string(v), "terminate");
            terminations.push_back(hpx::async<TreeConstructor::terminateSweep_action>(peer, v, winner.value, empty));
        }
    }

    hpx::lcos::wait_all(terminations);
#ifdef ENABLE_DEBUG_LOGGING
    LogInfo().tag(std::to_string(this->index)) << "Arc from " <<  (arc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (arc->extremum & VERTEX_INDEX_MASK)
                                              << " to "  << (arc->saddle >> BLOCK_INDEX_SHIFT) << "::" << (arc->saddle & VERTEX_INDEX_MASK);
#endif
    if (winner.location == this->get_id()){
        // check if we can continue sweep at this saddle
        if (this->checkSaddle(arc->saddle)) {
#ifdef ENABLE_DEBUG_LOGGING
            LogDebug().tag(std::to_string(this->index)) << "Checking saddle true" << (arc->saddle & VERTEX_INDEX_MASK);
#endif
            if (!done.occurred()){
                graph.taskStart("start"+std::to_string(v), "startsweep", "start"+std::to_string(arc->saddle), "startsweep");
                hpx::apply(this->executor_sweeps, TreeConstructor::startSweep_action(), winner.location, arc->saddle, false);
            } else {
                hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
                hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
            }
        } else {
#ifdef ENABLE_DEBUG_LOGGING
            LogDebug().tag(std::to_string(this->index)) << "Checking saddle false" << (arc->saddle & VERTEX_INDEX_MASK);
#endif
            if (options.trunkskip) this->countSweepsLocal(-1l);
        }
    } else {
        if (terminations[winnerLoc].get()){
            if (!done.occurred()){
                graph.taskStart("start"+std::to_string(v), "startsweep", "start"+std::to_string(arc->saddle), "startsweep");
                hpx::apply(this->executor_sweeps, TreeConstructor::startSweep_action(), winner.location, arc->saddle, false);
            } else {
                hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
                hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);;
            }
        } else {
            if (options.trunkskip) this->countSweepsLocal(-1l);
        }
    }
    graph.taskStop("start"+std::to_string(v), "startsweep");
    return;
}

/*
 *
 */
bool TreeConstructor::searchUF(uint64_t start, uint64_t goal)
{
    // goal is actively running sweep

    // TODO distributed:
    // - go through chain
    // - if "cache miss" (no entry exists):
    // if locality of required id is here: return false / and one which caused cache miss (which is still running) for others to path-compress
    // else: send query to locality of required id -> wait for answer -> cache result locally

    if (start == INVALID_VERTEX)
        return false;

    uint64_t c = start;
    uint64_t next = this->UF[c];
    if (c == goal || next == goal)
        return true;
    if (next == INVALID_VERTEX)
        return false;

    // includes path compression
    while (true) {
        if (this->UF[next] == goal) {
            this->UF[c] = this->UF[next];
            return true;
        }

        if (this->UF[next] == INVALID_VERTEX)
            return false;

        this->UF[c] = this->UF[next];
        next = this->UF[next];
    }
}

/*
 * Sweep reaches vertex and checks if it can be swept by going through *all* its smaller neighbors and check if they *all* have already been swept by us
 */
bool TreeConstructor::touch(uint64_t c, uint64_t v)
{
    const Value value = this->dataManager->getValue(c);

    uint64_t neighbors[6];
    uint32_t numNeighbors = this->dataManager->getNeighbors(c, neighbors);

    for (uint32_t i = 0; i < numNeighbors; ++i) {
        if (neighbors[i] != INVALID_VERTEX) {
            if (this->dataManager->getValue(neighbors[i]) <= value) {
                if (!this->searchUF(this->swept[neighbors[i]], v)) {
                    return false;
                }
            }
        }
    }
    return true;
}

/*
 * Checks in arc map if saddle is already represented: Either saddle's arc is constructed (we're first level set component to reach this saddle), or we take existing arc (if someone else already created it).
 */
void TreeConstructor::assignSaddle(Arc* finishedArc, Value s, std::vector<hpx::naming::id_type> peers)
{
    //Profiler profiler(__func__);
    this->mapLock.lock();
    Arc* saddleArc = arcMap[finishedArc->saddle]; // parent of finished arc
    if (saddleArc == nullptr) {
        saddleArc = new Arc(finishedArc->saddle, dataManager, &swept);
        arcMap[finishedArc->saddle] = saddleArc;
        saddleArc->activate(this->dataManager, &swept);
    }
    this->mapLock.unlock();

    saddleArc->body->lock.lock();
    saddleArc->body->children.push_back(finishedArc->extremum);
//#ifdef ENABLE_DEBUG_LOGGING
//    LogDebug().tag(std::to_string(this->index)) << "Next Arc " <<  (saddleArc->extremum >> BLOCK_INDEX_SHIFT) << "::" << (saddleArc->extremum & VERTEX_INDEX_MASK)
//                                                << " assignSaddle adds child: " << (finishedArc->extremum  >> BLOCK_INDEX_SHIFT) << "::" << (finishedArc->extremum  & VERTEX_INDEX_MASK);
//#endif

    Augmentation bob = finishedArc->body->augmentation.heritage(s);
    saddleArc->body->inheritedAugmentations.push_back(bob);
#ifdef VTIOUT
#ifndef FLATAUGMENTATION
    for (auto iter = finishedArc->body->augmentation.vertices.begin(); (iter != finishedArc->body->augmentation.vertices.end()); iter.operator ++()){
        swept[(*iter)->key.vertex] = finishedArc->extremum;
    }
#else
    for (Value c : finishedArc->body->augmentation.vertices){
        swept[c.vertex] = finishedArc->extremum;
    }
#endif
#endif
    //    for (auto iter = bob.vertices.begin(); (iter != bob.vertices.end()); iter.operator ++()){
    //    LogDebug().tag(std::to_string(this->index)) << "Inheriting "
    //                                                << ( (*iter)->key.vertex  >> BLOCK_INDEX_SHIFT) << "::" << ( (*iter)->key.vertex  & VERTEX_INDEX_MASK) << " from "
    //                                                << ( finishedArc->extremum  >> BLOCK_INDEX_SHIFT) << "::" << ( finishedArc->extremum  & VERTEX_INDEX_MASK) << " to "
    //                                                << ( saddleArc->extremum  >> BLOCK_INDEX_SHIFT) << "::" << ( saddleArc->extremum  & VERTEX_INDEX_MASK);
    //    }
    //    for (auto iter = finishedArc->body->augmentation.vertices.begin(); (iter != finishedArc->body->augmentation.vertices.end()); iter.operator ++()){
    //    LogDebug().tag(std::to_string(this->index)) << "Not Inheriting "
    //                                                << ( (*iter)->key.vertex  >> BLOCK_INDEX_SHIFT) << "::" << ( (*iter)->key.vertex  & VERTEX_INDEX_MASK) << " from "
    //                                                << ( finishedArc->extremum  >> BLOCK_INDEX_SHIFT) << "::" << ( finishedArc->extremum  & VERTEX_INDEX_MASK) << " to "
    //                                                << ( saddleArc->extremum  >> BLOCK_INDEX_SHIFT) << "::" << ( saddleArc->extremum  & VERTEX_INDEX_MASK);
    //    }

    for (hpx::id_type peer : peers)
        saddleArc->body->peers.insert(peer);

    // Update union find structure: finished arc belongs to saddle arc
    this->UF[finishedArc->extremum] = saddleArc->extremum; // TODO distributed: nothing changes here
    saddleArc->body->lock.unlock();
    //finishedArc->deactivate();
}

/*
 * Checks if we're the last child arc to reach saddle -> only the last one can continue sweeping at this saddle
 */
bool TreeConstructor::checkSaddle(uint64_t s)
{
    //return false;
    Profiler profiler(__func__);
    this->mapLock.lock();
    Arc* branch = this->arcMap[s];
    this->mapLock.unlock();

//    if (branch->body == nullptr)  checkSaddles can arrive wayyy after even the saddle has found its saddle
//        return false;

    // Make sure only one checkSaddle call continues sweep
    std::lock_guard<hpx::lcos::local::mutex> lock(branch->body->lock);

    if (branch->body->started){
        return false;
    }

    const Value value = this->dataManager->getValue(s);

    uint64_t neighbors[6];
    uint32_t numNeighbors = this->dataManager->getNeighbors(s, neighbors);

    // Go through all smaller neighbors
    for (uint32_t i = 0; i < numNeighbors; ++i) {
        if (neighbors[i] != INVALID_VERTEX) {
            if (this->dataManager->getValue(neighbors[i]) <= value) {
                // if not swept or swept by someone else
                if (this->swept[neighbors[i]] == INVALID_VERTEX || (!this->searchUF(this->swept[neighbors[i]], s))){
                    return false;
                }
            }
        }
    }

    branch->body->started = true;
    return true;
}
