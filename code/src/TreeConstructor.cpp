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

    numVertices = this->dataManager->getNumVerticesLocal(true);

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

    hpx::apply(this->executor_low, TreeConstructor::sendSweepCount_action(), this->get_id());
}

void TreeConstructor::destroy()
{
    if (this->dataManager){
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
        uint64_t* data = static_cast<uint64_t*>(imageData->GetScalarPointer());

        int j = 0;

        for (uint64_t i = 0; (i < this->dataManager->getNumVerticesLocal(true)); i++){

            uint64_t v = this->dataManager->getLocalVertex(i);
            if (this->dataManager->isGhost(v))
                continue;

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


    if (danglings.size() > 1){
        for (uint64_t i = 0; (i < danglings.size() - 2); i++){
            if (this->dataManager->isLocal(danglings[i].vertex)){
                Arc* newArc;
                fetchCreateArc(newArc, danglings[i].vertex);
                newArc->saddle = danglings[i+1].vertex;
            }
#ifndef FLATAUGMENTATION
                for (auto iter = arc->body->augmentation.vertices.begin(); (iter != arc->body->augmentation.vertices.end()); iter.operator ++()){
                    swept[(*iter)->key.vertex] = INVALID_VERTEX;
#endif
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
                swept[v] = INVALID_VERTEX;
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

#ifdef VTPOUT
        if (this->index == 0ul){
            vtkSmartPointer<vtkPoints> pointData =
                    vtkSmartPointer<vtkPoints>::New();

            vtkSmartPointer<vtkCellArray> cellData =
                    vtkSmartPointer<vtkCellArray>::New();

            int pointcount = 0;

            for (Arc* arc : this->arcMap.local){
                if (arc != nullptr && arc->saddle != INVALID_VERTEX){
                    glm::uvec3 coords = this->dataManager->getGlobalCoords(arc->extremum);
                    pointData->InsertNextPoint(coords.x, coords.y, coords.z);
                    glm::uvec3 coords2 = this->dataManager->getGlobalCoords(arc->saddle);
                    pointData->InsertNextPoint(coords2.x, coords2.y, coords2.z);
                    vtkSmartPointer<vtkLine> currentArc =
                            vtkSmartPointer<vtkLine>::New();
                    currentArc->GetPointIds()->SetId(0, pointcount++);
                    currentArc->GetPointIds()->SetId(1, pointcount++);
                    cellData->InsertNextCell(currentArc);
                }
            }

            for (auto entry : this->arcMap.remote){
                if (entry.second != nullptr){
                    Arc* arc = entry.second;
                    glm::uvec3 coords = this->dataManager->getGlobalCoords(this->dataManager->convertToGlobal(arc->extremum));
                    pointData->InsertNextPoint(coords.x, coords.y, coords.z);
                    coords = this->dataManager->getGlobalCoords(this->dataManager->convertToGlobal(arc->saddle));
                    pointData->InsertNextPoint(coords.x, coords.y, coords.z);
                    vtkSmartPointer<vtkLine> currentArc =
                            vtkSmartPointer<vtkLine>::New();
                    currentArc->GetPointIds()->SetId(0, pointcount++);
                    currentArc->GetPointIds()->SetId(1, pointcount++);
                    cellData->InsertNextCell(currentArc);
                }
            }


            vtkSmartPointer<vtkPolyData> arcData =
                    vtkSmartPointer<vtkPolyData>::New();
            arcData->SetPoints(pointData);
            arcData->SetLines(cellData);

            vtkSmartPointer<vtkXMLPolyDataWriter> arcwriter =
              vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            std::string filename = "output.vtp";
            arcwriter->SetFileName(filename.c_str());
          #if VTK_MAJOR_VERSION <= 5
            writer->SetInput(unstructuredGrid);
          #else
            arcwriter->SetInputData(arcData);
          #endif
            arcwriter->Write();
        }
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
}

void TreeConstructor::countSweeps(std::int64_t diff){
    if (diff == std::numeric_limits<std::int64_t>::max()){
        if (std::atomic_fetch_add(&this->remoteLeavesDone, 1) == (int)(this->treeConstructors.size()-1)){
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
        if (std::atomic_fetch_add(&this->remoteLeavesDone, 1) == (int)(this->treeConstructors.size()-1)){
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

uint64_t TreeConstructor::construct()
{
    Profiler profiler(__func__);
    hpx::util::high_resolution_timer timer;

    // Find minima in local data
    std::vector<uint64_t> localVertices = this->dataManager->getLocalVertices();

    for (uint64_t v : localVertices) {
        if (this->dataManager->isMinimum(v)) {
            if (!this->dataManager->isGhost(v)){
                hpx::apply(this->executor_start_sweeps, TreeConstructor::startSweep_action(), this->get_id(), v, true);
                ++numMinima;
            }
        }
    }

    if (numMinima == 0l)
        numMinima = std::numeric_limits<std::int64_t>::max();

    if (!options.trunkskip){
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], std::numeric_limits<std::int64_t>::max());
    } else {
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], numMinima);
    }

    this->termination.wait();

    if (this->index == 0){
        this->postProc();
        this->postProcessing.wait();

        int numArcs = 0;
        int numArcss[hpx::get_num_worker_threads()];
        std::vector<std::vector<Value>> danglingArcss;
        danglingArcss.resize(hpx::get_num_worker_threads());
        for (uint64_t i = 0; (i < hpx::get_num_worker_threads()); i++){
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

        for (uint64_t i = 0; (i < hpx::get_num_worker_threads()); i++){
            numArcs += numArcss[i];
            danglingArcs.insert(danglingArcs.end(), danglingArcss[i].begin(), danglingArcss[i].end());
        }
        numArcs += danglingArcs.size();
        hpx::parallel::sort(hpx::parallel::execution::parallel_unsequenced_policy(), this->danglingArcs.begin(), this->danglingArcs.end(),
                            [&](const Value a1, const Value a2) {
                                return a1 < a2;
                            });

        for (uint64_t i = 1ul; (i < this->treeConstructors.size()); i++){
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
        for (uint64_t i = 0; (i < hpx::get_num_worker_threads()); i++){
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

        for (uint64_t i = 0; (i < hpx::get_num_worker_threads()); i++){
            numArcs += numArcss[i];
            registered.push_back(hpx::async(this->executor_high, TreeConstructor::registerArc_action(), this->treeConstructors[0], danglingArcs[i]));
        }

        hpx::wait_all(registered);

        LogInfo().tag(std::to_string(this->index))<< "Arcs: " << numArcs;

        hpx::async<TreeConstructor::postProc_action>(this->treeConstructors[0]);

        this->postProcessing.wait();

        return numArcs;
    }
}

void TreeConstructor::sweepLoop(Arc* arc){
    SweepQueue& queue = arc->body->queue;
    Boundary& boundary = arc->body->boundary;
    uint64_t& v = arc->extremum;

    while (!queue.empty()) {
        uint64_t c = queue.pop();
        if (c == INVALID_VERTEX)
            break;

        const Value value = this->dataManager->getValue(c);

        if (this->dataManager->isGhost(c)){

            boundary.remove(value); // remove from boundary

            swept[c] = v;           //mark as swept


            uint64_t neighbors[6];
            uint32_t numNeighbors = this->dataManager->getNeighbors(c, neighbors);

            for (uint64_t i = 0; (i < numNeighbors); i++){
                if ((!this->dataManager->isGhost(neighbors[i])) && (this->swept[neighbors[i]] == INVALID_VERTEX)){
                    queue.push(neighbors[i]);
                }
            }
            continue;
        }

        // check if can be swept
        if (this->touch(c, v)) {
            boundary.remove(value); // remove from boundary
            this->swept[c] = v;     // put into our augmentation and mark as swept
            arc->body->augmentation.sweep(value);

            uint64_t neighbors[6];
            uint32_t numNeighbors = this->dataManager->getNeighbors(c, neighbors);

            for (uint64_t i = 0; (i < numNeighbors); i++){
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
            boundary.add(value);
        }
    }
}

std::vector<RemoteAnswer> TreeConstructor::continueSweep(uint64_t v, uint64_t c, uint64_t from, bool leaf)
{
    uint64_t sent = dataManager->convertToLocal(c);
    uint64_t swept_from = dataManager->convertToLocal(from);

    bool leaveItToContinuePeer = false;

    //Fetch Arc
    Arc* arc;
    this->mapLock.lock();
    if (fetchCreateArc(arc, v)){
        this->mapLock.unlock();
    } else {
        this->mapLock.unlock();
        arc->body->lock.lock();
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
        arc->body->lock.unlock();
        RemoteAnswer answer(this->get_id());
        std::vector<RemoteAnswer> answerMap;
        answerMap.push_back(answer);
        return answerMap;
    }


    arc->body->queue.push(sent);
    //No boundary or augmentation initialization needed

    //Inject in possibly running sweep
    if (arc->body->state == 1) {
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

    if ((arc->saddle == INVALID_VERTEX || this->dataManager->getValue(arc->saddle) > this->dataManager->getValue(min)) || swept[arc->saddle] == v)
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

    arc->body->children = children;

    //Sweep Startpoint
    this->swept[v] = v;

    //No Initial Fill except boundary intersections
    for (uint64_t child : children){
        this->swept[child] = v;
        this->UF[child] = v;
    }

    // Mutual intersections of child boundaries -> this is where sweeps must continue
    mergeBoundaries(arc);

    arc->body->augmentation.inherit(arc->body->inheritedAugmentations);

    //Inject in possibly running sweep
    arc->body->lock.lock();
    if (arc->body->state == 1){
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
        arc->body->boundary.remove(s);
        this->assignSaddle(arc, s, peers);

        // check if we can continue sweep at this saddle
        if (this->checkSaddle(arc->saddle)) {
            if (!done.occurred()){
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        std::vector<hpx::id_type> empty;
        arc->saddle = s.vertex;
        this->assignSaddle(arc, s, empty);
        return false;
    }
}

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

    return answerMap;
}

/**
 * @brief Starts a new sweep at the given local vertex.
 * @param v
 * @param leaf
 */
void TreeConstructor::startSweep(uint64_t v, bool leaf)
{
    //Fetch Arc
    Arc* arc;
    this->mapLock.lock();
    bool created = fetchCreateArc(arc, v);
    if (created) {
        arc->activate(this->dataManager, &swept);
        arc->body->leaf = true;
    }
    this->mapLock.unlock();

    Boundary& boundary = arc->body->boundary;
    SweepQueue& queue = arc->body->queue;

    //Sweep Startpoint
    this->swept[v] = v;

    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
        return;
    }

    // Mutual intersections of child boundaries -> this is where sweeps must continue
    mergeBoundaries(arc);

    arc->body->augmentation.inherit(arc->body->inheritedAugmentations); // gathered and merged here because no lock required here
    arc->body->augmentation.sweep(this->dataManager->getValue(v)); // also add saddle/local minimum to augmentation

    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
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

    AnswerMap answerMap = awaitRemoteAnswers(arc);
    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
        return;
    }

    if (myWinner.counter > answerMap[this->get_id()].counter){
        answerMap[this->get_id()] = myWinner;
    }

    RemoteAnswer winner(this->get_id());

    for (auto& entry : answerMap){
        RemoteAnswer a = entry.second;
        if ((a.counter > 0) && (a.value < winner.value)){
            winner = a;
        }
    }

    if (winner.minimum == INVALID_VERTEX){
        //Last sweep
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        this->countSweepsLocal(-1l);
        return;
    } else {
        arc->saddle = winner.minimum;
        boundary.remove(winner.value);
    }

    std::vector<hpx::id_type> peers;
    peers.reserve(answerMap.size());
    for (auto& entry : answerMap){
        peers.push_back(entry.first);
    }

    if (done.occurred()){
        hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
        hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
        return;
    }

    this->assignSaddle(arc, winner.value, peers);

    std::vector<hpx::id_type> empty;
    std::vector<hpx::future<bool>> terminations;
    int winnerLoc;
    for (hpx::id_type peer : peers){
        if (peer == this->get_id())
            continue;
        if (peer == winner.location){
            terminations.push_back(hpx::async<TreeConstructor::terminateSweep_action>(winner.location, v, winner.value, peers));
            winnerLoc = terminations.size()-1;
        } else {
            terminations.push_back(hpx::async<TreeConstructor::terminateSweep_action>(peer, v, winner.value, empty));
        }
    }

    hpx::lcos::wait_all(terminations);
    if (winner.location == this->get_id()){
        // check if we can continue sweep at this saddle
        if (this->checkSaddle(arc->saddle)) {
            if (!done.occurred()){
                hpx::apply(this->executor_sweeps, TreeConstructor::startSweep_action(), winner.location, arc->saddle, false);
            } else {
                hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
                hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);
            }
        } else {
            if (options.trunkskip) this->countSweepsLocal(-1l);
        }
    } else {
        if (terminations[winnerLoc].get()){
            if (!done.occurred()){
                hpx::apply(this->executor_sweeps, TreeConstructor::startSweep_action(), winner.location, arc->saddle, false);
            } else {
                hpx::apply(this->executor_high, TreeConstructor::startTrunk_action(), this->treeConstructors[0], this->dataManager->getValue(arc->extremum));
                hpx::apply(this->executor_high, TreeConstructor::countSweeps_action(), this->treeConstructors[0], -1l);;
            }
        } else {
            if (options.trunkskip) this->countSweepsLocal(-1l);
        }
    }
    return;
}

/*
 *
 */
bool TreeConstructor::searchUF(uint64_t start, uint64_t goal)
{
    // goal is actively running sweep

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

    Augmentation heritage = finishedArc->body->augmentation.heritage(s);
    saddleArc->body->inheritedAugmentations.push_back(heritage);
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

    for (hpx::id_type peer : peers)
        saddleArc->body->peers.insert(peer);

    // Update union find structure: finished arc belongs to saddle arc
    this->UF[finishedArc->extremum] = saddleArc->extremum;
    saddleArc->body->lock.unlock();
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
