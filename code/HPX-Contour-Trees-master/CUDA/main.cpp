#define HPX_HPX_MAIN_IMPL_HPP

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <cstring>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <boost/heap/fibonacci_heap.hpp>

#include <stack>

#include <glm/glm.hpp>

#include <omp.h>

#include "DataManager.h"
#include "Visualizer.h"
#include "TreeConstructor.h"

#include <cuda_profiler_api.h>

vtkStandardNewMacro(MouseInteractorStyle2) ;

struct DataValueComparator {

    static Datamanager* data;

    //    static void init(Datamanager* d)
    //    {
    //        data = d;
    //    }

    static bool compare(const int& n1, const int& n2)
    {
        return data->getValue(n1) < data->getValue(n2);
    }
};

Datamanager* DataValueComparator::data = nullptr;

struct compare_node
{
    bool operator()(const ctValue& n1, const ctValue& n2) const
    {
       return n1 > n2;
    }
};

bool searchUF(int start, int goal, int* swept){
    if (start < 0)
        return false;
    int c = start;
    int next = swept[c];
    if (c == goal || next == goal)
        return true;
    if (next < 0)
        return false;

    while (true){
        if (swept[next] == goal){
            //swept[c] = swept[next];
            return true;
        }

        if (swept[next] < 0 || swept[next] == next)
            return false;

        //swept[c] = swept[next];
        next = swept[next];
    }
}

int main(int argc, char* argv[])
{
    boost::program_options::options_description desc_commandline("Usage: [options]");

    // clang-format off
    desc_commandline.add_options()
            ("resample", boost::program_options::value<std::string>()->default_value("100"), "Resamples the input file with the given axis scale percentages (format: X,Y,Z or XxYxZ or N).")
            ("input", boost::program_options::value<std::string>()->default_value("ctBones.vti"), "Path to input data that should be used.")
            ("rounds", boost::program_options::value<int>()->default_value(3), "How many rounds of sweeps are performend on the gpu before passing back to cpu")
            ("minimarounds", boost::program_options::value<int>()->default_value(200), "How many steps are local minima permitted before interruption for next round")
            ("growthrounds", boost::program_options::value<int>()->default_value(100), "How many steps are sweeps permitted to do before interruption for next round")
            ("visualize", "Resulting tree embedded in input data will be rendered with vtk")
            ("no-trunkskip", "Perform explicit trunk computation and instead of collecting dangling saddles.");

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc_commandline), vm);
    boost::program_options::notify(vm);

    Options opt;
    opt.visualize = false;
    opt.trunkSkipping = true;

    if (vm.count("visualize"))
        opt.visualize = true;

    if (vm.count("no-trunkskip"))
        opt.trunkSkipping = false;

    opt.inputFile = vm["input"].as<std::string>();

    glm::ivec3 resample;

    std::vector<std::string> resampleComponents;
    boost::split(resampleComponents, vm["resample"].as<std::string>(), boost::is_any_of(",x"));

    if (resampleComponents.size() == 3)
        resample = glm::ivec3(std::stoi(resampleComponents[0]), std::stoi(resampleComponents[1]), std::stoi(resampleComponents[2]));
    else if (resampleComponents.size() == 1)
        resample = glm::ivec3(std::stoi(resampleComponents[0]));
    else
        throw std::runtime_error("Invalid resample size");

    opt.resampleX = resample.x;
    opt.resampleY = resample.y;
    opt.resampleZ = resample.z;

    TreeConstructor tree;
    Datamanager data = Datamanager(opt.resampleX, opt.resampleY, opt.resampleZ);
    Visualizer vis;
    vis.init(data);

    data.readFromFile(opt.inputFile);

    DataValueComparator::data = &data;



    //float scalars[data.getNumVertices()];
    //if (data.mesh != nullptr)
    //{
    //    for (int i = 0; (i < data.getNumVertices()); i++){
    //        scalars[i] = data.getValue(i).value;
    //    }
    //    tree.fvalues = scalars;
    //} else {
        tree.values = data.rawValues();
        tree.fvalues = data.fvalues;
    //}

    tree.numVertices = data.getNumVertices();
    tree.neighborsBuffer = data.rawNeighbors();
    tree.neighborsMap = data.neighborMap();
    tree.min_val = data.min_value;
    tree.conversionFactor = data.conversionFactor;
    tree.max_x = data.dimensions[0];
    tree.max_y = data.dimensions[1];
    tree.max_z = data.dimensions[2];

    std::vector<int*> minima;
    std::vector<int*> saddles;
    std::vector<int>  counts;
    int* swept;
    int* augmentation;
    //cudaProfilerStart();
    std::chrono::high_resolution_clock::time_point startTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_secs_total;

    doCudaStuff(tree, minima, saddles, counts, swept, augmentation, vm["rounds"].as<int>(), vm["minimarounds"].as<int>(), vm["growthrounds"].as<int>());

    //cudaProfilerStop();
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_secs;

    std::stack<int> tasks;
    std::set<int> dangling;
    std::vector<boost::heap::fibonacci_heap<ctValue, boost::heap::compare<compare_node>>*> heaps;
    std::vector<omp_lock_t> heapLocks;
    heaps.resize(tree.numVertices, nullptr);
    heapLocks.resize(tree.numVertices);
    omp_lock_t taskslock;
    omp_init_lock(&taskslock);

    for (int i = 0; (i < counts[counts.size()-1]); i++){
        int start = minima[minima.size()-1][i];
        tasks.push(start);
        //vis.createBoundaryActor(start, 0, 1, 0, 0.5);
        heaps[start] = new boost::heap::fibonacci_heap<ctValue, boost::heap::compare<compare_node>>();
        heaps[start]->push(data.getValue(start));
        //swept[start] = start;
    }

#pragma omp parallel for
    for (int j = 0; (j < tree.numVertices); j++){
        omp_init_lock(&heapLocks[j]);
    }

#pragma omp parallel for
    for (int j = 0; (j < tree.numVertices); j++){
        if (swept[j] == -1){

            ctValue cValue = data.getValue(j);

            int numNeighbors;
            const int* neighbors = data.getNeighbors(j, &numNeighbors);
            for (int i = 0; i < numNeighbors; ++i) {
                if (data.getValue(neighbors[i]) < cValue) {
                    if (swept[neighbors[i]] >= 0){
                        int loc = swept[neighbors[i]];
                        omp_set_lock(&heapLocks[loc]);
                        boost::heap::fibonacci_heap<ctValue, boost::heap::compare<compare_node>>* insertheap = heaps[loc];
                        if (insertheap != nullptr){
                            insertheap->push(data.getValue(j));
                            //omp_set_lock(&taskslock);
                            //vis.createBoundaryActor(j, 1, 1, 0, 0.5);
                            //omp_unset_lock(&taskslock);
                        } else {
                            omp_set_lock(&taskslock);
                            dangling.insert(loc);
                            //vis.createBoundaryActor(loc, 1, 1, 1, 0.5);
                            omp_unset_lock(&taskslock);
                            heaps[loc] = new boost::heap::fibonacci_heap<ctValue, boost::heap::compare<compare_node>>();
                            heaps[loc]->push(data.getValue(j));
                        }
                        omp_unset_lock(&heapLocks[loc]);
                    }
                }
            }
        }
    }

    elapsed_secs = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Extract Fibonacci Heaps Time: " << elapsed_secs.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();

/*
    for (int i = 0; (i < tree.numVertices); i++){
        if (swept[i] >= 0){
        int value = (data.getValue(swept[i]).value - data.min_value)/(data.max_value-data.min_value)*1000000;
        value = swept[i];
        vis.createBoundaryActor(i, (value % 5)/4.0, (value % 2)/1.0, (value % 11)/10.0, 1);
        }
    }
*/


    std::vector<bool> saddleStarted;
    saddleStarted.resize(tree.numVertices, false);

#pragma omp parallel
{
    omp_set_lock(&taskslock);
    while (tasks.size() > 1){


        int v = tasks.top();
        tasks.pop();

        swept[v] = v;

        omp_unset_lock(&taskslock);

        ctValue next;
        while (!heaps[v]->empty()){
            next = heaps[v]->top();
            heaps[v]->pop();

            if ((swept[next.key] == v) && (next.key != v))
                continue;

            bool mine = true;

            int numNeighbors;
            const int* neighbors = data.getNeighbors(next.key, &numNeighbors);
            for (int i = 0; i < numNeighbors; ++i) {
                if (data.getValue(neighbors[i]) < next) {
                    if (!searchUF(neighbors[i], v, swept))
                        mine = false;
                }
            }

            if (mine){
                swept[next.key] = v;
                for (int i = 0; i < numNeighbors; i++) {
                    if (swept[neighbors[i]] == -1){
                        heaps[v]->push(data.getValue(neighbors[i]));
                    }
                }
            } else {
                swept[next.key] = -1;
                swept[v] = next.key; //union
                if (heaps[next.key] == nullptr){
                    heaps[next.key] = heaps[v];
                } else {
                    heaps[next.key]->merge(*heaps[v]);
                }
                bool checksaddle = true;

                omp_set_lock(&taskslock);

                if (!saddleStarted[next.key]){
                    int numNeighbors;
                    const int* neighbors = data.getNeighbors(next.key, &numNeighbors);
                    for (int i = 0; i < numNeighbors; ++i) {
                        if (data.getValue(neighbors[i]) < next) {
                            if (!searchUF(neighbors[i], next.key, swept)){
                                checksaddle = false;
                            }
                        }
                    }
                    if (checksaddle){
                        auto it = dangling.find(next.key);
                        if (it != dangling.end())
                            dangling.erase(it);
                        saddleStarted[next.key] = true;
                        tasks.push(next.key);
                    } else {
                        dangling.insert(next.key);
                    }
                }
                omp_unset_lock(&taskslock);
                break;
            }
        }
        omp_set_lock(&taskslock);
    }
    omp_unset_lock(&taskslock);
}

    //std::cout << "Minima:" << counts[0] << std::endl;

#pragma omp parallel for
    for (int i = 0; (i < tree.numVertices); i++){
        if (augmentation[i] >= 0){
            swept[i] = augmentation[i];
        }
    }

    std::vector<int> danglingvect;
    //vis.createBoundaryActor(tasks.top(), 0, 1, 0, 1);
    for (int v : dangling){
            danglingvect.push_back(v);
            //vis.createBoundaryActor(v, 0, 0, 1, 1);
    }

    std::sort(danglingvect.begin(), danglingvect.end(), DataValueComparator::compare);

    swept[tasks.top()] = danglingvect[0];

#pragma omp parallel for
    for (int i = 0; (i < danglingvect.size()-1); i++){
        swept[danglingvect[i]] = danglingvect[i+1];
    }

    augmentTrunk(tree, tasks.top(), danglingvect, swept);

    elapsed_secs = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Host Time: " << elapsed_secs.count() << std::endl;

    if (opt.visualize){/*
        for (int i = 0; (i < tree.numVertices); i++){

            int value = (data.getValue(swept[i]).value - data.min_value)/(data.max_value-data.min_value)*10000;
            vis.createBoundaryActor(i, (value % 5)/5.0, (value % 2)/2.0, (value % 11)/11.0, 1);
        }*/

        for (int i = 0; (i < minima.size()-1); i++){
            std::cout << "Round " << i << " Size " << counts[i] << std::endl;
            for (int j = 0; (j < counts[i]); j++){
                if (minima[i][j] != saddles[i][j])
                vis.visualizeBranch(minima[i][j], saddles[i][j], (i == minima.size()-1));
            }
        }

        for (int i = 0; (i < counts[counts.size()-1]); i++){
            int start = minima[minima.size()-1][i];
            while (swept[start] >= 0){
                vis.visualizeBranch(start, swept[start], true);
                start = swept[start];
                if(start == swept[start])
                    break;
            }
        }

        vis.show();
    } else {
        for (int i = 0; (i < minima.size()); i++){
            std::cout << "Round " << i << " Size " << counts[i] << std::endl;
        }
    }

    elapsed_secs_total = std::chrono::high_resolution_clock::now() - startTotal;
    std::cout << "Total Time: " << elapsed_secs_total.count() << std::endl;

    return 0;
}
