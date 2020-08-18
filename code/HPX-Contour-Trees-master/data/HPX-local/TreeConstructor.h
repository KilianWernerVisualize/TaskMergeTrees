#ifndef TREECONSTRUCTOR_H
#define TREECONSTRUCTOR_H

#include "Arc.h"
#include "Visualizer.h"
#include <hpx/hpx.hpp>
#include <map>
#include <vector>

class Options {
public:
    bool visualize;
    bool trunkSkipping;
    std::string inputFile;
    int resampleX;
    int resampleY;
    int resampleZ;

private:
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & this->visualize & this->resampleX & this->resampleY & this->resampleZ;
    }
};

namespace server {

class TreeConstructor : public hpx::components::component_base<TreeConstructor> {

public:
    void init(const Options& opt);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, init);

    void construct();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, construct);

    void startSweep(int v, bool leaf);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, startSweep);

private:
    bool searchUF(int start, int goal);

    bool touch(int c, int v);

    void assignSaddle(Arc* finishedArc, bool trunk);

    bool checkSaddle(int s);

    bool isMinimum(int v) const;

    Arc* trunkStart;

public:
    Visualizer vis;
    Datamanager data;

    bool visualize;
    bool trunkSkipping;

    std::atomic<int> sweeps;
    hpx::lcos::local::event leavesDone;

    // Map for each vertex to ID of arc extremum
    std::vector<int> swept; // un-fixed augmentation (who has been swept by which saddle/local minimum)

    // Union-find-structure containing child-parent relations
    std::vector<int> UF; // TODO???

    // Map for each arc extremum to its arc
    std::set<ctValue> pendingArcList;

    std::vector<Arc*> arcMap; // -> hash map in distributed case because of global vertex indexing
    hpx::lcos::local::mutex mapLock; // locally: replace by vector<mutex> -> mutex for each arc (don't have to lock complete "map" because vector)

    hpx::lcos::local::event done;

    vtkIdType numVertices;
};
}

HPX_REGISTER_ACTION_DECLARATION(::server::TreeConstructor::init_action, treeConstructor_init_action);
HPX_REGISTER_ACTION_DECLARATION(::server::TreeConstructor::construct_action, treeConstructor_construct_action);
HPX_REGISTER_ACTION_DECLARATION(::server::TreeConstructor::startSweep_action, treeConstructor_startSweeps_action);

#endif
