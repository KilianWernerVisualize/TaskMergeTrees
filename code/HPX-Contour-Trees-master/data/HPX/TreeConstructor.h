#pragma once

#include "Arc.h"
#include "DataManager.h"

#include <hpx/hpx.hpp>

#include <string>
#include <vector>

class TreeConstructor : public hpx::components::component_base<TreeConstructor> {
    friend class SweepQueue;

public:
    void init(const std::vector<hpx::id_type>& treeConstructors, const std::string& input);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, init);

    void construct();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, construct);

    void destroy();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, destroy);

    void startSweep(uint64_t v, bool leaf);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, startSweep);

private:
    bool searchUF(uint64_t start, uint64_t goal);

    bool touch(uint64_t c, uint64_t v);

    void assignSaddle(Arc* finishedArc);

    bool checkSaddle(uint64_t s);

private:
    uint32_t index;
    std::vector<hpx::id_type> treeConstructors;

    DataManager* dataManager = nullptr;

    std::atomic<uint32_t> sweeps;
    hpx::lcos::local::event done;
    hpx::lcos::local::event leavesDone;

    // Map for each vertex to ID of arc extremum
    std::vector<uint64_t> swept; // who has been swept by which saddle/local minimum (inverse un-fixed augmentation)

    // Union-find-structure containing child-parent relations
    std::vector<uint64_t> UF; // TODO distributed:  dynamic

    std::vector<Arc*> arcMap; // -> hash map in distributed case because of global vertex indexing      (map: minima/saddle where started -> arc)
    hpx::lcos::local::mutex mapLock; // locally: replace by vector<mutex> -> mutex for each arc (don't have to lock complete "map" because vector)
};

HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::init_action, treeConstructor_init_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::construct_action, treeConstructor_construct_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::destroy_action, treeConstructor_destroy_action);
