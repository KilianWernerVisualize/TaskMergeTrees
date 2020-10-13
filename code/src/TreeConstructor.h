#pragma once

class TreeConstructor;

#include "Arc.h"
#include "DataManager.h"
#include "Log.h"

#include <hpx/hpx.hpp>
#include <hpx/include/serialization.hpp>
#include <hpx/parallel/executors.hpp>
#include <hpx/runtime/threads/executors/limiting_executor.hpp>

#include <string>
#include <vector>
#include <map>
#include <boost/container/flat_map.hpp>

#ifdef VTIOUT
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#endif

#ifdef VTPOUT
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#endif

class Options {
public:
    bool trunkskip;

private:
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & trunkskip;
    }

};

class TreeConstructor : public hpx::components::component_base<TreeConstructor> {
    friend class SweepQueue;

public:

    TreeConstructor() : executor_low("default", hpx::threads::thread_priority_low, hpx::threads::thread_stacksize_default),
    executor_start_sweeps(hpx::threads::executors::pool_executor("default", hpx::threads::thread_priority_low, hpx::threads::thread_stacksize_small), 1000, 2000),
    executor_sweeps("default", hpx::threads::thread_priority_default, hpx::threads::thread_stacksize_small),
    executor_high("default", hpx::threads::thread_priority_high, hpx::threads::thread_stacksize_small)
    {

    }

    void init(const std::vector<hpx::id_type>& treeConstructors, const std::string& input, const Options& options);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, init);

    uint64_t construct();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, construct);

    void destroy();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, destroy);

    void startSweep(uint64_t v, bool leaf);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, startSweep);

    bool terminateSweep(uint64_t v, Value s, std::vector<hpx::id_type> peers);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, terminateSweep);

    void continueLocalSweep(uint64_t v);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, continueLocalSweep);

    std::vector<RemoteAnswer> continueSweep(uint64_t v, uint64_t c, uint64_t from, bool leaf);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, continueSweep);

    std::vector<RemoteAnswer> continueSweepPeer(uint64_t v, const std::vector<std::uint64_t> children);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, continueSweepPeer);

    void sendSweepCount();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, sendSweepCount);

    void countSweeps(std::int64_t diff);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, countSweeps);

    void registerArc(std::vector<Value> extrema);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, registerArc);

    void finish();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, finish);

    void terminate();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, terminate);

    void postProc();
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, postProc);

    void stitch(std::vector<Value> danglings);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, stitch);

    void startTrunk(Value v);
    HPX_DEFINE_COMPONENT_ACTION(TreeConstructor, startTrunk);

private:
    void countSweepsLocal(int64_t diff);

    bool searchUF(uint64_t start, uint64_t goal);

    bool touch(uint64_t c, uint64_t v);

    void assignSaddle(Arc* finishedArc, Value s, std::vector<hpx::id_type> peers);

    void sweepLoop(Arc* arc);

    bool checkSaddle(uint64_t s);

    bool fetchCreateArc(Arc*& arc, uint64_t v);

    void mergeBoundaries(Arc*& arc);

    AnswerMap awaitRemoteAnswers(Arc*& arc);

private:
    uint32_t index;

    Options options;

    std::vector<hpx::id_type> treeConstructors;

    DataManager* dataManager = nullptr;

    std::atomic<int64_t> sweepsLocal;
    int64_t sweepsLocal_last;
    hpx::lcos::local::mutex countLock;

    std::atomic<int64_t> sweeps;
    std::atomic<int> remoteLeavesDone;

    hpx::lcos::local::event termination;
    hpx::lcos::local::event postProcessing;
    hpx::lcos::local::event done;
    hpx::lcos::local::event leavesDone;
    hpx::lcos::local::event trunkStarted;

    uint64_t scanCount;
    uint64_t numVertices;
    int64_t numMinima;
    hpx::lcos::local::mutex numMinimaLock;

    hpx::threads::executors::pool_executor executor_low;
    hpx::threads::executors::limiting_executor<hpx::threads::executors::pool_executor> executor_start_sweeps;
    hpx::threads::executors::pool_executor executor_sweeps;
    hpx::threads::executors::pool_executor executor_high;

    Value trunkStart;

    // Map for each vertex to ID of arc extremum
    DistVec<uint64_t> swept; // what vertex has been swept by which saddle/local minimum

    // Union-find-structure containing child-parent relations
    DistVec<uint64_t> UF;

    DistVec<Arc*> arcMap; // (map: starting minima/saddle -> arc)

    std::vector<Value> danglingArcs;

    hpx::lcos::local::mutex mapLock;
};

HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::init_action, treeConstructor_init_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::construct_action, treeConstructor_construct_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::destroy_action, treeConstructor_destroy_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::startSweep_action, treeConstructor_startSweep_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::continueLocalSweep_action, treeConstructor_continueLocalSweep_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::terminateSweep_action, treeConstructor_terminateSweep_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::continueSweep_action, treeConstructor_continueSweep_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::continueSweepPeer_action, treeConstructor_continueSweepPeer_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::sendSweepCount_action, treeConstructor_sendSweepCount_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::countSweeps_action, treeConstructor_countSweeps_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::finish_action, treeConstructor_finish_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::terminate_action, treeConstructor_terminate_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::registerArc_action, treeConstructor_registerArc_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::postProc_action, treeConstructor_postProc_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::stitch_action, treeConstructor_stitch_action);
HPX_REGISTER_ACTION_DECLARATION(TreeConstructor::startTrunk_action, treeConstructor_startTrunk_action);
