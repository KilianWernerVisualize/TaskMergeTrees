#pragma once

#include "Augmentation.h"
#include "Boundary.h"
#include "SweepQueue.h"
#include "DistVec.h"

#include <limits>
#include <vector>
#include <boost/container/flat_set.hpp>

#include <hpx/hpx.hpp>
#include <hpx/include/serialization.hpp>

#ifdef LAZYREMOTECALLS
struct RemoteCall {
    hpx::id_type peer;
    std::uint64_t neighbor;
    std::uint64_t origin;

    RemoteCall(std::uint64_t neighbor, std::uint64_t origin, hpx::id_type peer = hpx::naming::id_type(0ul, (hpx::naming::id_type::management_type)0)) : neighbor(neighbor), origin(origin), peer(peer){
    }
};
#endif

class RemoteAnswer {
public:
    RemoteAnswer() : minimum(INVALID_VERTEX), value(std::numeric_limits<Value>::max()), counter(0){

    }

    RemoteAnswer(hpx::id_type location) : location(location), minimum(INVALID_VERTEX), value(std::numeric_limits<Value>::max()), counter(0){

    }

    hpx::id_type location;
    uint64_t minimum;
    int counter;
    Value value;

private:
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & location & minimum & value & counter;
    }

};

typedef boost::container::flat_map<hpx::id_type, RemoteAnswer> AnswerMap;

class ArcBody {

public:

    ArcBody(uint64_t extremum, DataManager* data, DistVec<std::uint64_t>* swept) : queue(extremum, data, swept){

    }

    bool leaf = false;
    bool started = false;

    // 0 = not started yet
    // 1 = active
    // 2 = finalizing
    // 3 = inactive
    int state = 0;
    int counter = 0;


    Boundary boundary;
    //Augmentation augmentation;

#ifndef LAZYREMOTECALLS
    std::vector<hpx::future<std::vector<RemoteAnswer>>> remoteCalls = std::vector<hpx::future<std::vector<RemoteAnswer>>>();
#else
    std::vector<RemoteCall> remoteCalls = std::vector<RemoteCall>();
#endif
    std::vector<std::uint64_t> children = std::vector<std::uint64_t>();
    boost::container::flat_set<hpx::id_type> peers = boost::container::flat_set<hpx::id_type>();
    //std::vector<Augmentation> inheritedAugmentations = std::vector<Augmentation>();
    SweepQueue queue;

    hpx::lcos::local::mutex remoteCallLock;
    hpx::lcos::local::mutex lock;
};

class Arc {
public:

    Arc(uint64_t extremum, DataManager* data, DistVec<std::uint64_t>* swept) : extremum(extremum) {
        activate(data, swept);
    }

    void activate(DataManager* data, DistVec<std::uint64_t>* swept){
        if (body == nullptr)
            body = new ArcBody(extremum, data, swept);
    }

    void deactivate(){
        deactivated = true;
        if (body != nullptr)
            delete body;
        body = nullptr;
    }

public:

    bool deactivated = false;
    uint64_t extremum = INVALID_VERTEX;
    uint64_t saddle = INVALID_VERTEX;
    ArcBody* body = nullptr;
};
