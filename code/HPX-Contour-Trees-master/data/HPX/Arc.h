#pragma once

#include "Augmentation.h"
#include "Boundary.h"

#include <limits>
#include <vector>

#include <hpx/hpx.hpp>

class Arc {
public:
    Arc()
    {
    }

    Arc(uint64_t extremum)
        : extremum(extremum)
    {
    }

public:
    bool leaf = false;
    bool started = false;

    uint64_t extremum = INVALID_VERTEX;
    uint64_t saddle = INVALID_VERTEX;

    Boundary boundary;
    Augmentation augmentation;

    std::vector<Arc*> children;
    std::vector<Augmentation> inheritedAugmentations;

    hpx::lcos::local::mutex lock;
};

//struct ArcDataValueComparator {

//    static Datamanager* data;

//    //    static void init(Datamanager* d)
//    //    {
//    //        data = d;
//    //    }

//    static bool compare(Arc*& n1, Arc*& n2)
//    {
//        return (data->getValue(n1->extremum) < data->getValue(n2->extremum));
//    }
//};
