#ifndef ARC_H
#define ARC_H

#include "Augmentation.h"
#include "Boundary.h"
#include <list>
#include <stdlib.h>
#include <unordered_set>
#include <vector>

class Arc {
public:
    Arc()
    {
    }

    Arc(int extremum)
    {
        this->extremum = extremum;
    }

public:
    //    bool last = false;
    bool leaf = false;
    bool started = false;
    int extremum = -1;
    int saddle = -1;
    Boundary boundary;
    Augmentation augmentation;

    //struct ctBranch* parent = nullptr; Parent is always branchMap.at(saddle)
    std::vector<Arc*> children = std::vector<Arc*>(); //Contains only direct children
    std::vector<Augmentation> inherited_augmentations = std::vector<Augmentation>();
    //std::list<ctBranch*> descendants = std::list<ctBranch*>(); //Contains all the ancestors for fast iteration
    //std::unordered_set<int>* fast_descendants = new std::unordered_set<int>(); //Contains all the ancestors for fast lookup

    hpx::lcos::local::mutex lock;
};

struct ArcDataValueComparator {

    static Datamanager* data;

    //    static void init(Datamanager* d)
    //    {
    //        data = d;
    //    }

    static bool compare(Arc*& n1, Arc*& n2)
    {
        return (data->getValue(n1->extremum) < data->getValue(n2->extremum));
    }
};

#endif
