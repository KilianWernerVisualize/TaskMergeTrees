#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "DataManager.h"
#include <boost/container/flat_set.hpp>
#include <hpx/hpx.hpp>
#include <limits>
#include <set>
#include <stdlib.h>
#include <vector>

class Boundary {
public:
    void add(ctValue t)
    {
        vertices.emplace(t);
    }

    void remove(ctValue t)
    {
        vertices.erase(t);
    }

    void unite(Boundary& other)
    {
        if (other.vertices.size() > vertices.size())
            vertices.swap(other.vertices);
        for (ctValue c : other.vertices)
            vertices.emplace(c);
        other.vertices.clear();
    }

    std::vector<ctValue> intersect(Boundary& other)
    {
        std::vector<ctValue> result;
        std::set<ctValue>* big;
        std::set<ctValue>* small;
        big = &vertices;
        small = &other.vertices;
        if (other.vertices.size() > vertices.size()) {
            big = &other.vertices;
            small = &vertices;
        }

        for (ctValue c : *small) {
            if (big->count(c)) {
                result.push_back(c);
            }
        }

        for (ctValue c : result) {
            big->erase(c);
            small->erase(c);
        }

        return result;
    }

    int min(const Datamanager& d)
    {
        return (*vertices.begin()).key;
    }

    bool empty()
    {
        return vertices.empty();
    }

private:
    // boost::container::flat_set<std::uint64_t> vertices;
    std::set<ctValue> vertices;
    hpx::lcos::local::mutex lock;
};

#endif
