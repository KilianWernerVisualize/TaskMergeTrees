#pragma once

#include "DataManager.h"

#include <hpx/hpx.hpp>

#include <set>
#include <vector>

class Boundary {
public:
    void add(Value t)
    {
        vertices.emplace(t);
    }

    void remove(Value t)
    {
        vertices.erase(t);
    }

    void unite(Boundary& other)
    {
        if (other.vertices.size() > vertices.size())
            vertices.swap(other.vertices);
        for (const Value& c : other.vertices)
            vertices.emplace(c);
        other.vertices.clear();
    }

    std::vector<Value> intersect(Boundary& other)
    {
        std::vector<Value> result;
        std::set<Value>* big = (vertices.size() > other.vertices.size()) ? &vertices : &other.vertices;
        std::set<Value>* small = (vertices.size() > other.vertices.size()) ? &other.vertices : &vertices;

        for (Value c : *small) {
            if (big->count(c)) {
                result.push_back(c);
            }
        }

        for (Value c : result) {
            big->erase(c);
            small->erase(c);
        }

        return result;
    }

    uint64_t min() const
    {
        return (*vertices.begin()).vertex;
    }

    bool empty() const
    {
        return vertices.empty();
    }

private:
    std::set<Value> vertices;
    hpx::lcos::local::mutex lock;
};
