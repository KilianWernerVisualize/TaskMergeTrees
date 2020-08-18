#pragma once

#include "DataManager.h"
#include "TreeConstructor.h"

#include <boost/heap/fibonacci_heap.hpp>

class DataComparator {
public:
    DataComparator(DataManager* d)
        : data(d)
    {
    }

    bool operator()(uint64_t n1, uint64_t n2) const
    {
        return this->data->getValue(n1) > this->data->getValue(n2);
    }

private:
    DataManager* data;
};

class SweepQueue {
public:
    SweepQueue(uint64_t v, TreeConstructor* tree)
        : tree(tree)
        , queue(DataComparator(tree->dataManager))
        , id(v)
    {
    }

    void push(const std::vector<Value>& boundaryVertices)
    {
        for (const Value& c : boundaryVertices) {
            this->queue.push(c.vertex);
        }
    }

    void push(const uint64_t* neighbors, uint32_t n)
    {
        for (uint32_t i = 0; i < n; ++i) {
            if (neighbors[i] != INVALID_VERTEX) {
                if (this->tree->swept[neighbors[i]] == INVALID_VERTEX) {
                    this->queue.push(neighbors[i]);
                }
            }
        }
    }

    bool empty() const
    {
        return this->queue.empty();
    }

    uint64_t pop()
    {
        while (!this->queue.empty()) {
            uint64_t result = this->queue.top();
            this->queue.pop();
            if (this->tree->swept[result] == INVALID_VERTEX)
                return result;
        }
        return INVALID_VERTEX;
    }

private:
    TreeConstructor* tree;
    boost::heap::fibonacci_heap<uint64_t, boost::heap::compare<DataComparator>> queue;

    uint64_t id;
};
