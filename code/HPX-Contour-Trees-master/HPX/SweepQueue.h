#pragma once

class SweepQueue;

#include "DataManager.h"

#include "DistVec.h"

#include <boost/heap/fibonacci_heap.hpp>
#include <deque>

#include <hpx/hpx.hpp>

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
    SweepQueue(uint64_t v, DataManager* data, DistVec<uint64_t>* swept): swept(swept), id(v), queue(6){

    }

    void push(const std::vector<Value>& boundaryVertices){
        std::lock_guard<hpx::lcos::local::mutex> tmplock(lock);
        for (const Value& c : boundaryVertices) {
            this->queue.push_front(c.vertex);
        }
    }

    void push(uint64_t const& target){
        std::lock_guard<hpx::lcos::local::mutex> tmplock(lock);
        this->queue.push_front(target);
    }

    void push(const uint64_t* neighbors, uint32_t n) {
         std::lock_guard<hpx::lcos::local::mutex> tmplock(lock);
        for (uint32_t i = 0; i < n; ++i) {
            if (neighbors[i] != INVALID_VERTEX) {
                if (swept->operator [](neighbors[i]) == INVALID_VERTEX) {
                    //this->queue.push(neighbors[i]);
                    this->queue.push_front(neighbors[i]);
                }
            }
        }
    }

    bool empty() const
    {
        std::lock_guard<hpx::lcos::local::mutex> tmplock(lock);
        return this->queue.empty();
    }

    uint64_t pop() {
        std::lock_guard<hpx::lcos::local::mutex> tmplock(lock);
        while (!this->queue.empty()) {
            uint64_t result = this->queue.front();
            this->queue.pop_front();
            if (swept->operator [](result) == INVALID_VERTEX)
                return result;
        }
        return INVALID_VERTEX;
    }

private:
    mutable hpx::lcos::local::mutex lock;

    DistVec<uint64_t>* swept;
    //boost::heap::fibonacci_heap<uint64_t, boost::heap::compare<DataComparator>, boost::heap::constant_time_size<false>> queue;
    //std::stack<uint64_t> queue;
    std::deque<uint64_t> queue;

    uint64_t id;
};

