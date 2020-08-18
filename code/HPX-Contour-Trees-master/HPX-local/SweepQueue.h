#ifndef SWEEPQUEUE_H
#define SWEEPQUEUE_H

#include "Boundary.h"
#include "TreeConstructor.h"
#include <boost/container/flat_set.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <queue>
#include <stdlib.h>
#include <vector>

struct DataValueComparator {

    static Datamanager* data;

    //    static void init(Datamanager* d)
    //    {
    //        data = d;
    //    }

    bool operator()(const int& n1, const int& n2) const
    {
        return data->getValue(n1) > data->getValue(n2);
    }
};

class SweepQueue {
public:
    SweepQueue(int v, ::server::TreeConstructor* t)
    {
        this->queue = boost::heap::fibonacci_heap<int, boost::heap::compare<DataValueComparator>>();
        // this->queue = std::queue<int>();
        this->id = v;
        this->tree = t;
    }

    SweepQueue(int v, const std::vector<int>& neighbors)
    {
        this->queue = boost::heap::fibonacci_heap<int, boost::heap::compare<DataValueComparator>>();
        // this->queue = std::queue<int>();
        this->id = v;
        this->push(neighbors);
    }

    void push(const std::vector<ctValue>& boundaryVertices)
    {
        for (ctValue c : boundaryVertices) {
            this->queue.push(c.key);
        }
    }

    void push(const std::vector<int>& neighbors)
    {
        for (int c : neighbors) {
            if (this->tree->swept[c] == -1)
                this->queue.push(c);
        }
    }

    void push(const int* neighbors, int n)
    {
        for (int i = 0; i < n; ++i) {
            if (this->tree->swept[neighbors[i]] == -1)
                this->queue.push(neighbors[i]);
        }
    }

    void push(const boost::container::flat_set<int>& neighbors)
    {
        for (auto itr = neighbors.begin(); itr != neighbors.end(); ++itr) {
            if (this->tree->swept[*itr] == -1)
                this->queue.push(*itr);
        }
    }

    bool empty()
    {
        return this->queue.empty();
    }

    int pop()
    {
        while (!this->queue.empty()) {
            int result = this->queue.top();
            //int result = this->queue.front();
            this->queue.pop();
            if (this->tree->swept[result] == -1)
                return result;
        }
        return -1;
    }

private:
    ::server::TreeConstructor* tree;

    boost::heap::fibonacci_heap<int, boost::heap::compare<DataValueComparator>> queue;
    //std::queue<int> queue;
    int id;
};

#endif
