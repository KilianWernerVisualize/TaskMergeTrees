#pragma once

#include "DataManager.h"
#include "profiler.h"

#include <vector>
#include <algorithm>

struct SkipNode {
    Value key;
    std::vector<SkipNode*> forward;

    SkipNode(Value k, int level)
        : key(k)
    {
        forward.resize(level, nullptr);
    }
};

class SkipListSet {
public:
    SkipListSet()
    {
    }

    SkipListSet(float prob, int maxLevel)
    {
        this->probability = prob;
        this->maxLevel = maxLevel;
        head = new SkipNode(std::numeric_limits<Value>::min(), maxLevel);
        back = new SkipNode(std::numeric_limits<Value>::max(), maxLevel);
        for (int i = 0; (i < head->forward.size()); i++) {
            head->forward[i] = back;
            back->forward[i] = head;
        }
    }

    void clear(){
        for (int i = 0; (i < head->forward.size()); i++) {
            head->forward[i] = back;
            back->forward[i] = head;
        }
    }

    void reclaim()
    {
        delete head;
        delete back;
    }

    class Iterator {
        friend class SkipListSet;

    private:
        Iterator(SkipNode* node)
            : node(node)
        {
        }

    public:
        Iterator()
            : node(0)
        {
        }

        Iterator& operator++()
        {
            node = node->forward[0];
            return *this;
        }

        SkipNode* operator*()
        {
            return node;
        }

        bool operator!=(const Iterator& other)
        {
            return !(this->operator==(other));
        }

        bool operator==(const Iterator& other)
        {
            return node == other.node;
        }

    protected:
        SkipNode* node;
    };

    Iterator begin() { return Iterator(head->forward[0]); }
    Iterator end() { return Iterator(back); }

    SkipListSet split(Value splitKey)
    {
        SkipNode* x = head;
        SkipListSet result = SkipListSet(0.5f, 24);

        for (int i = maxLevel - 1; i >= 0; i--) {
            while (x->forward[i]->key < splitKey) {
                x = x->forward[i];
            }
            if (x->forward[i] != back) {
                result.injectRange(x->forward[i], back->forward[i], i);
                x->forward[i] = back;
                back->forward[i] = x;
            }
        }

        return result;
    }

    void injectLast(SkipNode* node, int level)
    {
        for (int i = level; i >= 0; i--) {
            back->forward[i]->forward[i] = node;
            node->forward[i] = back;
            back->forward[i] = node;
        }
    }

    void injectRange(SkipNode* front, SkipNode* last, int level)
    {
        head->forward[level] = front;
        last->forward[level] = back;
        back->forward[level] = last;
    }

    void insert(Value searchKey)
    {
        SkipNode* x = nullptr;
        x = find(searchKey);
        if (x) {
            return;
        }

        // vector of pointers that needs to be updated to account for the new node
        std::vector<SkipNode*> update(head->forward);
        int currentMaximum = nodeLevel(head->forward);
        x = head;

        // search the list
        for (int i = currentMaximum - 1; i >= 0; i--) {
            while (x->forward[i]->key < searchKey) {
                x = x->forward[i];
            }
            update[i] = x;
        }
        x = x->forward[0];

        // create new node
        int newNodeLevel = 1;
        if (x->key != searchKey) {

            newNodeLevel = randomLevel();
            int currentLevel = nodeLevel(update);

            if (newNodeLevel > currentLevel) {
                for (int i = newNodeLevel-1; i >= currentLevel-1; i--) {
                    update[i] = head;
                }
            }

            x = makeNode(searchKey, newNodeLevel);
        }

        // connect pointers of predecessors and new node to successors
        for (int i = 0; i < newNodeLevel; i++) {

            x->forward[i] = update[i]->forward[i];
            update[i]->forward[i] = x;
            if (x->forward[i] == back)
                back->forward[i] = x;
        }
    }

private:
    SkipNode* head;
    SkipNode* back;

    int randomLevel()
    {
        int v = 1;

        while ((((double) std::rand() / RAND_MAX)) < probability && v < maxLevel) {

            v += 1;
        }
        return v;
    }

    int nodeLevel(const std::vector<SkipNode*>& v)
    {
        int currentLevel = 0;
        // last element's key is the largest
        Value nilKey = std::numeric_limits<Value>::max();

        for (size_t i = 0; i < v.size(); i++) {
            currentLevel++;
            if (v[i]->key == nilKey) {
                break;
            }
        }
        return currentLevel;
    }

    SkipNode* makeNode(Value key, int level)
    {
        return new SkipNode(key, level);
    }

    SkipNode* find(Value searchKey)
    {
        SkipNode* x = head;
        int currentMaximum = nodeLevel(head->forward);

        for (int i = currentMaximum - 1; i >= 0; i--) {
            while (x->forward[i]->key < searchKey) {
                x = x->forward[i];
            }
        }
        x = x->forward[0];

        if (x->key == searchKey) {
            return x;
        } else {
            return nullptr;
        }
    }

    float probability;
    int maxLevel;
};

#ifndef FLATAUGMENTATION
class Augmentation {
public:
    Augmentation()
        : vertices(0.5f, 24)
    {
    }

    void sweep(Value v)
    {
        vertices.insert(v);
    }

    void clear(){
        vertices.clear();
    }

    void inherit(std::vector<Augmentation>& heritage)
    {
        Profiler profiler(__func__);
        if (heritage.size() == 0)
            return;
        if (heritage.size() == 1) {
            vertices = heritage.at(0).vertices;
            return;
        }
        std::vector<SkipListSet::Iterator> iters;
        for (Augmentation& current : heritage) {
            iters.push_back(current.vertices.begin());
        }
        vertices = SkipListSet(0.5f, 24);
        while (true) {
            int smallest = 0;
            for (int i = 1; i < iters.size(); i++) {
                if ((*iters[i])->key < (*iters[smallest])->key)
                    smallest = i;
            }
            if ((*iters[smallest])->key == std::numeric_limits<Value>::max())
                break;
            SkipNode* next = *iters[smallest];
            ++iters[smallest];
            vertices.injectLast(next, next->forward.size() - 1);
        }
        for (Augmentation& current : heritage) {
            current.vertices.reclaim();
        }
        heritage.clear();
    }

    // inherit from child to parent
    Augmentation heritage(Value saddle)
    {
        Profiler profiler(__func__);
        Augmentation result;
        result.vertices = vertices.split(saddle);
        return result;
    }

    SkipListSet vertices = SkipListSet(0.5f, 24);
};
#else
class Augmentation {
public:
    Augmentation(){
    }

    void sweep(Value v)
    {
        vertices.push_back(v);
    }

    void clear(){
        vertices.clear();
    }

    void inherit(std::vector<Augmentation>& heritage)
    {
        Profiler profiler(__func__);
        if (heritage.size() == 0)
            return;
        if (heritage.size() == 1) {
            vertices = heritage.at(0).vertices;
            return;
        }

        for (Augmentation& current : heritage){
            vertices.insert(vertices.end(), make_move_iterator(current.vertices.begin()), make_move_iterator(current.vertices.end()));
            current.vertices.clear();
        }

        heritage.clear();
    }

    // inherit from child to parent
    Augmentation heritage(Value saddle)
    {
        Profiler profiler(__func__);
        Augmentation result;

        auto splitIt = std::partition_point(vertices.begin(), vertices.end(), [saddle](const auto& em){ return em < saddle; });
        result.vertices.insert(result.vertices.end(), make_move_iterator(splitIt), make_move_iterator(vertices.end()));

        return result;
    }

    std::vector<Value> vertices;
};

#endif
