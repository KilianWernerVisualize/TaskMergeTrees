#pragma once

#include <vector>
#include <map>
#include "DataManager.h"
#include <hpx/hpx.hpp>

//Stores data associated with vertices, local vertices in array, remote vertices in map
template <typename T>
class DistVec {

public:
    DistVec(){
        data = nullptr;
    }

    DistVec(std::uint64_t size, const T& emptyValue = T(), DataManager* data = nullptr){
        local.resize(size, emptyValue);
        empty = emptyValue;
        this->data = data;
    }

    void init(std::uint64_t size, const T& emptyValue = T(), DataManager* data = nullptr){
        local.resize(size, emptyValue);
        empty = emptyValue;
        this->data = data;
    }

    T& operator [](uint64_t idx){
        if (data->isLocal(idx)){
            return local[idx & VERTEX_INDEX_MASK];
        } else {
            std::lock_guard<hpx::lcos::local::mutex> lock(mapLock);
            if (!remote.count(idx))
                remote[idx] = empty;
            return remote[idx];
        }
    }

    bool contains(uint64_t idx){
        if (data->isLocal(idx)){
            return true;
        } else {
            return remote.count(idx);
        }
    }

    auto begin(){
        return local.begin();
    }

    auto end(){
        return local.end();
    }

    DataManager* data;
    std::vector<T> local;
    std::map<std::uint64_t, T>  remote;
    T empty;

    hpx::lcos::local::mutex mapLock;
};
