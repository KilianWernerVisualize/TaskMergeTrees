#ifndef TREECONSTRUCTOR_H
#define TREECONSTRUCTOR_H

#include <vector>
#include <cstdint>
#include <iostream>

struct ctValue {
    double value; // TODO float faster?
    std::uint64_t key;

    inline bool operator<(const ctValue& other) const
    {
        if (value == other.value)
            return key < other.key;
        return value < other.value;
    }

    inline bool operator>(const ctValue& other) const
    {
        if (value == other.value)
            return key > other.key;
        return value > other.value;
    }

    inline bool operator<=(const ctValue& other) const
    {
        return !(this->operator>(other));
    }

    inline bool operator>=(const ctValue& other) const
    {
        return !(this->operator<(other));
    }

    inline bool operator==(const ctValue& other) const
    {
        return ((value == other.value) && (key == other.key));
    }

    inline bool operator!=(const ctValue& other) const
    {
        return !(this->operator==(other));
    }
};

struct NeighborsEntry {
    int offset;
    int numNeighbors;
};

class TreeConstructor {

public:
   int numVertices;
   double* values;
   float* fvalues;
   std::vector<NeighborsEntry>* neighborsMap;
   std::vector<int>* neighborsBuffer;
   double min_val;
   double conversionFactor;
   int max_x;
   int max_y;
   int max_z;

   const int* getNeighbors(int v, int* numNeighborsOut) const
   {
       *numNeighborsOut = this->neighborsMap->at(v).numNeighbors;
       return &this->neighborsBuffer->at(this->neighborsMap->at(v).offset);
   }
/*
   ctValue getValue(int v) const
   {
       ctValue result;
       result.value = values[v];
       result.key = v;
       return result;
   }*/
}; 

void doCudaStuff(TreeConstructor& tree, std::vector<int*>& minima, std::vector<int*>& saddles, std::vector<int> &counts, int*& swepto, int*& augmentationo, int rounds, int minimaroundsh, int growthroundsh);

void augmentTrunk(TreeConstructor &tree, int trunkStarter, std::vector<int> &danglingvect, int *&swept);
#endif
