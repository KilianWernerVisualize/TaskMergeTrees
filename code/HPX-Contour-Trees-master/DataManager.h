#pragma once

#include <limits>
#include <string>
#include <vector>
#include <glm/glm.hpp>

#include <hpx/hpx.hpp>

const uint64_t INVALID_VERTEX = std::numeric_limits<uint64_t>::max();
//10 MSB encode blockID --> 1024 blocks possible
const uint64_t BLOCK_INDEX_SHIFT = 54;
const uint64_t BLOCK_INDEX_MASK = 0xFFC0000000000000ull;
const uint64_t VERTEX_INDEX_MASK = ~BLOCK_INDEX_MASK;
const uint64_t INVALID_BLOCK = 0xFFC0000000000000ull;

class Value {
public:
    float value;
    uint64_t vertex;

    Value(){
    }

    Value(float value, uint64_t vertex)
        : value(value)
        , vertex(vertex)
    {
    }

    friend std::ostream & operator<< (std::ostream &out, Value const &t){
        out << "(" << t.value << "," << t.vertex << ")";
        return out;
    }

private:
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & value & vertex;
    }

public:
    inline bool operator<(const Value& other) const
    {
        if (this->value == other.value)
            return this->vertex < other.vertex;
        return this->value < other.value;
    }

    inline bool operator>(const Value& other) const
    {
        if (this->value == other.value)
            return this->vertex > other.vertex;
        return this->value > other.value;
    }

    inline bool operator<=(const Value& other) const
    {
        return !(this->operator>(other));
    }

    inline bool operator>=(const Value& other) const
    {
        return !(this->operator<(other));
    }

    inline bool operator==(const Value& other) const
    {
        return ((this->value == other.value) && (this->vertex == other.vertex));
    }

    inline bool operator!=(const Value& other) const
    {
        return !(this->operator==(other));
    }
};

namespace std {
template <>
class numeric_limits<Value> {
public:
    static Value min()
    {
        return Value(std::numeric_limits<float>::min(), std::numeric_limits<uint64_t>::min());
    }

    static Value max()
    {
        return Value(std::numeric_limits<float>::max(), std::numeric_limits<uint64_t>::max());
    }
};
}

/**
 * @brief Base class for distributed data loading and acccess.
 */
class DataManager {
public:
    static DataManager* load(const std::string& path, uint32_t blockIndex, uint32_t numBlocks);
    virtual ~DataManager() = default;

    virtual uint64_t getNumVertices() const = 0;
    virtual uint64_t getNumVerticesLocal(bool withGhost = false) const = 0;

    virtual uint64_t getLocalVertex(uint64_t idx) const = 0;
    virtual std::vector<uint64_t> getLocalVertices() const = 0;

    virtual Value getValue(uint64_t v) const = 0;
    virtual uint32_t getNeighbors(uint64_t v, uint64_t* neighborsOut) const = 0;
    virtual uint64_t getNeighbor(uint64_t v, int i) const = 0;

    virtual bool isMinimum(uint64_t v) const = 0;
    virtual bool isGhost(uint64_t v) const = 0;
    virtual bool isLocal(uint64_t v) const = 0;

    virtual uint32_t getBlockIndex(uint64_t v) const = 0;

    virtual glm::uvec3 getDimensions() const = 0;
    virtual glm::uvec3 getLocalDimensions() const = 0;

    virtual glm::uvec3 getLocalCoords(uint64_t v) const = 0;
    virtual glm::uvec3 getGlobalCoords(uint64_t v) const = 0;
    virtual uint64_t convertToGlobal(uint64_t v) const = 0;
    virtual uint64_t convertToLocal(uint64_t v) const = 0;

protected:
    DataManager()
        = default;

    virtual void init(uint32_t blockIndex, uint32_t numBlocks) = 0;

private:
    DataManager(const DataManager&) = delete;
};
