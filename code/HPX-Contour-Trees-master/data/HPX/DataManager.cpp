#include "DataManager.h"
#include "Log.h"

#include <boost/algorithm/string/predicate.hpp>
#include <glm/glm.hpp>

#include <teem/nrrd.h>

#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>

const uint64_t BLOCK_INDEX_SHIFT = 54;
const uint64_t BLOCK_INDEX_MASK = 0xFFC0000000000000ull;
const uint64_t VERTEX_INDEX_MASK = ~BLOCK_INDEX_MASK;

/**
 * @brief Base class for regular grids.
 */
class RegularGridManager : public DataManager {
public:
    virtual ~RegularGridManager()
    {
        if (this->blockData)
            delete[] this->blockData;
    }

    /**
     * @brief Returns the total number of vertices across all localities.
     * @return
     */
    uint64_t getNumVertices() const
    {
        return this->gridSize.x * this->gridSize.y * this->gridSize.z;
    }

    /**
     * @brief Returns the number of vertices on this local instance.
     * @param withGhost
     * @return
     */
    uint64_t getNumVerticesLocal(bool withGhost) const
    {
        if (withGhost)
            return this->blockSizeWithGhost.x * this->blockSizeWithGhost.y * this->blockSizeWithGhost.z;

        return this->blockSize.x * this->blockSize.y * this->blockSize.z;
    }

    /**
     * @brief Returns the (global) indices of the vertices on this locality.
     * @return
     */
    std::vector<uint64_t> getLocalVertices() const
    {
        std::vector<uint64_t> localVertices;
        localVertices.reserve(this->blockSize.x * this->blockSize.y * this->blockSize.z);

        for (uint64_t v = 0; v < this->getNumVerticesLocal(true); ++v) {
            if (!(this->blockMask[v] & 0x80)) // only non-ghost vertices
                localVertices.push_back(v | this->blockIndex);

            //            // Test: only vertices with full neighborhood
            //            const uint8_t mask = this->blockMask[v & VERTEX_INDEX_MASK];
            //            if ((mask & 0x20) && (mask & 0x10) && (mask & 0x8) && (mask & 0x4) && (mask & 0x2) && (mask & 0x1))
            //                localVertices.push_back(v | this->blockIndex);
        }

        return localVertices;
    }

    /**
     * @brief Returns the value of the vertex with given id.
     * @param v
     * @return
     */
    Value getValue(uint64_t v) const
    {
        assert((v & BLOCK_INDEX_MASK) == this->blockIndex);

        const uint64_t i = v & VERTEX_INDEX_MASK;
        return Value(this->blockData[i], i);
    }

    /**
     * @brief Returns the neighbors of the vertex with given id.
     * @param v We know this is only called for local non-ghost vertices.
     * @param neighborsOut Can contain invalid indices
     * @return (Maximum) Number of neighbors
     */
    uint32_t getNeighbors(uint64_t v, uint64_t* neighborsOut) const
    {
        assert((v & BLOCK_INDEX_MASK) == this->blockIndex);

        const uint8_t mask = this->blockMask[v & VERTEX_INDEX_MASK];

        neighborsOut[0] = (mask & 0x20) ? (v - 1) : INVALID_VERTEX;
        neighborsOut[1] = (mask & 0x10) ? (v + 1) : INVALID_VERTEX;
        neighborsOut[2] = (mask & 0x8) ? (v - this->blockSizeWithGhost.x) : INVALID_VERTEX;
        neighborsOut[3] = (mask & 0x4) ? (v + this->blockSizeWithGhost.x) : INVALID_VERTEX;
        neighborsOut[4] = (mask & 0x2) ? (v - this->blockSizeWithGhost.x * this->blockSizeWithGhost.y) : INVALID_VERTEX;
        neighborsOut[5] = (mask & 0x1) ? (v + this->blockSizeWithGhost.x * this->blockSizeWithGhost.y) : INVALID_VERTEX;

        return 6;
    }

    /**
     * @brief Checks if the given vertex is a local minimum.
     * @param v
     * @return
     */
    bool isMinimum(uint64_t v) const
    {
        uint64_t neighbors[6];
        this->getNeighbors(v, neighbors);

        const Value value = this->getValue(v);

        for (uint32_t i = 0; i < 6; ++i) {
            const uint64_t neighbor = neighbors[i];

            if (neighbor != INVALID_VERTEX && this->getValue(neighbor) < value)
                return false;
        }

        return true;
    }

    /**
     * @brief Returns the block index of the vertex with given global id.
     * @param v
     * @return
     */
    uint32_t getBlockIndex(uint64_t v) const
    {
        return static_cast<uint32_t>(v >> BLOCK_INDEX_SHIFT);
    }

    /**
     * @brief Converts a remotely produced ghost vertex id to the actual local vertex id.
     * @param v
     * @return
     */
    uint64_t convertToLocal(uint64_t v) const
    {
        // TODO
        return 0;
    }

protected:
    RegularGridManager() = default;

    /**
     * @brief Loads the requested sub block of the data.
     * @param blockIndex
     * @param numBlocks
     */
    virtual void init(uint32_t blockIndex, uint32_t numBlocks)
    {
        this->gridSize = this->getSize();

        if (numBlocks > 1) {
            // Compute block size based on prime factorization of block count
            std::vector<uint32_t> primeFactors;
            uint32_t z = 2;
            uint32_t n = numBlocks;
            while (z * z <= n) {
                if (n % z == 0) {
                    primeFactors.push_back(z);
                    n = (n / z);
                } else {
                    ++z;
                }
            }
            if (n > 1)
                primeFactors.push_back(n);

            this->baseBlockSize = this->gridSize;
            for (uint32_t f : primeFactors) {
                uint32_t maxDim = 0;
                for (uint32_t d = 1; d < 3; ++d)
                    if (this->baseBlockSize[d] > this->baseBlockSize[maxDim])
                        maxDim = d;
                this->baseBlockSize[maxDim] /= f;
            }

            this->numBlocks = this->gridSize / this->baseBlockSize;

            // Compute index of local block
            n = blockIndex;
            this->blockIndex3D.z = n / (this->numBlocks.x * this->numBlocks.y);
            n = n % (this->numBlocks.x * this->numBlocks.y);
            this->blockIndex3D.y = n / this->numBlocks.x;
            this->blockIndex3D.x = n % this->numBlocks.x;

            // Determine offset and size of this block
            this->blockOffset = this->blockIndex3D * this->baseBlockSize;
            this->blockSize = this->baseBlockSize;

            // Fix size of last block along each dimension
            for (uint32_t d = 0; d < 3; ++d)
                if (this->blockIndex3D[d] == this->numBlocks[d] - 1)
                    this->blockSize[d] = this->gridSize[d] - (this->numBlocks[d] - 1) * this->baseBlockSize[d];
        } else {
            this->blockOffset = glm::uvec3(0, 0, 0);
            this->blockSize = this->gridSize;
            this->numBlocks = glm::uvec3(1, 1, 1);
            this->baseBlockSize = this->gridSize;
            this->blockIndex3D = glm::uvec3(0, 0, 0);
        }

        // Compute ghost layer
        this->blockOffsetWithGhost = this->blockOffset;
        this->blockSizeWithGhost = this->blockSize;

        glm::uvec3 beginNonGhost = glm::uvec3(0, 0, 0);
        glm::uvec3 endNonGhost = this->blockSize;

        for (uint32_t d = 0; d < 3; ++d) {
            if (this->blockIndex3D[d] > 0) {
                --this->blockOffsetWithGhost[d];
                ++this->blockSizeWithGhost[d];

                ++beginNonGhost[d];
                ++endNonGhost[d];
            }

            if (this->blockIndex3D[d] < this->numBlocks[d] - 1) {
                ++this->blockSizeWithGhost[d];
            }
        }

        // Compute block mask
        uint64_t numVerticesWithGhost = this->getNumVerticesLocal(true);
        this->blockMask.resize(numVerticesWithGhost, 0);

        for (uint32_t z = 0; z < this->blockSizeWithGhost.z; ++z) {
            for (uint32_t y = 0; y < this->blockSizeWithGhost.y; ++y) {
                for (uint32_t x = 0; x < this->blockSizeWithGhost.x; ++x) {
                    uint8_t mask = 0;

                    // Check if ghost
                    if (x < beginNonGhost.x || x >= endNonGhost.x || y < beginNonGhost.y || y >= endNonGhost.y || z < beginNonGhost.z || z >= endNonGhost.z)
                        mask |= 0x80;

                    // Has -x neighbor
                    if (x > 0)
                        mask |= 0x20;

                    // Has +x neighbor
                    if (x < this->blockSizeWithGhost.x - 1)
                        mask |= 0x10;

                    // Has -y neighbor
                    if (y > 0)
                        mask |= 0x8;

                    // Has +y neighbor
                    if (y < this->blockSizeWithGhost.y - 1)
                        mask |= 0x4;

                    // Has -z neighbor
                    if (z > 0)
                        mask |= 0x2;

                    // Has +z neighbor
                    if (z < this->blockSizeWithGhost.z - 1)
                        mask |= 0x1;

                    this->blockMask[z * this->blockSizeWithGhost.x * this->blockSizeWithGhost.y + y * this->blockSizeWithGhost.x + x] = mask;
                }
            }
        }

        // Read data and release internal reader memory
        this->blockData
            = new float[numVerticesWithGhost];
        if (this->blockData == nullptr) {
            LogError().tag(std::to_string(blockIndex)) << "Failed to allocate block memory";
            return;
        }

        this->readBlock(this->blockOffsetWithGhost, this->blockSizeWithGhost, this->blockData);
        this->release();

        // Store block index in msb
        this->blockIndex = static_cast<uint64_t>(blockIndex) << BLOCK_INDEX_SHIFT;

        // Print info
        if (blockIndex == 0) {
            Log() << "Grid size: " << this->gridSize;
            Log() << "Num blocks: " << this->numBlocks;
            Log() << "Vertices: " << this->getNumVertices();
        }
        Log().tag(std::to_string(blockIndex)) << "Block index: " << this->blockIndex3D;
        Log().tag(std::to_string(blockIndex)) << "Block offset: " << this->blockOffset;
        Log().tag(std::to_string(blockIndex)) << "Block size: " << this->blockSize;
        Log().tag(std::to_string(blockIndex)) << "Vertices (local): " << this->getNumVerticesLocal(false);
        Log().tag(std::to_string(blockIndex)) << "Vertices (ghost): " << this->getNumVerticesLocal(true) - this->getNumVerticesLocal(false);

        Log().tag(std::to_string(blockIndex)) << "Values: " << byteString(numVerticesWithGhost * sizeof(float));
        Log().tag(std::to_string(blockIndex)) << "Mask: " << byteString(numVerticesWithGhost * sizeof(uint8_t));
    }

    virtual glm::uvec3 getSize() = 0;
    virtual void readBlock(const glm::uvec3& offset, const glm::uvec3& size, float* dataOut) = 0;
    virtual void release() = 0;

private:
    RegularGridManager(const RegularGridManager&) = delete;

private:
    uint64_t blockIndex;

    glm::uvec3 gridSize;
    glm::uvec3 numBlocks;
    glm::uvec3 blockIndex3D;
    glm::uvec3 baseBlockSize;

    // Exact block of this locality
    glm::uvec3 blockOffset;
    glm::uvec3 blockSize;

    // Block including ghost layer
    glm::uvec3 blockOffsetWithGhost;
    glm::uvec3 blockSizeWithGhost;

    float* blockData = nullptr;
    std::vector<uint8_t> blockMask; // msb to lsb: [ghost cell, unused, -x, +x, -y, +y, -z, +z]
};

/**
 * @brief Implementation for nrrd grids.
 */
class NrrdManager : public RegularGridManager {
public:
    NrrdManager(const std::string& path)
    {
        this->nrrd = nrrdNew();
        if (nrrdLoad(this->nrrd, path.c_str(), 0)) {
            LogError() << "Failed to read nrrd file: " << path;
            nrrdNuke(this->nrrd);
            return;
        }

        if (this->nrrd->type != nrrdTypeFloat) {
            Nrrd* tmp = nrrdNew();
            nrrdConvert(tmp, this->nrrd, nrrdTypeFloat);
            nrrdNuke(this->nrrd);
            this->nrrd = tmp;
        }
    }

    virtual ~NrrdManager() = default;

    glm::uvec3 getSize()
    {
        size_t size[4];
        nrrdAxisInfoGet_va(this->nrrd, nrrdAxisInfoSize, &size[0], &size[1], &size[2], &size[3]);

        return glm::uvec3(size[0], size[1], size[2]);
    }

    void readBlock(const glm::uvec3& offset, const glm::uvec3& size, float* blockOut)
    {
        const float* data = static_cast<float*>(this->nrrd->data);
        const glm::uvec3 dataSize = this->getSize();

        for (uint32_t z = 0; z < size.z; ++z)
            for (uint32_t y = 0; y < size.y; ++y)
                memcpy(blockOut + z * size.y * size.x + y * size.x, data + (offset.z + z) * dataSize.y * dataSize.x + (offset.y + y) * dataSize.x + offset.x, size.x * sizeof(float));
    }

    void release()
    {
        if (this->nrrd)
            nrrdNuke(this->nrrd);
    }

private:
    Nrrd* nrrd = nullptr;
};

/**
 * @brief Implementation for VTK's vti grids.
 */
class VTIManager : public RegularGridManager {
public:
    VTIManager(const std::string& path)
    {
        this->reader = vtkXMLImageDataReader::New();
        this->reader->SetFileName(path.c_str());
        this->reader->Update();
    }

    virtual ~VTIManager() = default;

    glm::uvec3 getSize()
    {
        int dims[3];
        this->reader->GetOutput()->GetDimensions(dims);
        return glm::uvec3(dims[0], dims[1], dims[2]);
    }

    void readBlock(const glm::uvec3& offset, const glm::uvec3& size, float* blockOut)
    {
        const glm::uvec3 dataSize = this->getSize();

        if (this->reader->GetOutput()->GetScalarType() == VTK_FLOAT) {
            const float* data = static_cast<float*>(this->reader->GetOutput()->GetScalarPointer());

            for (uint32_t z = 0; z < size.z; ++z)
                for (uint32_t y = 0; y < size.y; ++y)
                    memcpy(blockOut + z * size.y * size.x + y * size.x, data + (offset.z + z) * dataSize.y * dataSize.x + (offset.y + y) * dataSize.x + offset.x, size.x * sizeof(float));
        } else if (this->reader->GetOutput()->GetScalarType() == VTK_DOUBLE) {
            const double* data = static_cast<double*>(this->reader->GetOutput()->GetScalarPointer());

            for (uint32_t z = 0; z < size.z; ++z)
                for (uint32_t y = 0; y < size.y; ++y)
                    for (uint32_t x = 0; x < size.x; ++x)
                        blockOut[z * size.y * size.x + y * size.x + x] = static_cast<float>(data[(offset.z + z) * dataSize.y * dataSize.x + (offset.y + y) * dataSize.x + offset.x + x]);
        }
    }

    void release()
    {
        if (this->reader)
            this->reader->Delete();
    }

private:
    vtkXMLImageDataReader* reader = nullptr;
};

/**
 * @brief Constructs the correct instantiation for the given file format.
 * @param path
 * @param blockIndex
 * @param numBlocks
 * @return
 */
DataManager* DataManager::load(const std::string& path, uint32_t blockIndex, uint32_t numBlocks)
{
    DataManager* data = nullptr;

    if (boost::algorithm::ends_with(path, ".nrrd"))
        data = new NrrdManager(path);
    else if (boost::algorithm::ends_with(path, ".vti"))
        data = new VTIManager(path);

    if (data) {
        data->init(blockIndex, numBlocks);
        return data;
    } else {
        LogError() << "Unknown file format: " << path;
    }
}
