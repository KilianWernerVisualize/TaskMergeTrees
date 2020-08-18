#ifndef TREECONSTRUCTOR_CU
#define TREECONSTRUCTOR_CU
#include "TreeConstructor.h"
#include <chrono>
#include <stack>
//#include <cooperative_groups.h>

#include <stdexcept>

__device__ int minimaQueueCounter;
__device__ int doneCounter;
__device__ int minimaCounter;
__device__ double min_value;
__device__ double conversionFactor;

__device__ int minimarounds;
__device__ int growthrounds;

__device__ int max_x;
__device__ int max_y;
__device__ int max_z;

__device__ int boundarySize;
__device__ int* boundaryStart;

__device__ int trunkStart;

// -------------------------------------------------------------------------

/// throw exception if a CUDA error was generated in the wrapped call
inline void cuda_check( cudaError_t code )
{
    if( code != cudaSuccess )
        throw std::runtime_error( std::string( "CUDA error: " ) + cudaGetErrorString(code) );
}

/// throw exception if a CUDA error occurred in the past
inline void cuda_check()
{
    cuda_check( cudaGetLastError() );
}


// -------------------------------------------------------------------------

/// a simple wrapper for cudaArray3D
template<typename T>
struct array3D
{
    dim3 size() const;

    /// resize to (nx, ny, nz); array will be uninitialized
    void resize( dim3 size );

    /// copy host data (nx*ny*nz) elements into array
    void copy( const T* );

    /// return result of an texture lookup (interpolated read)
    __device__ T get(int idx) const;

    cudaArray*          m_array = 0;
    cudaExtent          m_extent;
    cudaTextureObject_t m_texture;
};

// -------------------------------------------------------------------------

template<typename T>
dim3 array3D<T>::size() const
{
    return dim3( m_extent.width, m_extent.height, m_extent.depth );
}
// -------------------------------------------------------------------------

template<typename T>
void array3D<T>::resize( dim3 size )
{
    if( m_array )
    {
        // free a poossibly previously allocated array
        // and the associated texture object
        cudaDestroyTextureObject( m_texture );
        cudaFreeArray( m_array );

        m_array = 0;
    }

    m_extent = make_cudaExtent( size.x, size.y, size.z );

    auto cdesc = cudaCreateChannelDesc<T>();

    cuda_check( cudaMalloc3DArray(
        &m_array,
        &cdesc,
        make_cudaExtent(
            m_extent.width*sizeof(T),
            m_extent.height,
            m_extent.depth
        ), 0 )
    );

    // set up texture
    cudaResourceDesc tr;
    memset( &tr, 0, sizeof(cudaResourceDesc) );
    tr.resType         = cudaResourceTypeArray;
    tr.res.array.array = m_array;

    cudaTextureDesc td;
    memset( &td, 0, sizeof(cudaTextureDesc) );

    td.filterMode       = cudaFilterModeLinear;
    td.addressMode[0]   = cudaAddressModeBorder;
    td.addressMode[1]   = cudaAddressModeBorder;
    td.addressMode[2]   = cudaAddressModeBorder;
    td.readMode         = cudaReadModeElementType;
    td.normalizedCoords = false;

    cuda_check( cudaCreateTextureObject(
        &m_texture,
        &tr,
        &td,
        NULL
    ) );
}

// -------------------------------------------------------------------------

template<typename T>
void array3D<T>::copy( T const *data  )
{
    cudaMemcpy3DParms copyParams = { 0 };
    copyParams.srcPtr = make_cudaPitchedPtr(
        (void*)data,
        m_extent.width*sizeof(T),
        m_extent.width,
        m_extent.height
    );

    copyParams.dstArray = m_array;
    copyParams.extent   = m_extent;
    copyParams.kind     = cudaMemcpyHostToDevice;

    cuda_check( cudaMemcpy3D( &copyParams ) );
}

// -------------------------------------------------------------------------

template<typename T>
__device__ T array3D<T>::get( int idx ) const
{
    return tex3D<T>( m_texture, idx % max_x, ( idx / max_x ) % max_y, ( idx / ( max_x*max_y )));
}

__global__ void min_search_preNeigh(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, int* augmentation, unsigned long long int* saddleMap){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    NeighborsEntry myNeigh;

    bool isMini;
    double neighVal;
    double myVal;
    int j;

    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        swept[idx] = -1;
        augmentation[idx] = -1;
        saddleMap[idx] = 0-1;
        isMini = true;
        myNeigh = neighborsMap[idx];
        myVal = values.get(idx);
        for (j = 0; (j < myNeigh.numNeighbors); j++){
            neighVal = values.get(neighborsBuffer[myNeigh.offset+j]);
            if (neighVal < myVal || ((neighVal == myVal) && (neighborsBuffer[myNeigh.offset+j] < idx))){
                isMini = false;
            }
        }
        if (isMini){
            j = atomicAdd(&minimaQueueCounter, 1);
            minimaQueue[j] = idx;
        }
    }
}

__global__ void min_growth_preNeigh(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap, int** sweepQueueStarters){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    int v;
    int touched;
    bool sweepable;

    double neighVal;
    double myVal;
    int j;
    int i;

    int stackCount = -1;
    int stack[1000];

    NeighborsEntry myNeigh;

    for (; idx < minimaQueueCounter; idx+=gridDim.x*blockDim.x){

        v = minimaQueue[idx];
        myNeigh = neighborsMap[v];
        stackCount = -1;
        for (j = 0; (j < myNeigh.numNeighbors); j++){
            stack[stackCount + j +1] = neighborsBuffer[myNeigh.offset + j];
        }
        stackCount += myNeigh.numNeighbors;
        swept[v] = v;

        i = 0;
        while (stackCount >= 0 && (i < minimarounds)){
            i++;
            sweepable = true;
            touched = stack[stackCount];
            myVal = values.get(touched);
            stackCount--;

            myNeigh = neighborsMap[touched];
            for (j = 0; (j < myNeigh.numNeighbors); j++){
                neighVal = values.get(neighborsBuffer[myNeigh.offset + j]);
                if ((neighVal < myVal || ((neighVal == myVal) && (neighborsBuffer[myNeigh.offset+j] < touched))) && ((swept[neighborsBuffer[myNeigh.offset + j]] >= -1) && (swept[neighborsBuffer[myNeigh.offset + j]] != v))){
                    sweepable = false;
                }
            }
            if (sweepable){
                swept[touched] = v;
                for (j = 0; (j < myNeigh.numNeighbors); j++){
                    if ((swept[neighborsBuffer[myNeigh.offset + j]] >= -1) && (swept[neighborsBuffer[myNeigh.offset + j]] != v)){
                        stackCount++;
                        if (stackCount >= 1000){
                            printf("Oh no!");
                            stackCount = -1;
                        } else {
                            stack[stackCount] = neighborsBuffer[myNeigh.offset + j];
                        }
                    }
                }
            }
        }
        if (stackCount >= 0){
            sweepQueueStarters[idx + (int)(numVertices/2)] = (int*)malloc((stackCount+1)*sizeof(int));
            for (i = 0; (i <= stackCount); i++){
                sweepQueueStarters[idx + (int)(numVertices/2)][i] = stack[i];
            }
            saddleMap[v] = (((unsigned long long int)idx) << 32) + stackCount+1 + 1;
            swept[v] = -2 - idx;
        } else {
            sweepQueueStarters[idx] = nullptr;
            swept[v] = v;
        }
    }
}

__global__ void saddle_growth(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap, int** sweepQueueStarters, int round){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    int v;
    int touched;
    bool sweepable;

    double neighVal;
    double myVal;
    int j;
    int i;

    int stackCount = -1;
    int stack[1000];

    NeighborsEntry myNeigh;

    for (; idx < minimaQueueCounter; idx+=gridDim.x*blockDim.x){
        v = minimaQueue[idx];
        myNeigh = neighborsMap[v];
        stackCount = -1;
        for (j = 0; (j < myNeigh.numNeighbors); j++){
            stack[stackCount + j +1] = neighborsBuffer[myNeigh.offset + j];
        }
        stackCount += myNeigh.numNeighbors;

        swept[v] = v;

        i = 0;
        while ((stackCount >= 0 || ((unsigned int)saddleMap[v]) > 1) && (i < growthrounds)){
            i++;
            sweepable = true;
            if (stackCount >= 0){
                touched = stack[stackCount];
                stackCount--;
            } else {
                saddleMap[v]--;
                touched = sweepQueueStarters[idx][((unsigned int)saddleMap[v])-1];
                //printf("Taking %d from %u at %u \n", touched, v, (((unsigned int)saddleMap[v])-1));
            }
            myVal = values.get(touched);

            myNeigh = neighborsMap[touched];
            for (j = 0; (j < myNeigh.numNeighbors); j++){
                neighVal = values.get(neighborsBuffer[myNeigh.offset + j]);
                if ((neighVal < myVal || ((neighVal == myVal) && (neighborsBuffer[myNeigh.offset+j] < touched))) && ((swept[neighborsBuffer[myNeigh.offset + j]] >= -1) && (swept[neighborsBuffer[myNeigh.offset + j]] != v))){
                    sweepable = false;
                }
            }
            if (sweepable){
                swept[touched] = v;
                for (j = 0; (j < myNeigh.numNeighbors); j++){
                    if ((swept[neighborsBuffer[myNeigh.offset + j]] >= -1) && (swept[neighborsBuffer[myNeigh.offset + j]] != v)){
                        stackCount++;
                        if (stackCount >= 1000){
                            printf("Oh no!");
                            stackCount = -1;
                        } else {
                            stack[stackCount] = neighborsBuffer[myNeigh.offset + j];
                        }
                    }
                }
            }
        }
        if (stackCount >= 0 || ((unsigned int)saddleMap[v]) > 1){
            j = (stackCount+((unsigned int)saddleMap[v]));
            int* tmp = (int*)malloc(j*sizeof(int));
            for (i = 0; (i <= stackCount); i++)
                tmp[i] = stack[i];
            for (i = stackCount+1; ( i < j); i++){
                tmp[i] = sweepQueueStarters[idx][i-stackCount-1];
            }
            free(sweepQueueStarters[idx]);
            sweepQueueStarters[idx + (int)(numVertices/2)] = tmp;
            saddleMap[v] = (((unsigned long long int)idx) << 32) + j + 1;
            swept[v] = -2 - idx;
        } else {
            free(sweepQueueStarters[idx]);
            sweepQueueStarters[idx] = nullptr;
            swept[v] = v;
        }
    }
}

__global__ void clean_saddleMap_to_infinity(int numVertices, int* swept, unsigned long long int* saddleMap){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] >= -1)
            saddleMap[idx] = 0-1;
    }
}

__global__ void clean_saddleMap_to_zero(int numVertices, int* swept, unsigned long long int* saddleMap){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (saddleMap[idx] != 0-2 && swept[idx] >= -1)
            saddleMap[idx] = 0;
    }
}

__global__ void saddle_search(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    NeighborsEntry myNeigh;
    int j;
    double neighVal;
    double myVal;
    unsigned int convertedValue;

    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] == -1){
            myNeigh = neighborsMap[idx];
            myVal = values.get(idx);
            for (j = 0; (j < myNeigh.numNeighbors); j++){
                neighVal = values.get(neighborsBuffer[myNeigh.offset+j]);
                if (neighVal < myVal || ((neighVal == myVal) && (neighborsBuffer[myNeigh.offset+j] < idx))){
                    if ((swept[neighborsBuffer[myNeigh.offset+j]] >= 0) && (swept[swept[neighborsBuffer[myNeigh.offset+j]]] >= -1)){
                        convertedValue = static_cast<unsigned int>((myVal - min_value)*conversionFactor);
                        atomicMin(&saddleMap[swept[neighborsBuffer[myNeigh.offset+j]]], (((unsigned long long int)convertedValue) << 32) + static_cast<unsigned long long int>(idx));
                    }
                }
            }
        }
    }
}

__global__ void compress_saddles(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap, int** sweepQueueStarters){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    NeighborsEntry myNeigh;
    int j;
    double neighVal;
    double myVal;
    bool saddle;
    bool checksaddle;


    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] == -1){
            saddle = false;
            checksaddle = true;
            myNeigh = neighborsMap[idx];
            myVal = values.get(idx);

            for (j = 0; (j < myNeigh.numNeighbors); j++){
                neighVal = values.get(neighborsBuffer[myNeigh.offset+j]);
                if (neighVal < myVal || ((neighVal == myVal) && (neighborsBuffer[myNeigh.offset+j] < idx))){
                    if ((swept[neighborsBuffer[myNeigh.offset+j]] >= 0) && (swept[swept[neighborsBuffer[myNeigh.offset+j]]] >= -1)){
                        if (((unsigned int)saddleMap[swept[neighborsBuffer[myNeigh.offset+j]]]) == idx){
                            //if (swept[neighborsBuffer[myNeigh.offset+j]] != idx)
                            //swept[swept[neighborsBuffer[myNeigh.offset+j]]] = idx;
                            saddle = true;
                        } else {
                            checksaddle = false;
                        }
                    } else {
                        checksaddle = false;
                    }
                }
            }
            if (!checksaddle)
                saddleMap[idx] = 0-2;
            if (saddle && checksaddle){
                j = atomicAdd(&minimaQueueCounter, 1);
                minimaQueue[j] = idx;
            }
        } else if (swept[idx] <= -2) {
            j = atomicAdd(&minimaQueueCounter, 1);
            minimaQueue[j] = idx;
            sweepQueueStarters[j] = sweepQueueStarters[-2-swept[idx] + (int)(numVertices/2)];
            sweepQueueStarters[-2-swept[idx] + (int)(numVertices/2)] = nullptr;
        }
    }
}

__global__ void size_sweepqueues(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap, int** sweepQueueStarters){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    NeighborsEntry myNeigh;
    int j;
    int i;
    double neighVal;
    double myVal;
    int mySaddle;
    int myMin[12];


    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] == -1){
            mySaddle = -1;
            myNeigh = neighborsMap[idx];
            myVal = values.get(idx);

            for (j = 0; (j < myNeigh.numNeighbors); j++){
                myMin[j] = -1;
                neighVal = values.get(neighborsBuffer[myNeigh.offset+j]);
                if (neighVal < myVal || ((neighVal == myVal) && (neighborsBuffer[myNeigh.offset+j] < idx))){
                    myMin[j] = swept[neighborsBuffer[myNeigh.offset+j]];
                }
            }
            for (j = 0; (j < myNeigh.numNeighbors-1); j++){
                for (i = j+1; (i < myNeigh.numNeighbors); i++){
                    if (myMin[i] >= 0 && myMin[j] >= 0 && myMin[i] != myMin[j] && swept[myMin[i]] == swept[myMin[j]])
                        mySaddle = swept[myMin[i]];
                }
            }
            if (mySaddle != -1 && saddleMap[mySaddle] != 0-2)
                atomicAdd(&saddleMap[mySaddle], 1);
        }
    }
}

__global__ void create_sweepqueues(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap, int** sweepQueueStarters){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    for(; idx<minimaQueueCounter; idx+=gridDim.x*blockDim.x){
        //printf("Idx: %d, MinimaQueue[idx]: %d \n", idx, minimaQueue[idx]);
        //if (saddleMap[minimaQueue[idx]] > 0)
            //printf("Idx: %d; sweepqueuestartersize: %llu \n",minimaQueue[idx], saddleMap[minimaQueue[idx]]);
        if (swept[minimaQueue[idx]] >= -1){
            sweepQueueStarters[idx] = (int*)malloc(saddleMap[minimaQueue[idx]]*sizeof(int));
            saddleMap[minimaQueue[idx]] = (((unsigned long long int)idx) << 32) + 1;
        }
    }
}

__global__ void fill_sweepqueues(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap, int** sweepQueueStarters){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    NeighborsEntry myNeigh;
    int j;
    int i;
    double neighVal;
    double myVal;
    int myMin[12];
    int mySaddle;
    unsigned long long int pos;


    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] == -1){
            mySaddle = -1;
            myNeigh = neighborsMap[idx];
            myVal = values.get(idx);

            for (j = 0; (j < myNeigh.numNeighbors); j++){
                myMin[j] = -1;
                neighVal = values.get(neighborsBuffer[myNeigh.offset+j]);
                if (neighVal < myVal || ((neighVal == myVal) && (neighborsBuffer[myNeigh.offset+j] < idx))){
                    myMin[j] = swept[neighborsBuffer[myNeigh.offset+j]];
                }
            }
            for (j = 0; (j < myNeigh.numNeighbors-1); j++){
                for (i = j+1; (i < myNeigh.numNeighbors); i++){
                    if (myMin[i] >= 0 && myMin[j] >= 0 && myMin[i] != myMin[j] && swept[myMin[i]] == swept[myMin[j]])
                        mySaddle = swept[myMin[i]];
                }
            }
            if (mySaddle != -1 && saddleMap[mySaddle] != 0-2){
                pos = atomicAdd(&saddleMap[mySaddle], 1);
                //printf("Adding %d to %u at %u \n", idx, minimaQueue[(unsigned int)(pos >> 32)], ((unsigned int)pos) -1);
                sweepQueueStarters[(unsigned int)(pos >> 32)][((unsigned int)pos) -1] = idx;
            }
        }
    }
}

__global__ void compress_UF(int numVertices, array3D<float> values, int* swept, int* augmentation){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] >= 0){
            if (swept[swept[idx]] >= 0){
                if (augmentation[idx] == -1){
                    if ((values.get(idx) < values.get(swept[swept[idx]])) || ((values.get(idx) == values.get(swept[swept[idx]])) && (idx < swept[swept[idx]])) ){
                        augmentation[idx] = swept[idx];
                    }
                }
                swept[idx] = swept[swept[idx]];
            }
        }
    }
}

__global__ void get_saddles(int numVertices, array3D<float> values, int* neighborsBuffer, NeighborsEntry* neighborsMap, int* minimaQueue, int* swept, unsigned long long int* saddleMap){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    unsigned long long int saddle;

    for (; idx < minimaCounter; idx+=gridDim.x*blockDim.x){
       if (swept[minimaQueue[idx]] < -1){
           minimaQueue[idx] = minimaQueue[idx];
       } else {
           saddle = saddleMap[minimaQueue[idx]];
           swept[minimaQueue[idx]] = saddle;
           minimaQueue[idx] = saddle;
       }
    }

}

__global__ void flatten_swept(int numVertices, array3D<float> values, int* swept){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] >= 0){
            if (values.get(swept[idx]) < values.get(idx) || ((values.get(swept[idx]) == values.get(idx)) && (swept[idx] < idx))){
                swept[idx] = -1;
            }
        }
    }
}

__global__ void trunk_augment(int numVertices, int dangcount, array3D<float> values, int* swept, int* dangling){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    for(; idx<numVertices; idx+=gridDim.x*blockDim.x){
        if (swept[idx] == -1){
            if (values.get(idx) < values.get(trunkStart)){
               swept[idx] = trunkStart;
            } else if (values.get(idx) > values.get(dangling[dangcount -1])){
                swept[idx] = dangling[dangcount-1];
            } else {
                int min = 0;
                int max = dangcount-1;
                while (true){
                    int j = (max+min)/2;
                    if (values.get(idx) > values.get(dangling[j])){
                        if (values.get(idx) < values.get(dangling[j+1])){
                            swept[idx] = dangling[j];
                            break;
                        } else {
                            min = j+1;
                            if (min > max){
                                swept[idx] = dangling[min];
                                break;
                            }
                        }
                    } else {
                        max = j-1;
                        if (min > max){
                            swept[idx] = dangling[min];
                            break;
                        }
                    }
                }
            }
        }
    }
}

    array3D<float> valuest;

void augmentTrunk(TreeConstructor &tree, int trunkStarter, std::vector<int> &danglingvect, int *&swept){

    int* sweptd;
    cudaMalloc(&sweptd, tree.numVertices*sizeof(int));
    cudaMemcpy(sweptd,  swept, tree.numVertices*sizeof(int), cudaMemcpyHostToDevice);

    int* dangling;
    cudaMalloc(&dangling, danglingvect.size()*sizeof(int));
    cudaMemcpy(dangling, &*danglingvect.begin(), danglingvect.size()*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(trunkStart, &trunkStarter, sizeof(int));

    trunk_augment<<<18,1024>>>(tree.numVertices, danglingvect.size(), valuest, sweptd, dangling);

    cudaMemcpy(swept, sweptd, tree.numVertices*sizeof(int), cudaMemcpyDeviceToHost);
}

void doCudaStuff(TreeConstructor& tree, std::vector<int*>& minima, std::vector<int*>& saddles, std::vector<int>& counts, int *&swepto, int*& augmentationo, int rounds, int minimaroundsh, int growthroundsh)
{
    int mincount = 0;
    int constexpr zero = 0;

    int*                    neighborsBuffer;
    NeighborsEntry*         neighborsMap;
    int*                    swept;
    int*                    augmentation;
    int*                    minimaQueue;
    unsigned long long int* saddleMap;
    int**                   sweepQueueStarters;

    valuest.resize(dim3(tree.max_x, tree.max_y, tree.max_z));
    valuest.copy(tree.fvalues);

    //Fill device neighbors
    cudaMalloc(             &neighborsBuffer,       tree.neighborsBuffer->size()*sizeof(int));
    cudaMemcpy(             neighborsBuffer,        &*tree.neighborsBuffer->begin(), tree.neighborsBuffer->size()*sizeof(int), cudaMemcpyHostToDevice);
    //Fill device neighborMap
    cudaMalloc(             &neighborsMap,          tree.numVertices*sizeof(NeighborsEntry));
    cudaMemcpy(             neighborsMap,           &*tree.neighborsMap->begin(), tree.numVertices*sizeof(NeighborsEntry), cudaMemcpyHostToDevice);
    //Reserve swept, minimaQueue and saddleMap
    cudaMalloc(             &swept,                 tree.numVertices*sizeof(int));
    cudaMalloc(             &augmentation,          tree.numVertices*sizeof(int));
    cudaMalloc(             &minimaQueue,           tree.numVertices*sizeof(int));
    cudaMalloc(             &saddleMap,             tree.numVertices*sizeof(unsigned long long int));
    cudaMalloc(             &sweepQueueStarters,    tree.numVertices*sizeof(int*));

    cudaMemcpyToSymbol(minimaQueueCounter, &mincount, sizeof(int));
    cudaMemcpyToSymbol(doneCounter, &zero, sizeof(int));
    cudaMemcpyToSymbol(min_value, &tree.min_val, sizeof(double));
    cudaMemcpyToSymbol(conversionFactor, &tree.conversionFactor, sizeof(double));
    cudaMemcpyToSymbol(max_x, &tree.max_x, sizeof(int));
    cudaMemcpyToSymbol(max_y, &tree.max_y, sizeof(int));
    cudaMemcpyToSymbol(max_z, &tree.max_z, sizeof(int));
    cudaMemcpyToSymbol(minimarounds, &minimaroundsh, sizeof(int));
    cudaMemcpyToSymbol(growthrounds, &growthroundsh, sizeof(int));



    //Start First Round
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_secs;

    //Search for start points
    min_search_preNeigh<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, augmentation, saddleMap);
    //min_search_<<<18,1024>>>(tree.numVertices, values, minimaQueue, swept, saddleMap);

    //Read number of found start points
    cudaMemcpyFromSymbol(&mincount, minimaQueueCounter, sizeof(int));

    //Read start points from device
    minima.push_back((int*)malloc(mincount*sizeof(int)));
    cudaMemcpy(minima.at(minima.size()-1), minimaQueue, mincount*sizeof(int), cudaMemcpyDeviceToHost);
    counts.push_back(mincount);

    //Remember number of start points on device
    cudaMemcpyToSymbol(minimaCounter, &mincount, sizeof(int));

    //Sweep from all starts
    //min_growth<<<18,1024>>>(tree.numVertices, values, minimaQueue, swept, saddleMap, sweepQueueStarters);
    min_growth_preNeigh<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap, sweepQueueStarters);

    //Clear minimaQueueCounter
    cudaMemcpyToSymbol(minimaQueueCounter, &zero, sizeof(int));

    //Search for end points
    saddle_search<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap);

    //prepare end points
    get_saddles<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap);

    //Read end points from device
    saddles.push_back((int*)malloc(mincount*sizeof(int)));
    cudaMemcpy(saddles.at(saddles.size()-1), minimaQueue, mincount*sizeof(int), cudaMemcpyDeviceToHost);

    //Search for start points
    compress_saddles<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap, sweepQueueStarters);

    //Read number of found start points
    cudaMemcpyFromSymbol(&mincount, minimaQueueCounter, sizeof(int));

    int i = 0;
    while (i < rounds){
        i++;
        //Read start points from device
        minima.push_back((int*)malloc(mincount*sizeof(int)));
        cudaMemcpy(minima.at(minima.size()-1), minimaQueue, mincount*sizeof(int), cudaMemcpyDeviceToHost);
        counts.push_back(mincount);

        //Remember number of start points on device
        cudaMemcpyToSymbol(minimaCounter, &mincount, sizeof(int));

        //Prepare sweepQueueStarters
        clean_saddleMap_to_zero<<<18,1024>>>(tree.numVertices, swept, saddleMap);

        size_sweepqueues<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap, sweepQueueStarters);
        create_sweepqueues<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap, sweepQueueStarters);
        fill_sweepqueues<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap, sweepQueueStarters);

        //Compress UF
        compress_UF<<<18,1024>>>(tree.numVertices, valuest, swept, augmentation);

        //cudaMemcpy(swepto, swept, tree.numVertices*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpyToSymbol(doneCounter, &zero, sizeof(int));
        //Sweep from all starts
        saddle_growth<<<36,512>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap, sweepQueueStarters, i);
        //min_growth<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept);
        //host_min_growth(tree.numVertices, mincount, tree.valuest, &*tree.neighborsBuffer->begin(), &*tree.neighborsMap->begin(), minima.at(i), swepto);

        //Clear minimaQueueCounter
        cudaMemcpyToSymbol(minimaQueueCounter, &zero, sizeof(int));

        //cudaMemcpy(swept, swepto, tree.numVertices*sizeof(int), cudaMemcpyHostToDevice);

        //Search for end points
        clean_saddleMap_to_infinity<<<18,1024>>>(tree.numVertices, swept, saddleMap);

        saddle_search<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap);

        //prepare end points
        get_saddles<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap);

        //Read end points from device
        saddles.push_back((int*)malloc(mincount*sizeof(int)));
        cudaMemcpy(saddles.at(saddles.size()-1), minimaQueue, mincount*sizeof(int), cudaMemcpyDeviceToHost);

        //Search for start points
        compress_saddles<<<18,1024>>>(tree.numVertices, valuest, neighborsBuffer, neighborsMap, minimaQueue, swept, saddleMap, sweepQueueStarters);

        //Read number of found start points
        cudaMemcpyFromSymbol(&mincount, minimaQueueCounter, sizeof(int));

    }

    minima.push_back((int*)malloc(mincount*sizeof(int)));
    cudaMemcpy(minima.at(minima.size()-1), minimaQueue, mincount*sizeof(int), cudaMemcpyDeviceToHost);
    counts.push_back(mincount);

    cudaMemcpyToSymbol(minimaCounter, &mincount, sizeof(int));

    compress_UF<<<18,1024>>>(tree.numVertices, valuest, swept, augmentation);

    flatten_swept<<<18,1024>>>(tree.numVertices, valuest, swept);

    swepto = (int*)malloc(tree.numVertices * sizeof(int));
    cudaMemcpy(swepto, swept, tree.numVertices * sizeof(int), cudaMemcpyDeviceToHost);

    augmentationo = (int*)malloc(tree.numVertices * sizeof(int));
    cudaMemcpy(augmentationo, augmentation, tree.numVertices * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(neighborsBuffer);
    cudaFree(neighborsMap);
    cudaFree(swept);
    cudaFree(augmentation);
    cudaFree(minimaQueue);
    cudaFree(saddleMap);
    cudaFree(sweepQueueStarters);

    elapsed_secs = std::chrono::high_resolution_clock::now() - start;
    //cudaMemcpy(swepto, swept, tree.numVertices*sizeof(int), cudaMemcpyDeviceToHost);

    std::cout << "Device Time: " << elapsed_secs.count() << std::endl;
}

#endif
