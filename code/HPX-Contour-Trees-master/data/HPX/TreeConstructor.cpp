#include "TreeConstructor.h"
#include "DataManager.h"
#include "Log.h"

#include "Arc.h"
#include "Augmentation.h"
#include "Boundary.h"
#include "SweepQueue.h"

HPX_REGISTER_COMPONENT_MODULE();

typedef hpx::components::simple_component<TreeConstructor> TreeConstructorComponent;

HPX_REGISTER_COMPONENT(TreeConstructorComponent, TreeConstructor);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::init_action, treeConstructor_init_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::construct_action, treeConstructor_construct_action);
HPX_REGISTER_ACTION(TreeConstructorComponent::wrapped_type::destroy_action, treeConstructor_destroy_action);

/**
 * @brief Initialize data.
 * @param treeConstructors
 * @param input
 */
void TreeConstructor::init(const std::vector<hpx::id_type>& treeConstructors, const std::string& input)
{
    // Store list of tree constructor components and determine index of this component
    this->treeConstructors = treeConstructors;
    for (uint32_t i = 0; i < this->treeConstructors.size(); ++i) {
        if (this->treeConstructors[i] == this->get_id()) {
            this->index = i;
            break;
        }
    }

    // Load data
    this->dataManager = DataManager::load(input, this->index, this->treeConstructors.size());

    // Init structures
    this->sweeps.store(0u);

    uint64_t numVertices = this->dataManager->getNumVerticesLocal(true); // TODO

    this->swept.resize(numVertices, INVALID_VERTEX);
    this->UF.resize(numVertices, INVALID_VERTEX);
    this->arcMap.resize(numVertices, nullptr);

    this->done.reset();
}

/**
 * @brief Shutdown.
 */
void TreeConstructor::destroy()
{
    if (this->dataManager)
        delete this->dataManager;
}

/**
 * @brief Construction main entry point.
 */
void TreeConstructor::construct()
{
    hpx::util::high_resolution_timer timer;

    // Find minima in local data
    uint32_t numMinima = 0;
    std::vector<uint64_t> localVertices = this->dataManager->getLocalVertices();

    for (uint64_t v : localVertices) {
        if (this->dataManager->isMinimum(v)) {
            std::atomic_fetch_add(&this->sweeps, 1u);
            hpx::apply<TreeConstructor::startSweep_action>(this->get_id(), v, true);
            ++numMinima;
        }
    }

    Log().tag(std::to_string(this->index)) << "Minima (local): " << numMinima;

    this->leavesDone.set();
    this->done.wait();

    /*
     * Construct trunk
     *
     * TODO:
     *  - parallelize
     *  - assign non-augmented vertices to correct arc/interval in trunk
     */
    std::vector<Arc*> arcList;
    for (Arc* arc : arcMap)
        if (arc != nullptr && arc->saddle == INVALID_VERTEX)
            arcList.push_back(arc);

    std::sort(arcList.begin(), arcList.end(),
        [&](const Arc* a1, const Arc* a2) {
            return this->dataManager->getValue(a1->extremum) < this->dataManager->getValue(a2->extremum);
        });

    for (uint64_t i = 0; i < arcList.size() - 1; ++i)
        arcList[i]->saddle = arcList[i + 1]->extremum;

    Log() << "Tree Constructed: " << timer.elapsed() << " s";

    // Print some stats
    // Number of arcs
    uint32_t numArcs = 0;
    for (Arc* arc : this->arcMap)
        if (arc != nullptr)
            ++numArcs;

    Log() << "Arcs: " << numArcs;
}

/**
 * @brief Starts a new sweep at the given local vertex.
 * @param v
 * @param leaf
 */
void TreeConstructor::startSweep(uint64_t v, bool leaf)
{
    Arc* arc;
    if (leaf) {
        arc = new Arc(v);
        this->mapLock.lock();
        this->arcMap[v] = arc;
        this->mapLock.unlock();
        arc->leaf = true;
    } else {
        this->mapLock.lock();
        arc = this->arcMap[v];
        this->mapLock.unlock();
    }

    SweepQueue queue(v, this); // v is not on sweep queue!
    // TODO distributed: move queue into arc

    Boundary& boundary = arc->boundary;

    this->swept[v] = v;

    uint64_t neighbors[6];
    uint32_t numNeighbors = this->dataManager->getNeighbors(v, neighbors);

    queue.push(neighbors, numNeighbors);

    std::vector<Arc*>& children = arc->children;

    // Mutual intersections of child boundaries -> this is where sweeps must continue
    for (uint32_t i = 0; i < children.size(); ++i) {
        Arc* child = children[i];
        for (uint32_t j = i + 1; j < children.size(); ++j) {
            Arc* child2 = children[j];
            queue.push(child->boundary.intersect(child2->boundary));

            // TODO distributed:
            // - maintain list of all localities that have contributed to at least one of the children
            // - send each list of all children's ids
            // - each locality: for each arc lookup boundary -> mutual intersection -> push to sweep queue -> start sweep with that queue (id: v)
        }
        boundary.unite(child->boundary); // intersections are removed in above intersect call; only non-queued vertices become boundary
    }

    arc->augmentation.inherit(arc->inheritedAugmentations); // gathered and merged here because no lock required here
    arc->augmentation.sweep(this->dataManager->getValue(v)); // also add saddle/local minimum to augmentation

    // actual sweep happens here
    while (!queue.empty()) {
        uint64_t c = queue.pop();
        if (c == INVALID_VERTEX) // empty (queue skips elements which have been swept in the mean time)  -> could move skip check to empty()
            break;

        const Value value = this->dataManager->getValue(c);

        // check if can be swept
        if (this->touch(c, v)) {
            boundary.remove(value); // remove from boundary

            // put into our augmentation and mark as swept
            this->swept[c] = v;
            arc->augmentation.sweep(value);

            numNeighbors = this->dataManager->getNeighbors(c, neighbors);

            // if ghost:  start or continue remote sweep
            // - parameter: v (where sweep started)
            // - return value: id + value of local minimum (which potentially is saddle of this sweep) + all involved localities
            // - check if that sweep already running (or already terminated) (check in arc map)
            // - if yes: push to that sweep queue (and return invalid local minimum)
            // - if terminated: start queue-not-empty loop again (continue its work) and return minimum after queue is empty
            // - if no: start new sweep queue and return minimum after queue is empty

            queue.push(neighbors, numNeighbors);
        } else {
            // not allowed to sweep yet -> add as potential boundary candidate
            boundary.add(value);
        }
    }

    // TODO distributed:
    // - when queue is empty and all remote calls/sweeps have returned
    // merge lists of remote calls / add self
    // find minimum over all -> saddle if v is local, else return as above

    // if boundary is empty after sweep, we are the last sweep
    if (boundary.empty()) {
        this->assignSaddle(arc);
        this->done.set(); // complete!
    } else {
        // non-empty boundary ->  minimum boundary vertex becomes saddle of this sweep
        arc->saddle = boundary.min();
        boundary.remove(this->dataManager->getValue(arc->saddle));

        this->assignSaddle(arc);

        // check if we can continue sweep at this saddle
        if (this->checkSaddle(arc->saddle)) {
            if (this->leavesDone.occurred() && this->sweeps.load() == 1u) {
                this->done.set();
            } else {
                hpx::apply<TreeConstructor::startSweep_action>(this->get_id(), arc->saddle, false);
            }
        } else {
            std::atomic_fetch_sub(&this->sweeps, 1u);
        }
    }
}

/*
 *
 */
bool TreeConstructor::searchUF(uint64_t start, uint64_t goal)
{
    // goal is actively running sweep

    // TODO distributed:
    // - go through chain
    // - if "cache miss" (no entry exists):
    // if locality of required id is here: return false / and one which caused cache miss (which is still running) for others to path-compress
    // else: send query to locality of required id -> wait for answer -> cache result locally

    if (start == INVALID_VERTEX)
        return false;

    uint64_t c = start;
    uint64_t next = this->UF[c];
    if (c == goal || next == goal)
        return true;
    if (next == INVALID_VERTEX)
        return false;

    // includes path compression
    while (true) {
        if (this->UF[next] == goal) {
            this->UF[c] = this->UF[next];
            return true;
        }

        if (this->UF[next] == INVALID_VERTEX)
            return false;

        this->UF[c] = this->UF[next];
        next = this->UF[next];
    }
}

/*
 * Sweep reaches vertex and checks if it can be swept by going through *all* its smaller neighbors and check if they *all* have already been swept by us
 */
bool TreeConstructor::touch(uint64_t c, uint64_t v)
{
    const Value value = this->dataManager->getValue(c);

    uint64_t neighbors[6];
    uint32_t numNeighbors = this->dataManager->getNeighbors(c, neighbors);

    for (uint32_t i = 0; i < numNeighbors; ++i) {
        if (neighbors[i] != INVALID_VERTEX) {
            if (this->dataManager->getValue(neighbors[i]) <= value) {
                if (!this->searchUF(this->swept[neighbors[i]], v)) {
                    return false;
                }
            }
        }
    }
    return true;
}

/*
 * Checks in arc map if saddle is already represented: Either saddle's arc is constructed (we're first level set component to reach this saddle), or we take existing arc (if someone else already created it).
 */
void TreeConstructor::assignSaddle(Arc* finishedArc)
{
    this->mapLock.lock();
    Arc* saddleArc = arcMap[finishedArc->saddle]; // parent of finished arc
    if (saddleArc == nullptr) {
        saddleArc = new Arc(finishedArc->saddle);
        arcMap[finishedArc->saddle] = saddleArc;
    }
    this->mapLock.unlock();

    saddleArc->lock.lock();
    saddleArc->children.push_back(finishedArc);
    saddleArc->inheritedAugmentations.push_back(finishedArc->augmentation.heritage(this->dataManager->getValue(finishedArc->saddle)));
    saddleArc->lock.unlock();

    // Update union find structure: finished arc belongs to saddle arc
    this->UF[finishedArc->extremum] = saddleArc->extremum; // TODO distributed: nothing changes here
}

/*
 * Checks if we're the last child arc to reach saddle -> only the last one can continue sweeping at this saddle
 */
bool TreeConstructor::checkSaddle(uint64_t s)
{
    this->mapLock.lock();
    Arc* branch = this->arcMap.at(s);
    this->mapLock.unlock();

    // Make sure only one checkSaddle call continues sweep
    std::lock_guard<hpx::lcos::local::mutex> lock(branch->lock);

    if (branch->started)
        return false;

    const Value value = this->dataManager->getValue(s);

    uint64_t neighbors[6];
    uint32_t numNeighbors = this->dataManager->getNeighbors(s, neighbors);

    // Go through all smaller neighbors
    for (uint32_t i = 0; i < numNeighbors; ++i) {
        if (neighbors[i] != INVALID_VERTEX) {
            if (this->dataManager->getValue(neighbors[i]) <= value) {
                // if not swept or swept by someone else
                if (this->swept[neighbors[i]] == INVALID_VERTEX || (!this->searchUF(this->swept[neighbors[i]], s)))
                    return false;
            }
        }
    }

    branch->started = true;
    return true;
}
