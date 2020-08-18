#include "TreeConstructor.h"
#include "SweepQueue.h"
#include <boost/range/irange.hpp>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/execution_policy.hpp>
#include <hpx/include/parallel_sort.hpp>

HPX_REGISTER_COMPONENT_MODULE();

typedef hpx::components::simple_component<server::TreeConstructor> TreeConst_type;

HPX_REGISTER_COMPONENT(TreeConst_type, TreeConstructor);
HPX_REGISTER_ACTION(TreeConst_type::wrapped_type::init_action, treeConstructor_init_action);
HPX_REGISTER_ACTION(TreeConst_type::wrapped_type::construct_action, treeConstructor_construct_action);
HPX_REGISTER_ACTION(TreeConst_type::wrapped_type::startSweep_action, treeConstructor_startSweeps_action);

namespace server {

void TreeConstructor::init(const Options& opt)
{
    data = Datamanager(opt.resampleX, opt.resampleY, opt.resampleZ);
    DataValueComparator::data = &data;
    ArcDataValueComparator::data = &data;

    vis.init(data);

    sweeps.store(0);
    this->visualize = opt.visualize;
    this->trunkSkipping = opt.trunkSkipping;

    data.readFromFile(opt.inputFile);
    numVertices = data.getNumVertices();

    swept.resize(numVertices, -1);
    UF.resize(numVertices, -1);

    arcMap.resize(numVertices, nullptr);
    done.reset();
}

void TreeConstructor::construct()
{
    hpx::util::high_resolution_timer t;
    t.restart();

    /*
     * Start sweeps at minima
     */
    int numMinima = 0;
    for (int v = 0; v < this->numVertices; ++v) {
        if (this->isMinimum(v)) {
            ++numMinima;
            std::atomic_fetch_add(&sweeps, 1);
            hpx::apply<TreeConstructor::startSweep_action>(this->get_id(), v, true);
        }
    }

    std::cout << "Minima: " << numMinima << std::endl;

    leavesDone.set();
    done.wait();


    std::vector<Arc*> arcList;
    for (int i = 0; (i < this->numVertices); i++) {
        if (arcMap.at(i) != nullptr && arcMap.at(i)->saddle == -1) {
            arcList.push_back(arcMap.at(i));
        }
    }

    std::sort(arcList.begin(), arcList.end(), ArcDataValueComparator::compare);


#pragma omp parallel for
    for (int i = 0; (i < arcList.size() - 1); i++) {
        arcList[i]->saddle = arcList[i + 1]->extremum;

        if (visualize) {
            vis.lock.lock();
            vis.visualizeBranch(arcList[i], false);
            vis.lock.unlock();
        }
    }

    for (int i = 0; (i < numVertices); i++){
        swept[i] = -1;
    }

#pragma omp parallel for
    for (int i = 0; (i < numVertices); i++){
        if (arcMap[i] != nullptr){
            Arc* arc = arcMap[i];
            for (auto iter = arc->augmentation.vertices.begin(); (iter != arc->augmentation.vertices.end()); iter.operator ++()){
                if (swept[(*iter)->key.key] >= 0){
                    if (arcMap[swept[(*iter)->key.key]]->extremum < arc->extremum)
                        swept[(*iter)->key.key] = arc->extremum;
                } else {
                    swept[(*iter)->key.key] = arc->extremum;
                }
            }
        }
    }

    hpx::parallel::for_each(
        hpx::parallel::execution::parallel_unsequenced_policy(),
        boost::irange(0, (int)this->numVertices).begin(), boost::irange(0, (int)this->numVertices).end(),
        [&](std::size_t i)
        {
            if (swept[i] == -1){
                if (data.getValue(i) < data.getValue(arcList[0]->extremum)){
                   swept[i] = trunkStart->extremum;
                } else if (data.getValue(i) > data.getValue(arcList[arcList.size()-1]->extremum)){
                    swept[i] = arcList[arcList.size()-1]->extremum;
                } else {
                    int min = 0;
                    int max = arcList.size()-1;
                    while (true){
                        int j = (max+min)/2;
                        if (data.getValue(i) > data.getValue(arcList[j]->extremum)){
                            if (data.getValue(i) < data.getValue(arcList[j+1]->extremum)){
                                swept[i] = arcList[j]->extremum;
                                break;
                            } else {
                                min = j+1;
                                if (min > max){
                                    swept[i] = arcList[min]->extremum;
                                    break;
                                }
                            }
                        } else {
                            max = j-1;
                            if (min > max){
                                swept[i] = arcList[min]->extremum;
                                break;
                            }
                        }
                    }
                }
            }
        }
    );
/*
    for (int i = 0; (i < numVertices); i++){

        int value = (data.getValue(swept[i]).value - data.min_value)/(data.max_value-data.min_value)*10000;
        vis.createBoundaryActor(i, (value % 5)/5.0, (value % 2)/2.0, (value % 11)/11.0, 1);
    }
*/

    std::cout << "Tree Constructed: " << t.elapsed() << " s" << std::endl;

    if (visualize)
        vis.show();

    // Print some stats
    // Number of arcs
    int numArcs = 0;
    for (Arc* arc : this->arcMap)
        if (arc != nullptr)
            ++numArcs;

    std::cout << "Arcs: " << numArcs << std::endl;
}

void TreeConstructor::startSweep(int v, bool leaf)
{
    Arc* arc;
    if (leaf) {
        arc = new Arc(v);
        mapLock.lock();
        arcMap.at(v) = arc;
        mapLock.unlock();
        arc->leaf = true;
    } else {
        mapLock.lock();
        arc = arcMap.at(v);
        mapLock.unlock();
    }

    SweepQueue queue = SweepQueue(v, this); // v is not on sweep queue!
    Boundary& boundary = arc->boundary;

    swept[v] = v;

    int numNeighbors;
    const int* neighbors = data.getNeighbors(v, &numNeighbors);
    queue.push(neighbors, numNeighbors);

    std::vector<Arc*>& children = arc->children;

    // Mutual intersections of child boundaries -> this is where sweeps must continue
    for (int i = 0; i < children.size(); i++) {
        Arc* child = children.at(i);
        for (int j = i + 1; j < children.size(); j++) {
            Arc* child2 = children.at(j);
            queue.push(child->boundary.intersect(child2->boundary));
        }
        boundary.unite(child->boundary); // intersections are removed in above intersect call; only non-queued vertices become boundary
    }

    arc->augmentation.inherit(arc->inherited_augmentations); // gathered and merged here because no lock required here
    arc->augmentation.sweep(data.getValue(v)); // also add saddle/local minimum to augmentation

    // actual sweep happens here
    while (!queue.empty()) {
        int c = queue.pop();
        if (c == -1) // empty (queue skips elements which have been swept in the mean time)  -> could move skip check to empty()
            break;

        // check if can be swept
        if (touch(c, v)) {
            boundary.remove(data.getValue(c)); // remove from boundary

            // put into our augmentation and mark as swept
            swept[c] = v;
            arc->augmentation.sweep(data.getValue(c));

            int numNeighbors;
            const int* neighbors = data.getNeighbors(c, &numNeighbors);
            queue.push(neighbors, numNeighbors);
        } else {
            // not allowed to sweep yet -> add as potential boundary candidate
            boundary.add(data.getValue(c));
        }
    }

    // if boundary is empty after sweep, we are the last sweep
    if (boundary.empty()) {
        assignSaddle(arc, false);
        done.set(); // complete!
    } else {
        // non-empty boundary ->  minimum boundary vertex becomes saddle of this sweep
        arc->saddle = boundary.min(data);
        boundary.remove(data.getValue(arc->saddle));

        assignSaddle(arc, leavesDone.occurred());

        // check if we can continue sweep at this saddle
        if (checkSaddle(arc->saddle)) {
            if (trunkSkipping && leavesDone.occurred() && sweeps.load() == 1) {
                trunkStart = arc;
                done.set();
            } else {
                hpx::apply<TreeConstructor::startSweep_action>(this->get_id(), arc->saddle, false);
            }
        } else {
            std::atomic_fetch_add(&sweeps, -1);
        }
    }
}

/*
 *
 */
bool TreeConstructor::searchUF(int start, int goal)
{
    if (start == -1)
        return false;
    int c = start;
    int next = UF[c];
    if (c == goal || next == goal)
        return true;
    if (next == -1)
        return false;

    // includes path compression
    while (true) {
        if (UF[next] == goal) {
            UF[c] = UF[next];
            return true;
        }

        if (UF[next] == -1)
            return false;

        UF[c] = UF[next];
        next = UF[next];
    }
}

/*
 * Sweep reaches vertex and checks if it can be swept by going through *all* its smaller neighbors and check if they *all* have already been swept by us
 */
bool TreeConstructor::touch(int c, int v)
{
    ctValue cValue = data.getValue(c);

    int numNeighbors;
    const int* neighbors = data.getNeighbors(c, &numNeighbors);
    for (int i = 0; i < numNeighbors; ++i) {
        if (data.getValue(neighbors[i]) <= cValue) {
            if (!searchUF(swept[neighbors[i]], v)) {
                return false;
            }
        }
    }
    return true;
}

/*
 * Checks if we're the last child arc to reach saddle -> only the last one can continue sweeping at this saddle
 */
bool TreeConstructor::checkSaddle(int s)
{
    mapLock.lock();
    Arc* branch = arcMap.at(s);
    mapLock.unlock();

    branch->lock.lock(); // only used by checkSaddle -> have to be synchronized , so that only one continues sweep

    if (branch->started) {
        branch->lock.unlock();
        return false;
    }
    ctValue sValue = data.getValue(s);

    int numNeighbors;
    const int* neighbors = data.getNeighbors(s, &numNeighbors);
    for (int i = 0; i < numNeighbors; ++i) {
        if (data.getValue(neighbors[i]) <= sValue) {
            // go through all smaller neighbors
            // if not swept or swept by someone else
            if (swept[neighbors[i]] == -1 || (!searchUF(swept[neighbors[i]], s))) {
                branch->lock.unlock();
                return false;
            }
        }
    }
    branch->started = true;
    branch->lock.unlock();
    return true;
}

bool TreeConstructor::isMinimum(int v) const
{
    ctValue vValue = data.getValue(v);
    int numNeighbors;
    const int* neighbors = data.getNeighbors(v, &numNeighbors);
    for (int i = 0; i < numNeighbors; ++i)
        if (data.getValue(neighbors[i]) < vValue)
            return false;

    return true;
}

/*
 * Checks in arc map if saddle is already represented: Either saddle's arc is constructed (we're first level set component to reach this saddle), or we take existing arc (if someone else already created it).
 */
void TreeConstructor::assignSaddle(Arc* finishedArc, bool trunk)
{
    mapLock.lock();

    Arc* saddleArc; // parent of finished arc
    if (nullptr == arcMap.at(finishedArc->saddle)) {
        saddleArc = new Arc(finishedArc->saddle);
        arcMap.at(finishedArc->saddle) = saddleArc;
    } else
        saddleArc = arcMap.at(finishedArc->saddle);

    saddleArc->children.push_back(finishedArc);

    saddleArc->inherited_augmentations.push_back(finishedArc->augmentation.heritage(data.getValue(finishedArc->saddle)));

    mapLock.unlock();

    // Update union find structure: finished arc belongs to saddle arc
    UF[finishedArc->extremum] = saddleArc->extremum;

    if (visualize) {
        vis.lock.lock();
        vis.visualizeBranch(finishedArc, trunk);
        vis.lock.unlock();
    }
}
}
