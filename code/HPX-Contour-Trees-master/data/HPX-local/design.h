#ifndef DESIGN_H
#define DESIGN_H

#include <hpx/hpx.hpp>


/**
 * @brief Non-AGAS, but globally unique identifier of a vertex in the complete Mesh
 * OR!!!
 * Local; Identifies a vertex within a Mesh (!!! Conversion of local IDs)
 */
struct ID{

};


/**
 * @brief Represents morse function value of vertices (simulation of simplicity!)
 */
struct ValType{

};


/**
 * @brief Component representing a Sweep in the AGAS
 * Holding AGAS-Adresses of sweeps leading to it (child branches)
 */
class /*component*/ Arc {
    //Threadsafe list of AGAS adresses holding all descendants (child branches)
    //Knows own AGAS adress
};

/**
 * @brief Local; List of vertex IDs on this locality that are saddle candidates for arc (or were saddle candidates for a child arc)
 * Allows for fast access to minimal element
 * Allows for fast lookup with vertex ID as key
 * No Thread-Safeness needed: Only one sweep pushes, only local calls, can be used for merging into father boundary only after locally done
 */
class Boundary {

};


/**
 * @brief Local; One-Per-Locality; Holds all Boundaries that are known on this locality for a certain arc
 * Allows for fast lookup with Arc-AGAS adress as key
 * Allows for fast merging of multiple child boundaries
 * Threadsafe Lookup vs. Insertion
 */
class BoundaryMap {

};


/**
 * @brief Local; List of vertex IDs on this locality that still need to be visited for arc
 * Thread safe Insertion vs. Deletion
 */
class SweepFront {

};


/**
 * @brief Local; One-Per-Locality; Holds all Sweep-Fronts that are known on this locality for a certain arc
 * Allows for fast lookup with Arc-AGAS adress as key
 * Threadsafe Lookup vs. Insertion
 */
class SweepFrontMap {

};

/**
 * @brief Local; One-Per-Locality; Holds data + ghost cells, can identify neighbors and morse value (simulation of simplicity!)
 * Needs to be thread-safe after creation. Is mainly read-only though.
 */
class /*VTK-Inside*/ Mesh {

};

/**
 * @brief Read from VTK-File OR represents VTK-File, holds complete input data and allows for construction of local meshes with ghost cells
 * Possibly allows for globally unique identification of vertices
 */
class /*VTK-Inside*/ GlobalMesh {

};


/**
 * @brief Extracts this localities partition + ghost cells from overall input mesh
 * @param localityCount Number of hpx localities to share inputData with
 * @param localityID HPX-ID of local hpx locality
 * @return Mesh object representing the local partition of the input Data + Ghost Cells
 */
Mesh* partition(int localityCount, int localityID, GlobalMesh* inputData);

/**
 * @brief Determines locally (possibly with ghost cells) if vertex id is a local minimum
 * @param mesh Mesh of my locality (must contain id OR id must be valid)
 * @param id Id of vertex to check, must never be a ghost cell
 * @return True: id is local extremum; False: id is not local extremum
 */
bool isExtremum(Mesh* mesh, ID id);

/**
 * @brief Creates an Arc with id as extremum or points to the Arc already existing for the id
 * @param id of the vertex that is the extremum of the Arc, must never be a ghost cell
 * @return Pointer to newly created or found Arc
 */
Arc* createOrFindArc(ID id);

/**
 * @brief
 * Szenario1: arc is colocated and newly created for local minimum
 * Creates empty boundary to hold boundary-candidates
 * Szenario2: arc already has coresponding boundary on this locality
 * Points to already known boundary
 * Szenario3: arc is colocated or copy of remote Arc and has no coresponding boundary on this locality yet
 * Create new boundary that is union of all boundaries on this locality for children of arc
 * @param arc reference to colocated Arc or copy of remote Arc that needs a boundary for this locality
 * @return Boundary containing all elements of locally known boundaries of child arcs of arc
 */
Boundary* createOrFindBoundary(Arc const& arc);

/**
 * @brief Used to create a SweepFront for an arc whose extremum is on this locality
 * or used to create a SweepFront after
 * @param arc
 * @return
 */
SweepFront* createSweepFront(Arc const& arc);

/**
 * @brief Called on locality 0 only, collects the arc and stores it for use as merge tree or for contour tree composition
 * @param arc copy of arc that should be stored
 */
/*HPX-Remote-Task*/
void recogniseArc(Arc arc);

/**
 * @brief Called on locality 0 only, collects the final arc of a merge tree and issues all localities to fix their augmentations for that tree
 * also checks wether the other trunk has arrived and contour tree composition can start
 * @param arc
 */
/*HPX-Remote-Task*/
void recogniseTrunk(Arc arc);

/*HPX-Remote-Task*/
/**
 * @brief Called on each locality from locality 0 after contour tree composition is done and triggers an event on which the main functions are waiting before termination
 */
void terminate();

/* All HPX-Remote-Tasks need to wait for event that signals local mesh is loaded */
/**
 * @brief Starts receiver-local sweep for arc which has a sender-located saddle as extremum and needs to continue at boundaries of all descendants (in this case remote)
 * @param arc sender-located arc spawned from a saddle holding receiver-located descendants whose boundaries need to start a sweep
 * @return minimal saddle candidate reached by local sweep (and recursively, the neighbors sweeps)
 */
/*HPX-Remote-Task*/
ValType continueSweep(Arc arc);


/**
 * @brief Starts receiver-local sweep for arc whose sweep has reached the senders ghost cells
 * Recursion to neighbor-localities' neighbor-localities! Localities might be visited multiple times from multiple neighbors!
 * --> if no sweep for this arc has started start one (variant of startSweep)
 * --> if a sweep for this arc is already in progress return infinity (but still add id to sweepfront)
 * --> if a sweep for this arc is already finished continue it from id and return new minimum
 * Updating neighbors ghost-cells.visited so he can see that no message is necessary (?)
 * @param arc arc for wich to continue the sweep on this locality. it's descendants are used to identify all visited vertices and left-over boundaries
 * @param id vertex which should be added to local SweepFront of arc
 * @return minimal saddle candidate reached by local sweep (and recursively, the neighbors sweeps)
 */
/*HPX-Remote-Task*/
ValType continueSweep(Arc arc, ID id);


/**
 * @brief finishSweep Either calls finishSweep on locality which reported the minimal saddle candidate towards this locality
 * or, if this locality found the minimal saddle candidate locally, calls linkToSaddle on that saddle
 * @param arc copy of arc whose sweep is entirely finished and whose minimal saddle candidate was reported here (or on recursively visited localities)
 * @return AGAS reference to created or found Arc starting at the saddle (to store in original arc.saddle)
 */
/*HPX-Remote-Task*/
hpx::id_type/*Arc*/ finishSweep(Arc arc);


/**
 * @brief Finds or Creates Arc for the minimal vertex in Boundary(arc)
 * writes arc and all of arc.descendants to the found/created saddle descendants
 * tests if saddle is ready for his own sweep and if so starts it (with async)
 * @param arc Reference to colocated Arc or copy of remote Arc that shall be linked to a saddle on this locality
 * @return pointer to created or found local Arc starting at the saddle of arc
 */
Arc* linkToSaddle(Arc const& arc);

/**
 * @brief
 * If arc has remote descendants call continueSweep(arc) on them
 * If arc has local descendants merge their boundaries to new Boundary(arc)
 * Move all elements from Boundary(arc) to new SweepFront(arc)
 * Add id to SweepFront(arc)
 * Start unordered BFS with SweepFront(arc) as queue through local mesh
 * If a boundary candidate is hit add it to Boundary(arc)
 * If a ghost cell is hit, call continueSweep on locality holding that vertex
 *
 * If SweepFront runs empty and boundary is empty and no continueSweep returns less than infinity
 *  send arc AGAS as trunk to contour tree compositing locality (it will issue augmentation fixing. it will start compositing the contour tree when both trunks arrived)
 *
 * Store minimal boundary candidate (from continueSweep-returns or boundary pushes)
 * call finishSweep on locality which returned minimal candidate or
 * call linkToSaddle on locally found minimal candidate
 * Store returned SaddleArc AGAS address in arc.saddle
 * Start Augmentation fix for each vertex visited (and on visited neighbor localities) for arc (?)
 * Send arc AGAS and SaddleArc AGAS to contour tree compositing locality
 * @param arc
 * Arcs do not change, once startSweep has been called on them, except for writing saddle index to original
 * Copies of remote Arcs are therefore valid
 * @param mesh
 * @param id
 * @return
 */
/*HPX-Local-Task*/
void startSweep(Arc* arc, Mesh* mesh, ID id);
/*
 *
*/


void mainOnEach(){
    //Load GlobalMesh (possibly no data loading here)
    //Mesh = partition() GlobalMesh
    //for each vertex v in Mesh (no ghost cells)
    //  if (isExtremum(v))
    //      arc = createOrFindArc(v) //must here always be create not find
    //      async startSweep() arc, v;
    //wait for signal of contour tree compositing locality that contour tree composition is done
}

#endif // DESIGN_H
