# TaskMergeTrees

Replicability repo for the paper

Unordered Task-Parallel Augmented Merge Tree Construction

Authors: Kilian Werner, Christoph Garth

To reproduce please run "runReplicabilityStamp.sh" on an Ubuntu 18.04.3 machine or newer.
Note that some dependency installations require super user privileges.

output consists of one VTK image data set file per locality containing augmentation ids for the respective domain as a scalar field 
A VTK polygon data set file contains the overall join tree.

The data ctBones.vti used in this repo corresponds to the foot data set used in the paper. 
Running the script runs the program on this data and thus reproduces the experiment of the first line of table 1 and (depending on the utilized hardware) the blue line (with squares) of the upper diagram in figure 5, as well as the orange bar (third from the right) in figure 7. Thus the vtp file should include around 0.54 million edges (tree arcs compare table 1) and the vti file should include around the same number of unique values (augmentation ids).


