# icp
A tutorial on iterative closest point using Python

The following has been implemented here:

* Basic point to plane matching has been done using a Least squares approach and a Gauss-Newton approach
* Point to point matching has been done using Gauss-Newton only

____________________________

All the important code snippets are in *basicICP.py*. The main functions are:

* icp_point_to_plane 
* icp_point_to_point_lm
* icp_point_to_plane_lm

*deformation.py* has been used to deform the point cloud, so that we may validate the ICP based registration. We had only one set of point cloud and their correspinding normal vectors as the input. That was deformed using deformation.py. And then it is being registered with basicICP.py. That gives us an easy way to validate the ICP results.


*transformations.py* has been taken from [the source code by Christoph Gohlke](http://www.lfd.uci.edu/~gohlke/code/transformations.py.html)

