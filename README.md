# Iterative Closest Point (ICP)
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

_____________________________

An example of how to use the code has been given in *basicICP.py*

```python
fileOriginal = '/icp/data/original.xyz'
deformed = '/icp/data/deformed.xyz'

source_points = read_file_original(fileOriginal)
dest_points_et_normal = read_file_deformed(deformed)

initial = np.array([[0.01], [0.05], [0.01], [0.001], [0.001], [0.001]])

# **********************************************************************
# Uncomment one of the following lines to trigger the respective module:
# **********************************************************************


#icp_point_to_plane(source_points,dest_points_et_normal,0)

#icp_point_to_point_lm(source_points,dest_points_et_normal,initial,0)

icp_point_to_plane_lm(source_points,dest_points_et_normal,initial,0)
```

This implementation does not use a CAD model for point-to-plane registration. You will need the two point clouds in .xyz format.

