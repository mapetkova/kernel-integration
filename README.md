# kernel-integration
This repository contains a C++ implementation of the code described in Petkova et al. 2018 (https://arxiv.org/abs/1710.07108).
The algorithm maps the density of SPH particles onto a Vornoi grid using the analytic integral of a cubic spline kernel function.
# Dependencies
The Voronoi grid is constructed using the free library Voro++ which can be downloaded from: http://math.lbl.gov/voro++/about.html
# How to use this code?
The best place to start is src/cell_density.cc. This simple program introduces the vertex integral function which is the final product of the derivation presented in Petkova et al.

The next program to review is src/single_particle_voronoi_density.cc which applies the vertex integral function to a Voronoi grid. This setup uses one SPH particle which is mapped onto 50 cells. The vertex integral is calculated for each vertex of each edge of each face of a cell and these calculations are added up to produce the mass of the cell. The density of each cell (cell mass divided by cell volume) can be viewed in 3D using plot/density_structure_3d.py, by plotting the Voronoi polyhedra with transparency set as the cell density.
