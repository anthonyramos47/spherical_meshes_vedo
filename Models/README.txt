Crumpled Developable - Keenan Crane (kmcrane@cs.cmu.edu)

Though this mesh looks like it has a lot of bumps and wrinkles, every vertex is
"intrinsically" flat: the angles around each around each interior vertex sum to
exactly 2π, which means the mesh can be flattened perfectly into the plane
without any distortion of lengths, areas, or interior angles.  (Such a
flattening is stored in the texture coordinates.)  The mesh may be useful for
testing the "isometry invariance" of geometry processing algorithms, i.e., the
degree to which their results are unchanged by motions of the vertex positions
that do not change edge lengths, areas, etc.

Naively, one sometimes calls such a mesh "discretely developable," though this
naive interpretation of developability neglects other important features of
smooth developable surfaces---such as the presence of straight lines passing
through each point of the surface.  This phenomenon is studied in the paper

   "Developability of Triangle Meshes"
   Stein, Grinspun, Crane
   (SIGGRAPH 2018)

The mesh is stored as a Wavefront OBJ file with vertex and texture coordinates.

License
--------------------------------------------

As the sole author of this data, I hereby release it into the public domain.

