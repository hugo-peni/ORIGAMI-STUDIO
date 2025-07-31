// Load the STEP geometry
Merge "origami_thickness.step";

// Use quadratic elements
Mesh.ElementOrder = 2;
Mesh.Optimize = 1;
Mesh.HighOrderOptimize = 2;
Mesh.CharacteristicLengthMax = 2;

// Recombine triangles into quads
Recombine Surface {*} ;

// Set meshing options for surfaces only
Mesh.SurfaceFaces = 1;
Mesh.SurfaceEdges = 1;
Mesh.VolumeFaces = 0;
Mesh.VolumeEdges = 0;

// Generate the surface mesh (2D)
Mesh 2;

// Save mesh groups and export to CalculiX format
Mesh.SaveGroupsOfNodes = 1;
Save "gmsh.inp";

// Viewer settings (optional)
General.Trackball = 0;
General.RotationX = -50;

