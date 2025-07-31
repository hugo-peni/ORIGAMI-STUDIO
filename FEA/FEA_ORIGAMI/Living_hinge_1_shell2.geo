Merge "origami_shell.step";

// Mesh settings
Mesh.ElementOrder = 2;
Mesh.Optimize = 1;
Mesh.HighOrderOptimize = 2;
Mesh.CharacteristicLengthMax = 4;

// Force all surfaces to be recombined into quads
Mesh.RecombineAll = 1;
Mesh.SecondOrderIncomplete = 1; // ⬅️ Reduces 9-node quads to 8-node (needed for S8R)

// This is safer than Recombine Surface {*}, especially with STEP
Physical Surface("all_surfs") = Surface{:}; // required to export to .inp

// Generate mesh
Mesh 2;

// Export
Mesh.SaveGroupsOfNodes = 1;
Save "gmsh.inp";

// Viewer
General.Trackball = 0;
General.RotationX = -50;