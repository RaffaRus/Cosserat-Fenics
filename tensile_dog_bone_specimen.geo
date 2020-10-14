Merge "tensile-test-specimen.step"; // read the step file
Mesh.CharacteristicLengthMax = 2.5; // set the max element size lc
Mesh.ElementOrder = 1;              // ask for second-order elements

// define physical groups for BCs and materials
// the name in the LHS has to appear in the Fino input
// the number in the RHS is the numerical id of the entity in the CAD file
//Physical Surface (1) =  {1};   // left face, to be fixed
//Physical Surface (2) =  {7};  // right face, will hold the load
//+
Physical Surface(37) = {7};
Physical Surface(38) = {1};
Physical Volume ("bulk") =  {1};    // bulk material elements
