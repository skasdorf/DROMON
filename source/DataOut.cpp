//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/19/21.
//

#include "DataOut.h"


unsigned int vtk_point_index_from_ijk(unsigned int i, unsigned int j, unsigned int order)
{
const bool ibdy = (i == 0 || i == order);
const bool jbdy = (j == 0 || j == order);
// How many boundaries do we lie on at once?
const int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0);

if (nbdy == 2) // Vertex DOF
{ // ijk is a corner node. Return the proper index (somewhere in [0,3]):
return (i ? (j ? 2 : 1) : (j ? 3 : 0));
}

int offset = 4;
if (nbdy == 1) // Edge DOF
{
if (!ibdy)
{ // On i axis
return (i - 1) + (j ? order - 1 + order - 1 : 0) + offset;
}

if (!jbdy)
{ // On j axis
return (j - 1) +
(i ? order - 1 : 2 * (order- 1) + order - 1) +
offset;
}
}

offset += 2 * (order - 1 + order- 1);
// nbdy == 0: Face DOF
return offset + (i - 1) + (order - 1) * ((j - 1));
}