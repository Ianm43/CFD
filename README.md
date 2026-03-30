# CFD

In this repo: inco.cpp & inco.h

inco, or incompressible, Is my first naive implementation of discretization techniques showcased in this video: https://youtu.be/dxt3zbVhL-k?si=XDydFeucaWwLc6IK

The basic flow of the sovler is as follows:

  1. Manipulate the properties of the continuum to allign with our governing equations. ( Div($$\vec{u}$$) = 0 )

  2. Propagate the properties throught the continuum along the velocity field. ( Advection )

Divergence:
  For this first atempt we are only taking into account the conservation of mass term, which for an incompressible fluid simply states that the total flux of the velocity field over some volume must be zero. This essentaily mean that for anny amount of fluid flowing into a volume that same amount of fluid must also be flowing out. To achive this We will just calculate the divergence of a cell and then split that divergence up and subtract it from each of the cell edges, if we iterate over all of the cells enough times this should manage to get our divergences to about zero.

Advection:
  to aproximate the time rate of change of our fluid we will utilize a rather simple scheme called semi-lagrangian advection. This scheme treats the tranport of continuum properties as the physical movement of a parcel to which these propeties belong. this is where the name comes from since lagrangiean simulations discritize the fluid domain as a bunch of parcels of fluid that move around ( carrying their properties with them ), whereas Eulerian simulations discritize the domain as a mesh of static* cells whose properites move around.  

