# CFD

In this repo: inco.cpp & inco.h

inco, or incompressible, Is my first naive implementation of discretization techniques showcased in this video: https://youtu.be/dxt3zbVhL-k?si=XDydFeucaWwLc6IK

The basic flow of the sovler is as follows:

  1. update cell veloceites due to external forces(gravity)

  2. Manipulate the properties of the continuum to allign with our governing equations. ( Div($$\vec{u}$$) = 0 )

  3. Propagate the properties throught the continuum. ( Advection )

Divergence:

  For this first atempt we are only taking into account the conservation of mass term, which for an incompressible fluid simply states that the total flux of the velocity field over some volume must be zero. This essentaily mean that for anny amount of fluid flowing into a volume that same amount of fluid must also be flowing out. To achive this We will just calculate the divergence of a cell and then split that divergence up and subtract it from each of the cell edges, if we iterate over all of the cells enough times this should manage to get our divergences to about zero. For ease of calculating the edge fluxes we will "store" the cell velocities at the edges, rather than at the cell centroids(as most other solvers do). Notice that this "solver" doesn't even atempt to implement any fluid behavior, and actually doesn't even touch on the real meat of the navier stokes equation. This solver only atempts to make a velocity field that obeys conservation of mass, which is one of the closure equations needed to solve the other terms of the navier stokes equations.

Advection:

  To approximate the time rate of change of our fluid we will utilize a rather simple scheme called semi-lagrangian advection. This scheme treats the tranport of continuum properties as the physical movement of a parcel to which those propeties belong. This is where the name "semi-Lagrangian" comes from, since lagrangiean simulations discritize the fluid domain as a bunch of parcels of fluid that move around ( carrying their properties with them ), whereas Eulerian simulations discritize the domain as a mesh of static* cells whose properites move around. As it turns out it tends to be much easier to describe the motion of coninuum properties from the lagrangian perspective. Semi-lagrangian advection takes advantage of this by essentially pretending that there is a fluid parcel that will move to the exact position of our eulerian cell carrying it's properties with it. To implement this idea in code we need an approximate velocity for each storage point in each cell( remember cell velocities are stored on the edges ), then we step backwards from the storage point using our approximated velocites to some arbitrary point in the mesh. This point will be where we draw our updated properites from but the chances of this point actually lying on one of our "storage points" is essentailly zero so we need a way to approximate our continuum properties at any arbitrary point within the mesh. This is where cell interpolation comes in.

Interpolation: 



Usage:

  To use this "Sovler" you will need to set the solver and mesh settings, and then compile. currently there is no settings file, that is the next big update.

  Make commands:

  1. make bin/INCO: compiles INCO.exe

  2. make all: compiles INCO.exe

  3. make run: complies and then executes INCO.exe


  Advection Notes: [text](https://twister.caps.ou.edu/CFD2007/Chapter6.pdf)