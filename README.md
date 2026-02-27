# CFD

In this repo: inco.cpp & inco.h

inco, or incompressible, Is my first naive implementation of discretization techniques showcased in this video: https://youtu.be/dxt3zbVhL-k?si=XDydFeucaWwLc6IK

The basic flow of the sovler is as follows:

  1. Manipulate the properties of the continuum to allign with our governing equations. ( Div($$\vec{u}$$) = 0 )

  2. Propagate the properties throught the continuum along the velocity field. ( Advection )

