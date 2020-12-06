# Q_SchrodingerPoisson1D_Dirichlet_demo
Schrodinger-Poisson solver with Dirichlet boundary condition in 1D demonstrator

This program computes the Schrodinger-Poisson equations in heterostructures with the Dirichlet boundary conditions. This means that the doping level on the left and right side of the structure are the boundary conditions since they are setting the Fermi level.
In order to be educatif and simple, few approximations are made:

-> The quantum structure is sandwitched between contacts with spacers. The Schrodinger solver is working ONLY in this domain.

-> The mass=meff is kept constant all over the structure. It means that meff should be set at the value of the well. 

-> The non-parabolicity of the bands is also NOT considered in the Schrodinger solver and the density of states.

-> Schottky contact can be simulated by setting the doping of the contact at zero and by setting the bandgap energy of the contact material.. 

As an input, the program only needs a set of layer with thickness, conduction band offset and doping.

Enjoy! If you like it, don t forget the star!
