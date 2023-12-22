# equipy

`equpy` is a Python package specialized in solving multiple chemical equilibria.
`equpy` was designed to provide a user-friendly experience for a class-based
implementations of the algorithm developed by Thomas Wayne Wall to handle the
solution of mixed linear/non-linear systems of equations describing the simultaneous
equilibration of multiple species reacting in a close system. `equpy` relies 
on numpy to handle data in matrix form and matplotlib to generate figures.

`equpy` was designed to be used by chemical, biological and physical/chemical
researchers studying systems comprising multiple species reacting/interacting.
This tool is particularly well-suited for researchers that (i) need to quickly 
calculate equilibria with a matrix-form input, and/or (ii) cannot rely on 
slower and more general solvers, for example when multiple kinetic traces have 
to be integrated as the species equilibrate and/or global fits on large sets 
of parameters have to be performed.

# Summary
Determining the distribution of multiple chemical species at equilibrium for a given system is a common problem that must be routinely addressed by scholars. While simple systems consisting of a few species and reactions can be solved manually, most of these problems require the definition and solution of higher-order equations and are intractable without reliable numerical methods, that can be slow and inefficient. In this work, we present straightforward Python and MATLAB implementations of the geometric-programming algorithm developed by Thomas Wayne Wall (1984) and we provide clear and easy-to-use scripts and examples for researchers approaching the problem. The performance and stability of the algorithm is tested versus out-of-the-box MATLAB numerical solver (vpasolve) and the solver available in chempy, showing a speed-up as high as two orders of magnitudes.

# Introduction
The mathematical treatment of multiple chemical equilibria is a problem that can be intimidating for inexperienced researchers and community members that are not familiar with computer science and linear algebra. While tutorials for the solution of simple systems are readily available on basic chemistry textbooks and online resources, the task of upscaling these is non-trivial.
While open source or commercial packages and standalone software are available to address this task, such as Cantera, EQS4WIN, TOMSYM, COMSOL and chempy, these often rely on black-box solvers and may require some background in coding and potentially inconvenient input requirements, needing to write the set of reactions and mass conservations in extended form that may not be well-suited for automated processes.
In this work, we focus on the implementation of an approach developed by Thomas Wayne Wall1,2 for linearization and solution of the system of equations describing a complex set of chemical reactions. The compact scripts provided here accept as inputs the reaction stoichiometry matrix, the associated array of equilibrium constants, mass conservation matrix and associated array of amounts of starting materials. All of these are conveniently input as matrices and arrays and in the examples provided here these and can be directly read through .csv files editable using common software such as Microsoft Excel, Apple Numbers, Apache OpenOffice Calc and LibreOffice Calc.

# Problem Definition
Let’s start with a simple example of interacting species:
$$A + 2B &harr; AB_2$$
$$AB_2 + C &harr; AB_2C$$

And their associated equilibrium constants, defined as the ratio between forward and backward reaction rates:
$$\frac{[AB_2]}{([A][B]^2)} = K_1$$
$$\frac{[AB_2C]}{([AB_2][C])} = K_2$$

And the associated mass conservations:

$$[A]_{tot} = [A] \+ [AB_2] \+ [AB_2C]$$

$$[B]_{tot} = [B] \+ 2\*[AB_2] \+ 2\*[AB_2C]$$

$$[C]_{tot} = [C] \+ [AB_2C]$$

We can define a system comprising these equations to be simultaneously solved.
In this system, we can see that chemical equilibria consist of nonlinear functions, meaning that they cannot be expressed as a sum of their variables each raised to the power of one.
`equpy` solves this problem by linearizing these equations to make them suitable to be solved employing linear algebra and an iterative numerical method equivalent to the Newton search of the logarithmic equations over the logarithm of the variables.

# Examples
## Example 1
Using equpy we can solve the example by simply expressing the reactions and mass conservation relationship in literal form.
In the example below, a simple implementation is shown:

```
from equpy import ChemicalReaction, EquationSystem
import numpy as np
import matplotlib.pyplot as plt

# we can define chemical reactions
reactions = ['A + 2B = AB2',
    'AB2 + C = AB2C']

# we can establish mass conservation for A, B and C
mass_conservation = ['A + AB2 + AB2C',
                'B + 2AB2 + 2AB2C',
                'C + AB2C']

K = [1, 10] # define equilibrium constants
total_masses = [1,2,3] # define total masses of A, B and C to be conserved

# use equpy to build an EquationSystem object
# starting from equations expressed in literal form
eq_system = EquationSystem.from_literal_equations(reactions, mass_conservation)

# set up the ChemicalReaction object comprising 
# the system of equation, equilibrium constants and total_masses
reaction = ChemicalReaction(eq_system, K, total_masses)

# solve and plot the results
x, delta = reaction.solve(20, 1e2, 0)
reaction.plotter()
```

## Example 2
Using equpy we can also solve the example by providing reactions and mass conservation relationship in matrix form.
In the example below, a simple implementation is shown:
```
from equpy import ChemicalReaction, EquationSystem
import numpy as np
import matplotlib.pyplot as plt

stoichiometry = np.array([[-1, -2, 0, 1, 0],[0, 0, -1, -1, 1]])
K = np.array([1, 10])
total_masses = np.array([1, 2, 3])
mass_conservation = np.array([[1, 0, 0, 1, 1],[0, 1, 0, 2, 2],[0, 0, 1, 0, 1]])
species_names = {'A':0, 'B':1, 'C':2, 'AB2':3, 'AB2C':4}

# use equpy to build an EquationSystem object
# starting from equations expressed in matrix form.
# since we have not provided reactions in literal form,
# a dictionary containing species_names is provided,
# with species listed in the same order as the 
# stoichiometry matrix
eq_system = EquationSystem(stoichiometry, mass_conservation, species_names)
reaction = ChemicalReaction(eq_system, K, total_masses)

x, delta = reaction.solve(20, 1e2, 0)
reaction.plotter()
```

## Example 3
Finally, using equpy we can also solve the example by reading reactions, mass conservation, equibrium constants and total masses in matrix form from .csv files using a simple utility function.
In the example below, a simple implementation is shown:
```
from equpy import ChemicalReaction, EquationSystem
import numpy as np
import matplotlib.pyplot as plt
from utils import csv_loader

filename_stoichiometry = '/examples/readme_example/st.csv'
filename_mass_conservations = '/examples/readme_example/mc.csv'

stoichiometry, K, mass_conservation, total_masses = csv_loader(filename_stoichiometry, filename_mass_conservations)

eq_system = EquationSystem(stoichiometry, mass_conservation, {'A':0, 'B':1, 'C':2, 'AB2':3, 'AB2C':4})
reaction = ChemicalReaction(eq_system, K, total_masses)

x, delta = reaction.solve(20, 1e2, 0)
reaction.plotter()
```