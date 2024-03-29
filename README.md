# equpy

`equpy` is a package specialized in solving multiple chemical equilibria.
`equpy` was designed to provide a user-friendly experience for a modern implementation of the algorithm developed by [Thomas Wayne Wall](https://repository.mines.edu/bitstream/handle/11124/13991/Wall_10782543.pdf?sequence=1to) to handle the solution of mixed linear/non-linear systems of equations describing the simultaneous equilibration of multiple species reacting in a close system. `equpy` relies on numpy to handle data in matrix form and matplotlib to generate figures.

`equpy` was designed to be used by chemical, biological and physical/chemical researchers studying systems comprising multiple species reacting/interacting. This tool is particularly well-suited for researchers that (i) need to quickly calculate equilibria with a matrix-form input, and/or (ii) cannot rely on slower and more general solvers, for example when multiple kinetic traces have to be integrated as the species equilibrate and/or global fits on large sets of parameters have to be performed.

# Summary
Determining the distribution of multiple chemical species at equilibrium for a given system is a common problem that must be routinely addressed by scholars. While simple systems consisting of a few species and reactions can be solved manually, most of these problems require the definition and solution of higher-order equations and are intractable without reliable numerical methods, that can be slow and inefficient. In this work, we present straightforward Python and MATLAB implementations of the geometric-programming algorithm developed by [Thomas Wayne Wall](https://repository.mines.edu/bitstream/handle/11124/13991/Wall_10782543.pdf?sequence=1to) (1984) and we provide clear and easy-to-use scripts and examples for researchers approaching the problem. The performance and stability of the algorithm is tested versus out-of-the-box MATLAB numerical solver (*vpasolve*) and the solver available in chempy - one of the most complete open source chemistry packages available to this date - showing an execution time reduced by as much as two orders of magnitudes.

# Installation
To incorporate `equpy` into your project, follow these simple steps:

1) Click on the green "Code" button and choose "Download ZIP" to download the package files.
2) Extract the downloaded ZIP file to a location of your choice.
3) Copy the extracted files into the local folder of your project.

Now, you're ready to import and seamlessly integrate `equpy` into your project and leverage its functionality. If you encounter any issues or have questions, refer to the documentation or feel free to reach out to us.

# Introduction
The mathematical treatment of multiple chemical equilibria is a problem that can be intimidating for inexperienced researchers and community members that are not familiar with computer science and linear algebra. While tutorials for the solution of simple systems are readily available on basic chemistry textbooks and online resources, the task of upscaling these is non-trivial.
While open source or commercial packages and standalone software are available to address this task, such as Cantera, EQS4WIN, TOMSYM, COMSOL and chempy, these often rely on black-box solvers and may require some background in coding and potentially inconvenient input requirements, needing to write the set of reactions and mass conservations in extended form that may not be well-suited for automated processes.
In this work, we focus on the implementation of an approach developed by [Thomas Wayne Wall](https://repository.mines.edu/bitstream/handle/11124/13991/Wall_10782543.pdf?sequence=1to) for linearization and solution of the system of equations describing a complex set of chemical reactions. The compact scripts provided here accept as inputs the reaction stoichiometry matrix, the associated array of equilibrium constants, mass conservation matrix and associated array of amounts of starting materials. All of these are conveniently input as matrices and arrays and in the examples provided here these and can be directly read through .csv files editable using common software such as Microsoft Excel, Apple Numbers, Apache OpenOffice Calc and LibreOffice Calc.

# Problem definition
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

# Solving the example using equpy
In the following section, we present different approaches on how to solve the previous example using `equpy`.
These are also available and ready to use in the [interactive Jupyter Notebook](equpy_test.ipynb).

In order to use equpy, download the files from the repository and place them in the same folder as your current project.

## Approach (1) - text input
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

# call solver: by default it is set to run with the following settings:
# iter = 1e2,
# x0 = np.array([1, 1, 1, 1, 1]),
# tolerance = 1e2, 
# w = 0
# all these can be passed as custom arguments by the user, for example:
# x, delta = reaction.solve(iter = 10, x0 = np.array([1,2,3,4,5]), tolerance = 1e5, w = 0.7)

x, delta = reaction.solve()

# plot the results
reaction.plotter()
```
<p align="center">
  <img src="Figures/example_result.png" width="1100" title="equpy plotter output">
</p>

## Approach (2) - matrix input
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

x, delta = reaction.solve()
reaction.plotter()
```

## Approach (3) - matrix input loaded from .csv
Finally, using equpy we can also solve the example by reading reactions, mass conservation, equibrium constants and total masses in matrix form from .csv files using a simple utility function.
In the example below, a simple implementation is shown:
```
from equpy import ChemicalReaction, EquationSystem
import numpy as np
import matplotlib.pyplot as plt
from utils import csv_loader

filename_stoichiometry = '/examples/readme_example/st.csv'
filename_mass_conservations = '/examples/readme_example/mc.csv'

stoichiometry, K, mass_conservation, total_masses = csv_loader(
    filename_stoichiometry, filename_mass_conservations)

eq_system = EquationSystem(
    stoichiometry, mass_conservation, {'A':0, 'B':1, 'C':2, 'AB2':3, 'AB2C':4})
reaction = ChemicalReaction(eq_system, K, total_masses)

x, delta = reaction.solve(iter = 20)
reaction.plotter()
```

# Performance benchmark
To test the performance of equpy, we focus on the most complete chemistry package available for Python users and the built-in MATLAB solver *vpasolve*. In order to test them on the same problem, we have ran the simple example provided in [chempy documentation](https://github.com/bjodah/chempy#chemical-equilibria) and timed uniquely the solver execution time with the "time" module, leaving out all overhead required for variables definitions and setting up the problem.

## chempy benchmark
```
from chempy import Equilibrium
from chempy.chemistry import Species
from chempy.equilibria import EqSystem
from collections import defaultdict
import time

water_autop = Equilibrium({'H2O'}, {'H+', 'OH-'}, 10**-14)  # unit "molar" assumed
ammonia_prot = Equilibrium({'NH4+'}, {'NH3', 'H+'}, 10**-9.24)  # same here
substances = [Species.from_formula(f) for f in 'H2O OH- H+ NH3 NH4+'.split()]
eqsys = EqSystem([water_autop, ammonia_prot], substances)
print('\n'.join(map(str, eqsys.rxns)))  # "rxns" short for "reactions"

init_conc = defaultdict(float, {'H2O': 1, 'NH3': 0.1})

start_time = time.time()
for i in range(1000):
    x_chempy, sol, sane = eqsys.root(init_conc)
print("")
print("execution time --- %s milliseconds ---" % ((time.time() - start_time)/i*1000))
print("")
```

Which provides a time of ~19 ms with the code executed in Visual Studio Code from a Jupyter Notebook.

## equpy benchmark
```
from equpy import ChemicalReaction, EquationSystem
import numpy as np
import time

reactions = ['OH + H = H2O',
    'NH3 + H = NH4']

mass_conservation = ['H2O + OH',
                     'H2O + H + NH4',
                     'NH3 + NH4']

K = [1e14, 10**(9.24)]
total_masses = [1, 1, 0.1]

eq_system = EquationSystem.from_literal_equations(reactions, mass_conservation)
reaction = ChemicalReaction(eq_system, K, total_masses)

start_time = time.time()
for j in range(1000):
    x, delta = reaction.solve(iter = 20, w = 0.5)
print("")
print("execution time --- %s milliseconds ---" % ((time.time() - start_time)/j*1000))
print("")
```

Which provides a time of ~1.9 ms with the code executed in Visual Studio Code from a Jupyter Notebook, for a 10x enhancement and approaching the same result in 20 iterations:

<p align="center">
  <img src="Figures/chempy_comparison.png" width="650" title="Comparison of chempy result and equpy result after 20 iterations">
</p>

## MATLAB *vpasolve* benchmark
We also provide the comparison with MATLAB built-in solver *[vpasolve](https://www.mathworks.com/help/symbolic/sym.vpasolve.html)*

```
x = solve_([1e14 10^(9.24)], 1, 0.1);

function values = solve_(K, H2Otot, NH3tot)
    guesses = zeros(5,2);
    guesses(:,2) = Inf;
    syms H2O OH H NH3 NH4 real
    
    eqns = [OH*H*K(1) == H2O;
        NH3*H*K(2) == NH4;
        OH + H2O == H2Otot;
        H + NH4 + H2O == H2Otot;
        NH3 + NH4 == NH3tot];
    
    tic
    for i = 1:100
    S = vpasolve(eqns, ...
        [H2O OH H NH3 NH4], guesses);
    end
    time = toc/i*1000
    
    fn = fieldnames(S);
    values = zeros(length(eqns),1);
    
    for k=1:numel(fn)
        values(k) = double(S.(fn{k}));
    end
```

Which provides a time of ~31 ms, with equpy giving more than a 15x enhancement.

# Final considerations on performance
While we focused here on benchmarking performance of the approach by Thomas Wayne Wall, it is worth noting how more important than the absolute time is the relative enhancement with respect to more general solvers. 

Regarding absolute execution time, a simple MATLAB implementation of equpy (see [MATLAB](MATLAB/example_code.m) folder) allows to solve the problem in 0.19 ms, for an extra 10x boost in absolute execution time. Running the python implementation here presented from terminal allows for a reduction of execution time both for equpy and chempy of ~30%. It would not be suprising for a simple implementation of the Python algorithm relying on numba or cython to execute in a time comparable or faster than MATLAB.