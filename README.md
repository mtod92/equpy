# equipy
---
title: 'equpy: Python and MATLAB implementations of a fast numerical algorithm for the solution of multiple binding and chemical equilibria'
tags:
  - Python
  - MATLAB
  - chemical equilibrium
  - binding
  - physical chemistry
authors:
  - name: Luca Citelli
      affiliation: 1
  - name: Marco Todisco
    corresponding: true
    orcid: 0000-0001-6627-5283
    affiliation: 2
affiliations:
 - name: Independent Researcher.
   index: 1
 - name: Howard Hughes Medical Institute, Department of Chemistry, The University of Chicago, Chicago, Illinois 60637, USA.
   index: 2
date: 11 December 2023
bibliography: paper.bib

# Summary

Determining the distribution of multiple chemical species at equilibrium for a given system is a common problem that must be routinely addressed by scholars. While simple systems consisting of a few species and reactions can be solved manually, most of these problems require the definition and solution of higher-order equations and are intractable without reliable numerical methods, that can be slow and inefficient. In this work, we present straightforward Python and MATLAB implementations of the geometric-programming algorithm developed by Thomas Wayne Wall (1984) and we provide clear and easy-to-use scripts and examples for researchers approaching the problem. The performance and stability of the algorithm is tested versus out-of-the-box MATLAB numerical solver (vpasolve) and the solver available in chempy, showing a speed-up as high as two orders of magnitudes.

# Statement of need

`equpy` is a Python package specialized in solving multiple chemical equilibria.
`equpy` was designed to provide a user-friendly experience for a class-based
implementations of the algorithm developed by Thomas Wayne Wall to handle the
solution of mixed linear/non-linear systems of equations describing the simultaneous
equilibration of multiple species reacting in a close system. `equpy` relies 
heavily on numpy to handle data in matrix form and matplotlib to generate figures.

`equpy` was designed to be used by chemical, biological and physical/chemical
researchers studying systems comprising multiple species reacting/interacting.
This tool is particularly well-suited for researchers that (i) need to quickly 
calculate equilibria with a matrix-form input, and/or (ii) cannot rely on 
slower and more general solvers, for example when multiple kinetic traces have 
to be integrated as the species equilibrate and/or global fits on large sets 
of parameters have to be performed.

# Introduction
The mathematical treatment of multiple chemical equilibria is a problem that can be intimidating for inexperienced researchers and community members that are not familiar with computer science and linear algebra. While tutorials for the solution of simple systems are readily available on basic chemistry textbooks and online resources, the task of upscaling these is non-trivial.
While open source or commercial packages and standalone software are available to address this task, such as Cantera, EQS4WIN, TOMSYM, COMSOL and chempy, these often rely on black-box solvers and may require some background in coding and potentially inconvenient input requirements, needing to write the set of reactions and mass conservations in extended form that may not be well-suited for automated processes.
In this work, we focus on the implementation of an approach developed by Thomas Wayne Wall1,2 for linearization and solution of the system of equations describing a complex set of chemical reactions. The compact scripts provided here accept as inputs the reaction stoichiometry matrix, the associated array of equilibrium constants, mass conservation matrix and associated array of amounts of starting materials. All of these are conveniently input as matrices and arrays and in the examples provided here these and can be directly read through .csv files editable using common software such as Microsoft Excel, Apple Numbers, Apache OpenOffice Calc and LibreOffice Calc.

In the first section of this work, we define the general problem, in the second section we present the general mathematical treatment and implementation of the algorithm, and finally in the last section we present practical examples, focused to readers that primarily need to apply this implementation to solve their own problem, and we discuss usage and performance of the algorithm, with attention to real case scenarios and applications for the unexperienced reader.

# Problem Definition
Let’s start with a simple example of interacting species:
$$A + 2B &harr; AB_2$$
$$AB_2 + C &harr; AB_2C$$

And their associated equilibrium constants, defined as the ratio between forward and backward reaction rates:
$$\frac{[AB_2]}{([A][B]^2)} = K_1 \label{Eq. 1}$$
$$\frac{[AB_2C]}{([AB_2][C])} = K_2 \label{Eq. 2}$$

And the associated mass conservations:
$$[A]_{tot} = [A] \+ [AB_2] \+ [AB_2C] \label{Eq. 3}$$
$$[B]_{tot} = [B] \+ 2\*[AB_2] \+ 2\*[AB_2C] \label{Eq. 4}$$
$$[C]_{tot} = [C] \+ [AB_2C] \label{Eq. 5}$$

We can define a system comprising these equations to be simultaneously solved.
In this system, we can see that chemical equilibria consist of nonlinear functions, meaning that they cannot be expressed as a sum of their variables each raised to the power of one. A common approach to solve this problem is to linearize these equations first to make them suitable to be solved employing numerical methods.

Common solutions in literature revolve around the so-called kinetic approach and thermodynamic approach, which involves finding the extent of reaction at the energy minimum of the system.

The algorithm presented here the advantage of operating on a user-friendly set of equations that is intuitively employed by any user with a basic chemistry knowledge. These equations are solved with an approach equivalent to the Newton search of the logarithmic equations over the logarithm of the variables.

# Mathematical Treatment
In a system with n different species $X_{1…n}$, the mass conservation relationship for the $i^{-th}$ species can be stated as the sum over all the species contributions with their relative stoichiometries (a). We can define the conservation of mass for species $X_i$ as:

$$a_1[X_1] + a_2[X_2] + ... + a_n[X_n] = [X_i]_{tot} \label{Eq. 6}$$

Or equivalently:
$$\sum_{j=1}^n a_j[X_j] = [X_i]_{tot} \label{Eq. 7}$$

In order to express such conservation of mass as a linear function of the logarithm of concentrations of the reactants, following the approach by Wall we must first transform the summations to products using the theory of the arithmetic-geometric mean inequality from Passy 3 as applied by Baker 4. We reorganize Eq. 7 so that the summation over all strictly positive terms *a* and *X* is rewritten as the following:

$$\frac{[X_i]_{tot}}{\sum_{j=1}^n a_j[X_j]} = 1 \label{Eq. 8}$$

Then we “condense” the sum in the denominator of Eq. 8 into a product:
$$\sum_{j=1}^n a_j[X_j] = \prod_{j=1}^n \bigg(\frac{a_j[X_j]}{W_j}\bigg)^{W_j} \label{Eq. 9}$$

With *W* for a given species *j* part of a mass conservation relationship being equal to:
$$W_j = \frac{a_j[X_j]}{\sum_{p=1}^n a_p[X_p]} \label{Eq. 10}$$

So that Eq. 8 becomes:
$$\frac{[X_i]_{tot}}{\bigg(\frac{a_1[X_1]}{W_1}\bigg)^{W_1} * \bigg(\frac{a_2[X_2]}{W_2}\bigg)^{W_2} * ... * \bigg(\frac{a_n[X_n]}{W_n}\bigg)^{W_n}} = 1 \label{Eq. 11}$$

We can then reorganize the fraction:
$$[X_1]^{-1}\*[X_2]^{-2}\*...\*[X_n]^{-n}$$

$$*[X_i]_{tot}*\bigg{\bigg(\frac{W_1}{a_1}\bigg)^{W_1}*\bigg(\frac{W_2}{a_2}\bigg)^{W_2}*...*\bigg(\frac{W_n}{a_n}\bigg)^{W_n}\bigg} \label{Eq. 12}$$

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

$$\frac{[AB_2]}{([A][B]^2)} = K_1$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

[AB_2\ ]/([A] [B]^2 )=K_1

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References