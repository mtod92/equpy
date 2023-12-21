import numpy as np
import matplotlib.pyplot as plt
import warnings
from typing import List, Tuple
from utils import eq_system_builder

class EquationSystem:
    """
    Assembles an object including the equations defining the systems reactions and mass conservations
    Parameters:
    - equations: set of reactions, either list of Str or list of numerical lists
    - mass_conservation: set of mass conservations, either list of Str or list of numerical lists
    - species: optional argument. If reactions in Str form are passed the species names are automatically generated.
      Custom species names can be passed when input is in matrix form as a dictionary of species, such as: 
      {'A':0, 'B':1, 'C':2} for a system involving three species (n=3) named A, B, C.
    """
    def __init__(self, equations, mass_conservation, species = None):
        self.stoichiometry, self.mass_conservation, self.species = equations, mass_conservation, species
        if isinstance(self.stoichiometry[0], str): #assuming the user is passing reactions in textual form
            self.stoichiometry, self.mass_conservation, self.species = eq_system_builder(self.stoichiometry, self.mass_conservation)

        if self.species == None:
            self.species = {i:i for i in range(np.shape(self.stoichiometry)[1])} #placeholder for species name if not specified

class ChemicalReaction:
    """
    Assembles an object with the equations defining the systems reactions and mass conservations 
    plus equilibrium constants and total masses.
    Parameters:
    - equation_system: EquationSystem object as defined in the respective Class
    - eq_constants: array of equilibrium constants, definde as the ratio between forward and reverse reactions.
    - initial_masses: array of initial masses.

    IMPORTANT: for equpy's result to be meaningful the unit of measurement have to be matching: if the equilibrium constants
    are given in (1/micromolar) units, the concentrations have to be given in micromolar units.
    """
    def __init__(
        self,
        equation_system: EquationSystem,
        eq_constants: np.ndarray,
        initial_masses: np.ndarray,
    ):
        initial_masses = np.array(initial_masses, dtype=float)
        if min(initial_masses) == 0:
            warnings.warn(
                "Species concentrations (S) should not be set to zero, eps has been set instead. The result may not be reliable.",
                stacklevel=2,
            )
            I = np.where(initial_masses == 0)[0]
            initial_masses[I] = 2.2e-16

        self.s = equation_system.species
        self.N = equation_system.stoichiometry
        self.K = eq_constants
        self.C = equation_system.mass_conservation
        self.S = initial_masses
        self.result = []
        self.residuals = []

    """
    solve method calls eqsolver to find equilibrium concentration of all species.
    Parameters:
    - iter: maximum number of iterations allowed. The algorithm typically converges in less then 10 steps, 
    but poor initial guesses for solution may require more steps.
    - tolerance: early stopping criteria. The algorithm is limited by numerical precision, so that the solution
    can (at best) be correct up to ~16 significant digits. A tolerance of 100 means that we are fine with two 
    orders of magnitude smaller precision, for example relying on ~14 significant digits.
    - w: weighted update for solutions. This can be kept as zero for well-behaving systems. Small values
    such as 0.5 or 1.0 grealy improve code stability for more difficult tasks, but will cause a slower
    approach to convergence.
    """
    def solve(
        self, iter: int, tolerance: float, w: float
    ) -> Tuple[np.ndarray, List[np.ndarray]]:
        iter = int(iter)
        x0 = np.ones(np.shape(self.N)[1])
        delta = []

        x, delta_ = eqsolver(self.N, self.K, self.C, self.S, x0, w)
        delta.append(delta_)

        for i in range(iter):
            x, delta_ = eqsolver(self.N, self.K, self.C, self.S, x, w)
            delta.append(delta_)

            if delta[-1] < tolerance * np.linalg.norm(x * 2.2e-16):
                self.result = np.exp(x)
                self.residuals = delta
                return np.exp(x), delta

        self.result = np.exp(x)
        self.residuals = delta
        warnings.warn(
            "Tolerance not reached. Manually check if the result is satisfying. Either change starting point, increase iterations and/or weight.",
            stacklevel=2,
        )
        return np.exp(x), delta

    """
    plotter method relies on matplotlib to generate a graph showing algorithm's progress to convergence and final results.
    """
    def plotter(self):
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=False, figsize=(15, 5))
        ax1.plot(np.arange(len(self.residuals)), self.residuals, lw=3, color="black")
        ax1.scatter(
            np.arange(len(self.residuals)), self.residuals, s=120, color="black"
        )

        ax1.set_xlabel("steps", fontsize=14, fontname = 'Arial')
        ax1.set_ylabel("||r||", fontsize=14, fontname = 'Arial')
        ax1.set_title('Algorithm Progress', fontsize=16, fontname = 'Arial')
        ax1.set_yscale("linear")

        ax2.bar(np.arange(0, len(self.result)), self.result)
        ax2.set_ylabel("concentration", fontsize=14, fontname = 'Arial')
        ax2.set_xticks(np.arange(0, len(self.result)))
        ax2.set_title('Equilibrium Concentrations', fontsize=16, fontname = 'Arial')
        ax2.set_xticklabels(
            sorted(self.s, key=lambda x: self.s[x]), fontsize=12, rotation=90, fontname = 'Arial'
        )


def eqsolver(
    N: np.ndarray, K: np.ndarray, C: np.ndarray, S: List[float], x: float, w: float
) -> Tuple[float, float]:
    """
    eqsolver is the core of equpy and handles the application of Thomas Wayne Wall algorithm on the input
    """
    Cx = C * np.exp(x)
    W = Cx / np.sum(Cx, axis=1)[:, None]
    Ks = np.prod((W / C) ** W, axis=1) * S

    M = np.vstack((N, W))
    y = np.log(np.hstack((K, Ks)))
    delta = np.linalg.norm(np.dot(M, x.T) - y[:, None].T)

    if np.shape(N)[1] == len(K) + len(S):
        x = (w * x + np.linalg.solve(M, y)) / (w + 1)
    else:  # solve overdetermined systems as M^(-1) * y
        x = (w * x + np.dot(np.linalg.pinv(M), y)) / (w + 1)
    return x, delta
