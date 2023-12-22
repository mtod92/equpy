import numpy as np
import matplotlib.pyplot as plt
import warnings
from typing import List, Tuple, Union, Optional, Dict
from utils import eq_system_builder


class EquationSystem:
    def __init__(
        self,
        equations: Union[List[List[float]], np.ndarray],
        mass_conservation: Union[List[List[float]], np.ndarray],
        species: Optional[Dict[str, int]] = None,
    ):
        """
        Assembles an object including the equations defining the systems reactions and mass conservations
        Parameters:
        - equations: set of reactions, either list of list of numbers or bidimensional numpy.ndarray
        - mass_conservation: set of mass conservations,  either list of list of numbers or bidimensional numpy.ndarray
        - species: optional argument. Custom species names can be passed when input is in matrix form as a dictionary of species, such as:
        {'A':0, 'B':1, 'C':2} for a system involving three species (n=3) named A, B, C.
        """
        self.stoichiometry = (
            equations if isinstance(equations, np.ndarray) else np.array(equations)
        )
        self.mass_conservation = (
            mass_conservation
            if isinstance(mass_conservation, np.ndarray)
            else np.array(mass_conservation)
        )
        self.species = (
            species
            if species is not None
            else {i: i for i in range(np.shape(self.stoichiometry)[1])}
        )  # placeholder for species name if not specified

    @classmethod
    def from_literal_equations(
        cls, equations: List[str], mass_conservation: List[str]
    ) -> None:
        """
        Alternative constructor that assembles an object including the equations defining the systems reactions and mass conservations,
        passed as text
        Parameters:
        - equations: set of reactions, expressed as a list of strings
        - mass_conservation: set of mass conservations, expressed as a list of strings
        """
        _stoichiometry, _mass_conservation, _species = eq_system_builder(
            equations, mass_conservation
        )
        return cls(_stoichiometry, _mass_conservation, _species)


class ChemicalReaction:
    """
    Assembles an object with the equations defining the systems reactions and mass conservations
    plus equilibrium constants and total masses.
    Parameters:
    - equation_system: EquationSystem object as defined in the respective Class
    - eq_constants: array of equilibrium constants, definde as the ratio between forward and reverse reactions.
    - total_masses: array of total masses.

    IMPORTANT: for equpy's result to be meaningful the unit of measurement have to be matching: if the equilibrium constants
    are given in (1/micromolar) units, the concentrations have to be given in micromolar units.
    """

    def __init__(
        self,
        equation_system: EquationSystem,
        eq_constants: np.ndarray,
        total_masses: np.ndarray,
    ):
        total_masses = np.array(total_masses, dtype=float)
        if min(total_masses) == 0:
            warnings.warn(
                "Species concentrations (S) should not be set to zero, eps has been set instead. The result may not be reliable.",
                stacklevel=2,
            )
            I = np.where(total_masses == 0)[0]
            total_masses[I] = 2.2e-16

        self.species = equation_system.species
        self.stoichiometry = equation_system.stoichiometry
        self.K = eq_constants
        self.mass_conservation = equation_system.mass_conservation
        self.total_masses = total_masses
        self.result = []
        self.residuals = []

    """
    solve method calls eqsolver to find equilibrium concentration of all species.
    Parameters:
    - iter: maximum number of iterations allowed. The algorithm typically converges in less then 10 steps, 
    but poor total guesses for solution may require more steps.
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
        x0 = np.ones(np.shape(self.stoichiometry)[1])
        delta = []

        x, delta_ = eqsolver(self.stoichiometry, self.K, self.mass_conservation, self.total_masses, x0, w)
        delta.append(delta_)

        for i in range(iter):
            x, delta_ = eqsolver(self.stoichiometry, self.K, self.mass_conservation, self.total_masses, x, w)
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

        ax1.set_xlabel("steps", fontsize=14, fontname="Arial")
        ax1.set_ylabel("||r||", fontsize=14, fontname="Arial")
        ax1.set_title("Algorithm Progress", fontsize=16, fontname="Arial")
        ax1.set_yscale("linear")

        ax2.bar(np.arange(0, len(self.result)), self.result)
        ax2.set_ylabel("concentration", fontsize=14, fontname="Arial")
        ax2.set_xticks(np.arange(0, len(self.result)))
        ax2.set_title("Equilibrium Concentrations", fontsize=16, fontname="Arial")
        ax2.set_xticklabels(
            sorted(self.species, key=lambda x: self.species[x]),
            fontsize=12,
            rotation=90,
            fontname="Arial",
        )


def eqsolver(
    stoichiometry: np.ndarray, K: np.ndarray, mass_conservation: np.ndarray, total_masses: List[float], x: float, w: float
) -> Tuple[float, float]:
    """
    eqsolver is the core of equpy and handles the application of Thomas Wayne Wall algorithm on the input
    """
    Cx = mass_conservation * np.exp(x)
    W = Cx / np.sum(Cx, axis=1)[:, None]
    Ks = np.prod((W / mass_conservation) ** W, axis=1) * total_masses

    M = np.vstack((stoichiometry, W))
    y = np.log(np.hstack((K, Ks)))
    delta = np.linalg.norm(np.dot(M, x.T) - y[:, None].T)

    if np.shape(stoichiometry)[1] == len(K) + len(total_masses):
        x = (w * x + np.linalg.solve(M, y)) / (w + 1)
    else:  # solve overdetermined systems as M^(-1) * y
        x = (w * x + np.dot(np.linalg.pinv(M), y)) / (w + 1)
    return x, delta
