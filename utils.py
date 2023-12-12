import re, csv
import numpy as np
from math import floor, log10
from typing import List, Set, Dict, Tuple


def define_species_set(eq: List[str]) -> Set[str]:
    """
    Extract species involved in the reaction from symbolic equations
    Parameters:
    - eq: list of equations in symbolic form
    """
    e_ = " ".join(eq)
    e_ = e_.split(" ")

    species_set = []

    for l in e_:
        if l != "+" and l != "=":
            species_set.append(l)

    for i, entry in enumerate(species_set):
        if entry[0].isdigit():
            entry_split = re.findall("\d+|\D+", entry)
            species_set[i] = "".join(entry_split[1:])

    return set(species_set)


def define_reactions(eq: List[str]) -> Tuple[np.ndarray, int, Dict[str, int]]:
    """
    Converts symbolic equations defining chemical reactions to matrix form suitable for the solver
    Parameters:
    - eq: list of equations in symbolic form
    """

    species_set = define_species_set(eq)
    n = len(species_set)
    species = {string: index for index, string in enumerate(sorted(species_set))}

    N = np.zeros((len(eq), n))

    for i, eq_ in enumerate(eq):
        eq_ = eq_.split("=")
        eq_left = eq_[0].split("+")

        for element in eq_left:
            element = element.strip()
            element_ = re.findall("\d+|\D+", element)  # obtain blocks for every species

            if len(element_) == 1:
                N[i, species[element_[0]]] = -1
            elif len(element_) > 1 and element_[0].isalpha():
                N[i, species["".join(element_)]] = -1
            else:
                N[i, species["".join(element_[1:])]] = -int(element_[0])

        eq_right = eq_[1].split("+")
        for element in eq_right:
            element = element.strip()
            element_ = re.findall("\d+|\D+", element)  # obtain blocks for every species

            if len(element_) == 1:
                N[i, species[element_[0]]] = 1
            elif len(element_) > 1 and element_[0].isalpha():
                N[i, species["".join(element_)]] = 1
            else:
                N[i, species["".join(element_[1:])]] = element_[0]

    return N, n, species


def define_conservations(eq: List[str], n: int, species: Dict[str, int]) -> np.ndarray:
    """
    Converts symbolic equations defining mass conservations to matrix form suitable for the solver
    Parameters:
    - eq: list of equations in symbolic form
    - n: number of reactions
    - species: dictionary containing {name: index} of the species involved in the reactions
    """

    pattern = r"(\d*\.?\d*)?([A-Za-z\d]+)"
    M = np.zeros((len(eq), n))
    eq = [eq_.replace("*", "") for eq_ in eq]

    for i, eq_ in enumerate(eq):
        eq_ = eq_.split("+")

        for element_ in eq_:
            element_ = re.findall(pattern, element_)
            element_ = [elem for elem in element_[0]]
            if len(element_[0]) == 0:
                M[i, species[element_[1]]] = 1
            else:
                M[i, species[element_[1]]] = element_[0]

    return M


def eq_system_builder(
    eq: List[str], mass_conservation: List[str]
) -> Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
    """
    Converts symbolic equations defining chemical reactions and mass conservations to matrix form suitable for the solver
    Parameters:
    - eq: list of reactions in symbolic form
    - mass_conservation: list of mass conservations in symbolic form
    """

    N, n, species = define_reactions(eq)
    C = define_conservations(mass_conservation, n, species)

    return N, C, species


def csv_loader(
    filename_N: str, filename_C: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Assembles matrices suitable for the solver starting from .csv files
    """

    N = []
    K = []
    C = []
    S = []

    with open(filename_N, newline="", encoding="utf-8-sig") as csvfile:
        filecontent = csv.reader(csvfile, delimiter=" ", quotechar="|")
        for i, row in enumerate(filecontent):
            if i == 0:
                lbls = row[0].split(",")[0:-1]
            else:
                N.append(row[0].split(",")[0:-1])
                K.append(row[0].split(",")[-1])

    with open(filename_C, newline="", encoding="utf-8-sig") as csvfile:
        filecontent = csv.reader(csvfile, delimiter=" ", quotechar="|")
        for i, row in enumerate(filecontent):
            if i > 0:
                C.append(row[0].split(",")[0:-1])
                S.append(row[0].split(",")[-1])

    N = np.array(N, dtype="float")
    K = np.array(K, dtype="float")
    C = np.array(C, dtype="float")
    S = np.array(S, dtype="float")

    return N, K, C, S
