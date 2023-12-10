import numpy as np
import matplotlib.pyplot as plt
import warnings
#github test branch_1

class ChemicalReaction:
    def __init__(self, species, stoichiometry, eq_constants, mass_conservation, initial_masses):
        initial_masses = np.array(initial_masses, dtype=float)
        if min(initial_masses) == 0:
            warnings.warn("Species concentrations (S) should not be set to zero, eps has been set instead. The result may not be reliable.", stacklevel=2)
            I = np.where(initial_masses == 0)[0]
            initial_masses[I] = 2.2e-16
    
        self.s = species
        self.N = stoichiometry
        self.K = eq_constants
        self.C = mass_conservation
        self.S = initial_masses
        self.result = []
        self.residuals = []

    def solve(self, iter, th, w):
        iter = int(iter)        
        x0 = np.ones(np.shape(self.N)[1])
        delta = []

        x, delta_ = eqsolver(self.N, self.K, self.C, self.S, x0, w)
        delta.append(delta_)

        for i in range(iter):
            x, delta_ = eqsolver(self.N, self.K, self.C, self.S, x, w)
            delta.append(delta_)

            if delta[-1] < th*np.linalg.norm(x*2.2e-16):
                self.result = np.exp(x)
                self.residuals = delta
                return np.exp(x), delta
        
        self.result = np.exp(x)
        self.residuals = delta
        warnings.warn("Tolerance not reached. Manually check if the result is satisfying. Either change starting point, increase iterations or weight.", stacklevel=2)
        return np.exp(x), delta
    
    def plotter(self):
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=False, figsize=(15,5))
        ax1.plot(np.arange(len(self.residuals)), self.residuals, lw=3, color='black')
        ax1.scatter(np.arange(len(self.residuals)), self.residuals, s=120, color='black')

        ax1.set_xlabel('steps', fontsize=12)
        ax1.set_ylabel('||r||', fontsize=12)
        ax1.set_yscale('linear')

        ax2.bar(np.arange(0, len(self.result)), self.result)
        ax2.set_ylabel('conc', fontsize=12)
        ax2.set_xticks(np.arange(0, len(self.result)))
        ax2.set_xticklabels(sorted(self.s, key=lambda x: self.s[x]), fontsize=12, rotation=90);

def eqsolver(N, K, C, S, x, w):
        Cx = C*np.exp(x)
        W = Cx/np.sum(Cx, axis=1)[:,None]
        Ks = np.prod((W/C)**W, axis=1)*S

        M = np.vstack((N, W))
        y = np.log(np.hstack((K, Ks)))
        delta = np.linalg.norm(np.dot(M, x.T) - y[:,None].T)

        if np.shape(N)[1] == len(K)+len(S):
            x = (w*x + np.linalg.solve(M, y))/(w+1)
        else: #solve overdetermined systems as M^(-1) * y
            x = (w*x + np.dot(np.linalg.pinv(M), y))/(w+1)
        return x, delta