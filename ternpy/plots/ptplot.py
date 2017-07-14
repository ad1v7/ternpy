import numpy as np
from scipy.interpolate import griddata
import math

from ternpy.configcreator import phaseconfig

# TODO: Creates a full P-T phase diagram for a specified phase from Gibbs free energies, including decomposition


def ptplot(phase):
    pass


# Returns the Gibbs free energy surface of a structure in a specified P-T range.
def ptinterp(phase, structure, Pmin, Pmax, Pstep, Tmin, Tmax, Tstep, method="cubic"):
    T = np.empty(0)
    P = np.empty(0)
    G = np.empty(0)

    with open("energies/" + phase + "-" + structure + ".dat", "r") as f:
        for l in f:
            t, p, g = l.split()
            T = np.append(T, [t])
            P = np.append(P, [p])
            G = np.append(G, [g])

    Tq = np.linspace(Tmin, Tmax, (Tmax - Tmin) / Tstep)
    Pq = np.linspace(Pmin, Pmax, (Pmax - Pmin) / Pstep)

    Gq = griddata((T, P), G, (Tq[None, :], Pq[:, None]), method)

    return Tq, Pq, Gq
