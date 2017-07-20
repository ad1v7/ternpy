import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from scipy.interpolate import griddata
import math

from ternpy.configcreator import phaseconfig

# TODO: Creates a full P-T phase diagram for a specified phase from Gibbs free energies, including decomposition


def ptplot(phase, Pmin, Pmax, Pstep, Tmin, Tmax, Tstep, method="linear"):
    def frange(x, y, jump):
        while x < y:
            yield x
            x += jump

    phasedict = phaseconfig.read_config("phases.conf")

    structlabels = phasedict[phase]["structures"]
    print(structlabels)

    # TODO: replace structlabels with another array, including decompositions; call it 'structanddecomps' or something
    coords = []
    imvals = np.zeros([int((Tmax-Tmin)/Tstep), int((Pmax-Pmin)/Pstep)], dtype=int)

    cmap, norm = mpl.colors.from_levels_and_colors([0, 1, 2, 3], ['red', 'green', 'blue', 'cyan'], extend="min")
    for T in frange(Tmin, Tmax, Tstep):
        for P in frange(Pmin, Pmax, Pstep):
            gibbsenergies = []  # Empty the list contents
            for struct in structlabels:
                gibbsenergies.append(gibbsinterp(phase, struct, T, P, method))

            idxlowest = np.argmin(np.array(gibbsenergies))
            #coords.append((T, P, structlabels[idxlowest]))
            imvals[int((T-Tmin)/Tstep), int((P-Pmin)/Pstep)] = idxlowest
            #plt.scatter(T, P, c=idxlowest, cmap=cmap, norm=norm)

    plt.figure(figsize=(4, 4))
    im = plt.imshow(imvals, cmap='jet', interpolation='none', origin='lower')

    plt.show()


# Returns the Gibbs free energy surface of a structure in a specified P-T range.
def gibbsinterp(phase, structure, Tval, Pval, method):
    T = np.empty(0)
    P = np.empty(0)
    G = np.empty(0)

    with open("energies/" + phase + "-" + structure + ".dat", "r") as f:
        for l in f:
            values = l.split()
            if not values[0].isalpha():
                p = values[0]
                t = values[1]
                g = values[5]
                T = np.append(T, [t])
                P = np.append(P, [p])
                G = np.append(G, [g])

    Gq = griddata((T, P), G, (Tval, Pval), method)

    return Gq



if __name__ == "__main__":
    ptplot("Quartz", 0, 160, 5, 0, 2000, 100, method="linear")