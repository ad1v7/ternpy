from functools import reduce


# Figures out the chemical formula of a VASP calculation using POTCAR and OUTCAR.
def vasp_chemformula(dir, OUTCAR, POSCAR):
    atoms = []
    numbers = []

    with open(dir + "/" + OUTCAR) as f:
        for line in f:
            if "VRHFIN" in line:
                # Append the atom's name from the line containing VRHFIN in POTCAR
                atoms.append(line.split("=")[1].split(":")[0])

    with open(dir + "/" + POSCAR) as f:
        lines = f.readlines()
        numbers = [int(number) for number in lines[5].split()]

    def gcd(a, b):
        if b == 0:
            return a
        else:
            return gcd(b, a % b)

    def gcd_list(list):
        return reduce(gcd, tuple(list))

    numbers = [int(num / gcd_list(numbers)) for num in numbers]

    chemform = ""

    for idx, atom in enumerate(atoms):
        chemform += atom
        if numbers[idx] != 1:
            chemform += str(numbers[idx])

    return chemform


# Returns the number of f.u.'s per unit cell.
def vasp_fu_per_ucell(dir):
    pass


# Obtains the enthalpy term from a VASP OUTCAR file.
def vasp_enthalpy(dir, OUTCAR):
    energies = []
    with open(dir + "/" + OUTCAR) as f:
        for line in f:
            if "enthalpy" in line:
                energies = line.split()
    return energies[4]


# Obtains the PV term from a VASP OUTCAR file.
def vasp_pv(dir, OUTCAR):
    energies = []
    with open(dir + "/" + OUTCAR) as f:
        for line in f:
            if "enthalpy" in line:
                energies = line.split()
    return energies[8]


def phasename(dir):
    pass


# Obtains a dict of free energies from a thermal_properties.yaml file.
def free_energy(dir):
    pass

if __name__ == "__main__":
    print(vasp_chemformula("joboutput/Quartz/alpha-quartz/80", "OUTCAR", "POSCAR"))
