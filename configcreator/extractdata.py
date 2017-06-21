from functools import reduce
import os


# Figures out the chemical formula of a VASP calculation using POTCAR and OUTCAR.
def _vasp_chemformula(dir, outcar="OUTCAR", poscar="POSCAR"):
    atoms = []
    numbers = []

    if not os.path.isfile(dir + "/" + outcar):
        return None
    if not os.path.isfile(dir + "/" + poscar):
        return None

    with open(dir + "/" + outcar) as f:
        for line in f:
            if "VRHFIN" in line:
                # Append the atom's name from the line containing VRHFIN in POTCAR
                atoms.append(line.split("=")[1].split(":")[0])

    with open(dir + "/" + poscar) as f:
        lines = f.readlines()
        if all(item.isdigit() for item in lines[5].split()):
            numbers = [int(number) for number in lines[5].split()]
        else:
            numbers = [int(number) for number in lines[6].split()]

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


def _vasp_atoms(dir, outcar="OUTCAR"):
    atoms = []
    if not os.path.isfile(dir + "/" + outcar):
        return None
    with open(dir + "/" + outcar) as f:
        for line in f:
            if "VRHFIN" in line:
                # Append the atom's name from the line containing VRHFIN in POTCAR
                atoms.append(line.split("=")[1].split(":")[0])
    return atoms


def _vasp_natoms(dir, poscar="POSCAR"):
    if not os.path.isfile(dir + "/" + poscar):
        return None
    with open(dir + "/" + poscar) as f:
        lines = f.readlines()
        if all(item.isdigit() for item in lines[5].split()):
            numbers = [int(number) for number in lines[5].split()]
        else:
            numbers = [int(number) for number in lines[6].split()]

    def gcd(a, b):
        if b == 0:
            return a
        else:
            return gcd(b, a % b)

    def gcd_list(list):
        return reduce(gcd, tuple(list))

    numbers = [int(num / gcd_list(numbers)) for num in numbers]

    return numbers


# TODO
# Returns the number of f.u.'s per unit cell.
def _vasp_fu_per_ucell(dir):
    pass


# Obtains the enthalpy term from a VASP OUTCAR file.
def _vasp_enthalpy(dir, outcar="OUTCAR"):
    energies = []
    if not os.path.isfile(dir + "/" + outcar):
        return None

    with open(dir + "/" + outcar) as f:
        for line in f:
            if "enthalpy" in line:
                energies = line.split()
    if energies:
        return energies[4]
    else:
        return None


def _vasp_internalenergy(dir, outcar="OUTCAR"):
    if _vasp_enthalpy(dir, outcar) is not None and _vasp_enthalpy(dir, outcar) is not None:
        return str(float(_vasp_enthalpy(dir, outcar)) - float(_vasp_pv(dir, outcar)))
    else:
        return None


# Obtains the PV term from a VASP OUTCAR file.
def _vasp_pv(dir, outcar="OUTCAR"):
    energies = []
    if not os.path.isfile(dir + "/" + outcar):
        return None
    with open(dir + "/" + outcar) as f:
        for line in f:
            if "enthalpy" in line:
                energies = line.split()
    if energies:
        return energies[8]
    else:
        return None


def _vasp_press(dir, outcar="OUTCAR"):
    press = -1
    if not os.path.isfile(dir + "/" + outcar):
        return None
    with open(dir + "/" + outcar) as f:
        for line in f:
            if "PSTRESS" in line:
                press = line.split("=")[1].split()[0]
    return press


# TODO
# Obtains a dict of free energies from a thermal_properties.yaml file.
def _free_energy(dir):
    pass


if __name__ == "__main__":
    print(_vasp_press("joboutput/Quartz/alpha-quartz/80", "OUTCAR"))
    #print(vasp_chemformula("joboutput/Quartz/alpha-quartz/80", "OUTCAR", "POSCAR"))
