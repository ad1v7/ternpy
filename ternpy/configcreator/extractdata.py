from functools import reduce
import os
import yaml


# Conversion factor for kJ/mol to eV
KJMOL_TO_EV = 0.01036410


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


# Returns the no. of formula units per unit cell, which is the gcd of the numbers of each atom present in a unit cell.
def _vasp_fu_per_ucell(dir, poscar="POSCAR"):
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

    return gcd_list(numbers)


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


# Obtains a dict of free energies from a thermal_properties.yaml file, mapping temperature to the corresponding free
# energy in units of eV per unit cell.
def _free_energy(dir, phonopyfile="thermal_properties.yaml"):
    energy = {}

    # Read yaml file for vibrational free energies
    with open(dir + "/" + phonopyfile, 'r') as stream:
        data_dict = yaml.load(stream)
        data_dict = data_dict['thermal_properties']

        for entry in data_dict:
            tempstr = entry.get('temperature')
            fenergystr = entry.get('free_energy')
            temperature = float(tempstr)
            fenergy = float(fenergystr) * KJMOL_TO_EV
            energy[temperature] = fenergy

    return energy


if __name__ == "__main__":
    # print(_vasp_press("joboutput/Quartz/alpha-quartz/80", "OUTCAR"))
    # print(vasp_chemformula("joboutput/Quartz/alpha-quartz/80", "OUTCAR", "POSCAR"))
    print(_free_energy("joboutput/Quartz/alpha-quartz/80"))
