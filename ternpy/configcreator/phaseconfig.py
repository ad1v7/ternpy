try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import math
import os

from ternpy.configcreator import extractdata
#import extractdata


# Creates a phase dictionary from a file written by the user
def _get_phasedirdict(phaselist):
    phasedict = {}
    with open(phaselist, 'r') as f:
        for line in f:
            if not (line.startswith("#") or line in ['\n', '\r\n']):
                phasedict[line.split(":")[0]] = [phase.strip("\n").strip() for phase in line.split(":")[1].split(",")]
    return phasedict


# Creates a config file by searching the folders given by the phase directory dictionary.
def create_config(phaselistfile, jobdir, confdir, outcar="OUTCAR", poscar="POSCAR"):
    def all_same(items):
        return all(x == items[0] for x in items)
    phasedirdict = _get_phasedirdict(phaselistfile)
    # Dictionary containing all found information; has one entry for all phases, which contains a list of all
    # attributes found (which should all be the same)
    infodict = {}

    for phase in phasedirdict:
        infodict[phase] = {}
        infodict[phase]["structurenames"] = phasedirdict[phase]
        infodict[phase]["atoms"] = []
        infodict[phase]["natoms"] = []
        infodict[phase]["chemforms"] = []

    # Go through all the directories in the job output folder
    for (dirpath, dirnames, filenames) in os.walk(jobdir):
        if extractdata._vasp_enthalpy(dirpath, outcar) is not None or extractdata._vasp_internalenergy(dirpath, poscar) is not None:
            # Check if the current folder matches any combination of the phase and the structure
            for phase in phasedirdict:
                if phase in dirpath:
                    for struct in phasedirdict[phase]:
                        if struct in dirpath:
                            print("Extracting from:")
                            print(dirpath)
                            infodict[phase]["atoms"].append(extractdata._vasp_atoms(dirpath, outcar))
                            infodict[phase]["natoms"].append(extractdata._vasp_natoms(dirpath, poscar))
                            infodict[phase]["chemforms"].append(extractdata._vasp_chemformula(dirpath, outcar, poscar))

    # Do consistency checks; all atoms, natoms etc should be the same for all structures within a phase
    for phase in phasedirdict:
        if not all_same(infodict[phase]["atoms"]):
            print("[WARNING] Not all structures in " + phase + " have the same atoms.")
        if not all_same(infodict[phase]["natoms"]):
            print("[WARNING] Not all structures in " + phase + " have the same number of atoms.")
        if not all_same(infodict[phase]["chemforms"]):
            print("[WARNING] Not all structures in " + phase + " have the same chemical formulae.")

    # Finally, write to the configuration file
    cfgfile = confdir + "/phases.conf"
    if not os.path.exists(confdir):
        os.makedirs(confdir)
    open(confdir + "/phases.conf", 'w').close()     # Delete file contents
    with open(cfgfile, "a") as f:
        # Write the path of the job output directory for future reference
        f.write("[Directories]\n")
        f.write("\tdftdir=" + jobdir + "\n")
        f.write("\tprojectdir=" + confdir + "\n")
        f.write("\n")
        for phase in phasedirdict:
            f.write("[" + phase + "]\n")
            f.write("\tname=" + phase + "\n")
            f.write("\tplotname=\n")
            f.write("\tatoms=" + " ".join(atom for atom in infodict[phase]["atoms"][0]) + "\n")
            f.write("\tnatoms=" + " ".join(str(x) for x in infodict[phase]["natoms"][0]) + "\n")
            f.write("\tchemform=" + str(infodict[phase]["chemforms"][0]) + "\n")
            f.write("\tplotchemform=\n")
            f.write("\n")
            f.write("\tstructures=" + ", ".join(str(struct) for struct in phasedirdict[phase]) + "\n")
            f.write("\tstructureplotnames=\n")
            f.write("\n\n")


# Reads a config file and returns a dictionary with all the information
def read_config(configfile):
    phasedict = {}

    cfg = configparser.ConfigParser()
    cfg.read(configfile)
    for phase in dict(cfg.items()):
        if phase != "DEFAULT" and phase != "Directories":
            phasedict[phase] = {}
            phasedict[phase]["name"] = cfg.get(phase, "name")
            phasedict[phase]["plotname"] = cfg.get(phase, "plotname")
            phasedict[phase]["atoms"] = cfg.get(phase, "atoms").split(' ')
            phasedict[phase]["natoms"] = [int(x) for x in cfg.get(phase, "natoms").split(' ')]
            phasedict[phase]["chemform"] = cfg.get(phase, "chemform")
            phasedict[phase]["plotchemform"] = cfg.get(phase, "plotchemform")
            phasedict[phase]["structures"] = cfg.get(phase, "structures").split(", ")
            phasedict[phase]["structureplotnames"] = cfg.get(phase, "structureplotnames").split(", ")

    return phasedict


def read_dftdir(configfile):
    cfg = configparser.ConfigParser()
    cfg.read(configfile)
    return cfg.get("Directories", "dftdir")


def read_projectdir(configfile):
    cfg = configparser.ConfigParser()
    cfg.read(configfile)
    return cfg.get("Directories", "projectdir")


# Creates mesh files for P, T data containing many types of energies.
# Uses a config file to search for the correct phase/structures; datafiles will be in the same directory as confdir.
def create_datafiles(configfile, jobdir, outcar="OUTCAR", poscar="POSCAR"):
    phasedict = read_config(configfile)

    # Dictionary to temporarily hold data
    energies = {}
    for phase in phasedict:
        energies[phase] = {}
        for struct in phasedict[phase]["structures"]:
            energies[phase][struct] = {}
            energies[phase][struct]["P"] = []
            energies[phase][struct]["T"] = []
            energies[phase][struct]["E"] = []
            energies[phase][struct]["H"] = []
            energies[phase][struct]["PV"] = []
            energies[phase][struct]["G"] = []
            energies[phase][struct]["F"] = []

    # Go through all the directories in the job output folder
    for (dirpath, dirnames, filenames) in os.walk(jobdir):
        if extractdata._vasp_enthalpy(dirpath, outcar) is not None or extractdata._vasp_internalenergy(dirpath, poscar) is not None:
            # Check if the current folder matches any combination of the phase and the structure
            # TODO: make it work for temperatures, and thus for G's and F's
            for phase in phasedict:
                if phase in dirpath:
                    for struct in phasedict[phase]["structures"]:
                        if struct in dirpath:
                            print("Extracting from:")
                            print(dirpath)

                            # Rescale everything in terms of energies per formula unit
                            fu_per_ucell = float(extractdata._vasp_fu_per_ucell(dirpath))

                            press = extractdata._vasp_press(dirpath, outcar)

                            enthalpy = float(extractdata._vasp_enthalpy(dirpath, outcar))
                            if enthalpy is not None:
                                enthalpy /= fu_per_ucell
                            intenergy = float(extractdata._vasp_internalenergy(dirpath, outcar))
                            if intenergy is not None:
                                intenergy /= fu_per_ucell
                            pvterm = float(extractdata._vasp_pv(dirpath, outcar))
                            if pvterm is not None:
                                pvterm /= fu_per_ucell

                            energies[phase][struct]["P"].append(float(press))
                            energies[phase][struct]["T"].append(float('inf'))
                            energies[phase][struct]["E"].append(float(intenergy))
                            energies[phase][struct]["H"].append(float(enthalpy))
                            energies[phase][struct]["PV"].append(float(pvterm))
                            energies[phase][struct]["G"].append(float('inf'))
                            energies[phase][struct]["F"].append(float('inf'))

    def sortedindices(list):
        return sorted(range(len(list)), key=lambda x: list[x])

    # Sort the values according to pressure values
    for phase in phasedict:
        for struct in phasedict[phase]["structures"]:
            sortedidx = sortedindices(energies[phase][struct]["P"])
            energies[phase][struct]["P"] = [energies[phase][struct]["P"][i] for i in sortedidx]
            energies[phase][struct]["T"] = [energies[phase][struct]["T"][i] for i in sortedidx]
            energies[phase][struct]["E"] = [energies[phase][struct]["E"][i] for i in sortedidx]
            energies[phase][struct]["H"] = [energies[phase][struct]["H"][i] for i in sortedidx]
            energies[phase][struct]["PV"] = [energies[phase][struct]["PV"][i] for i in sortedidx]
            energies[phase][struct]["G"] = [energies[phase][struct]["G"][i] for i in sortedidx]
            energies[phase][struct]["F"] = [energies[phase][struct]["F"][i] for i in sortedidx]

    # Finally, write energy values to file
    for phase in phasedict:
        if not os.path.exists("energies"):
            os.makedirs("energies")
        for struct in phasedict[phase]["structures"]:
            with open("energies/" + struct + ".dat", "w") as f:
                f.write("P\tT\tE\tH\tPV\tG\tF\n")
                for idx, press in enumerate(energies[phase][struct]["P"]):
                    line = "{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\n".format(
                                                                 float(energies[phase][struct]["P"][idx]),
                                                                 float(energies[phase][struct]["T"][idx]),
                                                                 float(energies[phase][struct]["E"][idx]),
                                                                 float(energies[phase][struct]["H"][idx]),
                                                                 float(energies[phase][struct]["PV"][idx]),
                                                                 float(energies[phase][struct]["G"][idx]),
                                                                 float(energies[phase][struct]["F"][idx]))

                    f.write(line)














def extract(phaselist, dftdir, poscar="POSCAR", outcar="OUTCAR", phonopyfile="thermal_properties.yaml"):
    def all_same(items):
        return all(x == items[0] for x in items)

    projdir = os.path.dirname(os.path.realpath(phaselist))

    phasedirdict = _get_phasedirdict(phaselist)

    print("\nSearching for the following folders:")
    for phase in phasedirdict:
        print(phase + ":")
        for struct in phasedirdict[phase]:
            print("\t" + struct)

    # Dictionary containing all found information; has one entry for all phases, which contains a list of all
    # attributes found (which should all be the same)
    infodict = {}

    # Dictionary to temporarily hold data
    energies = {}

    # Create entries in dictionaries
    for phase in phasedirdict:
        energies[phase] = {}

        infodict[phase] = {}
        infodict[phase]["structurenames"] = phasedirdict[phase]
        infodict[phase]["atoms"] = []
        infodict[phase]["natoms"] = []
        infodict[phase]["chemforms"] = []

        # Add entries for structure energies
        for struct in phasedirdict[phase]:
            energies[phase][struct] = {}
            energies[phase][struct]["P"] = []
            energies[phase][struct]["T"] = []
            energies[phase][struct]["E"] = []
            energies[phase][struct]["H"] = []
            energies[phase][struct]["PV"] = []
            energies[phase][struct]["G"] = []
            energies[phase][struct]["F"] = []

    print("\nStarting directory walk...")
    # Go through all the directories in the job output folder
    for (dirpath, dirnames, filenames) in os.walk(dftdir):
        if extractdata._vasp_enthalpy(dirpath, outcar) is not None or extractdata._vasp_internalenergy(dirpath, poscar) is not None:
            # Check if the current folder matches any combination of the phase and the structure
            for phase in phasedirdict:
                if phase in dirpath:
                    for struct in phasedirdict[phase]:
                        if struct in dirpath:
                            print("Extracting from:")
                            print(dirpath)

                            # Extract phase information
                            infodict[phase]["atoms"].append(extractdata._vasp_atoms(dirpath, outcar))
                            infodict[phase]["natoms"].append(extractdata._vasp_natoms(dirpath, poscar))
                            infodict[phase]["chemforms"].append(extractdata._vasp_chemformula(dirpath, outcar, poscar))

                            # Extract energies
                            # Rescale everything in terms of energies per formula unit
                            fu_per_ucell = float(extractdata._vasp_fu_per_ucell(dirpath))

                            press = extractdata._vasp_press(dirpath, outcar)

                            enthalpy = float(extractdata._vasp_enthalpy(dirpath, outcar))
                            if enthalpy is not None:
                                enthalpy /= fu_per_ucell
                            intenergy = float(extractdata._vasp_internalenergy(dirpath, outcar))
                            if intenergy is not None:
                                intenergy /= fu_per_ucell
                            pvterm = float(extractdata._vasp_pv(dirpath, outcar))
                            if pvterm is not None:
                                pvterm /= fu_per_ucell
                            free_energies = extractdata._free_energy(dirpath, phonopyfile)
                            # Rescale to eV per formula unit to avoid having to do it later
                            for temp in free_energies.keys():
                                free_energies[temp] /= fu_per_ucell

                            # Only write free energies if there are any
                            if free_energies:
                                for temp in free_energies.keys():
                                    energies[phase][struct]["P"].append(float(press))
                                    energies[phase][struct]["T"].append(temp)
                                    energies[phase][struct]["E"].append(float(intenergy))
                                    energies[phase][struct]["H"].append(float(enthalpy))
                                    energies[phase][struct]["PV"].append(float(pvterm))
                                    energies[phase][struct]["G"].append(float(enthalpy) + free_energies[temp])
                                    energies[phase][struct]["F"].append(free_energies[temp])
                            else:
                                energies[phase][struct]["P"].append(float(press))
                                energies[phase][struct]["T"].append(float('inf'))
                                energies[phase][struct]["E"].append(float(intenergy))
                                energies[phase][struct]["H"].append(float(enthalpy))
                                energies[phase][struct]["PV"].append(float(pvterm))
                                energies[phase][struct]["G"].append(float('inf'))
                                energies[phase][struct]["F"].append(float('inf'))

                    # Check for duplicates
                    dup = [x for x in energies[phase][struct]["P"] if energies[phase][struct]["P"].count(x) >= 2]
                    if dup:
                        print("\n\n[WARNING] Found duplicate pressure value: ")
                        print("Structure: " + struct + ", press=" + str(dup[0]) + "\n\n")

    # Do consistency checks; all atoms, natoms etc should be the same for all structures within a phase
    for phase in phasedirdict:
        if not all_same(infodict[phase]["atoms"]):
            print("[WARNING] Not all structures in " + phase + " have the same atoms.")
        if not all_same(infodict[phase]["natoms"]):
            print("[WARNING] Not all structures in " + phase + " have the same number of atoms.")
        if not all_same(infodict[phase]["chemforms"]):
            print("[WARNING] Not all structures in " + phase + " have the same chemical formulae.")

    # Finally, write to the configuration file
    print("\nCreating phases configuration file...\n")

    cfgfile = projdir + "/phases.conf"
    if not os.path.exists(projdir):
        os.makedirs(projdir)
    open(projdir + "/phases.conf", 'w').close()     # Delete file contents
    with open(cfgfile, "a") as f:
        # Write the path of the job output directory for future reference
        f.write("[Directories]\n")
        f.write("\tdftdir=" + os.path.abspath(dftdir) + "\n")
        f.write("\tprojectdir=" + projdir + "\n")
        f.write("\n")
        for phase in phasedirdict:
            f.write("[" + phase + "]\n")
            f.write("\tname=" + phase + "\n")
            f.write("\tplotname=\n")
            f.write("\tatoms=" + " ".join(atom for atom in infodict[phase]["atoms"][0]) + "\n")
            f.write("\tnatoms=" + " ".join(str(x) for x in infodict[phase]["natoms"][0]) + "\n")
            f.write("\tchemform=" + str(infodict[phase]["chemforms"][0]) + "\n")
            f.write("\tplotchemform=\n")
            f.write("\n")
            f.write("\tstructures=" + ", ".join(str(struct) for struct in phasedirdict[phase]) + "\n")
            f.write("\tstructureplotnames=\n")
            f.write("\n\n")

    print("\nWriting energies to files...\n")

    def sortedindices(list):
        return sorted(range(len(list)), key=lambda x: list[x])

    # Sort the values according to pressure values
    for phase in phasedirdict:
        for struct in phasedirdict[phase]:
            sortedidx = sortedindices(energies[phase][struct]["P"])
            energies[phase][struct]["P"] = [energies[phase][struct]["P"][i] for i in sortedidx]
            energies[phase][struct]["T"] = [energies[phase][struct]["T"][i] for i in sortedidx]
            energies[phase][struct]["E"] = [energies[phase][struct]["E"][i] for i in sortedidx]
            energies[phase][struct]["H"] = [energies[phase][struct]["H"][i] for i in sortedidx]
            energies[phase][struct]["PV"] = [energies[phase][struct]["PV"][i] for i in sortedidx]
            energies[phase][struct]["G"] = [energies[phase][struct]["G"][i] for i in sortedidx]
            energies[phase][struct]["F"] = [energies[phase][struct]["F"][i] for i in sortedidx]

    # Finally, write energy values to file
    for phase in phasedirdict:
        if not os.path.exists("energies"):
            os.makedirs("energies")
        for struct in phasedirdict[phase]:
            with open("energies/" + phase + "-" + struct + ".dat", "w") as f:
                f.write("P\tT\tE\tH\tPV\tG\tF\n")
                for idx, press in enumerate(energies[phase][struct]["P"]):
                    line = "{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\n".format(
                                                                 float(energies[phase][struct]["P"][idx]),
                                                                 float(energies[phase][struct]["T"][idx]),
                                                                 float(energies[phase][struct]["E"][idx]),
                                                                 float(energies[phase][struct]["H"][idx]),
                                                                 float(energies[phase][struct]["PV"][idx]),
                                                                 float(energies[phase][struct]["G"][idx]),
                                                                 float(energies[phase][struct]["F"][idx]))

                    f.write(line)


    print("\nFinished.\n")

















if __name__ == "__main__":
    # create_config("phaselist.conf", "joboutput", "configs")
    # print(read_config("configs"))
    # create_datafiles("configs/phases.conf", "joboutput")
    extract("phaselist.conf", "joboutput")
