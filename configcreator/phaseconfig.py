from configcreator import extractdata
import os
import configparser


# Creates a phase dictionary from a file written by the user
def __get_phasedirdict(phaselist):
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
    phasedirdict = __get_phasedirdict(phaselistfile)
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
        if extractdata.__vasp_enthalpy(dirpath, outcar) is not None or extractdata.__vasp_internalenergy(dirpath, poscar) is not None:
            # Check if the current folder matches any combination of the phase and the structure
            for phase in phasedirdict:
                if phase in dirpath:
                    for struct in phasedirdict[phase]:
                        if struct in dirpath:
                            print("Extracting from:")
                            print(dirpath)
                            infodict[phase]["atoms"].append(extractdata.__vasp_atoms(dirpath, outcar))
                            infodict[phase]["natoms"].append(extractdata.__vasp_natoms(dirpath, poscar))
                            infodict[phase]["chemforms"].append(extractdata.__vasp_chemformula(dirpath, outcar, poscar))

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
    open(confdir + "/phases.conf", 'w').close()     # Delete file contents
    with open(cfgfile, "a") as f:
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
def read_config(confdir):
    phasedict = {}

    cfg = configparser.ConfigParser()
    cfg.read(confdir + "/phases.conf")
    for phase in dict(cfg.items()):
        if phase != "DEFAULT":
            phasedict[phase] = {}
            phasedict[phase]["name"] = cfg.get(phase, "name")
            phasedict[phase]["plotname"] = cfg.get(phase, "plotname")
            phasedict[phase]["atoms"] = cfg.get(phase, "atoms")
            phasedict[phase]["natoms"] = cfg.get(phase, "natoms")
            phasedict[phase]["chemform"] = cfg.get(phase, "chemform")
            phasedict[phase]["plotchemform"] = cfg.get(phase, "plotchemform")
            phasedict[phase]["structures"] = cfg.get(phase, "structures").split(", ")
            phasedict[phase]["structureplotnames"] = cfg.get(phase, "structureplotnames").split(", ")

    return phasedict


# Creates mesh files for P, T data containing many types of energies.
# Searches for any data for the compound 'name', and puts the output data file in 'outdir'.
# 'namedepth' defines how many folders back will be used to name the structure, which by default is the folder the
# OUTCAR, POSCAR etc is found in.
# TODO: make it scan pressure directories
def create_datafile(name, outdir, filename, foldernameoffset=0):
    with open(filename, "w") as f:
        f.write("P\tT\tH\tPV\tG\tF")


if __name__ == "__main__":
    # create_datafile("joboutput/Quartz/alpha-quartz/80", "./", "test.dat")
    create_config("phaselist.conf", "joboutput", "configs")
    read_config("configs")
