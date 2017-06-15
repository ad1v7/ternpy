phase_names = ['Magnesia',
               'Quartz',
               'Ice',
               'Brucite']
phase_atoms = [['Mg', 'O'],
               ['Si', 'O'],
               ['H', 'O'],
               ['Mg', 'O', 'H']]
phase_natoms = [[1, 1],
                [1, 2],
                [2, 1],
                [1, 2, 2]]


def gen_simple_chem_formula(atoms, natoms):
    formula = ''
    for atom, natom in zip(atoms, natoms):
        if natom == 1:
            formula += atom
        else:
            formula += atom+str(natom)
    return formula


# this is very limited and currently only consider case
# when number of O match number of H
# otherwise it will gen simpllified formula
def gen_proper_chem_formula(atoms, natoms):
    new_atoms = atoms
    new_natoms = natoms
    if 'O' in atoms and 'H' in atoms:
        #find their indexes
        o_idx = atoms.index('O')
        h_idx = atoms.index('H')
        # find numer of atoms:
        o_atoms = natoms[o_idx]
        h_atoms = natoms[h_idx]
        if o_atoms == h_atoms:
            new_atoms.remove('O')
            new_atoms.remove('H')
            # carefully remove corresponding number of atoms
            if o_idx > h_idx:
                del new_natoms[o_idx]
                del new_natoms[h_idx]
            else:
                del new_natoms[h_idx]
                del new_natoms[o_idx]
            new_atoms.append('(OH)')
            new_natoms.append(o_atoms)
            gen_simple_chem_formula(new_atoms, new_natoms)
        else:
            gen_simple_chem_formula(atoms, natoms)
    else:
        gen_simple_chem_formula(atoms, natoms)


if __name__ == '__main__':
    gen_simple_chem_formula(phase_atoms[3], phase_natoms[3])
    gen_proper_chem_formula(phase_atoms[3], phase_natoms[3])
