'''
Module contains class to represent phase
And a bunch of usefull functions
such as
-------------------------------------
balance stoichimetry
find unique atoms in two arrays
flatten array
gcd and gcd for a list of numbers
-------------------------------------
it would make sense to move them somewhere else
but temporarly they resides here
sympy.Matrix is imported only to compute
reduced row echelon form of the matrix
I wish I could do it with numpy without back and forth convertion
from numpy matrix to sympy matrix
but it seems that algebraic method as employed in sympy is the way to go.
From scipy can use LU factorization but it is prone to errors
'''
from random import randint
from somefunctions import gen_simple_chem_formula
import numpy as np
from sympy import Matrix
import numpy.linalg as l
from fractions import Fraction
import functools


# name can come from directory name
# formula from CONTCAR
# atoms from CONTCAR
# x, y coords can be generated given numbers of constituents
# but how to deal with case when we want to use phase in
# multiple number of ternary diagrams
# it is probably better not to include xy coords in this class
# convenient class to represent phase
class Phase:
    def __init__(self, **kwargs):
        r = str(randint(0, 99))
        self.name = kwargs.get('name', 'phasename'+r)
        self.atoms = kwargs.get('atoms', 'atom types'+r)
        self.natoms = kwargs.get('natoms', [1, 2, 3])
        self.formula = kwargs.get('formula', 'formula'+r)
        self.niceformula = kwargs.get('niceformula', 'niceformula'+r)
        if 'infile' in kwargs:
            self.get_phase_data(kwargs['infile'])

    # suppose to load data from the config file
    # config file is generated during data extraction process
    # from VASP run. Config file can be modified by the user
    # for example to provide nice chemical formula
    def get_phase_data(self, infile):
        name, atoms, natoms, formula, niceformula = self.load_config(infile)
        self.name = name
        self.atoms = atoms
        self.natoms = natoms
        self.formula = gen_simple_chem_formula(atoms, natoms)
        self.niceformula = niceformula

    def load_config(self, infile):
        return 'Brucite', ['Mg', 'O', 'H'], [1, 2, 2], 'MgO2H2', 'Mg(OH)2'


# given two lists of atom symbols
# returns one with all atoms found
# without duplicates
def find_unique_atoms(p, r):
    return list(set(flatten(p) + flatten(r)))


# flatten nested lists of lists
def flatten(l):
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out


# this is currently redundant as better way
# is implemented inside balance()
# returns array of integers given array of fractions
# it is asumed that given floats have common multiplier
# which casts them to integers
# if this is not the case input array is returned
def find_int(arr):
    test = False
    epsilon = 1e-18
    maxint = 1000
    for i in range(2, maxint, 1):
        for item in arr:
            if abs(i*item-int(i*item)) < epsilon:
                test = True
            else:
                test = False
                break
        if test:
            print i
            return [int(i*item) for item in arr]
    print "Could not find one"
    return arr


# balance stoichiometry function as presented in
# https://arxiv.org/pdf/1110.4321.pdf
# the input is a numpy matrix, e.g.
# H2 + O2 <-> H2O
#
#                H2 O2 H2O
#           ---------------
#           H |  2  0   2
#           O |  0  2   1
#
# so the input matrix is
#
#           2 0 2
#           0 2 1
#
# the output is an 1D list
# [2, 1, 2]
# which gives
# 2*H2 + 1*O2 <-> 2*H2O
def balance(mtx):
    rows, cols = mtx.shape
    if rows == cols:
        # compute reduced row echelon form
        mtx = np.matrix(Matrix(mtx).rref()[0])
        counter = 0
        # find number of rows containing only zeroes
        for row in np.absolute(mtx):
            if row.sum() == 0:
                counter += 1
        # bottom right counter x counter partition matrix
        # must be converted to identity matrix
        for i in range(rows-1, rows-1-counter, -1):
            mtx[i, i] = 1
    else:
        rank = l.matrix_rank(mtx)
        nullity = cols - rank
        # add nullity number of extra rows to the bottom
        # the RHS nullity x nullity dimension partition matrix
        # must be idenity, set rest to zeros
        for i in range(nullity, 0, -1):
            newrow = np.zeros(cols)
            newrow[-i] = 1.
            newrow = np.matrix(newrow)
            mtx = np.vstack([mtx, newrow])
    # need to convert matrix dtype back to float
    mtx = mtx.astype(float)
    # compute invers
    mtx = l.inv(mtx)
    # extract rightmost column and convert to list
    solution = (flatten(mtx[..., cols-1].tolist()))
    # divide each element by abs of min value in the column
    solution = [abs(item) for item in solution]
    print 'This is a non integer solution:'
    print solution
    # convert floats to fractions and list all denominators
    denoms = [Fraction(x).limit_denominator().denominator for x in solution]
    # find the least common multiple from the denoms list
    factor = functools.reduce(lambda a, b: a*b//gcd(a, b), denoms)
    print 'This is an integer solution:'
    print [int(round(factor*item)) for item in solution]


# prepare matrix input for balance stoichiometry function
# input should be a list of Phase objects
# BUT for now it will take two lists
# one with sublist containing atom names found in a given phase
# second with sublist containg number of atoms for a given phase
# the ordering of items in lists and sublists must correspond to each other
# but this is guaranteed given data is extracted from CONTCAR file
def get_balance_matrix(phases, phases_atoms):
    atomslist = []
    for phase in phases:
        atomslist = find_unique_atoms(phase, atomslist)
    mtx = []
    for atom in atomslist:
        row = []
        for phase, phase_atoms in zip(phases, phases_atoms):
            if atom in phase:
                idx = phase.index(atom)
                row.append(phase_atoms[idx])
            else:
                row.append(0)
        mtx.append(row)
    #print 'your matrix\n', np.matrix(mtx), '\nend'
    return np.matrix(mtx)


# return minimum value from
# 1D arr of positive numbers
# negative numbers are converted to positive
# before doing test
def absminvalue(arr):
    low = abs(arr[0])
    for i in arr:
        if abs(i) < low:
            low = abs(i)
    return low


# returns greatest common divisor for a, b
def gcd(a, b):
    if b == 0:
        return a
    else:
        return gcd(b, a % b)


# computes greatest common divisor for the list
# of non-integers... gcd? non-integers? yeah I know, I explain later ;)
# currenty redundand
def gcd_list(l):
    # rounding to arbitrary number of decimals
    # to avoid numerical issues
    # 14 decimals is plenty in this case
    # and should round small error on 64 bit float
    l = np.around(l, decimals=14)
    res = l[0]
    for c in l[1:]:
        res = gcd(res, c)
    return res

if __name__ == '__main__':
    Brucite = Phase(name='Brucite',  infile='somefilename')
    products = [['Mg', 'O'], ['H', 'O'], ['Si', 'F']]
    reactants = [['Mg', 'O', 'H'], ['K', 'L']]
    phaseA = ['Mg', 'Si', 'O', 'H']
    phaseA_a = [7, 2, 14, 6]
    MgO = ['Mg', 'O']
    MgO_a = [1, 1]
    SHyB = ['Mg', 'Si', 'O', 'H']
    SHyB_a = [10, 3, 18, 4]
    phase365 = ['Mg', 'Si', 'O', 'H']
    phase365_a = [1, 1, 6, 6]
    phases = [phaseA, MgO, SHyB, phase365]
    phases_atoms = [phaseA_a, MgO_a, SHyB_a, phase365_a]
    #phases = [['Mg', 'O', 'H'], ['Mg', 'O'], ['H', 'O']]
    #phases_atoms = [[1, 2, 2], [1, 1], [2, 1]]
    #phases = [['Mg', 'O', 'H', 'Si'], ['Mg', 'O'], ['H', 'O'], ['Si', 'O']]
    mtx = get_balance_matrix(phases, phases_atoms)
    print '###############'
    balance(mtx)
    print '###############'

    print '------------------------'
    #a= np.array([[1, 0, -1], [1, 1, -2], [0, 2, -2]])
    #b = np.array([0, 0, 0])
    matrix = np.matrix([[1, 1, 0], [2, 1, 1], [2, 0, 2]])
    print matrix
    print l.matrix_rank(matrix)
    mat = np.matrix([[1, 1, 0, 0, 0, 1],
                    [1, 0, 0, 2, 0, 0],
                    [0, 3, 0, 0, 1, 0],
                    [0, 0, 1, 0, 2, 0],
                    [0, 1, 1, 0, 0, 1]])
    #matrix = np.matrix(symmatrix.rref()[0])
    #matrix = np.matrix([[1,0,1],[0,2,3]])
    #matrix = np.matrix([[2, 0, 2], [0, 2, 1]])
    matrix = np.matrix([[7, 0, 1, 0], [6, 0, 0, 2], [2, 2, 2, 1]])
    #balance(matrix)
