from random import randint
from somefunctions import gen_simple_chem_formula
from sympy.solvers import solve
from sympy import Symbol
from numpy.linalg import solve
from numpy.linalg import lstsq
import numpy as np
from scipy.optimize import nnls
from sympy import Matrix
# this does not work
# do I need to write my function
# to balance stoichiometry?


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


def AAbalance_stoichiometry(products, reactants):
    all_atoms = find_unique_atoms(products, reactants)
    prod_sym = [Symbol(item) for item in all_atoms]
    #print solve(2*prod_sym[0]+1, prod_sym[0])
    return all_atoms


def find_unique_atoms(p, r):
    return list(set(flatten(p) + flatten(r)))


def flatten(l):
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out


Brucite = Phase(name='Brucite',  infile='somefilename')
print Brucite.name
print Brucite.atoms
print Brucite.natoms
print Brucite.formula
print Brucite.niceformula
#reac, prod = balance_stoichiometry({'NH4ClO4', 'Al'}, {'Al2O3', 'HCl', 'H2O', 'N2'})
#print reac
#print prod

products = [['Mg', 'O'], ['H', 'O'], ['Si', 'F']]
reactants = [['Mg', 'O', 'H'], ['K', 'L']]

print flatten(products)
print find_unique_atoms(products, reactants)
#print balance_stoichiometry(products, reactants)

# https://arxiv.org/pdf/1110.4321.pdf
def balance(matrix):
    rows, cols = matrix.shape
    if rows == cols:
        matrix = np.matrix(Matrix(matrix).rref()[0])
        counter = 0
        for row in np.absolute(matrix):
            if row.sum() == 0:
                counter += 1
        for i in range(rows-1, rows-1-counter, -1):
            matrix[i, i] = 1
    else:
        rank = l.matrix_rank(matrix)
        nullity = cols - rank
        # add nullity number of extra rows to the bottom
        # the RHS nullity x nullity dimension submatrix
        # must be idenity, set rest to zeros
        for i in range(nullity, 0, -1):
            newrow = np.zeros(cols)
            newrow[-i] = 1.
            newrow = np.matrix(newrow)
            matrix = np.vstack([matrix, newrow])
    # need to convert matrix dtype back to float
    matrix = matrix.astype('float')
    print matrix
    # compute invers
    matrix = l.inv(matrix)
    print matrix
    # extract rightmost column
    solution = (flatten(matrix[...,cols-1].tolist()))
    # divide each element by abs of min value in the column
    minvalue = absminvalue(solution)
    solution = [abs((item/minvalue)) for item in solution]
    print solution
    # we want integer result??
    #num = gcd_list(solution)
    #print num
    #print [(item/num) for item in solution]

def absminvalue(arr):
    low = abs(arr[0])
    for i in arr:
        if abs(i) < low:
            low = abs(i)
    return low

def gcd (a,b):
    if (b == 0):
        return a
    else:
        return gcd (b, a % b)
def gcd_list(l):
    res = l[0]
    for c in l[1::]:
        res = gcd(res , c)
    return res

#a= np.array([[1, 0, -1], [1, 1, -2], [0, 2, -2]])
#b = np.array([0, 0, 0])
from sympy import Matrix
matrix = np.matrix([[1, 1, 0], [2, 1, 1], [2, 0, 2]])
print matrix
import numpy.linalg as l
print l.matrix_rank(matrix)
mat = np.matrix([[1, 1, 0, 0, 0, 1],
                [1, 0, 0, 2, 0, 0],
                [0, 3, 0, 0, 1, 0],
                [0, 0, 1, 0, 2, 0],
                [0, 1, 1, 0, 0, 1]])
print l.matrix_rank(mat)
print 'shape', mat.shape
symmatrix = Matrix(matrix)
print 'aaaaaaaaaaaaaaaaa'
#matrix = np.matrix(symmatrix.rref()[0])
#matrix = np.matrix([[1,0,1],[0,2,3]])
matrix = np.matrix([[7,0,1,0], [6,0,0,2], [2,2,2,1]])
balance(matrix)
