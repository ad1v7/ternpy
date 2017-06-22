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
import numpy as np
from sympy import Matrix
import numpy.linalg as l
from fractions import Fraction
import functools
from collections import OrderedDict


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
            print(i)
            return [int(i*item) for item in arr]
    print("Could not find one")
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
    #solution = [abs(item) for item in solution]
    # convert floats to fractions and list all denominators
    denoms = [Fraction(x).limit_denominator().denominator for x in solution]
    # find the least common multiple from the denoms list
    factor = functools.reduce(lambda a, b: a*b//gcd(a, b), denoms)
    solution = [int(round(factor*item)) for item in solution]
    #
    # currently return solution which contains negative numbers
    # use balance_interface() instead
    #
    #
    return solution


# prepare matrix input for balance stoichiometry function
# input should be a list of Phase objects
# BUT for now it will take two lists
# one with sublist containing atom names found in a given phase
# second with sublist containg number of atoms for a given phase
# the ordering of items in lists and sublists must correspond to each other
# but this is guaranteed given data is extracted from CONTCAR file

def get_balance_matrix(phases_atoms, phases_natoms):
    atomslist = []
    for atoms in phases_atoms:
        atomslist = find_unique_atoms(atoms, atomslist)
    mtx = []
    for atom in atomslist:
        row = []
        for atoms, phase_atoms in zip(phases_atoms, phases_natoms):
            if atom in atoms:
                idx = atoms.index(atom)
                row.append(phase_atoms[idx])
            else:
                row.append(0)
        mtx.append(row)
    return np.matrix(mtx)


# interface to balance()
# takes dictionary of phases
def balance_interface(dictdata, *args):
    lhs = 0
    for arg in args:
        lhs = arg
    if type(dictdata) != type(OrderedDict()):
        print(type(dictdata))
        print('Unordered dictionary supplied. Result might be wrong!')
    phases_atoms = [dictdata[phase]['atoms'] for phase in dictdata]
    phases_natoms = [dictdata[phase]['natoms'] for phase in dictdata]
    mtx = get_balance_matrix(phases_atoms, phases_natoms)
    solution = balance(mtx)
    # test if lhs and rhs of the equation have same signs
    # if not then the equation can not be balanced
    # THIS MAY NOT ALWAYS WORK: AD HOC SOLUTION
    if ((abs(sum(solution[:lhs])) == sum([abs(i) for i in solution[:lhs]])) and
        (abs(sum(solution[lhs:])) == sum([abs(i) for i in solution[lhs:]]))):
        return True, [abs(i) for i in solution]
    else:
        return False, []


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
    #balance(mtx)

    #a= np.array([[1, 0, -1], [1, 1, -2], [0, 2, -2]])
    #b = np.array([0, 0, 0])
    matrix = np.matrix([[1, 1, 0], [2, 1, 1], [2, 0, 2]])
    mat = np.matrix([[1, 1, 0, 0, 0, 1],
                    [1, 0, 0, 2, 0, 0],
                    [0, 3, 0, 0, 1, 0],
                    [0, 0, 1, 0, 2, 0],
                    [0, 1, 1, 0, 0, 1]])
    def get_corners_data(self, phase_info, corners):
        corners_dict = {}
        for corner in corners:
            try:
                corners_dict[corner] = phase_info[corner]
            except KeyError:
                print('No key in main dict: '+corner)
        return corners_dict


    #matrix = np.matrix(symmatrix.rref()[0])
    #matrix = np.matrix([[1,0,1],[0,2,3]])
    matrix = np.matrix([[2, 0, 2], [0, 2, 1]])
    #matrix = np.matrix([[7, 0, 1, 0], [6, 0, 0, 2], [2, 2, 2, 1]])
    print(balance(matrix))
