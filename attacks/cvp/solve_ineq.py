from dataclasses import dataclass
from sage.all    import *
from itertools   import chain
from CVP         import CVP

"""
    DISCLAIMER:
    This code is not mine.
    Idea is taken from "https://org.anize.rs/SECCON-2021/crypto/signwars"
    Re-coded for understanding purpose only.
"""

"""
    //class// Inequation:

    "" Purpose:
        Define a multivariate linear function and its constraints.

        The corresponding formula is:
            lower_bound <= sum(coefficients[var] * var, for all var) <= upper_bound
"""
@dataclass
class Inequation:
    function:     dict
    lower_bound:  int
    upper_bound:  int

    def __str__(self) -> str:
        formula = ' + '.join(f'{coefficient}*{varName}' for varName, coefficient in self.function.items() if coefficient != 0)
        return f'{self.lower_bound} <= {formula} <= {self.upper_bound}'

"""
    /function/ inequations_to_matrix():

    "" Purpose:
        

    "" Args:
        inequations:   List of inequations.
        verbose=False: Decide if we're going to display more info.
"""
def inequations_to_matrix(inequations: list[Inequation], verbose=False):
    variables = sorted(list(set(chain.from_iterable(
        inequation.function.keys() for inequation in inequations
    ))))

    basis = [[0] * len(inequations) for _ in range(len(variables))]
    for i, inequation in enumerate(inequations):
        for varName, coefficient in inequation.function.items():
            basis[variables.index(varName)][i] = coefficient

    if verbose:
        print(f'[ i ] Inequations info:')
        print(f'      [ + ] Variables: {variables}')
        print(f'      [ + ] Matrix rows: {len(variables)} variables.')
        print(f'      [ + ] Matrix cols: {len(inequations)} inequations.')
        print(f'[ i ] Result matrix:')
        for varIndex, row in enumerate(basis):
            print(''.join('1' if col == 1 else '*' if col else '.' for col in row) + f' -> var: {variables[varIndex]}')

    return matrix(basis), variables

"""
    /function/ linear_functions_to_matrix():

    "" Purpose:
        

    "" Args:
        ln_fns:        List of linear functions.
        lower_bounds:  Lower bound of the functions.
        upper_bounds:  Upper bound of the functions.
        verbose=False: Decide if we're going to display more info.
"""
def linear_functions_to_inequations(ln_fns: list, lower_bounds: list, upper_bounds: list, verbose=False):
    assert len(ln_fns) == len(lower_bounds) == len(upper_bounds), ValueError(f"List of functions should be the same length ({len(ln_fns)}) as lower bounds ({len(lower_bounds)}) and upper bounds ({len(upper_bounds)}).")

    # Check if one of them is polynomial in ZZ?
    parent = None
    for f in ln_fns:
        if 'Polynomial Ring' in str(f.parent()):        # works in Sage version 9.5
            if 'Integer Ring' in str(f.parent()):       # works in Sage version 9.5
                parent = f.parent()
                break
    assert parent != None, ValueError("There should be a function in ln_fns that is a polynomial in ZZ.")

    # Can all of them be forced into type PolynomialRing in ZZ?
    _ln_fns = deepcopy(ln_fns)
    try:
        for i in range(len(ln_fns)):
            _ln_fns[i] = parent(ln_fns[i])
    except:
        raise ValueError(f"All functions in ln_fns should be able to be converted to type \"{parent}\", but \"{ln_fns[i]}\" cannot.")
    
    # Check if they're all linear functions?
    parent_gens = parent.gens()
    for f in _ln_fns:
        for f_mono in f.monomials():
            if f_mono == 1:
                continue
            elif f_mono not in parent_gens:
                raise ValueError(f"\"{f}\" is not a linear function.")

    # Now we're building the dictionaries!
    parent_gens_map_str = {
        str(var): var
            for var in parent_gens
    }

    inequations = []
    for i, f in enumerate(_ln_fns):
        if f.is_constant():
            assert lower_bounds[i] <= f <= upper_bounds[i], ValueError(f"Constant function \"{f}\" does not satisfy {lower_bounds[i]} <= {f} <= {upper_bounds[i]}!")
            continue

        inequation_terms = {
            varname: f[parent_gens_map_str[varname]]
                for varname in parent_gens_map_str
        }

        inequation_lower_bound = lower_bounds[i] - f[parent(1)]
        inequation_upper_bound = upper_bounds[i] - f[parent(1)]

        inequations.append(Inequation(
            inequation_terms,
            inequation_lower_bound,
            inequation_upper_bound
        ))

    return inequations

"""
    /function/ solve():

    "" Purpose:
        

    "" Args:
        
"""
def solve(ln_fns: list, lower_bounds: list, upper_bounds: list, verbose=False, weight=None):
    # Convert to inequations format
    inequations = linear_functions_to_inequations(ln_fns, lower_bounds, upper_bounds, verbose)

    if verbose:
        print(f'[ i ] Start solving these inequations...')
        print(f'[ * ] List of inequations:')
        print('\n'.join(f'  L___ {inequation}' for inequation in inequations))

    # Convert to basis matrix to put into CVP
    basis, variables = inequations_to_matrix(inequations, verbose)

    # Get result CVP vector
    if verbose:
        print("[ i ] Running rkm0959's close vector algorithm...")
    close_vec = CVP(
        basis,
        [inequation.lower_bound for inequation in inequations],
        [inequation.upper_bound for inequation in inequations],
        weight,
        verbose=False
    )

    if not close_vec:
        raise Exception('[ ! ] Bad luck... No result!')

    # Solve matrix equations to get solution...
    if basis.nrows() != basis.ncols():
        weighted_lattice = matrix(basis)
        H, U             = weighted_lattice.hermite_form(transformation=True)
        result_vec       = H.solve_left(close_vec).change_ring(ZZ) * U
    else:
        result_vec = basis.transpose().inverse() * close_vec
    
    # Yey, finally!
    solution = dict(zip(variables, result_vec))
    if verbose:
        print('[ i ] Solution:')
        for variable in solution:
            print(f'      [ + ] {variable} = {solution[variable]}')

    # Sanity check for solution
    check = True
    for inequation in inequations:
        function, lower_bound, upper_bound = inequation.function, inequation.lower_bound, inequation.upper_bound
        function_value                     = sum(coefficient * solution[varName] for varName, coefficient in inequation.function.items())

        if not lower_bound <= function_value <= upper_bound:
            check = False

    if check:
        return solution
    else:
        raise Exception(f'[ ! ] Failed constraint checking!')
