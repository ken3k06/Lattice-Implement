from sage.all import *
from sage.modules.free_module_integer import IntegerLattice
import itertools, subprocess
"""
    DISCLAIMER:
    This code is not mine.
    Idea is taken from "https://github.com/rkm0959/Inequality_Solving_with_CVP/blob/main/solver.sage"
    Re-coded for understanding purpose only.
"""

"""
    /function/ BabaiCVP():

    "" Purpose:
        Find the closest vector to a target vector.
        Original algorithm created by Babai.

    "" Args:
        basis:  The matrix composed of basis vectors.
        target: The target vector.
"""
def dump(M):
    return "[{}]".format("\n".join("[{}]".format(" ".join(map(str, r))) for r in M))

def parse_row(r):
    return map(ZZ, r.split())

def parse(x):
    x = x.replace("\n", "")
    assert x[0] == "["
    assert x[-1] == "]"
    return Matrix(ZZ, map(parse_row, x[2:-2].split("][")))
def flatter(M):
    try:
        PATH = FLATTER_PATH
    except NameError:
        print("Warning: using `flatter` if it can be found on $PATH. Set `FLATTER_PATH` to the command if this fails.")
        PATH = "flatter"

    try:
        ARGS = FLATTER_ARGS
    except:
        ARGS = []

    return parse(subprocess.check_output([PATH, *ARGS], input=dump(M).encode()).decode())
# Or
def flatter_(M):
    import os
    import re
    from subprocess import check_output
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(os.cpu_count())
    ret = subprocess.check_output(["flatter"], input=z.encode(), env=env)
    return matrix(M.nrows(), M.ncols(), map(int, findall(rb"-?\d+", ret)))

def BabaiCVP(basis, target):
    basis = basis.LLL()
    # basis = flatter(basis)
    # print(basis)
    # basis   = IntegerLattice(basis, lll_reduce=True).reduced_basis # slow
    reduced = basis.gram_schmidt()[0]
    diff    = target

    for i in reversed(range(basis.nrows())):
        diff -= basis[i] * ((diff * reduced[i]) / (reduced[i] * reduced[i])).round()

    return target - diff


"""
    /function/ CVP_rkm0959():

    "" Purpose:
        Approximate a vector in the range [l, r].
        Seems to have been optimised with scaling and stuff (?)

    "" Author:
        rkm0959

    "" Args:
        basis:   The matrix composed of basis vectors.
        l:       The lower bound target vector.
        r:       The upper bound target vector.
        weight:  For scaling up vectors.
        verbose: To display more information.
"""
def CVP_rkm0959(basis, l, r, weight=None, verbose=False):
    nvec = basis.nrows()
    dims = basis.ncols()


    # Sanity checks
    if len(l) != dims:
        raise ValueError(f'[ ! ] The dimension of the lower-bound vector != dimension of basis vectors: {len(l)} != {dims}')
    if len(r) != dims:
        raise ValueError(f'[ ! ] The dimension of the upper-bound vector != dimension of basis vectors: {len(l)} != {dims}')
    for d in range(dims):
        if l[d] > r[d]:
            raise ValueError(f'[ ! ] Lower-bound vector should have all the elements smaller than upper-bound vectors!\nl={l}\nr={r}')


    # Estimate the number of possible solutions (verbose=True).
    # However, only calculable if the basis matrix is square one.
    if verbose:
        if nvec == dims:
            DET     = abs(basis.det())
            num_sol = 1
            for i in range(dims):
                num_sol *= (r[i] - l[i])

            if DET == 0:
                print(f"[ i ] Could not estimate the number of solutions, since the determinant of basis matrix is 0.")
            else:
                num_sol //= DET
                # + 1 added in for the sake of not making it zero...
                print(f"[ i ] Expected number of solutions: {num_sol + 1}")
        else:
            print(f"[ i ] Could not estimate the number of solutions, since the size of basis matrix is {nvec}x{dims}, not square one.")


    # Set weight... Why?
    maxElement = max([max([abs(basis[i, j]) for i in range(nvec)]) for j in range(dims)])
    if not weight:
        weight = dims * maxElement


    # Scaling the vectors (to form some sort of high-dimensional sphere?)
    maxDiff        = max([r[d] - l[d] for d in range(dims)])
    appliedWeights = []


    # Scaling some of the dimensions of the bound vectors & the matrix.
    for d in range(dims):
        scaleFactor = weight if l[d] == r[d] else maxDiff // (r[d] - l[d])
        appliedWeights.append(scaleFactor)
        for i in range(nvec):
            basis[i, d] *= scaleFactor
        l[d] *= scaleFactor
        r[d] *= scaleFactor


    # Find closest vector to middle vector of (l, r)
    # = close to (l, r)
    m      = vector([(l[d] + r[d]) // 2 for d in range(dims)])
    result = BabaiCVP(basis, m)


    # Sanity checking after this...
    for d in range(dims):
        if not (l[d] <= result[d] <= r[d]):
            if verbose:
                print('[ i ] The result target is not in between lower-bound & upper-bound vector!')


    # Recover input.
    for d in range(dims):
        result[d] /= appliedWeights[d]
        l[d]      /= appliedWeights[d]
        r[d]      /= appliedWeights[d]
        for i in range(nvec):
            basis[i, d] /= appliedWeights[d]

    return result

"""
    /function/ CVP():

    "" Purpose:
        Approximate a lattice vector in [l, r].

    "" Args:
        basis:              The matrix composed of basis vectors.
        l:                  The lower bound target vector.
        r:                  The upper bound target vector.

        weight=None:        For scaling up vectors.

        verbose=False:      Print more outputs.
        
"""
def CVP(basis, l, r, weight=None, verbose=False):
    if verbose:
        print(f'================= DEBUG CVP =======================')

        print(f'[ i ] Target range vectors: {l} -> {r}')
        print(f'[ i ] Basis:')
        print(basis)
        print()

    result = CVP_rkm0959(basis, l, r, weight, verbose)

    if verbose:
        def _dist(v, w):
            return float(sqrt(sum([(a-b)**2 for a, b in zip(v, w)])))

        print(f'[ i ] Result vector: {result if result else "-NaV-"}')
        print(f'[ i ] Lattice coordinate: {basis.transpose().inverse() * result}')
        print(f'[ i ] Distance:')
        print(f'      [ + ] Result -> Lower: {_dist(l, result)}')
        print(f'      [ + ] Result -> Upper: {_dist(r, result)}')
        print()

    return result
