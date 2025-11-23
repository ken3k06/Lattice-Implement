from sage.all import *
from sage.modules.free_module_integer import IntegerLattice
def flatter(M):
    import re
    from subprocess import check_output
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), list(map(int, re.findall(b"-?\\d+", ret))))
def Babai_CVP(mat,target):
    M = flatter(mat)
    G = M.gram_schmidt()[0]
    diff = target 
    for i in reversed(range(G.nrows())):
        diff -= M[i] * ((diff*G[i]) / (G[i]*G[i])).round()
    return target - diff 
def small_roots(self, X=None, beta=1.0, epsilon=None, **kwds):
    from sage.misc.verbose import verbose
    from sage.matrix.constructor import Matrix
    from sage.rings.real_mpfr import RR
 
    N = self.parent().characteristic()
 
    if not self.is_monic():
        raise ArithmeticError("Polynomial must be monic.")
 
    beta = RR(beta)
    if beta <= 0.0 or beta > 1.0:
        raise ValueError("0.0 < beta <= 1.0 not satisfied.")
 
    f = self.change_ring(ZZ)
 
    P, (x,) = f.parent().objgens()
 
    delta = f.degree()
 
    if epsilon is None:
        epsilon = beta / 8
    verbose("epsilon = %f" % epsilon, level=2)
 
    m = max(beta**2 / (delta * epsilon), 7 * beta / delta).ceil()
    verbose("m = %d" % m, level=2)
 
    t = int((delta * m * (1 / beta - 1)).floor())
    verbose("t = %d" % t, level=2)
 
    if X is None:
        X = (0.5 * N ** (beta**2 / delta - epsilon)).ceil()
    verbose("X = %s" % X, level=2)
 
    # we could do this much faster, but this is a cheap step
    # compared to LLL
    g = [x**j * N ** (m - i) * f**i for i in range(m) for j in range(delta)]
    g.extend([x**i * f**m for i in range(t)])  # h
 
    B = Matrix(ZZ, len(g), delta * m + max(delta, t))
    for i in range(B.nrows()):
        for j in range(g[i].degree() + 1):
            B[i, j] = g[i][j] * X**j
 
    B = flatter(B)
 
    f = sum([ZZ(B[0, i] // X**i) * x**i for i in range(B.ncols())])
    R = f.roots()
 
    ZmodN = self.base_ring()
    roots = set([ZmodN(r) for r, m in R if abs(r) <= X])
    Nbeta = N**beta
    return [root for root in roots if N.gcd(ZZ(self(root))) >= Nbeta]
