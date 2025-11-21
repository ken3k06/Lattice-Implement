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
    M = flatter(mat) # or .LLL() if you like 
    G = M.gram_schmidt()[0]
    diff = target 
    for i in reversed(range(G.nrows())):
        diff -= M[i] * ((diff*G[i]) / (G[i]*G[i])).round()
    return target - diff 
