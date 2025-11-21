from sage.all import * 
import random 
q = 37
assert is_prime(q)
F = GF(q)
n = 5
m = n**2
A = Matrix(GF(q),m,n,lambda i,j: randint(0,q-1))
e = vector((0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0),GF(q))
V = VectorSpace(GF(q),n)
S = V.random_element()
print(f"secret key={S}")
b = (A*S) + e
def flatter(M):
    import re
    from subprocess import check_output
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), list(map(int, re.findall(b"-?\\d+", ret))))
def primal_attack(A,b,fll = False):
    m = A.nrows()
    n = A.ncols()
    A_ = matrix(ZZ, A.T.rref().submatrix(0,n,n,m-n))
    b_zz = [int(x) for x in b]
    I_n = identity_matrix(ZZ,n)
    qAry = q*identity_matrix(ZZ,m-n)
    zero_mn_n = zero_matrix(ZZ, m-n, n)    
    zero_n_1  = zero_matrix(ZZ, n, 1)  
    zero_mn_1 = zero_matrix(ZZ, m-n, 1)    
    c1 = matrix(ZZ, 1, n,    b_zz[:n])    
    c2 = matrix(ZZ, 1, m-n,  b_zz[n:]) 
    t  = matrix(ZZ, 1, 1, [1])
    B = block_matrix([
    [I_n,        A_,   zero_n_1],
    [zero_mn_n,  qAry,   zero_mn_1],
    [c1,c2,t]
])
    if fll == True:
        L = flatter(B)
    else:
        L = B.LLL()
    return L 
L = primal_attack(A,b, fll = True)
print(L[0][:-1] == e)
S_ = A.solve_right(b-e)
print(f"recover secret key = {S_}")
