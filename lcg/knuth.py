import random 
import math 
from sage.all import *
from itertools import combinations
class LCG:
    m = 10734367385013619889
    a = 9807963723765715717
    b = 7226300108115682840
    def __init__(self,seed):
        self.seed = seed
    def next(self):
        self.seed = (self.seed*self.a + self.b) % self.m 
        return self.seed 
seed = 2877244225168654778
lcg = LCG(seed)
m = LCG.m 
a = LCG.a 
b = LCG.b
# print(m.bit_length())
# generate truncated xs 
xs = []
ys = []
zs = []
mod = 2**32
alpha = 1/2
k = 64
s = 32
for _ in range(17):
    x = lcg.next()
    xs.append(x)
    ys.append(x >> 32)
    zs.append(x % mod) 
assert all(y*(2**32)+z==x for x,y,z in zip(xs,ys,zs))
# print(ys) # được biết ys
t = 3
n = 14
Vs = []
for i in range(n):
    V = []
    for j in range(t):
        V.append(ys[i+j+1]-ys[i+j])
    Vs.append(V)
K = 356131

Ws = []
for i in range(n):
    W = [] 
    for j in range(t):
        W.append(xs[i+j+1]-xs[i+j])
    Ws.append(W)

Vs = [vector(ZZ, v) for v in Vs]   
B = matrix(ZZ, [K * v for v in Vs])  
I = identity_matrix(ZZ, n)
M = B.augment(I)
L = M.LLL()
Ws = [vector(ZZ, w) for w in Ws]
print(Ws)
for v in L:
    v_ = v[3:]
    s_  = 0 
    for i in range(len(Ws)):
        s_ += v_[i] * Ws[i]
    print(s_)
        



for r in L:
    print(r)
bound = 2**(k*(1-1/t))
b_norm = bound/(sqrt(n*t)*2**((1-alpha)*k))
print(b_norm)
# build polynomial 
R = PolynomialRing(ZZ,'x')
x = R.gen()
lamb = []
ps = []
for vec in L:
    P = sum(vec[3:][i]*x**i for i in range(n))
    print(P)
    print(P(a)%m) # for testing
    ps.append(P)
ms = []
# recover m 
for comb in combinations(ps,3):
    p0 = comb[0]
    p1 = comb[1]
    p2 = comb[2]
    m_ = math.gcd(p0.resultant(p1), p1.resultant(p2), p0.resultant(p2))
    print(m_)
    if (m_.bit_length() > 20):
        ms.append(m_)
ms = list(set(ms))
print(ms)


def HGCD(a, b):
    if 2 * b.degree() <= a.degree() or a.degree() == 1:
        return 1, 0, 0, 1
    m = a.degree() // 2
    a_top, a_bot = a.quo_rem(x**m)
    b_top, b_bot = b.quo_rem(x**m)
    R00, R01, R10, R11 = HGCD(a_top, b_top)
    c = R00 * a + R01 * b
    d = R10 * a + R11 * b
    q, e = c.quo_rem(d)
    d_top, d_bot = d.quo_rem(x**(m // 2))
    e_top, e_bot = e.quo_rem(x**(m // 2))
    S00, S01, S10, S11 = HGCD(d_top, e_top)
    RET00 = S01 * R00 + (S00 - q * S01) * R10
    RET01 = S01 * R01 + (S00 - q * S01) * R11
    RET10 = S11 * R00 + (S10 - q * S11) * R10
    RET11 = S11 * R01 + (S10 - q * S11) * R11
    return RET00, RET01, RET10, RET11
    
def GCD(a, b):
    print(a.degree(), b.degree())
    q, r = a.quo_rem(b)
    if r == 0:
        return b
    R00, R01, R10, R11 = HGCD(a, b)
    c = R00 * a + R01 * b
    d = R10 * a + R11 * b
    if d == 0:
        return c.monic()
    q, r = c.quo_rem(d)
    if r == 0:
        return d
    return GCD(d, r)

def gcd_zmod(f, g):
    while g:
        f, g = g, f % g
    return f
def polynomial_gcd_crt(a, b, factors):
    assert a.base_ring() == b.base_ring() == ZZ
    gs = []
    ps = []
    for p, _ in factors:
        zmodp = Zmod(p)
        gs.append(GCD(a.change_ring(zmodp), b.change_ring(zmodp)).change_ring(ZZ))
        ps.append(p)

    g, _ = fast_crt(gs, ps)
    return g
def fast_crt(X, M, segment_size=8):
    assert len(X) == len(M)
    assert len(X) > 0
    while len(X) > 1:
        X_ = []
        M_ = []
        for i in range(0, len(X), segment_size):
            if i == len(X) - 1:
                X_.append(X[i])
                M_.append(M[i])
            else:
                X_.append(crt(X[i:i + segment_size], M[i:i + segment_size]))
                M_.append(lcm(*M[i:i + segment_size]))
        X = X_
        M = M_

    return X[0], M[0]
# recover a 
a_lst = []

for comb in combinations(ps,3):
    p0 = comb[0]
    p1 = comb[1]
    p2 = comb[2]
    for m_ in ms: 
        fac = factor(m_)
        try:
            g = polynomial_gcd_crt(p0,polynomial_gcd_crt(p1,p2,fac),fac)
            for a_ in g.change_ring(Zmod(m_)).roots(multiplicities = False):
                    a_lst.append(a_)
                    print(a_)
        except Exception as e:
            continue 
print(set(a_lst))



