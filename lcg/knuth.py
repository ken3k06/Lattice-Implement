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
print(m.bit_length())
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
print(ys) # được biết ys
t = 3
n = 14
Vs = []
for i in range(n):
    V = []
    for j in range(t):
        V.append(ys[i+j+1]-ys[i+j])
    Vs.append(V)
K = 356131
Vs = [vector(ZZ, v) for v in Vs]   
B = matrix(ZZ, [K * v for v in Vs])  
I = identity_matrix(ZZ, n)
M = B.augment(I)
L = M.LLL()
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
    if (m_.bit_length() > 2**20):
        ms.append(m_)
    print(m_)
print(ms)
# recover a 
TODO
