"""Units mod 1800 = 8*9*25: which small tournament counts generate (Z/1800)^x ?  (A set generates (Z/2^a 3^b 5^c)^x
for all a, b, c iff it generates (Z/1800)^x: the kernel of reduction is inside the Frattini subgroup.)"""
from math import gcd
M = 1800
units = [x for x in range(M) if gcd(x, M) == 1]
def closure(gens):
    S = {1}; frontier = [1]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = x * g % M
                if y not in S: S.add(y); new.append(y)
        frontier = new
    return S
print("|units mod 1800| =", len(units))
counts = [x for x in range(1, 200, 2) if x not in (7, 21) and gcd(x, 30) == 1]
best = None
import itertools
for k in range(1, 5):
    for c in itertools.combinations(counts[:20], k):
        if len(closure(c)) == len(units):
            best = c; break
    if best: break
print("smallest generating set of counts coprime to 30:", best)
# Frattini check: kernel of (Z/1800)^x -> (Z/8)^x x (Z/9)^x x (Z/25)^x is trivial here (1800 = 8*9*25 exactly)
