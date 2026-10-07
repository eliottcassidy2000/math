# Independent decomposition analysis for L = Q(i, sqrt p : p | m), via Frobenius/inertia on the
# F_2-vector space of quadratic characters.  For a rational prime q:
#   L^I = compositum of quadratic subfields unramified at q;  L^D = compositum of those where q splits.
#   c in D <=> c fixes L^D <=> every quadratic subfield in which q splits is real.
#   c in I <=> every quadratic subfield unramified at q is real.
from itertools import product
from sympy import factorint, primerange, isprime
def sqf(n):
    s = -1 if n < 0 else 1
    for p, e in factorint(abs(n)).items():
        if e % 2: s *= p
    return s
def kron_quad(e, q):   # splitting of q in Q(sqrt e), e squarefree != 1: 'R','S','I'
    D = e if e % 4 == 1 else 4*e
    if D % q == 0: return 'R'
    if q == 2:
        return 'S' if e % 8 == 1 else 'I'
    return 'S' if pow(e % q, (q-1)//2, q) == 1 else 'I'
def analyse(m, q):
    gens = [-1] + (list(factorint(m)) if m > 1 else [])
    subs = []
    for S in product((0,1), repeat=len(gens)):
        if not any(S): continue
        e = 1
        for s, g in zip(S, gens):
            if s: e *= g
        e = sqf(e)
        subs.append((e, kron_quad(e, q)))
    cD = all(e > 0 for e, b in subs if b == 'S')
    cI = all(e > 0 for e, b in subs if b != 'R')
    if not cD: return None
    return 'ram' if cI else 'inert'
ms = [1,5,13,21,33,37,57,85,93,105,133,165,177,253,273,345,357,385,1365]
for m in ms:
    cs = []
    for q in primerange(2, 200):
        a = analyse(m, q)
        if a: cs.append(f"{q}:{a}")
    print(f"m={m:5d}  c-stable primes q<200: {cs}")
