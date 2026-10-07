# Lemma R table for the 19 "Hilbert class field planes" F^2, F = Q(sqrt p : p | m), m squarefree idoneal = 1 mod 4.
# L = F(i) multiquadratic; Galois group (Z/2)^r; complex conjugation c negates exactly the generators with negative radicand.
import itertools, math
from sympy import factorint, isprime, legendre_symbol
KAPPA = {2:(4,4),3:(3,3),4:(4,4),5:(4,4),7:(4,4),8:(4,4),9:(3,3),11:(5,5),13:(5,6),16:(4,4),17:(5,6),19:(5,5)}  # (lower, upper) for kappa(q)
def disc(e):          # discriminant of Q(sqrt e), e squarefree != 1
    return e if e % 4 == 1 else 4 * e
def sqfree(n):
    s = 1 if n > 0 else -1
    for p, k in factorint(abs(n)).items():
        if k % 2: s *= p
    return s
def behaviour(e, q):  # 'R','S','I' for q in Q(sqrt e)
    D = disc(e)
    if D % q == 0: return 'R'
    if q == 2: return 'S' if e % 8 == 1 else 'I'
    return 'S' if legendre_symbol(e % q, q) == 1 else 'I'
def bound(gens, q):
    """Lemma R bound at the primes over q for L = Q(sqrt g : g in gens): returns (kind, chi_upper) or None if no c-stable prime."""
    r = len(gens)
    subs = []
    for S in itertools.product((0, 1), repeat=r):
        if not any(S): continue
        e = 1
        for s, g in zip(S, gens):
            if s: e *= g
        e = sqfree(e)
        if e == 1: continue
        subs.append((S, e, behaviour(e, q)))
    # c in D  <=>  every split quadratic subfield is real; c in I <=> every unramified quadratic subfield is real
    cD = all(e > 0 for S, e, b in subs if b == 'S')
    if not cD: return None
    cI = all(e > 0 for S, e, b in subs if b != 'R')
    if cI: return ('ramified', 2 if q == 2 else 3)
    # residue degree f of L at q = order of Frobenius = 2 if some unramified subfield is inert, else 1; c in D\I acts nontrivially
    # on the residue field, so the residue field of L+ is F_(q^(f/2)) with f = 2 here (multiquadratic: f <= 2)
    return ('inert', KAPPA.get(q, (None, None))[1], q)
ms = [1,5,13,21,33,37,57,85,93,105,133,165,177,253,273,345,357,385,1365]
for m in ms:
    ps = [p for p in factorint(m)] if m > 1 else []
    gens = [-1] + ps
    best = None; how = None
    for q in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]:
        b = bound(gens, q)
        if b is None: continue
        ub = b[1]
        if ub is not None and (best is None or ub < best):
            best, how = ub, (q, b[0])
    has = lambda *ds: all(d in ps for d in ds)
    lo = 2
    if has(3) or has(7): lo = 3                               # triangle (sqrt3) or Madore's 9-cycle (sqrt7)
    if has(3, 11) or (has(3) and any(sqfree(12*n-1) in [sqfree(p1*p2) for p1 in ps for p2 in ps] for n in (3,7,9,13))): lo = 4
    if has(3, 5, 11) or has(3, 11, 13, 19): lo = 5           # Polymath field / Exoo-Ismailescu field
    print(f"m={m:5d} F=Q(sqrt {ps})  lower >= {lo}   Lemma R upper <= {best} via q={how}")
