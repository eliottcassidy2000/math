"""grassmann_chart_average.py -- mac-mini-2026-10-27-potentials

Generalisation of the concurrent three-chart identity (dyadic-module-spectrum R1: four child bins, three pairings,
mean intermediate energy (E + 2S)/3).  Children are indexed by F_2^r with charges z_x in a real inner-product space
(here the Eisenstein form Q(a, b) = a^2 - ab + b^2 on Z^2); for a k-dimensional subspace W the intermediate partition
is the set of cosets of W and M_W = sum over cosets of Q(sum of the charges in the coset).  Claim (exact, rational):
    mean over all k-dimensional W of M_W = S + theta_k (E - S),   theta_k = (2^k - 1)/(2^r - 1),
with E = Q(sum of all charges), S = sum Q(z_x).  A pair x != y shares a coset iff x + y lies in W, which happens for the
fraction theta_k of the k-subspaces.  Checked exhaustively over all subspaces for r <= 4 on random integer charges."""
import itertools, random
from fractions import Fraction

def Q(v):
    a, b = v
    return a * a - a * b + b * b

def subspaces(r, k):
    vecs = list(range(1, 1 << r))
    seen = set()
    for basis in itertools.combinations(vecs, k):
        span = {0}
        for b in basis:
            span |= {s ^ b for s in span}
        if len(span) == (1 << k):
            seen.add(frozenset(span))
    return seen

rng = random.Random(27)
bad = 0
checks = 0
for r in range(1, 5):
    for k in range(0, r + 1):
        subs = subspaces(r, k)
        for trial in range(20):
            z = [(rng.randrange(-5, 6), rng.randrange(-5, 6)) for _ in range(1 << r)]
            E = Q((sum(a for a, _ in z), sum(b for _, b in z)))
            S = sum(Q(v) for v in z)
            tot = Fraction(0)
            for W in subs:
                cosets = {}
                for x in range(1 << r):
                    rep = min(x ^ w for w in W)
                    cosets.setdefault(rep, []).append(x)
                M = sum(Q((sum(z[x][0] for x in c), sum(z[x][1] for x in c))) for c in cosets.values())
                tot += M
            mean = tot / len(subs)
            theta = Fraction((1 << k) - 1, (1 << r) - 1)
            checks += 1
            if mean != S + theta * (E - S):
                bad += 1
print(f"Grassmann chart averages: {checks} (r, k, charge) checks for r <= 4, all subspaces: failures {bad}")
print("r = 2, k = 1: theta = 1/3, mean = (E + 2S)/3 (the concurrent three-chart identity)")
