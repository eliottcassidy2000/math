# AUDIT task 3, follow-up to aud_task3.out:
#  (B') scale of the odd-#Shat entries of the reduced weight tensors W2(z) (are they zero core by core?) and of the
#       odd-#Shat homomorphism numbers Gamma2.
#  (E)  discriminating power: exact INTEGER injective sums  I(G) = sum over injective f: V(G) -> [23] of prod_e S[f(t_e), f(h_e)]
#       for 6-vertex patterns, on D23_1 and D23_2 (identical 5-census):
#         - 9-edge cubic simple graphs K_33 and prism (odd #S = 9, the labelling killed by parity in step 2);
#         - 10-edge 6-vertex patterns with all degrees >= 3 (even #S = 10: the excess-4 case of the 'Hierarchy');
#       and the same on P_23 and on a random REGULAR (not doubly regular) 23-tournament, to see what is DRT-specific.
import sys, json, itertools
from math import factorial
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE); sys.path.insert(0, HERE + '/audit2')
import loop_engine4 as LE4
from drt_switch import code_to_S
from common import rot_A, skew

cls = json.load(open(f'{HERE}/drt23_classes.json'))
n = 23
cores = LE4.load_cores(3)
six = [i for i, (V, e, a) in enumerate(cores) if V == 6]
print("(B') reduced weight tensors at z = 1: max |W2[t]| over odd-#Shat t vs over all t; loops/multi-edges of the core")
C1 = LE4.drt_coeffs(n, 1.0)
for i in six:
    V, edges, aut = cores[i]
    W2 = LE4.reduce_W(LE4.weight_tensor(3, i), C1)
    par = np.indices((2,) * 9).reshape(9, -1).sum(axis=0).reshape((2,) * 9) % 2
    loops = sum(1 for a, b in edges if a == b); multi = len(edges) - len(set(map(tuple, map(sorted, edges))))
    print(f"   core {i}: loops {loops}, repeated edges {multi}: max|W2 odd| = {np.abs(W2[par == 1]).max():.1e}, max|W2| = {np.abs(W2).max():.3e}")
S0 = code_to_S(cls[0]['code']); P = np.eye(n) - np.ones((n, n)) / n
G2 = {i: LE4.network(3, i).contract(np.array([P, S0 / n])) for i in six}
par = np.indices((2,) * 9).reshape(9, -1).sum(axis=0).reshape((2,) * 9) % 2
print("   Gamma2 (D23_0, 6-vertex cores): max |Gamma2 odd| / max |Gamma2| = "
      + ', '.join(f"{np.abs(G2[i][par == 1]).max() / np.abs(G2[i]).max():.1e}" for i in six))

def set_partitions(items):
    if not items: yield []; return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for j in range(len(part)):
            yield part[:j] + [[first] + part[j]] + part[j + 1:]
        yield [[first]] + part
LET = 'abcdefgh'
def hom_int(V, edges, S):
    """sum over all maps f: [V] -> [n] of prod_e S[f(t), f(h)] (exact: int64 / object)"""
    if any(a == b for a, b in edges): return 0          # S has zero diagonal
    spec = ','.join(LET[a] + LET[b] for a, b in edges) + '->'
    return int(np.einsum(spec, *([S.astype(np.int64)] * len(edges)), optimize='greedy'))
def inj_int(V, edges, S):
    tot = 0
    for part in set_partitions(list(range(V))):
        mu = 1
        for b in part: mu *= (-1) ** (len(b) - 1) * factorial(len(b) - 1)
        lab = {v: bi for bi, b in enumerate(part) for v in b}
        tot += mu * hom_int(len(part), [(lab[a], lab[b]) for a, b in edges], S)
    return tot

K33 = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
PRISM = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (0, 3), (1, 4), (2, 5)]
TEN = {'K33+(0,1)': K33 + [(0, 1)], 'prism+(0,4)': PRISM + [(0, 4)], 'K33+(0,1) reoriented': [(3, 0), (0, 4), (5, 0), (1, 3), (4, 1), (1, 5), (3, 2), (2, 4), (5, 2), (1, 0)]}
rng = np.random.default_rng(99)
A = rot_A(23).copy()
for _ in range(20000):
    u, v, w = rng.choice(23, 3, replace=False)
    if A[u, v] and A[v, w] and A[w, u]: A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
Sreg = skew(A); assert not Sreg.sum(axis=1).any()
print(f"   random regular 23-tournament (3-cycle chain, seed 99): doubly regular? {np.array_equal(Sreg @ Sreg.T, 23 * np.eye(23, dtype=int) - 1)}")
tests = [('D23_1', code_to_S(cls[1]['code'])), ('D23_2', code_to_S(cls[2]['code'])), ('P_23', code_to_S(cls[12]['code'])),
         ('regular non-DRT', Sreg), ('its converse', -Sreg)]
print("\n(E) exact integer injective sums I(G) (divide by 23^#edges for the Shat version):")
for gname, edges in [('K33 (9 edges)', K33), ('prism (9 edges)', PRISM)] + list(TEN.items()):
    vals = [(name, inj_int(6, edges, S)) for name, S in tests]
    print(f"   {gname:22s}: " + ', '.join(f"{name} {v}" for name, v in vals), flush=True)
