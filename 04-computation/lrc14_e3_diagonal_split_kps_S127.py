# kind-pasteur-2026-07-09-S127
# The E3 diagonal/off-diagonal split as the lens on the endgame's final rung.
#
# E3(S) = #{(a,b) in S^2 : a+b in S}  (ordered, incl. diagonal)
#       = D2(S)           [DIAGONAL  (a,a): 2a in S -- the DOUBLINGS]
#       + 2*T(S)          [OFF-DIAG {a,b} a!=b: a+b in S -- the SCHUR TRIPLES]
#
# THM-682(d) (monad): only DOUBLINGS carry the exact-load W0 (line wt 0.0382 each,
# W0>0.08 needs >=3); Schur triples are nearly free (0.0027). So the endgame's final
# rung lives on the DIAGONAL of E3 -- the 2-adic axis -- not the whole E3.
# This script confirms the split identity and the doubling<->2-adic-content correlation.
from itertools import combinations
from math import comb, gcd
from functools import reduce

def E3(S):
    Sset = set(S); return sum(1 for a in S for b in S if a + b in Sset)
def D2(S):  # doublings = diagonal of E3
    Sset = set(S); return sum(1 for a in S if 2 * a in Sset)
def T(S):   # unordered off-diagonal Schur triples
    Sset = set(S); return sum(1 for a, b in combinations(S, 2) if a + b in Sset)
def v2(n):  # 2-adic valuation
    k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k

print("=== (1) split identity  E3 = D2 + 2T  over all k-sets ===")
for k in [4, 5, 6, 7]:
    N = 2 * k + 3; n = ok = 0
    for S in combinations(range(1, N + 1), k):
        S = list(S); n += 1
        if E3(S) == D2(S) + 2 * T(S): ok += 1
    print(f"  k={k}: {'OK ' if ok==n else 'FAIL'} {ok}/{n}")

print("\n=== (2) doublings track 2-adic content; Schur triples track AP-ness ===")
# dilated interval {d,..,kd}: D2 = floor(k/2), T = rest; AP (not dilated) {a,a+d,..}: fewer doublings
def is_dilated(S):
    S = sorted(S); d = S[0]; return all(S[i] == (i + 1) * d for i in range(len(S)))
for k in [13]:
    dil = list(range(1, k + 1))                      # {1,..,13} dilated (d=1)
    ap  = list(range(3, 3 + 2 * k, 2))               # {3,5,..,27} AP diff 2, not dilated
    generic = [1,3,4,9,10,12,25,27,30,31,40,44,50]   # a scattered set
    for name, S in [("dilated {1..13}", dil), ("AP diff2 {3..27}", ap), ("generic", generic)]:
        print(f"  {name:20s}: E3={E3(S):3d}  D2(doublings)={D2(S):2d}  T(schur)={T(S):3d}  sum2adic={sum(v2(x) for x in S)}")

print("\n=== (3) doubling-rich <=> even-rich (the W0>0.08 corner is 2-adic) ===")
# sample sets, correlate D2>=3 with evenness / 2-adic sum
import random
random.seed(1)  # note: seeded so this is reproducible; Math.random-free
buckets = {}
for _ in range(200000):
    S = sorted(random.sample(range(1, 60), 13))
    if reduce(gcd, S) != 1: continue
    d = D2(S)
    ev = sum(1 for x in S if x % 2 == 0)
    buckets.setdefault(d, []).append(ev)
print("  D2 (doublings) -> mean #even speeds (of 13), sample count:")
for d in sorted(buckets):
    evs = buckets[d]
    if len(evs) >= 20:
        print(f"    D2={d:2d}: mean_even={sum(evs)/len(evs):5.2f}  n={len(evs)}")

print("\n=== (4) the doubling GRAPH is a forest of 2-adic chains (geometric ratio 2) ===")
# doublings a->2a form chains a,2a,4a,..; D2 = # edges = |S| - (#chains). Max doublings from long chains.
def doubling_chain_edges(S):
    Sset = set(S); return sum(1 for a in S if 2 * a in Sset)  # = D2
def num_chains(S):  # connected components under a~2a
    Sset = set(S); roots = sum(1 for a in S if a % 2 == 1 or (a // 2) not in Sset)
    return roots
for name, S in [("dilated {1..13}", list(range(1,14))),
                ("geom-rich [22,39,42,44,45,84,88,90,168,176,180,336,352]",
                 [22,39,42,44,45,84,88,90,168,176,180,336,352])]:
    print(f"  {name}: D2={D2(S)}  #chains={num_chains(S)}  (edges = 13 - chains = {13-num_chains(S)})")
