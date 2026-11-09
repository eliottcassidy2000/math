#!/usr/bin/env python3
"""lines_formula.py n -- the component-class line counts of the order-n metagraph from the cut-set formula
(THM-4953 (A)), to be compared with lines_census.c:
  N_sup(S)  = number of closed tilings t whose open end T_{t^c} has cut set containing S
            = 2^(m-1) if S is empty, else 2^(e(S)), e(S) = sum over the blocks B of S of C(|B|-1, 2);
  N(S)      = sum_{S' superset S} (-1)^|S' - S| N_sup(S')           (exact cut set S);
  blue:  N_sup^blue(S) = 2^((m+f)/2 - 1) if S empty, else 2^((e(Sbar) + f(Sbar))/2) with Sbar = S u sigma(S),
         f(Sbar) = number of sigma-fixed tiles inside the blocks of Sbar (tiles (a,b) with a + b = n + 1).
Prints, per k = |S| + 1 (number of strong components of the open end), the all-lines and blue-lines counts, and the
exact cut-set table in lines_census.c's encoding (position i -> bit i-1)."""
import sys, itertools
from math import comb
n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
m = comb(n - 1, 2); f = (n - 1) // 2
positions = list(range(1, n))
def blocks(S):
    cuts = sorted(S); out = []; start = 1
    for c in cuts: out.append((start, c)); start = c + 1
    out.append((start, n)); return out
def e_of(S): return sum(comb(b - a + 1 - 1, 2) for a, b in blocks(S))   # C(|B|-1, 2)
def f_of(S): return sum(1 for a, b in blocks(S) for x in range(a, b + 1) for y in range(a, b + 1) if x >= y + 2 and x + y == n + 1)
def sigma(S): return frozenset(n - i for i in S)
def Nsup(S): return 2 ** (m - 1) if not S else 2 ** e_of(S)
def Nsup_blue(S):
    Sb = S | sigma(S)
    return 2 ** ((m + f) // 2 - 1) if not Sb else 2 ** ((e_of(Sb) + f_of(Sb)) // 2)
allS = [frozenset(c) for r in range(0, n) for c in itertools.combinations(positions, r)]
sup = {S: Nsup(S) for S in allS}; supb = {S: Nsup_blue(S) for S in allS}
def exact(table, S): return sum((-1) ** (len(S2) - len(S)) * table[S2] for S2 in allS if S <= S2)
N = {S: exact(sup, S) for S in allS}; Nb = {S: exact(supb, S) for S in allS}
assert sum(N.values()) == 2 ** (m - 1) and sum(Nb.values()) == 2 ** ((m + f) // 2 - 1)
assert all(v >= 0 for v in N.values()) and all(v >= 0 for v in Nb.values())
print(f"# n={n} m={m} f={f}: lines {2**(m-1)}, blue {2**((m+f)//2-1)}")
print("#k")
for k in range(1, n + 1):
    print(k, sum(v for S, v in N.items() if len(S) == k - 1), sum(v for S, v in Nb.items() if len(S) == k - 1))
print("#cuts")
for S in sorted(allS, key=lambda S: sum(1 << (i - 1) for i in S)):
    code = sum(1 << (i - 1) for i in S)
    if N[S]: print(code, N[S], Nb[S])
