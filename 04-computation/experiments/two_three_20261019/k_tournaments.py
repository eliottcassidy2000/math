"""k-valued tournaments: a function on unordered pairs with values in Z/k, the swap acting by x -> 1 - x (tournaments at k = 2).
Isomorphism classes by Burnside: T_k(n) = sum over cycle types lambda of k^{o_+(lambda)} c^{o_-(lambda)} / z_lambda, where
o_+ = sum_i floor((L_i - 1)/2) + sum_{i<j} gcd(L_i, L_j) (orientation-preserving pair orbits), o_- = number of even cycles
(orientation-flipping orbits), and c = 1 if k is odd (the unique tie state x = 1/2 mod k), c = 0 if k is even (no tie state:
only permutations with all cycles odd contribute, exactly as for tournaments)."""
from fractions import Fraction as F
from math import gcd, factorial
from sympy.utilities.iterables import partitions
import sys, json
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
def zlam(parts):
    z = 1
    from collections import Counter
    for L, m in Counter(parts).items(): z *= L**m * factorial(m)
    return z
def T(k, n):
    c = 1 if k % 2 else 0
    tot = F(0)
    for p in partitions(n):
        parts = [L for L, m in p.items() for _ in range(m)]
        oplus = sum((L-1)//2 for L in parts) + sum(gcd(parts[i], parts[j]) for i in range(len(parts)) for j in range(i+1, len(parts)))
        ominus = sum(1 for L in parts if L % 2 == 0)
        if ominus and c == 0: continue
        tot += F(k**oplus, zlam(parts))
    assert tot.denominator == 1
    return int(tot)
seqs = {}
for k in range(2, 9):
    seqs[k] = [T(k, n) for n in range(1, NMAX+1)]
    print(f"k={k}: {seqs[k]}")
json.dump(seqs, open("k_tournaments.json", "w"))
# sanity: k=2 tournaments A000568, k=3 oriented graphs A001174
assert seqs[2][:8] == [1, 1, 2, 4, 12, 56, 456, 6880], seqs[2][:8]
assert seqs[3][:7] == [1, 2, 7, 42, 582, 21480, 2142288], seqs[3][:7]
print("k=2 = A000568 (tournaments), k=3 = A001174 (oriented graphs) confirmed")
# ratio to the naive count k^C(n,2)/n! (almost all structures rigid)
for k in (2, 3, 4, 5):
    print(f"k={k}: T_k(n) n!/k^C(n,2) for n=8..{NMAX}:", [round(seqs[k][n-1]*factorial(n)/k**(n*(n-1)//2), 6) for n in range(8, NMAX+1)])
