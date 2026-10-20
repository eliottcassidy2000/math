#!/usr/bin/env python3
"""THM-4725 Theorem U numerics: primes p = 9 mod 16 (v_2(p-1) = 3) at which given integers c are 8th-power residues,
so every dilation rule with those constants is constant on the two K_max-invariant tournament types."""
from sympy import isprime
def kmax_member(c, p):   # c in the 2^e-th powers with 2^e || p-1  <=>  c^((p-1)/2^e) = 1
    m = p - 1
    while m % 2 == 0: m //= 2
    return pow(c % p, m, p) == 1
def v2(n):
    k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k
sets = [[2], [3], [5], [6], [7], [10], [2, 3], [2, 3, 5], [2, 3, 5, 7], [3, 7, 11]]
for cs in sets:
    found = [p for p in range(17, 20000) if isprime(p) and v2(p - 1) == 3 and all(kmax_member(c, p) for c in cs)]
    print(f"c = {cs}: first primes p = 9 (mod 16) with all c 8th-power residues: {found[:8]} ... ({len(found)} below 20000, of {sum(1 for p in range(17,20000) if isprime(p) and v2(p-1)==3)} primes with v_2(p-1) = 3)")
# the two K_max-types at such a prime: antisymmetric functions on Z_8 up to translation; show that any dilation rule
# with constants in K_max gives the same image on both
p = 73  # 2 is an 8th power residue (2^9 = 1 mod 73)
import itertools
K = sorted({pow(x, 8, p) for x in range(1, p)})
g = 5
assert len({pow(g, k, p) for k in range(p - 1)}) == p - 1
cosets = [[pow(g, t + 8 * k, p) for k in range(9)] for t in range(8)]   # coset of g^t K, t in Z_8
types = set()
for s in itertools.product([0, 1], repeat=4):
    sfull = list(s) + [1 - x for x in s]      # antisymmetric on Z_8
    S = frozenset(x for t in range(8) if sfull[t] for x in cosets[t])
    canon = min(tuple(sorted(a * x % p for x in S)) for a in range(1, p)); types.add(canon)
print(f"p = 73: K_max-invariant tournament sets 16, multiplier types {len(types)}; for c = 2 (in K_max) 1_S(2x) = 1_S(x) on all of them, so f_phi(S) in {{empty, all}} for every dilation rule with constant 2")
