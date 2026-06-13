#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S706b
Probe: WHEN is the LRC worry-set field Q(zeta_C), C=2n-1, "rosette-maximal" in
its own dimension N=phi(C)? I.e. when is its root-of-unity count w=2C equal to
the MAXIMAL achievable count M(N)=max{m:phi(m)|N} in dimension N (THM-416)?

Equivalently: is 2C = max{ m : phi(m) | phi(C) } ?

Spot-checks suggested: YES iff 3 | C. And 3|C=2n-1 is EXACTLY the THM-413
order-3 silent-flip / most-degenerate-sign-orbit condition. If the equivalence
holds, THM-416 (totient density-quantum ladder) and THM-413 (order-3 degeneracy)
pick out the SAME arithmetic locus on the LRC modulus -- the prime-3 thread.

Test over all odd C in [3, CMAX], exact arithmetic.
"""
from math import gcd

CAP = 20000                       # sieve bound; M(N) for N<=54 needs m well below this
_phi = list(range(CAP+1))
for p in range(2, CAP+1):
    if _phi[p] == p:              # p prime
        for k in range(p, CAP+1, p):
            _phi[k] -= _phi[k] // p
def euler_phi(m):
    return _phi[m]

def M_of(N):
    """max{m : phi(m) | N}, scanned over the sieve (CAP >> any needed m)."""
    best = 1
    for m in range(1, CAP+1):
        if _phi[m] != 0 and N % _phi[m] == 0:
            best = m
    return best

CMAX = 81
print("Probe: is Q(zeta_C) rosette-maximal in dim phi(C)  <=>  3 | C ?")
print(" C  | n=(C+1)/2 | phi(C) | 2C | M(phi(C)) | maximal? | 3|C ? | match")
mismatches = 0
for C in range(3, CMAX+1, 2):
    n = (C+1)//2
    N = euler_phi(C)
    w = 2*C                      # #roots of unity in Q(zeta_C), C odd
    MN = M_of(N)
    maximal = (w == MN)
    div3 = (C % 3 == 0)
    match = (maximal == div3)
    if not match:
        mismatches += 1
    flag = "" if match else "  <<< MISMATCH"
    print(f" {C:2d} |    {n:3d}    |  {N:3d}   | {w:3d}|   {MN:4d}    | "
          f"{'YES' if maximal else 'no ':3s}     | {'YES' if div3 else 'no ':3s}  | "
          f"{'ok' if match else 'NO'}{flag}")

print()
print(f"mismatches (maximal != (3|C)): {mismatches}")
print("CONJECTURE HOLDS on this range" if mismatches == 0 else "CONJECTURE FAILS")

# Why: heuristic. M(N) is realized by m a product of prime powers q with phi(q)|N.
# Adjoining the prime 3 (phi(3)=2) is "cheap" (2 | phi(C) always, since C odd>1
# => phi(C) even). For C=2n-1 odd: phi(C) even. If 3 not already dividing C, then
# 3*C has phi(3C)=2*phi(C) which may EXCEED phi(C)'s divisibility... the precise
# mechanism is below.
print()
print("Mechanism check: is M(phi(C)) typically 2*C*(extra 3-factor) when 3 does NOT divide C?")
for C in [5,7,11,13,25,35]:
    N = euler_phi(C); MN = M_of(N)
    print(f"  C={C:2d}: phi={N:2d}, 2C={2*C:3d}, M(phi)={MN:3d}, ratio M/(2C)={MN/(2*C):.3f}, M/C={MN/C:.3f}")
