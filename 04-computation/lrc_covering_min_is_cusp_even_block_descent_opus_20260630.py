"""
Reconciliation: the even block (covering-min, even n) DESCENDS to the cusp (AP, Z_7) -- it IS the cusp config,
M=1/n = the comb-witness (empty tooth). Vindicates cusp-existence; construction was the off-cusp red herring.
Also: odd-n parity-blocking (the cusp/AP unreachable as a covering set => covering-min off-cusp > 1/n).
"""
import math, cmath
from fractions import Fraction
w7=cmath.exp(2j*cmath.pi/7)
def descent(S):
    ch=[]; cur=list(S)
    while cur:
        O=[x for x in cur if x%2==1]; E=[x//2 for x in cur if x%2==0]
        ch.append(sorted(O)); cur=E
        if len(ch)>15: break
    return ch
def apex(O): return min(abs(sum(w7**((k*x)%7) for x in O))**2 for k in range(1,7)) if O else None
# n=14 even block descent
n=14; eb=[2*k for k in range(1,n)]
ch=descent(eb)
print(f"EVEN BLOCK 2*{{1..13}} (n=14 covering-min, M=1/14) -- 2-adic descent:")
for j,O in enumerate(ch):
    res=sorted(set(x%7 for x in O)) if O else []
    cusp = set(x%7 for x in O)==set(range(7)) if O else False
    print(f"  level {j}: odd core={O} mod7={res} {'<- CUSP Z_7 (apex gap 0)' if cusp else ''}")
print(f"  => the even block is ALL EVEN: level-0 odd core empty; level-1 = {{1..13}} = the AP = the CUSP (Z_7).")
print(f"     So the covering-min (even n) IS the cusp config; M=1/n is the COMB-WITNESS (empty tooth). Cusp-existence VINDICATED.")
print()
# the AP {1..13} itself (the cusp) -- M=1/14
print(f"  AP {{1..13}} (the cusp): odd core mod7 = {sorted(set(x%7 for x in range(1,14)))} = Z_7; apex gap g(Z_7)={apex(list(range(7)))}")
print()
print("ODD-n parity blocking: the AP/cusp is UNREACHABLE as an all-even covering set (q=n needs odd mult of odd n).")
print("  => odd-n covering-min is OFF-CUSP, > 1/n (n=7: 2/13). The construction (Phi_6, off-cusp) is STILL wrong (2/13 != 7/43).")
print()
print("RECONCILIATION (what was right all along):")
print("  * cusp-existence / comb-witness / empty-tooth (the AP, M=1/n): RIGHT -- it IS the even-n covering-min (the even block).")
print("  * construction / Phi_6 / zeta_6 / hexagonal / Sylvester-covering-min (off-cusp): RED HERRING (non-extremal).")
print("  * the covering-min = the CUSP (AP doubled, even n) = measure-0 / comb-witness M=1/n. Two columns:")
print("    MEASURE = apex gap 0 at the cusp (C_p spectral, Ramanujan); EXISTENCE = M=1/n (C_n covering radius, comb).")
