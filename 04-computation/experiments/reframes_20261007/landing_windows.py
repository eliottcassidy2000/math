#!/usr/bin/env python3
"""Is the landing integer N_D of the +1-barrier child stationary in D?  Tail of N_D by D-window, against the stationary
level-recursion constant: b_(m-1) = 1/2 + (3/2) 2^(-v) b_m with v ~ Geom(1/2) gives P(b > x) ~ C/x, C = E[B]/E[A ln A]
= 0.5/0.1744 = 2.87 (Kesten-Goldie, kappa = 1).  Also the landing time ratio s_0/D and the number of levels spent with
b < 0 (the negative phase near -1)."""
import math
def landing(D):
    a, m, s, neg = 2 - 3**D, D, 0, 0
    while m > 0:
        if a & 1: a = (a + 3**(m-1)) // 2; m -= 1
        else: a //= 2
        s += 1
    return a, s
C = 0.5 / (math.log(1.5) - math.log(2) / 3)
print(f"stationary constant C = {C:.3f}")
for lo, hi in ((3, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 4500)):
    Ns = []; ss = []
    for D in range(lo, hi):
        N, s = landing(D); Ns.append(N); ss.append(s / D)
    n = len(Ns)
    row = " ".join(f"P(>{x})={sum(1 for v in Ns if v > x)/n:.4f}[{C/x:.4f}]" for x in (4, 16, 64, 256))
    print(f"D in [{lo},{hi}): n={n} max {max(Ns)} median {sorted(Ns)[n//2]}  {row}  mean s0/D {sum(ss)/n:.3f}", flush=True)
