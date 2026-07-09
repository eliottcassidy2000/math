"""
mac-mini-2026-07-09-S61 -- decompose the capstone ratio r_N (opus-S165) to test whether it is
a-priori bounded by the density-floor tail rate.

r_N = |Corr_N| / (N (6/7)^k),  Corr_N = S_N - N(6/7)^k,  S_N = sum_{j=1}^N W(j/Vmax),
N = ceil(7(k-1)/6).  Fourier (LEM-011): Corr_N = sum_{n!=0} What(n) D_N(n.e/Vmax).
Split by the resonance n.e mod Vmax:
  EXACT-resonance (n.e ≡ 0 mod Vmax): D_N = N  => contributes N * sum_{n!=0, Vmax|n.e} What(n)
      = N * (E_grid[W] - (6/7)^k)   [the SAME grid-decorrelation of THM-664/R2, a-priori].
  NEAR/off-resonance (n.e !≡ 0): the partial-sum correction (bounded by min(N,1/(2||.||))).
TEST: is Corr_N/N ~ (E_grid[W] - (6/7)^k)? If the decorrelation part dominates, r_N is a-priori
bounded by |E_grid[W]-(6/7)^k|/(6/7)^k -- the density-floor tail rate (R2), CLOSING the capstone.
"""
from fractions import Fraction as F
from math import ceil, gcd
from functools import reduce
import random
random.seed(61)

def W_at(E, j, V):
    ph = sorted({(e*j) % V for e in E}); m = len(ph)
    if m == 1: return F(6,7)
    Wv = F(0)
    for i in range(m):
        g = F(ph[i+1]-ph[i], V) if i < m-1 else F(ph[0]+V-ph[-1], V)
        if g > F(1,7): Wv += g - F(1,7)
    return Wv
def egrid_W(E, V):                 # (1/V) sum_{j=0}^{V-1} W(j/V)  = (6/7)^k + grid-decorr
    return sum((W_at(E,j,V) for j in range(V)), F(0)) / V
def SN(E, V, N):                    # sum_{j=1}^N W(j/V)
    return sum((W_at(E,j,V) for j in range(1, N+1)), F(0))
def prim(E):
    E = sorted(E); return len(E) >= 2 and reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

print("r_N decomposition: Corr_N/N  vs  (E_grid[W] - (6/7)^k)  [decorrelation part]\n")
for k in (8, 11, 13):
    iid = F(6,7)**k; N = ceil(7*(k-1)/6)
    print(f"k={k}  N=ceil(7(k-1)/6)={N}  iid=(6/7)^k={float(iid):.5f}:")
    print(f"   {'Vmax':>5} {'r_N':>7} {'Corr_N/N':>10} {'Egrid-iid(decorr)':>18} {'near-res=diff':>14}")
    for _ in range(6):
        # hard spread-dense cluster (j=1 fails)
        V = random.randint(200, 2500); lo = 6*V//7 + 1
        if lo >= V: continue
        s = random.randint(lo, V-1)
        mid = sorted(random.sample(range(1, s), k-2)); E = [0]+mid+[s]
        if len(set(E)) != k or not prim(E): continue
        sn = SN(E, V, N); corrN = sn - N*iid
        decorr = egrid_W(E, V) - iid          # the exact-resonance (decorrelation) part of Corr_N/N
        rN = abs(float(corrN)) / (N*float(iid))
        nearres = float(corrN)/N - float(decorr)
        print(f"   {V:>5} {rN:>7.4f} {float(corrN)/N:>10.5f} {float(decorr):>18.5f} {nearres:>14.5f}")
    print()
print("=> if Corr_N/N ~ decorr (near-res small), then r_N ~ |Egrid[W]-(6/7)^k|/(6/7)^k")
print("   = the density-floor tail rate (R2, a-priori via LEM-011) => capstone a-priori.")
