"""
mac-mini-2026-07-08-S59 -- the pairwise Weyl correction, exactly, for THM-527-A large-spread.

E_grid[W] = (6/7)^k + corrections;  the PAIRWISE correction is
  C_2 = sum_{i<j} ( T(V'_{ij}) - 1/49 ),   V'_{ij} = Vmax/gcd(e_i-e_j, Vmax),
  T(V') = (1/V') sum_{r=0}^{V'-1} (1/7 - ||r/V'||)_+   (grid-avg of the arc-overlap tent).
Exact closed form (R = #positive grid pts within 1/7 of 0):
  T(V') = (1+2R)/(7V') - R(R+1)/V'^2,   R = floor((V'-1)/7)  [verified below].
This script: (a) EXACT T(V') and max |T(V')-1/49| + its decay; (b) the FULL correction
E_grid[W]-(6/7)^k decomposed into pairwise vs higher-order, for structured clusters.
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce

def T_exact(Vp):
    """exact grid-avg of the tent (1/7-||t||)_+ over {r/Vp}."""
    tot = F(0)
    for r in range(Vp):
        d = min(r, Vp-r)                 # ||r/Vp|| = d/Vp
        val = F(1,7) - F(d, Vp)
        if val > 0: tot += val
    return tot / Vp

def T_formula(Vp):
    R = (Vp-1)//7
    return F(1+2*R, 7*Vp) - F(R*(R+1), Vp*Vp)

print("=== (a) exact T(V') vs 1/49, and the closed form ===")
oneq = F(1,49)
maxdev = 0; maxdevV = 0
for Vp in range(1, 400):
    t = T_exact(Vp)
    assert t == T_formula(Vp), f"formula mismatch at {Vp}: {t} vs {T_formula(Vp)}"
    dev = abs(float(t) - float(oneq))
    if dev > maxdev: maxdev, maxdevV = dev, Vp
print(f"  closed form T(V')=(1+2R)/(7V')-R(R+1)/V'^2, R=floor((V'-1)/7): VERIFIED V'=1..399")
print(f"  max |T(V')-1/49| over V'=1..399 = {maxdev:.5f} at V'={maxdevV}")
# decay: bound |T(V')-1/49| <= c/V'
c = max(abs(float(T_exact(Vp))-float(oneq))*Vp for Vp in range(1,400))
print(f"  max |T(V')-1/49|*V' over V'=1..399 = {c:.4f}  => |T(V')-1/49| <= {c:.3f}/V'")
print(f"  small V': T(1)={float(T_exact(1)):.4f} T(2)={float(T_exact(2)):.4f} T(7)={float(T_exact(7)):.4f} "
      f"T(14)={float(T_exact(14)):.4f}  (1/49={float(oneq):.4f})")

# === (b) full correction decomposition for structured clusters ===
def W_at(E, j, Vmax):
    ph = sorted({(e*j) % Vmax for e in E}); m = len(ph)
    if m == 1: return F(6,7)
    W = F(0)
    for i in range(m):
        g = F(ph[i+1]-ph[i], Vmax) if i < m-1 else F(ph[0]+Vmax-ph[-1], Vmax)
        if g > F(1,7): W += g - F(1,7)
    return W
def egrid_W(E, Vmax):
    return sum((W_at(E,j,Vmax) for j in range(Vmax)), F(0)) / Vmax
def prim(E): E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1

def pairwise_corr(E, Vmax):
    C2 = F(0)
    for a in range(len(E)):
        for b in range(a+1, len(E)):
            d = abs(E[a]-E[b]); Vp = Vmax // gcd(d, Vmax)
            C2 += T_exact(Vp) - F(1,49)
    return C2

print("\n=== (b) E_grid[W] - (6/7)^k = pairwise C_2 + higher-order, structured clusters ===")
def twoblock(k,s):
    a=k//2; b=k-a; return sorted(set(list(range(a))+[s-b+1+i for i in range(b)]))
for k in (11,13):
    iid = F(6,7)**k
    print(f"k={k} (iid (6/7)^k={float(iid):.5f}):")
    for s in [30,60,120]:
        E=twoblock(k,s)
        if len(E)!=k or not prim(E): continue
        for Vmax in [s+1, s+7]:
            if Vmax<=max(E): continue
            eg = egrid_W(E,Vmax); C2 = pairwise_corr(E,Vmax)
            full = eg - iid; ho = full - C2
            print(f"   2-block s={s:3d} Vmax={Vmax:3d}: E_grid[W]={float(eg):.5f} "
                  f"full corr={float(full):+.5f} = pairwise {float(C2):+.5f} + higher-order {float(ho):+.5f}")
