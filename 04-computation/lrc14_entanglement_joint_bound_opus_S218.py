"""
opus-2026-07-11-S218: finding the exact provable JOINT inequality for the ungapped Plat<->Delta entanglement.

THM-700 (kps, PROVED): for E=E' u {w}, w=max E, |p0(E)-Plat(E')| <= V(E')/(6w),
  Plat(E')=p0(E')+(1/7)p1(E'),  V(E')=sum_s V(f_s),  f_s=1{E' misses exactly sector s}.
So p0(E) <= Plat(E') + V(E')/(6w).  The ungapped case (w~max E') has V/(6w)=O(1) BUT Plat(E') small (wide).

CANDIDATE JOINT LEMMAS to test (find the tightest that holds with margin):
  (J1)  Plat(E') + V(E')/(6w) <= cap_k            [one-shot, crude V]
  (J2)  Plat(E') + Delta_w    <= cap_k            [exact, = p0(E) <= cap_k, the target]
  (J3)  the ENTANGLEMENT: does small Plat compensate large V/w?  plot Plat(E') vs V(E')/(6w).
  (J4)  full recursive peel: p0(E) = p0(core) + sum_i[(1/7)p1(E_i) + Delta_{w_i}]; does it telescope < cap?
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

def sectors_hit(E, x):
    occ = set()
    for e in E:
        occ.add(int(((e*x) % 1)*7))
    return occ

def breakpoints(E):
    Es = [abs(e) for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Es:
        for m in range(0, 7*e+1):
            bps.add(F(m, 7*e))
    return sorted(b for b in bps if 0 <= b < 1)

def measures(E):
    """Return p0 (all 7 hit), p1 (exactly 6 hit = exactly 1 missed), and V (sum_s TV of f_s)."""
    pts = breakpoints(E); pts2 = pts + [F(1)]
    p0 = F(0); p1 = F(0)
    # f_s value on each cell, for s=0..6 (missed-exactly-{s})
    cell_missone = []   # for each cell: the missed sector if exactly one missed, else None
    for i in range(len(pts2)-1):
        a, b = pts2[i], pts2[i+1]
        if b <= a: continue
        mid = (a+b)/2
        occ = sectors_hit(E, mid)
        miss = set(range(7)) - occ
        length = b - a
        if len(miss) == 0:
            p0 += length; cell_missone.append((None, length))
        elif len(miss) == 1:
            p1 += length; cell_missone.append((next(iter(miss)), length))
        else:
            cell_missone.append(("multi", length))
    # total variation of f_s (s=1..6): count transitions of the indicator 1{missed=={s}} across cells (cyclic)
    labels = [c[0] for c in cell_missone]   # None / int s / "multi"
    V = 0
    n = len(labels)
    for s in range(7):
        prev = (labels[-1] == s)
        for lab in labels:
            cur = (lab == s)
            if cur != prev: V += 1
            prev = cur
    return p0, p1, V

def cap(k):
    caps = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}
    return caps[k]

def peel_once(E):
    w = max(E); Ep = [e for e in E if e != w]
    p0E,_,_ = measures(E)
    p0p, p1p, Vp = measures(Ep)
    Plat = p0p + p1p/7
    Delta = p0E - Plat
    return w, Ep, p0E, p0p, p1p, Vp, Plat, Delta

# wide k=9 families (0 in E), spanning ungapped-multiscale to gapped
fams = {
    "consec_9 (bounded)   ": [0,1,2,3,4,5,6,7,8],
    "wide 3-cluster       ": [0,1,2,30,31,32,60,61,62],
    "wide 3-cluster x10   ": [0,1,2,100,101,102,200,201,202],
    "multiscale 2-2-2-2-1 ": [0,1,2,30,31,32,60,61,150],
    "odd-struct+far       ": [0,1,3,5,7,9,10,11,90],
    "gapped {0..7,500}    ": [0,1,2,3,4,5,6,7,500],
    "2-scale {0..3,300..} ": [0,1,2,3,300,301,302,303,600],
}
c9 = cap(9); Q8 = F(621,1715)   # Phi(consec_8)=0.36210
print(f"cap_9 = {float(c9):.5f}   Q(8)=Phi(consec_8) = {float(Q8):.5f}   margin_9 = cap_9-Q8 = {float(c9-Q8):.5f}")
print(f"Delta_cap = cap_9 - cap_8 = {float(c9-cap(8)):.5f}\n")
print(f"{'family':>22} {'p0(E)':>8} {'Plat(E)':>8} {'Delta_w':>9} {'p1(E)':>7} {'V(E)':>5} {'V/(6w)':>8} {'J1:Plat+V/6w':>13} {'<=cap9?':>7}")
for name, E in fams.items():
    w, Ep, p0E, p0p, p1p, Vp, Plat, Delta = peel_once(E)
    j1 = Plat + F(Vp, 6*w)
    ok = "OK" if j1 <= c9 else "FAIL"
    print(f"{name:>22} {float(p0E):>8.5f} {float(Plat):>8.5f} {float(Delta):>+9.5f} {float(p1p):>7.4f} {Vp:>5} {float(F(Vp,6*w)):>8.4f} {float(j1):>13.5f} {ok:>7}")

print("\n=== the ENTANGLEMENT: is Plat(E') + V(E')/(6w) <= cap always? and how does Plat trade vs V/w? ===")
print("(if J1 holds with margin, the ONE-SHOT crude-V joint bound closes the wide half -- no accumulation needed)\n")

print("=== full recursive peel of the worst multiscale set (accumulation) ===")
E = [0,1,2,30,31,32,60,61,62]
depth = 0; acc_tax = F(0)
p0_start,_,_ = measures(E)
print(f"start p0({E}) = {float(p0_start):.5f}, cap_9={float(c9):.5f}")
cur = E[:]
while max(cur) > 8 and len(cur) > 5:
    w, Ep, p0E, p0p, p1p, Vp, Plat, Delta = peel_once(cur)
    tax = p1p/7 + Delta
    acc_tax += tax
    print(f"  peel w={w:>4}: p0(E)={float(p0E):.5f} -> p0(E')={float(p0p):.5f} + tax[(1/7)p1+Delta]={float(tax):+.5f}  (V/6w={float(F(Vp,6*w)):.4f})")
    cur = Ep; depth += 1
    if depth > 6: break
print(f"  bounded core {cur}: p0={float(measures(cur)[0]):.5f}; accumulated tax={float(acc_tax):+.5f}; core+tax={float(measures(cur)[0]+acc_tax):.5f} <= cap_9={float(c9):.5f}?")
