"""
opus-2026-07-11-S219: the SHARP far-element boundary term -- the missed-sector-phase cancellation, PROVED structure.

THM-699/700 give Delta_b = meas(S7(E'+b)) - Plat(E') = -(1/b) sum_i W_i(b x_i),
  W_i(y) = sum_{T subset {0..6}} (-1)^{|T|} J_i^T H_T(y),  the signed boundary weight at jump x_i.
The crude bound sums |W_i| termwise over the 127 atoms (=> C(E')/b, ~8x loose).

THE CLAIM (proved by H_T linearity): H_T = -sum_{j in T} G_j  (G_j = antiderivative of g_j=1_sectorj-1/7),
and J_i^T = 1{T cap occ = 0}(1{s in T}-1{s+1 in T})  [runner e crosses sector s->s+1; occ=other runners' sectors].
=> W_i(y) = -sum_{j in A} G_j(y) * c_j,  c_j = sum_{T subset A, j in T}(-1)^{|T|}(1{s in T}-1{s+1 in T}),  A={0..6}\occ.
The subset-sum c_j VANISHES for |A|>=3 (a (1-1)^{|A|-2}=0 telescoping). So:
  W_i(y) = 0 unless |A|<=2; for |A|=2 (A={s,s+1}): W_i(y) = G_s(y) - G_{s+1}(y).
Hence Delta_b = -(1/b) sum_{i in S}[G_{s_i}(b x_i) - G_{s_i+1}(b x_i)],  S={crossings where the OTHER runners
occupy exactly the 5 sectors != {s,s+1}} -- the near-full-coverage adjacent-pivot crossings. |S| << V.
This script VERIFIES the identity exactly and measures |S| vs V.
"""
from fractions import Fraction as F
from itertools import combinations

def Gj(j, y):
    """G_j(y) = int_0^y (1{sector j} - 1/7) dt = clamp(y - j/7, 0, 1/7) - y/7,  G_j(0)=G_j(1)=0."""
    lo = F(j,7)
    ov = min(max(y - lo, F(0)), F(1,7))   # overlap of [0,y] with sector j
    return ov - y/7

def sector(e, x):
    return int(((e*x) % 1) * 7)   # which sector frac(ex) is in

def brk(E):
    Es=[abs(e) for e in E if e!=0]; bps={F(0),F(1)}
    for e in Es:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    return sorted(b for b in bps if 0<=b<1)

def measS7(E):
    pts=brk(E)+[F(1)]; tot=F(0)
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2
        if len(set(sector(e,mid) for e in E))==7: tot+=(b-a)
    return tot
def p1(E):
    pts=brk(E)+[F(1)]; tot=F(0)
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2
        if len(set(sector(e,mid) for e in E))==6: tot+=(b-a)
    return tot

def delta_exact(Ep, b):
    E=Ep+[b]
    return measS7(E) - (measS7(Ep) + p1(Ep)/7)

def Afree(Ep, mid):
    # include the stationary runner e=0 (always in sector 0) -> sector 0 never free, A subset {1..6}
    return set(range(7)) - set(sector(e, mid) for e in Ep)

def P(A, y):
    """P(A,y) = sum_{T subset A}(-1)^{|T|} H_T(y) = G_j(y) if A={j} (exactly one free sector), else 0."""
    if len(A) == 1:
        (j,) = tuple(A)
        return Gj(j, y)
    return F(0)

def delta_boundary(Ep, b):
    """Delta_b = -(1/b) sum_x [P(A_after, b x) - P(A_before, b x)],  supported on the p1-region boundary (|A|=1)."""
    pts = brk(Ep); cells = pts + [F(1)]; M = len(pts)
    mids = [(cells[i] + cells[i+1]) / 2 for i in range(M)]   # midpoint of cell i = [pts[i], pts[i+1])
    Afr = [Afree(Ep, m) for m in mids]
    acc = F(0); S_count = 0; V_total = M
    for i in range(M):
        x = pts[i]                       # breakpoint separating cell (i-1)%M (before) and cell i (after)
        A_after = Afr[i]; A_before = Afr[(i-1) % M]
        y = (b * x) % 1
        w = P(A_after, y) - P(A_before, y)
        if w != 0 or len(A_after) == 1 or len(A_before) == 1:
            if len(A_after) == 1 or len(A_before) == 1: S_count += 1
        acc += w
    return F(-1, b) * acc, S_count, V_total

fams = {
 "consec_8 core     ":[0,1,2,3,4,5,6,7],
 "consec_5          ":[0,1,2,3,4],
 "wide 3-cluster k7 ":[0,1,2,30,31,32,33],
 "spread            ":[0,1,3,6,10,15,21],
}
print(f"{'E-prime':>20} {'b':>5} {'Delta_exact':>13} {'Delta_boundary':>15} {'match?':>7} {'|S|':>5} {'V_total':>8} {'V/|S|':>6}")
for name, Ep in fams.items():
    for b in [50, 200]:
        de = delta_exact(Ep, b)
        db, Sc, Vt = delta_boundary(Ep, b)
        match = "YES" if de==db else "NO"
        ratio = f"{Vt/Sc:.1f}" if Sc>0 else "inf"
        print(f"{name:>20} {b:>5} {float(de):>+13.7f} {float(db):>+15.7f} {match:>7} {Sc:>5} {Vt:>8} {ratio:>6}")
