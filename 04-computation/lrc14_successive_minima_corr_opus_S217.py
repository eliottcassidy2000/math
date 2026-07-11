"""
opus-2026-07-11-S217: the successive-minima / correction experiment for the LRC(14) last gap.

THE OBJECT (THM-538/546/HYP-2606):  p0(E) = meas(S7(E)) = meas{ x in [0,1) : the runners {e_i} of E
miss at least one of the 7 sectors [j/7,(j+1)/7) }  (0 in E, |E|=k; the stationary runner fixes sector 0,
so it is really about the 6 inner sectors and the k-1 moving runners).  Exact signed expansion:
    p0(E) = M7(k) + corr(E),   corr(E) = Sum_{0 != n in Lambda-o(E)} K(n),  K(n)=D7(n mod 7)/prod n_j,
Lambda-o(E) = { n in Z^{k-1} : Sum_j n_j e_j = 0 }  (the offset relation lattice, rank k-2),
M7(k) = Sum_{t=0}^6 (-1)^t C(6,t)(1-t/7)^{k-1}  (the iid / dissociated limit).

The correction is a reciprocal-weighted lattice sum, DOMINATED BY SHORT RELATIONS (small prod|n_j|).
CLAIM UNDER TEST (the successive-minima bound): |corr(E)| is controlled by lambda_1(Lambda-o(E)) =
length of the SHORTEST nonzero relation.  AP has short relations (e.g. 1+1-2=0 -> (1,-2,1)) => large corr;
dissociated (Sidon) has only LONG relations => corr -> 0 => p0 -> M7(k) < cap_k.  If |corr| decays with
lambda_1, the ungapped-wide dissociated regime closes by geometry of numbers.
"""
from fractions import Fraction as F
from itertools import product as iproduct
import sympy

def p0_exact(E):
    """Exact meas{ x in [0,1) : occupied sectors {floor(7 frac(e_i x))} != {0..6} }, breakpoint sweep.
    Breakpoints: x with 7 e_i x = integer, i.e. x = m/(7 e_i)."""
    Es = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Es:
        e = abs(e)
        for m in range(0, 7*e + 1):
            bps.add(F(m, 7*e))
    pts = sorted(b for b in bps if 0 <= b < 1)
    pts.append(F(1))
    total = F(0)
    for i in range(len(pts)-1):
        a, b = pts[i], pts[i+1]
        if b <= a: continue
        mid = (a+b)/2
        occ = set()
        occ.add(0)  # stationary runner e=0 -> sector 0
        for e in Es:
            fr = (e*mid) % 1
            occ.add(int(fr*7))    # which sector [j/7,(j+1)/7)
        if len(occ) != 7:         # at least one sector missed
            total += (b-a)
    return total

def M7(k):
    from math import comb
    return sum(F((-1)**t)*comb(6,t)*F(7-t,7)**(k-1) for t in range(7))

def lambda1(E, cap=6):
    """Shortest nonzero relation Sum n_j e_j = 0 over the moving offsets, by ell1 norm Sum|n_j|.
    Brute force over |n_j| <= cap.  Returns (min ell1, the relation)."""
    Es = [e for e in E if e != 0]
    r = len(Es)
    best = None
    for n in iproduct(range(-cap, cap+1), repeat=r):
        if all(x==0 for x in n): continue
        if sum(nj*ej for nj,ej in zip(n,Es)) == 0:
            l1 = sum(abs(x) for x in n)
            if best is None or l1 < best[0]:
                best = (l1, n)
    return best if best else (None, None)

# families: AP (short relations) -> dissociated/Sidon (long-only relations), all k=8 (0 + 7 movers)
fams = {
    "consec {0..7}      ": [0,1,2,3,4,5,6,7],
    "AP dilated 3*      ": [0,3,6,9,12,15,18,21],
    "near-AP (1 detune) ": [0,1,2,3,4,5,6,10],
    "mild spread        ": [0,1,2,4,7,11,16,22],
    "Sidon-ish          ": [0,1,3,7,12,20,30,44],
    "Sidon (B2) small   ": [0,1,3,7,12,20],           # k=6 true Sidon
    "wide dissociated   ": [0,1,5,11,23,47,95,191],   # near-geometric, long relations
    "wide prime-ish     ": [0,7,17,31,53,83,127,179],
}

k_ref = 8
print(f"M7(8) = {float(M7(8)):.6f}   M7(6) = {float(M7(6)):.6f}   (the iid/dissociated limit; corr -> 0 there)\n")
print(f"{'family':>22} {'k':>3} {'p0(E)':>10} {'M7(k)':>9} {'corr':>10} {'lambda1':>8} {'shortest reln':>22}")
for name, E in fams.items():
    k = len(E)
    p0 = p0_exact(E)
    m = M7(k)
    corr = p0 - m
    l1, reln = lambda1(E, cap=6)
    rs = str(reln) if reln and len(reln)<=8 else "(>cap)"
    print(f"{name:>22} {k:>3} {float(p0):>10.6f} {float(m):>9.6f} {float(corr):>+10.6f} {str(l1):>8} {rs:>22}")

print("\n=== far-element peel (THM-546): Delta_w = p0(E'+w) - [p0(E') + (1/7)p1(E')] vs (6/49)V/w ===")
def p1_exact(E):
    """meas{ x : E misses EXACTLY one sector }."""
    Es = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Es:
        e=abs(e)
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    pts = sorted(b for b in bps if 0<=b<1); pts.append(F(1))
    tot=F(0)
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2; occ={0}
        for e in Es: occ.add(int(((e*mid)%1)*7))
        if len(occ)==6:  # exactly one of 7 missed
            tot+=(b-a)
    return tot
core = [0,1,2,3,4,5,6]        # consec_7 core
p0c, p1c = p0_exact(core), p1_exact(core)
Phi = p0c + p1c/7
V_actual = None  # (skip exact arc count; use the measured Delta rate)
print(f"core consec_7: p0={float(p0c):.6f} p1={float(p1c):.6f}  Phi=p0+p1/7={float(Phi):.6f}")
for w in [10, 20, 50, 100, 200, 400]:
    E = core + [w]
    dw = p0_exact(E) - Phi
    print(f"  w={w:>4}: Delta_w={float(dw):+.6f}   |Delta_w|*w={float(abs(dw)*w):>8.4f}   (6/49)V bound ~ const?")
