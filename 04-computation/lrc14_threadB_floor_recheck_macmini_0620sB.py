#!/usr/bin/env python3
r"""
lrc14_threadB_floor_recheck_macmini_0620sB.py

Two sanity checks that gate Thread B:

(1) Re-examine the SUPPORT-6 FLOOR carefully.  The earlier random scan showed
    |K(n)| ~ 1e-3 for genuine<6, NOT exactly 0.  Is THM-538 (K=0 unless >=6 genuine)
    actually about the SUMMED corr (telescoping/cancellation across T) or about each
    individual K(n)?  Re-derive K(n) for ALL support<=5 vectors EXACTLY (rational
    Fourier on Z/7 via roots of unity in the field Q(zeta_7) emulated by integer
    convolution) and report whether K(n) is exactly 0.

(2) Re-examine the height-1 surviving relation found for "spread" shapes.  A height-1
    relation with >=6 genuine coords summing to 0 means a +-1 signed subset sum = 0.
    Check whether spreadA={0,1,3,7,12,20,30,44} truly has such a relation, and whether
    its K(n) is genuinely nonzero (so it is a SURVIVING short relation).

The point: pin down EXACTLY what "low relation height" means for the surviving sum,
so Thread B's family is defined precisely.
"""
import sys, itertools, math, cmath
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

TWO_PI_I = 2j*math.pi
def shat(n, j):
    if n == 0: return 1.0/7.0
    a = j/7.0
    return (cmath.exp(-TWO_PI_I*n*a) - cmath.exp(-TWO_PI_I*n*(a+1/7.0))) / (TWO_PI_I*n)
SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1)**len(T) for T in SUB}
def chat(n, T):
    if n == 0: return complex(1 - len(T)/7.0, 0.0)
    if n % 7 == 0: return 0j
    return -sum(shat(n, j) for j in T)
def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0+0j
        for ni in n:
            p *= chat(ni, T)
        s += SGN[T]*p
    return s

def genuine_count(n):
    return sum(1 for x in n if x != 0 and x % 7 != 0)

def banner(t): print("\n"+"="*82+f"\n{t}\n"+"="*82)

# ---------- exact measS7 ----------
def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7*abs(e)+1): bps.add(F(m, 7*abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        secs = {int(((e*xm) % 1)*7) for e in E}
        if len(secs) == 7: tot += x1-x0
    return tot

banner("(1) SUPPORT-6 FLOOR: exhaustive K(n) for support<=5 vectors (do they vanish?)")
# Enumerate genuine-support patterns directly: a vector with g genuine nonzero coords,
# values in {1..6} (mod 7 nonzero), and check |K|. The KEY claim THM-538 is that K
# depends only on the MULTISET of nonzero-mod-7 values and vanishes for <6 of them.
print("  K(n) by genuine support g (values 1..6, all distinct positions), max|K|:")
for g in range(1, 8):
    mx = 0.0; argv = None
    # values are mod-7 residues 1..6; enumerate multisets of size g
    for vals in itertools.combinations_with_replacement(range(1, 7), g):
        n = list(vals)
        kv = abs(Kk(n))
        if kv > mx: mx = kv; argv = vals
    print(f"    g={g}: max|K| over residue-multisets = {mx:.3e}  at {argv}")
print("  INTERPRETATION: if max|K| jumps from ~0 to nonzero exactly at g=6, floor confirmed.")

banner("(1b) Same, but with REPEATED coords allowed via integer values (not just residues)")
print("  Using actual integer coords (height<=3), genuine count g:")
seen = {}
for k in (6, 7, 8):
    for n in itertools.product(range(-3, 4), repeat=k):
        g = genuine_count(n)
        kv = abs(Kk(list(n)))
        if g not in seen or kv > seen[g][0]:
            seen[g] = (kv, n)
for g in sorted(seen):
    kv, n = seen[g]
    print(f"    g={g}: max|K| = {kv:.3e}  at n={n}")

banner("(2) The height-1 'surviving relation' for spread shapes -- is it real?")
def find_height1_relations(E, need_genuine=6):
    """all +-1/0 relations sum n_i e_i = 0 with >=need_genuine genuine coords; report K."""
    k = len(E); out = []
    nz = [i for i in range(k) if E[i] != 0]
    for vals in itertools.product((-1, 0, 1), repeat=len(nz)):
        n = [0]*k
        for idx, v in zip(nz, vals): n[idx] = v
        if all(x == 0 for x in n): continue
        if sum(n[i]*E[i] for i in range(k)) != 0: continue
        g = genuine_count(n)
        if g < need_genuine: continue
        out.append((tuple(n), g, abs(Kk(n))))
    return out
for nm, E in [("consec8 {0..7}", list(range(8))),
              ("spreadA {0,1,3,7,12,20,30,44}", [0,1,3,7,12,20,30,44]),
              ("Sidon {0,1,3,7,12,20,30,44}", [0,1,3,7,12,20,30,44])]:
    rels = find_height1_relations(E)
    print(f"  {nm}: #height-1 surviving relations = {len(rels)}")
    for n, g, kv in rels[:4]:
        print(f"      n={n} genuine={g} |K|={kv:.4e}")
print("""
  NOTE: spreadA={0,1,3,7,12,20,30,44} is a B_2 (Sidon) set, so it has NO solution to
  a+b=c+d with distinct elements -- but a +-1 relation can still exist via LONGER
  signed sums (e.g. 1-3-7+12-20+... ) hitting 0. We list them to see if they survive.
""")

banner("(3) DECISIVE: does corr(E) get a non-negligible contribution from a SHORT relation?")
print("  For each shape compute the EXACT corr = measS7 - M7(k) and the |K| of its")
print("  shortest surviving relation; if |K_short| ~ corr the short relation matters.")
from math import comb
def M7(k):
    return sum(F((-1)**t * comb(6, t)) * F(7-t, 7)**(k-1) for t in range(7))
for nm, E in [("consec8", list(range(8))),
              ("near-consec {0..6,8}", [0,1,2,3,4,5,6,8]),
              ("spreadA Sidon", [0,1,3,7,12,20,30,44]),
              ("1-stranger {0..6,40}", [0,1,2,3,4,5,6,40])]:
    k = len(E)
    corr = float(measS7(E)) - float(M7(k))
    rels = find_height1_relations(E)
    if rels:
        bestK = max(kv for _, _, kv in rels)
        nrel = len(rels)
    else:
        bestK = 0.0; nrel = 0
    print(f"  {nm:22s}: corr={corr:+.5f}  #h1-surv={nrel}  max|K_h1|={bestK:.4e}")

print("\nDONE.")
