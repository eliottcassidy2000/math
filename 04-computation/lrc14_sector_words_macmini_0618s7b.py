#!/usr/bin/env python3
"""
lrc14_sector_words_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B: cutting sequences)

GOAL: re-express meas(S7(E)) in a word/character form and probe AP-extremality.

SETUP. sigma_e(x) = floor(7 e x) mod 7 in Z/7 is the CUTTING SEQUENCE of the slope-7e line.
S7(E) = {x in [0,1): {sigma_e(x): e in E} = Z/7 (surjective)}.  e=0 pins sigma_0(x)=0 always.

THREE-DISTANCE / WORD STRUCTURE.  For a single e>=1, the breakpoints of sigma_e are at
x = m/(7e), m=0..7e, and sigma_e cycles 0,1,2,3,4,5,6,0,1,...  over the 7e equal cells of
width 1/(7e). So over [0,1), sigma_e is a perfectly periodic word (0123456)^e on a grid of
7e cells.  The JOINT word over E is the superposition (common refinement) of these grids.

INCLUSION-EXCLUSION (codex's M7+corr).  meas(S7) = sum_{T subset Z/7, 0 not forced... }
Actually with 0 in E forcing residue 0 hit:
   meas(surjective) = sum_{T subset {1..6}} (-1)^{|T|} meas(B_T),
   B_T = {x: NO e in E has sigma_e(x) in T} = {x: for all e, sigma_e(x) notin T}.
For nonzero e, P[sigma_e(x) in B_T-allowed] over a full period = (7-|T|)/7 EXACTLY, and
across DIFFERENT e these events are INDEPENDENT iff the e are "dissociated enough".
M7(k) = sum_T (-1)^|T| (1-|T|/7)^{k-1}  is the independent (main) term.
corr(E) = meas(S7) - M7(k) = the DEPENDENCE between cutting words = relation-lattice tail.

THIS SCRIPT:
  (A) verify the per-e period structure and the independence main term exactly;
  (B) express meas(B_T) as the measure of a JOINT cutting word avoiding a sector-block T,
      and decompose corr by Fourier mode n in Z (the relation lattice);
  (C) test: does the AP maximize meas(B_T) for EACH single T, or only in aggregate?
      (If AP maximizes each meas(B_T) with the IE sign, extremality is structural.)
  (D) the character-word identity: meas(B_T) = sum_{n} prod_e c_T-hat(...) -- check 7|n vanish.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact meas of {x: for all e in E, frac(7 e x) avoids the union of sectors in T} ----
# sigma_e(x) = floor(7 e x) mod 7.  Avoid set T subset {0..6}: keep cells where sigma_e notin T.
def measB_T(E, T):
    """meas{x in [0,1): for all e in E, (floor(7 e x) mod 7) notin T}. e=0 -> sigma=0, so if 0 in T -> 0."""
    E = sorted(set(E))
    Tset = set(T)
    # zero handling
    if 0 in E and 0 in Tset:
        return F(0)
    Enz = [e for e in E if e != 0]
    if not Enz:
        return F(1)  # only e=0, and 0 not in T
    # breakpoints: x = m/(7e)
    bps = set([F(0), F(1)])
    for e in Enz:
        for m in range(0, 7*e+1):
            bps.add(F(m, 7*e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        ok = True
        for e in Enz:
            s = int((7*e*xm)) % 7   # floor(7 e xm) mod 7
            if s in Tset:
                ok = False; break
        if ok:
            total += x1 - x0
    return total

def measS7(E):
    """inclusion-exclusion over T subset {1..6} (sector 0 auto-hit by e=0)."""
    E = sorted(set(E))
    # general: surjective onto Z/7. Use IE over which residues are MISSED.
    # meas(surj) = sum_{M subset Z/7} (-1)^|M| meas(all e avoid M).  But residue 0 always hit if 0 in E.
    s = F(0)
    for r in range(0, 8):
        for M in itertools.combinations(range(7), r):
            s += F((-1)**r) * measB_T(E, M)
    return s

def M7(k):
    s = F(0)
    for t in range(0,7):
        s += F((-1)**t * comb(6,t)) * F(7-t,7)**(k-1)
    return s

# ---------------- (A) sanity: per-e period & independence main term ----------------
print("="*92)
print("(A) Single-e period structure and the independence main term M7(k)")
print("="*92)
for e in [1,2,3]:
    # meas{x: sigma_e(x) avoids T} for T of size t should be exactly (7-t)/7
    for t in [1,2,3]:
        T = list(range(1, 1+t))
        v = measB_T([e], T)
        print(f"  e={e}, |T|={t}: meas(avoid T) = {v} = {float(v):.5f}  (expect {(7-t)}/7 = {(7-t)/7:.5f})")
print()
for k in [8,9,10,11,12,13]:
    print(f"  M7({k}) = {M7(k)} = {float(M7(k)):.6f}")

# ---------------- (B) per-T extremality test: does AP maximize meas(B_T)? -----------
print()
print("="*92)
print("(B) Per-block extremality: for each sector-block T, is AP_k a max of (-1)^|T| meas(B_T)?")
print("    meas(S7) = sum_T (-1)^|T| meas(B_T).  If AP maximizes each IE-signed term, structural.")
print("="*92)

def normalize(E):
    E = sorted(set(E)); g = 0
    for e in E: g = gcd(g, e)
    if g: E = [e//g for e in E]
    return tuple(E)

def gen_shapes(k, maxE):
    """primitive normalized 0-based k-subsets of {0..maxE} with 0, gcd 1."""
    out = []
    for rest in itertools.combinations(range(1, maxE+1), k-1):
        E = (0,) + rest
        g = 0
        for e in E: g = gcd(g, e)
        if g != 1: continue
        out.append(E)
    return out

for k in [8, 9]:
    maxE = k + 4
    AP = tuple(range(k))
    shapes = gen_shapes(k, maxE)
    # For each T (by size, take representative blocks), find argmax of signed term
    print(f"\n  k={k}, box maxE<={maxE}, #shapes={len(shapes)}")
    # aggregate AP-extremality already known; here test per-T
    # We test all T grouped; report which T have a NON-AP maximizer of |meas(B_T)|
    bad_T = []
    # iterate over T up to symmetry: meas(B_T) depends on T only via... actually via the gaps of T on the circle Z/7
    seen_class = {}
    def Tclass(T):
        # circular gap multiset of complement positions? Use sorted cyclic structure of T
        T = sorted(T); n=len(T)
        if n==0: return ('empty',)
        # rotations
        best=None
        for sh in range(7):
            rot = tuple(sorted((t+sh)%7 for t in T))
            if best is None or rot<best: best=rot
        return (n, best)
    Treps = {}
    for r in range(1,7):
        for T in itertools.combinations(range(7), r):
            c = Tclass(T)
            if c not in Treps: Treps[c]=T
    for c,T in Treps.items():
        vals = [(measB_T(E,T), E) for E in shapes]
        mx = max(vals, key=lambda p: p[0])
        mn = min(vals, key=lambda p: p[0])
        ap_val = measB_T(AP, T)
        # IE sign favors: even |T| -> want max meas(B_T); odd |T| -> want MIN meas(B_T) to maximize -meas
        sign = (-1)**len(T)
        if sign>0:
            extremal = mx; want="MAX"
        else:
            extremal = mn; want="MIN(neg term)"
        is_ap = (ap_val == extremal[0])
        tag = "AP-extremal" if is_ap else f"NOT-AP (extremal E={extremal[1]}, val={float(extremal[0]):.4f} vs AP {float(ap_val):.4f})"
        if not is_ap:
            bad_T.append((T, want, extremal[1]))
    if not bad_T:
        print(f"    ALL T-blocks: AP is the IE-correct extremiser. (per-T structural extremality holds)")
    else:
        print(f"    {len(bad_T)} T-blocks where AP is NOT the per-T IE-extremiser:")
        for T,want,E in bad_T[:12]:
            print(f"       T={T} ({want}): extremal E={E}")
print()
print("DONE.")
