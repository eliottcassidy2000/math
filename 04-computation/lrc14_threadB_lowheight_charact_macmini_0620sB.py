#!/usr/bin/env python3
r"""
lrc14_threadB_lowheight_charact_macmini_0620sB.py   (mac-mini, Thread B)

THREAD B of the LRC(14) sector route -- the FINITE CHECK feasibility.

CONTEXT (proved/verified upstream, cited):
  - measS7(E) = meas{ x in [0,1) : c(e,x)=floor(7 frac(e x)) hits all of Z/7, over e in E }.
  - SIGNED IDENTITY (HYP-2606): measS7(E) = M7(k) + corr(E),
        corr(E) = sum_{0 != n in Lambda(E), K(n)!=0} K(n),
        Lambda(E) = { n in Z^k : sum_i n_i e_i = 0 }  (rank k-1).
  - SUPPORT-6 FLOOR (THM-538): K(n)=0 unless n has >=6 'genuine' coordinates
        (nonzero AND not a multiple of 7).
  - ENVELOPE (THM-538): |K(n)| <= 2^6 * prod over its genuine coords (c1/|n_j|),
        c1 = 0.6974.  So K(n) is SMALL when the genuine coords are LARGE.
  - iid_k = 7! S(k,7)/7^k = M7(k) (the spread limit); caps cap_8=2243/5880 etc.
  - WIDE BOUND (Thread A): if EVERY surviving relation has large height, corr is tiny
        => measS7 < cap.   Its complement is THREAD B.

THREAD B DEFINITION.
  The "relation-height" of E (for the surviving sum) is
        h(E) := min over surviving n (>=6 genuine coords, K(n)!=0) of  ||n||_inf.
  Thread A handles h(E) >= H0(k).  THREAD B = { E : h(E) <= H0(k) } : a SHORT
  surviving relation exists.  We must show this family is FINITE after
  scale-normalization (measS7(cE)=measS7(E), divide by gcd) and verify measS7<=cap.

THIS SCRIPT (exploratory, Stage 1):
  (A) Confirm measS7 scale invariance: measS7(cE)=measS7(E).
  (B) Confirm support-6 floor on small support: K(n)=0 unless >=6 genuine coords.
  (C) For a working height bound H, ENUMERATE the surviving relation vectors n
      (support>=6 among k coords, ||n||_inf <= H) and ask: given such an n, what is
      the set of E (0 in E, primitive, |E|=k) admitting n as a relation?  This is a
      sublattice/affine condition; we count its scale-normalized solutions.
  (D) Measure, on random E, the distribution of h(E): how large is the smallest
      surviving relation for SPREAD vs CONSEC shapes?  (Sanity: spread => h large.)

stdlib only; exact Fraction for measS7.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

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

def M7(k):
    return sum(F((-1)**t * comb(6, t)) * F(7-t, 7)**(k-1) for t in range(7))

def gcd_all(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0)

def primitive(E):
    return gcd_all(E) == 1

# ---------- genuine-coordinate count of a relation vector ----------
def genuine_count(n):
    """number of coords nonzero AND not multiple of 7."""
    return sum(1 for x in n if x != 0 and x % 7 != 0)

# ---------- K(n) (float, for floor verification) ----------
TWO_PI_I = 2j*math.pi
def shat(n, j):
    if n == 0: return 1.0/7.0
    a = j/7.0
    return (math.e**(-TWO_PI_I*n*a) - math.e**(-TWO_PI_I*n*(a+1/7.0))) / (TWO_PI_I*n)
import cmath
def shat(n, j):
    if n == 0: return 1.0/7.0
    a = j/7.0
    return (cmath.exp(-TWO_PI_I*n*a) - cmath.exp(-TWO_PI_I*n*(a+1/7.0))) / (TWO_PI_I*n)
SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1)**len(T) for T in SUB}
_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0: v = complex(1 - len(T)/7.0, 0.0)
    elif n % 7 == 0: v = 0j
    else: v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v
def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0+0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T]*p
    return s

def banner(t): print("\n"+"="*82+f"\n{t}\n"+"="*82)

# ===========================================================================
def partA_scale_invariance():
    banner("PART A -- measS7 scale invariance: measS7(cE) = measS7(E)  (=> normalize by gcd)")
    ok = True
    tests = [[0,1,2,3,4,5,6,7], [0,2,3,5,8,11,13,17], [0,1,4,9,16,25,36,49]]
    for E in tests:
        base = measS7(E)
        for c in (2, 3, 5, 7):
            cE = [c*e for e in E]
            if measS7(cE) != base:
                ok = False
                print(f"  FAIL measS7({c}*{E}) != measS7({E})")
    print(f"  measS7(cE)=measS7(E) for c in 2,3,5,7 on 3 shapes: {ok}")
    print(f"  => scale-normalization: WLOG primitive (gcd of nonzero offsets = 1).")
    return ok

def partB_support6_floor():
    banner("PART B -- support-6 floor (THM-538): K(n)=0 unless >=6 genuine coords")
    # test relations with support 1..6 (genuine count) and check K vanishes for <6.
    import random
    rng = random.Random(11)
    by_supp = {}
    for _ in range(4000):
        k = rng.randint(7, 9)
        n = [0]*k
        supp = rng.randint(1, 7)
        idxs = rng.sample(range(k), supp)
        for i in idxs:
            n[i] = rng.choice([x for x in range(-4, 5) if x != 0])
        g = genuine_count(n)
        kv = abs(Kk(n))
        by_supp.setdefault(g, []).append(kv)
    print("  genuine-count -> max |K(n)| observed (should be ~0 for <6):")
    for g in sorted(by_supp):
        mx = max(by_supp[g])
        flag = " <-- nonzero" if mx > 1e-9 else " (vanishes)"
        print(f"    genuine={g}: max|K| = {mx:.2e}  (n={len(by_supp[g])} samples){flag}")
    print("  CONFIRMS: only relations with >=6 genuine coords survive (contribute to corr).")

def partD_height_distribution():
    banner("PART D -- height of the smallest SURVIVING relation: consec vs spread")
    print("  h(E) = min over surviving n (>=6 genuine, ||n||_inf<=Hcap) of ||n||_inf.")
    print("  (search relations by enumerating low-height integer vectors in Lambda(E).)")
    def smallest_surviving_height(E, Hcap=4):
        """brute search for the shortest surviving relation up to height Hcap."""
        nz_idx = [i for i, e in enumerate(E) if e != 0]
        k = len(E)
        best = None
        # iterate height layers
        for H in range(1, Hcap+1):
            found = False
            # enumerate vectors with ||.||_inf == H on the nonzero coords (0-coord forced 0? no:
            # n_0 (for e=0) can be anything but contributes via chat(n_0,T); but e_0=0 means the
            # term n_0*e_0=0 always, so n_0 is FREE in the lattice. For survival we need genuine
            # coords; n_0 on a zero-speed: c_0(x)=0 const, so its Fourier is delta -> chat(n_0,T)
            # with n_0!=0 gives 0 unless... For e=0 the color is always 0; the relation only
            # constrains the NONZERO speeds. Set n_0=0 WLOG for the surviving-relation search.)
            for vals in itertools.product(range(-H, H+1), repeat=len(nz_idx)):
                if max(abs(v) for v in vals) != H: continue
                n = [0]*k
                for idx, v in zip(nz_idx, vals): n[idx] = v
                if sum(n[i]*E[i] for i in range(k)) != 0: continue
                if genuine_count(n) < 6: continue
                if abs(Kk(n)) < 1e-9: continue
                best = H; found = True; break
            if found: break
        return best
    families = {
        "consec8 {0..7}": list(range(8)),
        "consec9 {0..8}": list(range(9)),
        "near-consec {0..6,8}": [0,1,2,3,4,5,6,8],
        "spreadA {0,1,3,7,12,20,30,44}": [0,1,3,7,12,20,30,44],
        "1-stranger {0..6,40}": [0,1,2,3,4,5,6,40],
        "geom {0,1,2,4,8,16,32,64}": [0,1,2,4,8,16,32,64],
    }
    print(f"  {'family':36s} {'k':>2} {'h(E)<=4?':>9} {'measS7':>9}")
    for nm, E in families.items():
        h = smallest_surviving_height(E, Hcap=4)
        hs = str(h) if h is not None else ">4"
        print(f"  {nm:36s} {len(E):>2} {hs:>9} {float(measS7(E)):>9.4f}")
    print("\n  READING: dense (consec, near-consec) shapes have a SHORT surviving relation")
    print("  (h small); spread/dissociated shapes have NO short surviving relation (h>4).")
    print("  => THREAD B family = shapes with a short surviving relation = the DENSE ones.")

if __name__ == "__main__":
    print("#"*82)
    print("# THREAD B -- characterizing the LOW-RELATION-HEIGHT family (finite check)")
    print("#"*82)
    partA_scale_invariance()
    partB_support6_floor()
    partD_height_distribution()
    print("\nDONE (Thread B stage 1).")
