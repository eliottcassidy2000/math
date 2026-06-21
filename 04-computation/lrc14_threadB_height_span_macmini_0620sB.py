#!/usr/bin/env python3
r"""
lrc14_threadB_height_span_macmini_0620sB.py    (mac-mini, Thread B stage 2)

THE CRUX OF THREAD B.  Reframed correctly after THM-538-CORRECTED / HYP-2646:
short relations DO contribute to corr (no support-6 floor on the zero-padded kernel).
So the relevant notion is the height of the SHORTEST relation in the FULL lattice
   Lambda(E) = { n in Z^{k-1} : sum_j n_j e_j = 0 }   (over the NONZERO offsets; e_0=0 is free).

DEFINE (l_infty relation height of E):
   h(E) := min { ||n||_inf : 0 != n in Lambda(E) }      (the first successive minimum, l_inf).
THREAD A (wide bound) wants h(E) >= H0.   THREAD B = { E primitive, 0 in E, |E|=k, h(E) <= H0 }.

THE KEY QUESTION (finiteness):  is the family
   B_{k,H0} := { E : 0 in E, |E|=k, primitive, h(E) <= H0 }
FINITE (i.e. span bounded)?    A short relation sum n_j e_j = 0 with small coeffs ties the
offsets together.  Intuition: a primitive E with a single short relation can still have ONE
free large offset (e.g. {0,1,2,...,k-2, N}: the bounded part has a short relation, N is free).
=> B_{k,H0} is INFINITE in general.   The finiteness must come from a STRONGER condition:
the WHOLE lattice must be low-height (all k-2 successive minima small), i.e. E is
"totally low-height" / sub-AP, which DOES bound span.   We test BOTH:

  (Q1) h(E) (first minimum) small  -> span unbounded?  (expect YES: one free stranger.)
  (Q2) lambda_{k-2}(E) (LAST minimum) small -> span bounded?  (the real finiteness lever.)
  (Q3) The DICHOTOMY: either E has a SHORT FULL-RANK low-height lattice (=> bounded span,
       finite family, finite check) OR E has a HIGH last-minimum (=> a dissociated DIRECTION,
       wide bound applies).   We want the threshold and the count.

THE GROUNDED FACT (prompt): within span<=13 all shapes are relheight<=6; dissociation needs
large span.  We make this PRECISE: relate span to the successive minima, and bound the family.

stdlib only; exact Fraction for measS7; integer linear algebra for the lattice minima.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

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

def gcd_all(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0)
def banner(t): print("\n"+"="*82+f"\n{t}\n"+"="*82)

# ---------- successive-minima (l_infty) of the relation lattice ----------
def first_minimum_linf(speeds, Hcap=8):
    """h(E) = min ||n||_inf over nonzero n with sum n_j speeds_j = 0.
       speeds = the NONZERO offsets (the e_0=0 coordinate is free, omitted)."""
    d = len(speeds)
    for H in range(1, Hcap+1):
        # vectors with ||.||_inf == H
        for n in itertools.product(range(-H, H+1), repeat=d):
            if max(abs(x) for x in n) != H: continue
            if sum(n[i]*speeds[i] for i in range(d)) == 0:
                return H, n
    return None, None

def rank_of_lowheight_relations(speeds, H):
    """rank over Q of the set of relations with ||n||_inf <= H.  If = d-1 (full corank-1
       lattice rank for a 1-dim speed vector => relation lattice rank = d-1), the whole
       lattice is generated below height H => E is 'totally low-height'."""
    d = len(speeds)
    rels = []
    for n in itertools.product(range(-H, H+1), repeat=d):
        if all(x == 0 for x in n): continue
        if sum(n[i]*speeds[i] for i in range(d)) == 0:
            rels.append(list(n))
    if not rels: return 0
    # rank over Q via fraction gaussian elimination
    return mat_rank_Q(rels)

def mat_rank_Q(rows):
    M = [[F(x) for x in r] for r in rows]
    rank = 0; ncol = len(M[0]); r0 = 0
    for c in range(ncol):
        piv = None
        for r in range(r0, len(M)):
            if M[r][c] != 0: piv = r; break
        if piv is None: continue
        M[r0], M[piv] = M[piv], M[r0]
        pv = M[r0][c]
        for r in range(len(M)):
            if r != r0 and M[r][c] != 0:
                f = M[r][c]/pv
                M[r] = [M[r][j]-f*M[r0][j] for j in range(ncol)]
        r0 += 1; rank += 1
        if r0 == len(M): break
    return rank

# ===========================================================================
def Q1_first_min_unbounded():
    banner("Q1 -- first minimum h(E) small does NOT bound span (one free stranger)")
    print("  Family {0,1,2,...,k-2, N}: the bounded part {0..k-2} has h=1 (e.g. 1+1-2=0),")
    print("  but N is a FREE large offset => span unbounded with h(E)=1.  Confirm:")
    print(f"  {'k':>2} {'N':>6} {'h(E)':>5} {'span':>6} {'measS7':>9}")
    for k in (8,):
        for N in (10, 50, 200, 1000):
            E = list(range(k-1)) + [N]
            sp = sorted(e for e in E if e != 0)
            h, n = first_minimum_linf(sp, Hcap=3)
            print(f"  {k:>2} {N:>6} {str(h):>5} {max(E):>6} {float(measS7(E)):>9.4f}")
    print("  => first-minimum smallness alone gives an INFINITE family. Need a stronger notion.")

def Q2_last_min_bounds_span():
    banner("Q2 -- the LAST minimum / FULL low-height lattice DOES bound span")
    print("  If the relation lattice Lambda(E) is GENERATED by vectors of height <= H (i.e.")
    print("  the rank-(d-1) lattice has all successive minima <= H), then the speeds satisfy")
    print("  d-1 independent low-height relations, pinning the speed RATIOS to a bounded set.")
    print("  Test: for primitive E, if rank of {||n||<=H relations} = d-1, span is bounded by f(H,k).\n")
    print(f"  {'shape':34s} {'d':>2} {'rank(H=2)':>9} {'d-1':>4} {'span':>6} {'full-lowH?':>10}")
    shapes = {
        "consec8 {0..7}": list(range(8)),
        "{0,1,2,3,4,5,6,8}": [0,1,2,3,4,5,6,8],
        "{0,1,2,...,6,40} stranger": [0,1,2,3,4,5,6,40],
        "{0,2,4,6,8,10,12,14} AP*2": [0,2,4,6,8,10,12,14],
        "{0,1,3,7,12,20,30,44} Sidon": [0,1,3,7,12,20,30,44],
        "{0,3,6,9,12,15,18,21} AP*3": [0,3,6,9,12,15,18,21],
        "{0,1,2,4,8,16,32,64} geom": [0,1,2,4,8,16,32,64],
    }
    for nm, E in shapes.items():
        sp = sorted(e for e in E if e != 0)
        d = len(sp)
        rk = rank_of_lowheight_relations(sp, H=2)
        full = (rk == d-1)
        print(f"  {nm:34s} {d:>2} {rk:>9} {d-1:>4} {max(E):>6} {str(full):>10}")
    print("\n  => 'full-low-height' (rank of height<=H relations = d-1) is the FINITENESS criterion.")
    print("     These are the sub-AP / Freiman-structured shapes.  Strangers/Sidon are NOT full-low-H.")

def Q3_span_bound_from_full_lowheight():
    banner("Q3 -- EXPLICIT span bound for full-low-height (sub-AP) primitive shapes")
    print("  CLAIM: if the relation lattice of the nonzero speeds is generated by height<=H")
    print("  vectors (full rank d-1), then primitivity forces span <= S(d,H).  We MEASURE the")
    print("  largest span achievable by a primitive full-low-height shape, by enumeration.\n")
    print("  Enumerate primitive E (0 in E, |E|=k) with span<=SP and check full-low-height at H;")
    print("  report the MAX span among full-low-height ones (should saturate => finite family).")
    for k in (8,):
        for H in (1, 2, 3):
            maxspan = 0; cnt = 0; argmax = None
            SP = {1: 14, 2: 16, 3: 16}[H]
            for combo in itertools.combinations(range(1, SP+1), k-1):
                E = [0] + list(combo)
                if gcd_all(E) != 1: continue
                sp = list(combo)
                rk = rank_of_lowheight_relations(sp, H=H)
                if rk == len(sp)-1:
                    cnt += 1
                    if max(E) > maxspan: maxspan = max(E); argmax = E
            print(f"  k={k} H={H}: #full-low-H primitive with span<={SP}: {cnt}; "
                  f"MAX span among them = {maxspan}  e.g. {argmax}")
    print("\n  If MAX span < SP (the search bound), the full-low-height family is span-bounded => FINITE.")

if __name__ == "__main__":
    print("#"*82)
    print("# THREAD B stage 2 -- relation height vs span: the finiteness dichotomy")
    print("#"*82)
    Q1_first_min_unbounded()
    Q2_last_min_bounds_span()
    Q3_span_bound_from_full_lowheight()
    print("\nDONE (Thread B stage 2).")
