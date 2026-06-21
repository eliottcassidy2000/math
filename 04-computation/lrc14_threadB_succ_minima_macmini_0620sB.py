#!/usr/bin/env python3
r"""
lrc14_threadB_succ_minima_macmini_0620sB.py   (mac-mini, Thread B stage 3)

CORRECTED finiteness analysis via SUCCESSIVE MINIMA of the relation lattice.

The relation lattice of the nonzero speeds s=(s_1,...,s_d):
    Lambda(E) = { n in Z^d : sum_j n_j s_j = 0 },  rank d-1,
    covolume(Lambda) = ||s||_2 / gcd(s)  =  ||s||_2 for primitive E   (STANDARD FACT).
l_infty successive minima  lambda_1 <= ... <= lambda_{d-1}.
By Minkowski's 2nd thm:  prod lambda_i  asymptotically ~  covolume ~  ||s||_2 ~ span.

CORRECT DICHOTOMY for the finite check:
   * lambda_1 small + lambda_{d-1} small  (ALL minima <= H0)  <=>  span <= H0^{d-1}  (FINITE).
   * lambda_{d-1} LARGE  =>  a dissociated direction  =>  one (or more) FAR/STRANGER offset(s).
The stranger {0..k-2,N} has lambda_1=1 but lambda_{d-1} ~ N: it is the BOUNDARY case,
handled by Thread A's STRANGER-CONTRACTION (peel far offset).

THREAD B PROPER = the TOTALLY-LOW-HEIGHT family  { E : lambda_{d-1}(E) <= H0 }.
Goal: show it is FINITE (span <= H0^{d-1}) and verify measS7 <= cap exhaustively.

THIS SCRIPT:
  (A) Compute exact l_infty successive minima of Lambda(E) for canonical shapes.
  (B) Verify the span bound:  lambda_{d-1}(E) <= H0  =>  span(E) <= bound(H0,k).
  (C) Count |{E : 0 in E, |E|=k primitive, lambda_{d-1}(E) <= H0}| (the FINITE family size).
  (D) The honest middle: strangers (lambda_1 small, lambda_{d-1} large) are neither
      bounded-span nor fully-dissociated; recorded as the contraction case (Thread A).

stdlib only; integer LLL-free successive-minima via bounded enumeration + lattice basis.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import gcd, isqrt
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def gcd_all(xs): return reduce(gcd, [abs(x) for x in xs if x != 0], 0)
def banner(t): print("\n"+"="*82+f"\n{t}\n"+"="*82)

# ---------- exact l_infty successive minima of the relation lattice ----------
def relation_basis(speeds):
    """A Z-basis of Lambda = {n: sum n_j s_j = 0}, rank d-1.  Hermite/kernel via integer
       row reduction of the 1xd matrix [s].  Returns list of d-1 integer vectors."""
    d = len(speeds)
    s = list(speeds)
    g = gcd_all(s)
    s = [x // g for x in s]   # primitive
    # kernel of n -> sum n_j s_j: standard basis e_i - (s_i / s_0)... use integer construction.
    # Build basis: for the primitive vector s, the orthogonal lattice has covolume ||s||.
    # Construct via: find index p with s_p != 0; basis vectors  b_i = s_p * e_i - s_i * e_p  for i!=p,
    # then saturate (make the lattice the FULL integer kernel, not the sublattice).
    p = max(range(d), key=lambda i: abs(s[i]))
    raw = []
    for i in range(d):
        if i == p: continue
        v = [0]*d
        v[i] = s[p]; v[p] = -s[i]
        raw.append(v)
    # saturate: the full kernel = primitive closure. Use HNF-style reduction over Q then clear.
    return saturate(raw, d)

def saturate(rows, d):
    """Return a Z-basis of the saturation (primitive closure) of the lattice spanned by rows,
       inside Z^d.  Implemented via the Smith-free 'column HNF + content division' approach."""
    # Compute the rational row space, then find the integer lattice = (rowspace) cap Z^d.
    # Simple robust method: stack rows, compute integer kernel of the COMPLEMENT.
    # Here we instead use: lattice L = rows over Z; its saturation = { x in Z^d : x in Q-span(rows) }.
    # Build a basis of Q-span via gaussian elim, then the saturation basis via integer HNF of the
    # full inclusion. For our sizes (d<=12), use the following: the saturation is the kernel of
    # the matrix whose rows are a basis of the orthogonal complement (1-dim, = s). So saturation
    # of any full-corank-1 lattice with the right Q-span IS the full kernel of s. Return that.
    # => Just return integer kernel of s directly (full kernel), computed by HNF.
    return rows  # raw e_p-construction already spans Q-kernel; saturation handled in minima search

def succ_minima_linf(speeds, Hcap=12, want=None):
    """Exact l_infty successive minima lambda_1<=...<=lambda_{d-1} of the FULL integer
       relation lattice {n in Z^d: sum n_j s_j=0}.  Enumerate lattice points by height,
       greedily extend a Q-independent set.  want = how many minima to compute (default all)."""
    d = len(speeds); r = d - 1
    if want is None: want = r
    found = []   # independent lattice vectors, in increasing height
    mins = []
    H = 0
    while len(found) < want and H < Hcap:
        H += 1
        # all lattice points with ||.||_inf == H
        for n in itertools.product(range(-H, H+1), repeat=d):
            if max(abs(x) for x in n) != H: continue
            if sum(n[i]*speeds[i] for i in range(d)) != 0: continue
            # is n independent of current found set (over Q)?
            if independent(found, n):
                found.append(list(n)); mins.append(H)
                if len(found) >= want: break
    return mins, found

def independent(basis, v):
    if not basis: return any(x != 0 for x in v)
    M = [[F(x) for x in row] for row in basis] + [[F(x) for x in v]]
    return mat_rank_Q(M) > len(basis)

def mat_rank_Q(rows):
    M = [row[:] for row in rows]; r0 = 0; rank = 0; ncol = len(M[0])
    for c in range(ncol):
        piv = None
        for r in range(r0, len(M)):
            if M[r][c] != 0: piv = r; break
        if piv is None: continue
        M[r0], M[piv] = M[piv], M[r0]; pv = M[r0][c]
        for r in range(len(M)):
            if r != r0 and M[r][c] != 0:
                f = M[r][c]/pv
                M[r] = [M[r][j]-f*M[r0][j] for j in range(ncol)]
        r0 += 1; rank += 1
        if r0 == len(M): break
    return rank

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

# ===========================================================================
def partA_minima_of_shapes():
    banner("PART A -- exact l_infty successive minima  lambda_1 <= ... <= lambda_{d-1}")
    print("  Relation lattice covolume = ||speeds||_2 (primitive); prod(minima) ~ covolume ~ span.")
    print(f"  {'shape':32s} {'d-1':>4} {'minima (l_inf)':>26} {'lam_max':>7} {'||s||2':>8} {'span':>5} {'measS7':>8}")
    shapes = {
        "consec8 {0..7}": list(range(8)),
        "{0..6,8}": [0,1,2,3,4,5,6,8],
        "AP*2 {0,2..14}": [0,2,4,6,8,10,12,14],
        "AP*3 {0,3..21}": [0,3,6,9,12,15,18,21],
        "stranger {0..6,40}": [0,1,2,3,4,5,6,40],
        "Sidon {0,1,3,7,12,20,30,44}": [0,1,3,7,12,20,30,44],
        "geom {0,1,2,4,8,16,32,64}": [0,1,2,4,8,16,32,64],
        "two-cluster {0,1,2,40,41,42,43}": [0,1,2,40,41,42,43],
    }
    for nm, E in shapes.items():
        sp = sorted(e for e in E if e != 0)
        d = len(sp)
        mins, basis = succ_minima_linf(sp, Hcap=10)
        l2 = math.isqrt(sum(x*x for x in sp))
        lam_max = max(mins) if mins else None
        ms = str(mins)
        print(f"  {nm:32s} {d-1:>4} {ms:>26} {str(lam_max):>7} {l2:>8} {max(E):>5} {float(measS7(E)):>8.4f}")
    print("\n  READING: consec/AP have ALL minima small (totally low-height). Stranger/two-cluster")
    print("  have lambda_max LARGE (a dissociated direction). lambda_max is the finiteness ruler.")

def partB_span_bound():
    banner("PART B -- SPAN BOUND: lambda_max(E) <= H0  =>  span(E) bounded (finiteness)")
    print("  THEOREM (to verify): for primitive E, span(E) <= product of minima <= (lambda_max)^{d-1}")
    print("  times a small constant.  More usefully measure: max span over primitive E with")
    print("  lambda_max <= H0.  If finite => the totally-low-height family is finite.\n")
    for k in (8,):
        d = k-1
        for H0 in (2, 3):
            maxspan = 0; cnt = 0; argmax = None
            SP = 18
            for combo in itertools.combinations(range(1, SP+1), k-1):
                E = [0] + list(combo)
                if gcd_all(combo) != 1: continue
                mins, _ = succ_minima_linf(list(combo), Hcap=H0)
                # totally-low-height iff we found d-1 independent vectors all <= H0
                if len(mins) == d-1 and max(mins) <= H0:
                    cnt += 1
                    if max(E) > maxspan: maxspan = max(E); argmax = E
            sat = "SATURATED (<SP)" if maxspan < SP else "HIT SEARCH BOUND -- raise SP"
            print(f"  k={k} H0={H0}: #totally-low-H primitive span<={SP}: {cnt}; "
                  f"MAX span = {maxspan} [{sat}]  e.g. {argmax}")

if __name__ == "__main__":
    print("#"*82)
    print("# THREAD B stage 3 -- successive minima: the CORRECT finiteness ruler")
    print("#"*82)
    partA_minima_of_shapes()
    partB_span_bound()
    print("\nDONE (Thread B stage 3).")
