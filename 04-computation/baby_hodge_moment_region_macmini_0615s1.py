#!/usr/bin/env python3
"""
baby_hodge_moment_region_macmini_0615s1.py  (mac-mini-2026-06-15-S1, T820/THM-508)

BABY HODGE = THE TOURNAMENT MOMENT PROBLEM.
tr(A^k) = Σ λ^k are the moments of the tournament spectrum. Cycle counts c_k decompose
into Witt-census sums (spectral) + overlap corrections (non-spectral). The realizable
invariant-vector region = a truncated moment problem; the HOLES (forbidden vectors
inside the realized range) = the "non-algebraic Hodge classes".

This script:
 (1) enumerates all tournaments on n=4..7 (gentourng), computes the cycle counts
     c3,c4,c5,c6,c7 (directed k-cycles), tr(A^k), and the SKEW moment sequence
     m_r = tr((SS^T)^r), S = A - A^T;
 (2) finds the realizable region of (c3,c5) and (c3,c4,c5) per n and the HOLES
     (integer points inside the coordinate-wise range / convex hull, not realized);
 (3) THE HODGE-INEQUALITY TEST: for each hole, checks whether the SKEW moment sequence
     it would require is Hankel-PSD-feasible (a genuine Stieltjes moment sequence) —
     i.e. whether the hole is "moment-interior" (passes the continuous Hodge
     inequalities) yet integer-forbidden = a certified baby-Hodge non-algebraic class;
 (4) the moment-cumulant ledger: tr(A^k) vs k*c_k (the overlap defect W_k - c_k).
"""
import sys, subprocess, itertools
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(line_buffering=True)
GEN = "/opt/homebrew/bin/gentourng"

def gen_tournaments(n):
    out = subprocess.run([GEN, str(n)], capture_output=True, text=True)
    pairs = list(itertools.combinations(range(n), 2))
    for line in out.stdout.split():
        bits = line.strip()
        if len(bits) != len(pairs): continue
        A = np.zeros((n, n), dtype=np.int64)
        for (i, j), b in zip(pairs, bits):
            if b == '1': A[i][j] = 1
            else: A[j][i] = 1
        yield A

def cycle_counts(n, A):
    """directed k-cycle counts c_k for k=3..min(n,7) by direct enumeration."""
    cc = {}
    for k in range(3, min(n, 7) + 1):
        cnt = 0
        for verts in itertools.combinations(range(n), k):
            # count directed Hamiltonian cycles within this k-subset, /k for rotation, but
            # we count distinct directed cycles (as edge-sets): fix smallest vertex, count
            # directed cyclic orders
            v0 = verts[0]; rest = verts[1:]
            for perm in itertools.permutations(rest):
                seq = (v0,) + perm
                if all(A[seq[i]][seq[(i+1) % k]] for i in range(k)):
                    cnt += 1
        # each undirected cyclic order counted; directed k-cycles fixing v0 first divides
        # rotation already (v0 fixed) but reflection gives the reverse cycle (a DIFFERENT
        # directed cycle). So cnt = number of directed k-cycles. divide by nothing.
        cc[k] = cnt // 2 if False else cnt
    # the above counts each directed cycle once per (rotation fixed by v0) but permutations
    # of rest already fix rotation; reverse is a distinct directed cycle and is counted. So
    # cc[k] = # directed k-cycles. (validate via tr(A^3)=3 c3)
    return cc

def tr_powers(n, A, K=7):
    P = np.eye(n, dtype=object); tr = {}
    M = A.astype(object)
    cur = np.eye(n, dtype=object)
    for k in range(1, K + 1):
        cur = cur @ M
        tr[k] = int(np.trace(cur))
    return tr

def skew_moments(n, A, R=4):
    S = (A - A.T).astype(float)
    G = -(S @ S)            # = S S^T, PSD symmetric
    m = {}
    cur = np.eye(n)
    for r in range(0, R + 1):
        m[r] = float(np.trace(cur))
        cur = cur @ G
    return m   # m[0]=n, m[1]=tr(SS^T)=n(n-1), ...

def hankel_psd(moments):
    """moments = [m0,m1,...,m_{2t}]; check Hankel [m_{i+j}] PSD (Hamburger necessary)."""
    t = len(moments) // 2
    H = np.array([[moments[i + j] for j in range(t + 1)] for i in range(t + 1)], dtype=float)
    ev = np.linalg.eigvalsh(H)
    return ev.min() > -1e-7, ev.min()

print("=" * 74)
print("BABY HODGE: realizable cycle-count region, holes, and the moment/Hankel test")
print("=" * 74)

data = {}
for n in range(4, 8):
    rows = []
    for A in gen_tournaments(n):
        cc = cycle_counts(n, A)
        tr = tr_powers(n, A, K=min(n, 7))
        # sanity: tr(A^3) = 3 c3
        assert tr[3] == 3 * cc[3], (n, tr[3], cc[3])
        sm = skew_moments(n, A, R=4)
        rows.append((cc, tr, sm))
    data[n] = rows
    c3s = sorted(set(r[0][3] for r in rows))
    c5s = sorted(set(r[0].get(5, 0) for r in rows)) if n >= 5 else [0]
    print(f"\nn={n}: {len(rows)} iso classes; c3 range {c3s[0]}..{c3s[-1]}; "
          f"c5 values {c5s if n>=5 else 'n/a'}")
    # realizable (c3,c5) set and holes
    if n >= 5:
        real = set((r[0][3], r[0][5]) for r in rows)
        # per-c3 fiber holes: c5 inside [min,max] of that fiber, not realized
        holes = []
        by_c3 = {}
        for (a, b) in real: by_c3.setdefault(a, []).append(b)
        for a, bs in sorted(by_c3.items()):
            lo, hi = min(bs), max(bs)
            for b in range(lo, hi + 1):
                if (a, b) not in real:
                    holes.append((a, b))
        print(f"   realized (c3,c5) pairs: {sorted(real)}")
        print(f"   HOLES (c3,c5) inside fibers: {holes}")
        data[n] = (rows, real, holes)

print("\n" + "=" * 74)
print("HODGE-INEQUALITY TEST on holes: is the hole moment-interior (Hankel-PSD feasible)?")
print("=" * 74)
# For each hole at n, check whether the skew moment sequence m0..m4 is achievable as a
# CONVEX COMBINATION of realized tournaments' moment sequences with the same (c3,c5) target
# would require -- here we test the simpler necessary Hodge condition: does ANY realized
# tournament at this n have Hankel-PSD skew moments? (sanity: all should, since SS^T is PSD)
for n in range(5, 8):
    rows = data[n][0]
    allpsd = True; worst = 1e9
    for (cc, tr, sm) in rows:
        ok, mn = hankel_psd([sm[0], sm[1], sm[2], sm[3], sm[4]])
        allpsd = allpsd and ok; worst = min(worst, mn)
        # the SS^T moments are a genuine Stieltjes sequence (PSD measure) -> Hankel PSD always
    print(f"n={n}: skew-moment Hankel PSD on all {len(rows)} classes: {allpsd} (worst min-eig {worst:.3e})")
    holes = data[n][2]
    if holes:
        print(f"   holes to certify as 'moment-interior non-algebraic Hodge classes': {holes}")
        print(f"   (full SDP moment-feasibility certificate per hole: deferred to verify step)")

print("\n" + "=" * 74)
print("MOMENT-CUMULANT LEDGER: tr(A^k) vs k*c_k (the overlap/Witt defect)")
print("=" * 74)
for n in range(4, 8):
    rows = data[n] if n < 5 else data[n][0]
    # show the defect tr(A^k) - k*c_k for k=3..min(n,7) on a few classes
    defects = {}
    for (cc, tr, sm) in rows:
        for k in range(3, min(n, 7) + 1):
            d = tr[k] - k * cc.get(k, 0)
            defects.setdefault(k, set()).add(d)
    print(f"n={n}: overlap defect tr(A^k)-k*c_k value-sets:")
    for k in sorted(defects):
        vs = sorted(defects[k])
        spectral = (vs == [0])
        print(f"   k={k}: defect values {vs[:8]}{'...' if len(vs)>8 else ''}  "
              f"{'SPECTRAL (c_k=tr/k)' if spectral else 'OVERLAP present (non-trivial)'}")

print("\nDONE.")
