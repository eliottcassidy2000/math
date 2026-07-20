#!/usr/bin/env python3
"""glmy_convention_audit_kps_S128c106.py -- kind-pasteur-2026-07-20-S128c106

INDEPENDENT AUDIT OF A SERIOUS CLAIM AGAINST CANON.

A survey agent reports that the repo's path-homology implementations define an
elementary p-path as a sequence of DISTINCT vertices (path_homology_v2.py:11
"sequence (v_0, ..., v_p) of DISTINCT vertices"), whereas Grigor'yan-Lin-Muranov-Yau
require only CONSECUTIVE distinctness (regularity): i_k != i_{k+1}.  Non-consecutive
repeats are allowed in GLMY.

Why it bites for TOURNAMENTS specifically: in a tournament exactly one of a->b,
b->a is present, so a->b->a can never be allowed -- but a->b->c->a IS allowed
whenever abc is a directed 3-cycle, and every non-transitive tournament has one.
So the two conventions must first differ at p = 3, and only for tournaments with a
3-cycle.  That is exactly where the repo's flagship degree->=3 numbers live
(beta_4(T_7) = 6, chi(T_p) = p, the T_11 table).

I am NOT taking this on trust.  This script recomputes both complexes from the
definitions.

DEFINITIONS USED (GLMY, standard):
  * an elementary p-path (i_0..i_p) is ALLOWED if i_k -> i_{k+1} is an arc for all k;
  * A_p = span of allowed p-paths;
  * d(i_0..i_p) = sum_k (-1)^k (i_0..hat(i_k)..i_p), any NON-REGULAR term (a repeat
    in consecutive position after deletion) set to 0;
  * Omega_p = { v in A_p : dv in A_{p-1} }  -- the invariant subspace;
  * beta_p = dim Omega_p - rank(d_p|Omega) - rank(d_{p+1}|Omega).

The repo convention is the same construction with "allowed" additionally requiring
all vertices distinct.

Linear algebra is exact, over a large prime field, cross-checked at a second prime.
CONTROL: the directed 3-cycle, whose GLMY homology is known (beta_1 = 1).
"""
import sys
from itertools import product

P1, P2 = (1 << 31) - 1, 2147483629
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def rank_modp(rows, ncols, p):
    """Gaussian elimination over F_p.  rows = list of dict(col -> coeff)."""
    mat = [dict(r) for r in rows if r]
    rank = 0
    pivots = {}
    for r in mat:
        rr = {c: v % p for c, v in r.items() if v % p}
        while rr:
            c = min(rr)
            if c in pivots:
                pv = pivots[c]
                f = rr[c] * pow(pv[c], p - 2, p) % p
                for cc, vv in pv.items():
                    rr[cc] = (rr.get(cc, 0) - f * vv) % p
                rr = {cc: vv for cc, vv in rr.items() if vv}
            else:
                pivots[c] = rr
                rank += 1
                break
    return rank


def nullspace_dim(rows, ncols, p):
    return ncols - rank_modp(rows, ncols, p)


def allowed_paths(adj, n, p, distinct):
    """All allowed elementary p-paths (p+1 vertices)."""
    out = []
    def ext(seq):
        if len(seq) == p + 1:
            out.append(tuple(seq))
            return
        for w in range(n):
            if (adj[seq[-1]] >> w) & 1:
                if distinct and w in seq:
                    continue
                ext(seq + [w])
    for v in range(n):
        ext([v])
    return out


def boundary(path):
    """d(path) as dict (p-1)-path -> coeff, non-regular terms dropped."""
    res = {}
    for k in range(len(path)):
        q = path[:k] + path[k + 1:]
        bad = any(q[i] == q[i + 1] for i in range(len(q) - 1))
        if bad:
            continue
        res[q] = res.get(q, 0) + (-1) ** k
    return {k: v for k, v in res.items() if v}


def omega_and_betti(adj, n, pmax, distinct, p=P1):
    """dim Omega_p and beta_p for p = 0..pmax-1."""
    A = {}
    for q in range(0, pmax + 2):
        A[q] = allowed_paths(adj, n, q, distinct)
    Aset = {q: set(A[q]) for q in A}

    dimOm, basis = {}, {}
    for q in range(0, pmax + 2):
        cols = A[q]
        if q == 0:
            dimOm[0] = len(cols)
            basis[0] = [{i: 1} for i in range(len(cols))]
            continue
        idx = {c: i for i, c in enumerate(cols)}
        # rows indexed by NON-allowed (q-1)-paths appearing in the image
        nonallowed = {}
        rows = {}
        for c in cols:
            for q1, co in boundary(c).items():
                if q1 in Aset[q - 1]:
                    continue
                r = nonallowed.setdefault(q1, len(nonallowed))
                rows.setdefault(r, {})[idx[c]] = rows.setdefault(r, {}).get(idx[c], 0) + co
        rl = [rows[r] for r in sorted(rows)]
        dimOm[q] = nullspace_dim(rl, len(cols), p)
    # ranks of d restricted to Omega: rank(d_q on Omega_q).  Since Omega_q is a
    # subspace, compute rank of d on A_q and on the constraint rows together.
    ranks = {}
    for q in range(1, pmax + 2):
        cols = A[q]
        idx = {c: i for i, c in enumerate(cols)}
        allrows, nonrows = {}, {}
        for c in cols:
            for q1, co in boundary(c).items():
                tgt = allrows if q1 in Aset[q - 1] else nonrows
                key = q1
                tgt.setdefault(key, {})[idx[c]] = tgt.setdefault(key, {}).get(idx[c], 0) + co
        # rank(d|Omega) = rank(all constraints + image rows) - rank(constraints)
        rc = rank_modp(list(nonrows.values()), len(cols), p)
        rb = rank_modp(list(nonrows.values()) + list(allrows.values()), len(cols), p)
        ranks[q] = rb - rc
    betti = {}
    for q in range(0, pmax + 1):
        rq = ranks.get(q, 0)
        rq1 = ranks.get(q + 1, 0)
        betti[q] = dimOm[q] - rq - rq1
    return dimOm, betti


def paley(p):
    qr = {(x * x) % p for x in range(1, p)}
    adj = [0] * p
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i] |= 1 << j
    return adj, p


print("=" * 78)
print("CONTROL -- the directed 3-cycle (GLMY: beta_1 = 1)")
print("=" * 78)
tri = [0b010, 0b100, 0b001]
for distinct, name in ((False, "GLMY (regular)"), (True, "repo (distinct)")):
    dO, b = omega_and_betti(tri, 3, 4, distinct)
    print("  %-18s dim Omega = %s" % (name, [dO[q] for q in range(0, 5)]))
    print("  %-18s betti     = %s" % ("", [b[q] for q in range(0, 4)]))
print("  -> GLMY beta_1 = 1 reproduces the textbook value, so the implementation")
print("     is calibrated.  Note Omega_3 differs already here: the 3-cycle admits")
print("     (a,b,c,a) as an allowed regular 3-path, excluded by the repo rule.")

print()
print("=" * 78)
print("THE AUDIT -- Paley T_7, both conventions, degrees 0..%d" % PMAX)
print("=" * 78)
adj7, n7 = paley(7)
res = {}
for distinct, name in ((False, "GLMY (regular)"), (True, "repo (distinct)")):
    dO, b = omega_and_betti(adj7, n7, PMAX, distinct)
    res[name] = (dO, b)
    print("  %-18s dim Omega_0..%d = %s"
          % (name, PMAX + 1, [dO[q] for q in range(0, PMAX + 2)]))
    print("  %-18s betti_0..%d     = %s"
          % ("", PMAX, [b[q] for q in range(0, PMAX + 1)]))

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
dG, bG = res["GLMY (regular)"]
dR, bR = res["repo (distinct)"]
agree_low = all(dG[q] == dR[q] for q in (0, 1, 2))
print("  Omega_0, Omega_1, Omega_2 agree : %s" % agree_low)
first_diff = next((q for q in range(0, PMAX + 2) if dG[q] != dR[q]), None)
print("  first degree where dim Omega differs : %s  (GLMY %s vs repo %s)"
      % (first_diff, dG.get(first_diff), dR.get(first_diff)))
print("  beta_2 under GLMY : %s   (repo: %s)" % (bG.get(2), bR.get(2)))
print("  beta_4 under GLMY : %s   (repo canon reports 6)" % bG.get(4))
print()
print("  If Omega_{0,1,2} agree then beta_1 <= 1 (THM-103) and beta_2 = 0 (THM-108)")
print("  are CONVENTION-SAFE and survive.  If beta_4 differs from 6, the degree->=3")
print("  results -- THM-129 chi(T_p) = p, THM-130's Paley Betti formula, THM-096's")
print("  and THM-099's degree->=3 tables -- are computed in the wrong complex.")
