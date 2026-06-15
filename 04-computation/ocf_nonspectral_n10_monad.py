#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-505 : THE GROWTH LAW OF dim_nonspec(H)  -- INTRINSIC (basis-independent)
          non-spectral dimension via carrier-delta RANK, n = 8, 9, 10.
          (monad-explorer-2026-06-15-S3)

Builds on my n=9 result.  Two advances:

(A) INTRINSIC dimension via RANK, robust to small cospectral classes.
    Within cospectral classes (fixed trace vector tr A^3..A^n), collect the
    delta vectors of the OCF-carrier set and take their rank over Q.  Pooling
    deltas across MANY small classes makes the rank well-defined even when each
    class has only 2-3 members (the regime forced by sampling at n=10, where the
    8-component trace signature is very discriminating).  This sidesteps the
    bucket-chain's small-class artifact (singleton (sig,carrier) buckets falsely
    read as "determined").

(B) The OCF CARRIER BASIS (basis-independent).  H = I(Omega,2) = 1+2a1+4a2+8a3
    is built from ODD-cycle packings:
        a1 = sum of single odd cycles c_k (k odd)             [c3,c5 spectral]
        a2 = sum of DISJOINT odd-cycle PAIR counts D_{l,l'}   [l<=l' odd, l+l'<=n]
        a3 = sum of DISJOINT odd-cycle TRIPLE counts          [sum<=n]
    The NON-SPECTRAL OCF carriers (per n) are exactly:
        n=8 : {c7, D33, D35}                       -> predicted dim 3
        n=9 : {c7, c9, D33, D35, T333}             -> predicted dim 5
        n=10: {c7, c9, D33, D35, D37, D55, T333}   -> predicted dim 7
    KEY: in this basis the n=9 dim is 5, not the 6 of the trace basis
    {c6,c7,c8,c9,Q44,T333} -- because c8 and Q44 enter H ONLY via the SUM
    c8+Q44 = D35 - c3*c5 + W8 (one carrier).  We verify this directly.

GROWTH-LAW CONJECTURE (the headline):
    dim_nonspec(H)(n) = #{ partitions of s into odd parts >=3, 3<=s<=n }  - 2
                      = [ sum_{s<=n} p_{odd>=3}(s) ]  - 3
    (the -3 removes the spectral configs: empty packing, single c3, single c5).
    Predicts 1,2,3,5,7,9,12,... for n=6,7,8,9,10,11,12,...
    Generating function: P(x) = prod_{j>=1} 1/(1 - x^(2j+1)).

Usage: python3 ocf_nonspectral_n10_monad.py [N10_PHASE1] [N10_CAP] [N89_SAMPLES]
"""

import sys
import numpy as np
from collections import defaultdict
from math import comb
from fractions import Fraction

# ---------- phase-1 numpy trace-signature sampling (general n) ----------

def sample_sigs(n, nsamples, batch, seed):
    """Yield (sig_tuple, code_int) for random T_n.  code = upper-triangle bits."""
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(n, 1)
    m = len(iu[0])
    done = 0
    while done < nsamples:
        b = min(batch, nsamples - done)
        bits = rng.integers(0, 2, size=(b, m), dtype=np.int8)
        A = np.zeros((b, n, n), dtype=np.int64)
        A[:, iu[0], iu[1]] = bits
        A[:, iu[1], iu[0]] = 1 - bits
        P = np.matmul(A, A)            # A^2
        powtr = {}
        for k in range(3, n + 1):
            P = np.matmul(P, A)
            powtr[k] = np.einsum('bii->b', P)
        for i in range(b):
            sig = tuple(int(powtr[k][i]) for k in range(3, n + 1))
            code = int(np.dot(bits[i].astype(object), 1 << np.arange(m, dtype=object)))
            yield sig, code
        done += b


def code_to_adj(code, n):
    A = [[0] * n for _ in range(n)]
    t = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> t) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            t += 1
    return A

# ---------- pure-python carriers ----------

def matmul(X, Y, n):
    Z = [[0] * n for _ in range(n)]
    for i in range(n):
        Xi = X[i]; Zi = Z[i]
        for k in range(n):
            x = Xi[k]
            if x:
                Yk = Y[k]
                for j in range(n):
                    Zi[j] += x * Yk[j]
    return Z


def traces(A, n, kmax):
    out = {}
    P = [row[:] for row in A]
    out[1] = sum(P[i][i] for i in range(n))
    for k in range(2, kmax + 1):
        P = matmul(P, A, n)
        out[k] = sum(P[i][i] for i in range(n))
    return out


def cycles_by_len(A, n, maxL):
    adj = [[j for j in range(n) if A[i][j]] for i in range(n)]
    out = defaultdict(list)
    for start in range(n):
        smask = 1 << start
        stack = [(start, smask, 1)]
        while stack:
            v, vis, ln = stack.pop()
            for w in adj[v]:
                if w == start and ln >= 3:
                    out[ln].append(vis)
                elif w > start and not (vis >> w) & 1 and ln < maxL:
                    stack.append((w, vis | (1 << w), ln + 1))
    return out


def overlap_pairs(la, lb, same):
    cnt = 0
    if same:
        L = len(la)
        for i in range(L):
            ai = la[i]
            for j in range(i + 1, L):
                if ai & la[j]:
                    cnt += 1
    else:
        for x in la:
            for y in lb:
                if x & y:
                    cnt += 1
    return cnt


def disjoint_pairs(la, lb, same):
    cnt = 0
    if same:
        L = len(la)
        for i in range(L):
            ai = la[i]
            for j in range(i + 1, L):
                if not (ai & la[j]):
                    cnt += 1
    else:
        for x in la:
            for y in lb:
                if not (x & y):
                    cnt += 1
    return cnt


def disjoint_triples_same(la):
    L = len(la); cnt = 0
    for i in range(L):
        ai = la[i]
        for j in range(i + 1, L):
            aj = la[j]
            if ai & aj:
                continue
            aij = ai | aj
            for k in range(j + 1, L):
                if not (aij & la[k]):
                    cnt += 1
    return cnt


def disjoint_T335(tri, pent):
    """#(unordered triangle pair, pentagon) all pairwise vertex-disjoint."""
    L = len(tri); cnt = 0
    for i in range(L):
        ai = tri[i]
        for j in range(i + 1, L):
            aj = tri[j]
            if ai & aj:
                continue
            u = ai | aj
            for p in pent:
                if not (u & p):
                    cnt += 1
    return cnt


def H_dp(A, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            cnt = row[v]
            if cnt == 0 or not (mask & (1 << v)):
                continue
            Av = A[v]
            for u in range(n):
                if (mask >> u) & 1:
                    continue
                if Av[u]:
                    dp[mask | (1 << u)][u] += cnt
    return sum(dp[full][v] for v in range(n))


def analyze(A, n, sig):
    cyc = cycles_by_len(A, n, n)
    c = {k: len(cyc[k]) for k in range(3, n + 1)}
    rec = dict(sig=sig, H=H_dp(A, n))
    for k in (6, 7, 8, 9, 10, 11):
        rec[f"c{k}"] = c.get(k, 0)
    # disjoint odd-cycle pair counts (a2 carriers), present at this n
    def D(l, lp):
        if l + lp > n:
            return 0
        return disjoint_pairs(cyc[l], cyc[lp] if l != lp else None, l == lp)
    rec["D33"] = D(3, 3); rec["D35"] = D(3, 5); rec["D37"] = D(3, 7)
    rec["D55"] = D(5, 5); rec["D57"] = D(5, 7); rec["D39"] = D(3, 9)
    # a3 triples
    rec["T333"] = disjoint_triples_same(cyc[3]) if n >= 9 else 0
    rec["T335"] = disjoint_T335(cyc[3], cyc[5]) if n >= 11 else 0
    # trace-basis overlap carriers (for the n=9 reconciliation)
    rec["TF"] = overlap_pairs(cyc[3], cyc[5], False)
    tr = sig  # sig[k-3] = tr A^k
    def trk(k): return tr[k - 3]
    W8 = (trk(8) - trk(4)) // 8 if n >= 8 else 0
    rec["Q44"] = W8 - rec["c8"] - rec["TF"] if n >= 8 else 0
    # OCF check
    a1 = c.get(3, 0) + c.get(5, 0) + c.get(7, 0) + c.get(9, 0) + c.get(11, 0)
    a2 = rec["D33"] + rec["D35"] + rec["D37"] + rec["D55"] + rec["D57"] + rec["D39"]
    a3 = rec["T333"] + rec["T335"]
    rec["a1"], rec["a2"], rec["a3"] = a1, a2, a3
    rec["ok_ocf"] = (rec["H"] == 1 + 2 * a1 + 4 * a2 + 8 * a3)
    return rec

# ---------- exact rank over Q ----------

def matrix_rank_Q(rows, cols):
    M = [[Fraction(int(r[c])) for c in cols] for r in rows]
    R = len(M); C = len(cols)
    pr = 0
    for col in range(C):
        piv = None
        for r in range(pr, R):
            if M[r][col] != 0:
                piv = r; break
        if piv is None:
            continue
        M[pr], M[piv] = M[piv], M[pr]
        pv = M[pr][col]
        M[pr] = [x / pv for x in M[pr]]
        for r in range(R):
            if r != pr and M[r][col] != 0:
                f = M[r][col]
                M[r] = [M[r][k] - f * M[pr][k] for k in range(C)]
        pr += 1
        if pr == R:
            break
    return pr


def within_class_deltas(classes, cols):
    """All within-cospectral-class delta vectors over the given columns + H."""
    rows = []
    for sig, recs in classes.items():
        if len(recs) < 2:
            continue
        base = recs[0]
        for r in recs[1:]:
            d = {c: r[c] - base[c] for c in cols}
            d["H"] = r["H"] - base["H"]
            if any(d[c] != 0 for c in cols) or d["H"] != 0:
                rows.append(d)
    return rows

# ---------- partition growth law ----------

def partitions_odd_ge3_upto(N):
    """sum_{s<=N} #partitions of s into odd parts >=3."""
    # dp[s] = #partitions of s into parts from {3,5,7,...}
    dp = [0] * (N + 1)
    dp[0] = 1
    k = 3
    while k <= N:
        for s in range(k, N + 1):
            dp[s] += dp[s - k]
        k += 2
    return dp, sum(dp[:N + 1])

# ---------- per-n driver ----------

def run_n(n, samples, p2cap, seed):
    print(f"\n{'='*70}\n  n = {n}\n{'='*70}", flush=True)
    # OCF non-spectral carrier set at this n
    ocf = []
    if n >= 7: ocf.append("c7")
    if n >= 9: ocf.append("c9")
    if n >= 11: ocf.append("c11")
    if n >= 6: ocf.append("D33")
    if n >= 8: ocf.append("D35")
    if n >= 10: ocf += ["D37", "D55"]
    if n >= 12: ocf += ["D39", "D57"]
    if n >= 9: ocf.append("T333")
    if n >= 11: ocf.append("T335")
    print(f"  OCF non-spectral carriers: {ocf}  (predicted dim = {len(ocf)})", flush=True)

    # PHASE 1: bucket sigs
    sig_codes = defaultdict(list)
    cnt = 0
    for sig, code in sample_sigs(n, samples, 4000, seed):
        sig_codes[sig].append(code)
        cnt += 1
    colliding = {s: list(dict.fromkeys(cs)) for s, cs in sig_codes.items()
                 if len(set(cs)) >= 2}
    members = sum(len(cs) for cs in colliding.values())
    print(f"  phase1: {cnt} samples; colliding cospectral classes {len(colliding)}; "
          f"distinct members {members}", flush=True)

    # PHASE 2: carriers on colliding members (capped)
    classes = defaultdict(list)
    processed = bad = 0
    for sig in sorted(colliding, key=lambda s: -len(colliding[s])):
        for code in colliding[sig]:
            if processed >= p2cap:
                break
            r = analyze(code_to_adj(code, n), n, sig)
            if not r["ok_ocf"]:
                bad += 1
            classes[sig].append(r)
            processed += 1
        if processed >= p2cap:
            break
    classes = {s: rs for s, rs in classes.items() if len(rs) >= 2}
    print(f"  phase2: processed {processed}; OCF holds {processed-bad}/{processed}; "
          f"usable classes {len(classes)}", flush=True)
    nsplit = sum(1 for rs in classes.values() if len({r['H'] for r in rs}) >= 2)
    print(f"  cospectral classes with SPLIT H (non-spectral witnesses): {nsplit}",
          flush=True)
    if not classes:
        print("  no collisions -- increase samples.", flush=True)
        return None

    # INTRINSIC dimension = rank of OCF-carrier delta matrix
    drows = within_class_deltas(classes, ocf)
    print(f"  within-class delta vectors: {len(drows)}", flush=True)
    rk = matrix_rank_Q(drows, ocf)
    rkH = matrix_rank_Q(drows, ocf + ["H"])
    print(f"  RANK of OCF-carrier deltas             : {rk}", flush=True)
    print(f"  RANK including H column                : {rkH}  "
          f"({'H IN carrier span (linear)' if rkH == rk else 'H NEEDS extra carrier!'})",
          flush=True)
    print(f"  ==> intrinsic dim_nonspec(H) at n={n}  : {rk}", flush=True)

    # which OCF carriers are independent (drop-one rank drop)?
    for c in ocf:
        sub = [x for x in ocf if x != c]
        rsub = matrix_rank_Q(drows, sub)
        print(f"     drop {c:5s}: rank {rsub:2d}  ({'independent' if rsub < rk else 'REDUNDANT'})",
              flush=True)

    # n=9 reconciliation: does c8 enter only via c8+Q44 ?
    if n == 9:
        rk_trace6 = matrix_rank_Q(drows := within_class_deltas(classes,
                       ["c6", "c7", "c8", "c9", "Q44", "T333"]),
                       ["c6", "c7", "c8", "c9", "Q44", "T333"])
        # is the SUM (c8+Q44) the only combination of c8,Q44 that H sees?
        for r in classes:
            pass
        # build augmented col c8+Q44
        d2 = within_class_deltas(classes, ["c6", "c7", "c8", "c9", "Q44", "T333"])
        for r in d2:
            r["c8Q44"] = r["c8"] + r["Q44"]
        rk_sum = matrix_rank_Q(d2, ["c6", "c7", "c8Q44", "c9", "T333"])
        rkH_sum = matrix_rank_Q(d2, ["c6", "c7", "c8Q44", "c9", "T333", "H"])
        print(f"\n  RECONCILE TRACE BASIS vs OCF BASIS at n=9:", flush=True)
        print(f"    rank of trace-basis carriers {{c6,c7,c8,c9,Q44,T333}} = {rk_trace6}",
              flush=True)
        print(f"    rank of {{c6,c7,(c8+Q44),c9,T333}}                    = {rk_sum}",
              flush=True)
        print(f"    rank including H of latter                          = {rkH_sum}  "
              f"({'H in span -> c8,Q44 enter ONLY via their SUM' if rkH_sum==rk_sum else 'no'})",
              flush=True)
        print(f"    => trace-basis chain counts 6; INTRINSIC rank is {rk} "
              f"(c8 & Q44 are ONE non-spectral d.o.f. = D35)", flush=True)
    return rk


def main():
    N10_P1 = int(sys.argv[1]) if len(sys.argv) > 1 else 1500000
    N10_CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    N89 = int(sys.argv[3]) if len(sys.argv) > 3 else 120000
    seed = 20260615

    print("THM-505 : intrinsic non-spectral dimension of H, the GROWTH LAW", flush=True)

    N11_P1 = int(sys.argv[4]) if len(sys.argv) > 4 else 0   # 0 => skip n=11
    N11_CAP = int(sys.argv[5]) if len(sys.argv) > 5 else 4000

    dims = {}
    dims[8] = run_n(8, N89, 6000, seed)
    dims[9] = run_n(9, N89, 6000, seed + 1)
    dims[10] = run_n(10, N10_P1, N10_CAP, seed + 2)
    if N11_P1 > 0:
        dims[11] = run_n(11, N11_P1, N11_CAP, seed + 3)

    # growth-law check
    print(f"\n{'='*70}\n  GROWTH LAW\n{'='*70}", flush=True)
    dp, _ = partitions_odd_ge3_upto(14)
    print("  partitions of s into odd parts >=3:", flush=True)
    for s in range(0, 15):
        print(f"    s={s:2d}: {dp[s]}", flush=True)
    print("\n  n :  intrinsic dim   |  predicted = (sum_{s<=n} p_odd>=3(s)) - 3", flush=True)
    for nn in range(6, 15):
        _, cum = partitions_odd_ge3_upto(nn)
        pred = cum - 3
        meas = dims.get(nn, None)
        ms = str(meas) if meas is not None else " ? "
        print(f"  {nn:2d}:      {ms:>3s}          |   {pred}", flush=True)


if __name__ == "__main__":
    main()
