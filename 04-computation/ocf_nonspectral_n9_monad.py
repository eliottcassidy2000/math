#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-505 PUSHED TO n=9 + THE LINEARITY QUESTION  (monad-explorer-2026-06-15)

Builds on THM-505 (the OCF non-spectral defect). Three experiments:

  lin8 : Is H LINEAR in the simple-cycle carriers at n=8?  THM-505 gives
         H = skel + 2c7 + 4c6 + 4c8 + 4*Q44, and Q44 is DETERMINED by
         (spectrum, c6,c7,c8) -- but possibly NONLINEARLY.  Linearity of H in
         (c6,c7,c8) <=> exists UNIVERSAL integers (b6,b7,b8) with the
         within-cospectral-class law  dQ44 = b6*dc6 + b7*dc7 + b8*dc8.
         (Within a class all spectral invariants are fixed, so this kills the
         spectral part automatically.)  We collect every within-class delta and
         solve for a single (b6,b7,b8); exactness on ALL deltas => H linear.

  form9: Verify the n=9 closed form derived by substitution from the OCF
         (H = 1 + 2a1 + 4a2 + 8a3, truncates: a4=0 needs >=12 vertices) and
         THM-502 census defects.  At n=9 the TRIPLE level a3 first switches on
         (3 disjoint triangles = 9 vertices).  Derivation:
            a1 = c3+c5+c7+c9
            a2 = D33 + D35       (disjoint odd pairs, total<=9: (3,3),(3,5))
                 D33 = C(c3,2) - p33,   p33 = W6 - c6
                 D35 = c3*c5    - TF ,   TF  = W8 - c8 - Q44
            a3 = T333            (disjoint triangle TRIPLES, the new carrier)
         =>  H = [1 + 2c3 + 2c5 + 4C(c3,2) + 4c3c5 - 4W6 - 4W8]   (SKELETON)
                 + 2c7 + 2c9 + 4c6 + 4c8 + 4Q44 + 8*T333.
         The coefficient 8 = 2^3 on T333 is the fugacity-2 weight of the TRIPLE
         independent-set level -- the next rung of the "2^level" rule.

  dim9 : Non-spectral dimension at n=9.  Conjecture (HYP-2513): dim = n-5 = 4,
         carriers (c6,c7,c8,c9), with Q44 and T333 SPECTRALLY DEPENDENT on them.
         Probe by trace-vector buckets:
           - does (spectrum,c6,c7,c8) STILL split H?  (=> dim >= 4, c9 needed)
           - does (spectrum,c6,c7,c8,c9) determine H? (=> dim <= 4)
           - are Q44, T333 determined by (spectrum,c6,c7,c8,c9)?
         Also: within-class linearity test dH = lin(dc6,dc7,dc8,dc9).
"""

import sys, random, itertools
from collections import defaultdict
from math import comb
from fractions import Fraction

sys.setrecursionlimit(100000)


def random_tournament(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.getrandbits(1):
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def matmul(X, Y, n):
    Z = [[0] * n for _ in range(n)]
    for i in range(n):
        Xi = X[i]; Zi = Z[i]
        for k in range(n):
            a = Xi[k]
            if a:
                Yk = Y[k]
                for j in range(n):
                    Zi[j] += a * Yk[j]
    return Z


def traces(A, n, kmax):
    """Return dict k -> tr(A^k) for k=1..kmax."""
    out = {}
    P = [row[:] for row in A]
    out[1] = sum(P[i][i] for i in range(n))
    for k in range(2, kmax + 1):
        P = matmul(P, A, n)
        out[k] = sum(P[i][i] for i in range(n))
    return out


def cycles_by_len(A, n, maxL):
    """All simple directed cycles, grouped by length; each as a vertex bitmask."""
    adj = [[j for j in range(n) if A[i][j]] for i in range(n)]
    out = defaultdict(list)
    for start in range(n):
        smask = 1 << start
        stack = [(start, smask, 1)]   # (vertex, visited-mask, length)
        while stack:
            v, vis, ln = stack.pop()
            for w in adj[v]:
                if w == start and ln >= 3:
                    out[ln].append(vis)
                elif w > start and not (vis >> w) & 1 and ln < maxL:
                    stack.append((w, vis | (1 << w), ln + 1))
    return out


def overlap_pairs(la, lb, same):
    """count pairs sharing >=1 vertex (bitmask AND != 0)."""
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
    """count unordered triples of PAIRWISE-disjoint sets (bitmask)."""
    L = len(la)
    cnt = 0
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


# ---------- per-tournament analyzers ----------

def analyze8(A):
    n = 8
    cyc = cycles_by_len(A, n, 8)
    c3 = len(cyc[3]); c4 = len(cyc[4]); c5 = len(cyc[5])
    c6 = len(cyc[6]); c7 = len(cyc[7]); c8 = len(cyc[8])
    Q44 = overlap_pairs(cyc[4], None, True)
    TF = overlap_pairs(cyc[3], cyc[5], False)
    D33 = disjoint_pairs(cyc[3], None, True)
    D35 = disjoint_pairs(cyc[3], cyc[5], False)
    a2 = D33 + D35
    tr = traces(A, n, 8)
    W6 = (tr[6] - tr[3]) // 6
    W8 = (tr[8] - tr[4]) // 8
    H = H_dp(A, n)
    skel = 1 + 2 * c3 + 2 * c5 + 4 * comb(c3, 2) + 4 * c3 * c5 - 4 * W6 - 4 * W8
    ok_ocf = (H == 1 + 2 * (c3 + c5 + c7) + 4 * a2)
    ok_closed = (H == skel + 2 * c7 + 4 * c6 + 4 * c8 + 4 * Q44)
    sig = tuple(tr[k] for k in range(3, 9))
    return dict(c6=c6, c7=c7, c8=c8, Q44=Q44, TF=TF, H=H, sig=sig,
                ok=(ok_ocf, ok_closed))


def analyze9(A):
    n = 9
    cyc = cycles_by_len(A, n, 9)
    c3 = len(cyc[3]); c4 = len(cyc[4]); c5 = len(cyc[5])
    c6 = len(cyc[6]); c7 = len(cyc[7]); c8 = len(cyc[8]); c9 = len(cyc[9])
    p33 = overlap_pairs(cyc[3], None, True)
    Q44 = overlap_pairs(cyc[4], None, True)
    TF = overlap_pairs(cyc[3], cyc[5], False)
    D33 = disjoint_pairs(cyc[3], None, True)
    D35 = disjoint_pairs(cyc[3], cyc[5], False)
    a2 = D33 + D35
    T333 = disjoint_triples_same(cyc[3])   # alpha_3 at n=9
    a3 = T333
    tr = traces(A, n, 9)
    W6 = (tr[6] - tr[3]) // 6
    W8 = (tr[8] - tr[4]) // 8
    H = H_dp(A, n)
    a1 = c3 + c5 + c7 + c9
    skel = 1 + 2 * c3 + 2 * c5 + 4 * comb(c3, 2) + 4 * c3 * c5 - 4 * W6 - 4 * W8
    closed = skel + 2 * c7 + 2 * c9 + 4 * c6 + 4 * c8 + 4 * Q44 + 8 * T333
    ok_ocf = (H == 1 + 2 * a1 + 4 * a2 + 8 * a3)
    ok_p33 = (p33 == W6 - c6)
    ok_TF = (TF == W8 - c8 - Q44)
    ok_closed = (H == closed)
    sig = tuple(tr[k] for k in range(3, 10))
    return dict(c6=c6, c7=c7, c8=c8, c9=c9, Q44=Q44, TF=TF, T333=T333,
                H=H, a1=a1, a2=a2, a3=a3, sig=sig,
                ok=(ok_ocf, ok_p33, ok_TF, ok_closed))


# ---------- exact integer linear fit ----------

def exact_linear_fit(rows, xcols, ycol):
    """Least-squares-then-verify-exact.  rows: list of dict.  Solve y = sum a_i x_i
    (no intercept; deltas are already centered).  Returns (sol, exact_bool)."""
    K = len(xcols)
    M = [[Fraction(int(d[c])) for c in xcols] for d in rows]
    y = [Fraction(int(d[ycol])) for d in rows]
    N = len(rows)
    # normal equations
    Mt = [[sum(M[r][a] * M[r][b] for r in range(N)) for b in range(K)] for a in range(K)]
    Mty = [sum(M[r][a] * y[r] for r in range(N)) for a in range(K)]
    sol = gauss(Mt, Mty)
    if sol is None:
        return None, False
    exact = all(sum(M[r][a] * sol[a] for a in range(K)) == y[r] for r in range(N))
    return sol, exact


def gauss(Am, bv):
    n = len(Am)
    M = [[Am[i][j] for j in range(n)] + [bv[i]] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] != 0), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                fac = M[r][col]
                M[r] = [M[r][j] - fac * M[col][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]


def collect_within_class_deltas(classes, carriers, target):
    """classes: sig -> list of dicts. For each class, form deltas vs first rep.
    Returns list of delta-dicts (carriers + target)."""
    deltas = []
    for sig, recs in classes.items():
        # dedup identical carrier+target tuples to avoid degenerate rows
        seen = {}
        for r in recs:
            key = tuple(r[c] for c in carriers) + (r[target],)
            seen[key] = r
        reps = list(seen.values())
        if len(reps) < 2:
            continue
        base = reps[0]
        for r in reps[1:]:
            d = {c: r[c] - base[c] for c in carriers}
            d[target] = r[target] - base[target]
            deltas.append(d)
    return deltas


# ---------- main experiments ----------

def add_quadratic_monomials(rows, base_cols):
    """augment each row with squares and cross-products of base_cols."""
    for r in rows:
        for i, a in enumerate(base_cols):
            r[a + "^2"] = r[a] * r[a]
            for b in base_cols[i + 1:]:
                r[a + "*" + b] = r[a] * r[b]
    quad = []
    for i, a in enumerate(base_cols):
        quad.append(a + "^2")
        for b in base_cols[i + 1:]:
            quad.append(a + "*" + b)
    return base_cols + quad


def run_lin8(NS, rng):
    print(f"=== n=8 LINEARITY TEST  ({NS} random tournaments) ===", flush=True)
    classes = defaultdict(list)
    bad = [0, 0]
    for _ in range(NS):
        r = analyze8(random_tournament(8, rng))
        for i, b in enumerate(r["ok"]):
            if not b:
                bad[i] += 1
        classes[r["sig"]].append(r)
    print(f"  [OCF]    H=1+2a1+4a2          : {NS-bad[0]}/{NS}", flush=True)
    print(f"  [closed] H=skel+2c7+4c6+4c8+4Q44: {NS-bad[1]}/{NS}", flush=True)
    split = {s: rs for s, rs in classes.items()
             if len({r["H"] for r in rs}) >= 2}
    print(f"  cospectral classes: {len(classes)}; with split H: {len(split)}", flush=True)

    # (0) confirm THM-505: does (sig,c6,c7,c8) DETERMINE H?
    det_buckets = defaultdict(set)
    for sig, rs in classes.items():
        for r in rs:
            det_buckets[(sig, r["c6"], r["c7"], r["c8"])].add(r["H"])
    free = sum(1 for v in det_buckets.values() if len(v) >= 2)
    print(f"\n  (sig,c6,c7,c8) buckets still splitting H: {free}  "
          f"(0 => H is a FUNCTION of (spectrum,c6,c7,c8))", flush=True)

    # (1) LINEAR within-class fit:  dH = w6 dc6 + w7 dc7 + w8 dc8 ?
    dH = collect_within_class_deltas(classes, ["c6", "c7", "c8"], "H")
    solH, exactH = exact_linear_fit(dH, ["c6", "c7", "c8"], "H")
    print(f"\n  (1) LINEAR  dH = w6 dc6 + w7 dc7 + w8 dc8 :", flush=True)
    print(f"      weights = {tuple(solH) if solH else None}   EXACT = {exactH}", flush=True)
    if not exactH and solH:
        worst = max(abs(d["H"] - sum(solH[i] * d[c]
                    for i, c in enumerate(["c6", "c7", "c8"]))) for d in dH)
        print(f"      NONLINEAR: max |dH - linear pred| = {worst}", flush=True)

    # (2) QUADRATIC within-class fit (add squares + cross terms)
    base = ["c6", "c7", "c8"]
    qcols = base + [a + "^2" for a in base] + \
        [base[i] + "*" + base[j] for i in range(len(base)) for j in range(i + 1, len(base))]
    dHq = collect_within_class_deltas_poly(classes, base, "H")
    solQ, exactQ = exact_linear_fit(dHq, qcols, "H")
    print(f"\n  (2) QUADRATIC  dH = lin(c6,c7,c8) + quad(c6^2,..,c7*c8) :", flush=True)
    print(f"      EXACT (universal degree-2 polynomial in carriers) = {exactQ}", flush=True)

    # (3) PER-CLASS linear fit (coefficients allowed to be spectral/class-dependent)
    nclass = 0; nlin = 0; nrank = defaultdict(int)
    for sig, recs in classes.items():
        seen = {}
        for r in recs:
            seen[(r["c6"], r["c7"], r["c8"], r["Q44"])] = r
        reps = list(seen.values())
        if len(reps) < 2:
            continue
        nclass += 1
        base_r = reps[0]
        ds = [{c: r[c] - base_r[c] for c in ["c6", "c7", "c8", "Q44"]} for r in reps[1:]]
        # is Q44-delta exactly linear in (c6,c7,c8)-deltas, per this class?
        if len(ds) >= 1:
            solc, exc = exact_linear_fit(ds, ["c6", "c7", "c8"], "Q44") if len(ds) >= 3 \
                else (None, None)
            if exc:
                nlin += 1
    print(f"\n  (3) PER-CLASS Q44~lin(c6,c7,c8): classes with >=4 distinct carrier "
          f"tuples fit exactly: {nlin}/{nclass}", flush=True)

    # (4) THE UNIVERSAL-LINEAR forms (THM-505) DO hold -- verify weights universal
    dHt = collect_within_class_deltas(classes, ["c6", "c7", "TF"], "H")
    solT, exactT = exact_linear_fit(dHt, ["c6", "c7", "TF"], "H")
    print(f"\n  (4) UNIVERSAL-LINEAR in (c6,c7,TF):  dH = a*dc6+b*dc7+d*dTF", flush=True)
    print(f"      weights = {tuple(solT) if solT else None}  EXACT = {exactT}  "
          f"(expect (4,2,-4))", flush=True)
    dHq2 = collect_within_class_deltas(classes, ["c6", "c7", "c8", "Q44"], "H")
    solQ2, exactQ2 = exact_linear_fit(dHq2, ["c6", "c7", "c8", "Q44"], "H")
    print(f"      and in (c6,c7,c8,Q44): weights = {tuple(solQ2) if solQ2 else None}  "
          f"EXACT = {exactQ2}  (expect (4,2,4,4))", flush=True)


def collect_within_class_deltas_poly(classes, base_cols, target):
    """like collect_within_class_deltas but augments each rec with quadratic
    monomials BEFORE differencing (so deltas of c6^2 etc. are exact)."""
    deltas = []
    qcols = None
    for sig, recs in classes.items():
        aug = []
        for r in recs:
            d = {c: r[c] for c in base_cols}
            d[target] = r[target]
            for i, a in enumerate(base_cols):
                d[a + "^2"] = r[a] * r[a]
                for b in base_cols[i + 1:]:
                    d[a + "*" + b] = r[a] * r[b]
            aug.append(d)
        # dedup by base-carrier tuple + target
        seen = {}
        for d in aug:
            key = tuple(d[c] for c in base_cols) + (d[target],)
            seen[key] = d
        reps = list(seen.values())
        if len(reps) < 2:
            continue
        base = reps[0]
        for d in reps[1:]:
            delta = {k: d[k] - base[k] for k in d}
            deltas.append(delta)
    return deltas


def run_form9(NS, rng):
    print(f"=== n=9 CLOSED-FORM VERIFICATION  ({NS} random tournaments) ===", flush=True)
    print("  H = [1+2c3+2c5+4C(c3,2)+4c3c5-4W6-4W8] + 2c7+2c9+4c6+4c8+4Q44 + 8*T333", flush=True)
    bad = [0, 0, 0, 0]
    names = ["OCF H=1+2a1+4a2+8a3", "p33=W6-c6", "TF=W8-c8-Q44", "CLOSED FORM (n=9)"]
    a3_pos = 0
    for _ in range(NS):
        r = analyze9(random_tournament(9, rng))
        for i, b in enumerate(r["ok"]):
            if not b:
                bad[i] += 1
        if r["a3"] > 0:
            a3_pos += 1
    for nm, b in zip(names, bad):
        print(f"  [{'OK ' if b == 0 else 'FAIL'}] {nm:24s}: {NS-b}/{NS}", flush=True)
    print(f"  tournaments with a3>0 (triple-overlap level active): {a3_pos}/{NS}", flush=True)


def run_dim9(NS, rng):
    print(f"=== n=9 NON-SPECTRAL DIMENSION PROBE  ({NS} random tournaments) ===", flush=True)
    classes = defaultdict(list)
    bad = 0
    for _ in range(NS):
        r = analyze9(random_tournament(9, rng))
        if not r["ok"][3]:
            bad += 1
        classes[r["sig"]].append(r)
    print(f"  closed form holds: {NS-bad}/{NS}", flush=True)
    split = {s: rs for s, rs in classes.items()
             if len({r["H"] for r in rs}) >= 2}
    print(f"  cospectral classes sampled: {len(classes)}; with split H: {len(split)}", flush=True)
    if not split:
        print("  (no cospectral H-splits found at this sample size; increase NS)", flush=True)
        return

    # bucket tests
    def bucket_split(carrier_keys, collect_examples=False):
        buckets = defaultdict(list)
        for sig, rs in classes.items():
            for r in rs:
                key = (sig,) + tuple(r[c] for c in carrier_keys)
                buckets[key].append(r)
        nsplit = 0
        examples = []
        for key, rs in buckets.items():
            if len({r["H"] for r in rs}) >= 2:
                nsplit += 1
                if collect_examples and len(examples) < 6:
                    examples.append((key, sorted({(r["H"], r["c9"], r["Q44"], r["T333"])
                                                  for r in rs})))
        return nsplit, examples

    s678, ex678 = bucket_split(["c6", "c7", "c8"], collect_examples=True)
    s6789, _ = bucket_split(["c6", "c7", "c8", "c9"])
    print(f"\n  (sig,c6,c7,c8)    buckets still splitting H : {s678}  "
          f"(>0 => dim>=4, c9 genuinely needed)", flush=True)
    print(f"  (sig,c6,c7,c8,c9) buckets still splitting H : {s6789}  "
          f"(0 => dim<=4, carriers (c6..c9) determine H)", flush=True)
    print(f"  --- the (sig,c6,c7,c8)-splitting buckets, showing (H,c9,Q44,T333): ---", flush=True)
    for key, vals in ex678:
        sig, c6, c7, c8 = key
        print(f"     c6={c6} c7={c7} c8={c8}: {vals}", flush=True)

    # are Q44, T333 determined by (sig,c6,c7,c8,c9)?
    for tgt in ["Q44", "T333"]:
        b = defaultdict(set)
        for sig, rs in classes.items():
            for r in rs:
                b[(sig, r["c6"], r["c7"], r["c8"], r["c9"])].add(r[tgt])
        free = sum(1 for v in b.values() if len(v) >= 2)
        print(f"  (sig,c6,c7,c8,c9) buckets where {tgt:5s} still varies: {free}  "
              f"(0 => {tgt} spectrally dependent, NOT a new axis)", flush=True)

    # within-class linearity: dH linear in (dc6,dc7,dc8,dc9)?
    dH = collect_within_class_deltas(classes, ["c6", "c7", "c8", "c9"], "H")
    print(f"\n  within-class delta-vectors: {len(dH)}", flush=True)
    if len(dH) >= 4:
        solH, exactH = exact_linear_fit(dH, ["c6", "c7", "c8", "c9"], "H")
        print(f"  dH = w6 dc6 + w7 dc7 + w8 dc8 + w9 dc9 :", flush=True)
        print(f"       weights = {tuple(solH) if solH else None}  EXACT = {exactH}", flush=True)


def run_chain(n, NS, rng):
    """Nested-carrier dimension probe at general n.  Determines the MINIMAL
    carrier set (added to the full spectrum) that determines H, by checking
    which nested prefix of carriers reduces the H-split count to 0."""
    analyze = analyze8 if n == 8 else analyze9
    print(f"=== n={n} NESTED-CARRIER DIMENSION CHAIN  ({NS} random tournaments) ===",
          flush=True)
    classes = defaultdict(list)
    bad = 0
    for _ in range(NS):
        r = analyze(random_tournament(n, rng))
        if not r["ok"][-1]:
            bad += 1
        classes[r["sig"]].append(r)
    print(f"  closed form holds: {NS-bad}/{NS}", flush=True)
    nclass = len(classes)
    nsplit = sum(1 for rs in classes.values() if len({r['H'] for r in rs}) >= 2)
    print(f"  cospectral classes: {nclass}; with split H: {nsplit}", flush=True)

    # the carrier chain (simple cycles first, then overlap/triple configs)
    if n == 8:
        chain = [[], ["c6"], ["c6", "c7"], ["c6", "c7", "c8"],
                 ["c6", "c7", "c8", "Q44"], ["c6", "c7", "c8", "Q44", "TF"]]
    else:
        chain = [[], ["c6", "c7", "c8"], ["c6", "c7", "c8", "c9"],
                 ["c6", "c7", "c8", "c9", "Q44"],
                 ["c6", "c7", "c8", "c9", "Q44", "T333"]]

    print(f"\n  carrier set                          -> #buckets splitting H", flush=True)
    prev = None
    for carriers in chain:
        buckets = defaultdict(set)
        for sig, rs in classes.items():
            for r in rs:
                buckets[(sig,) + tuple(r[c] for c in carriers)].add(r["H"])
        nsp = sum(1 for v in buckets.values() if len(v) >= 2)
        lbl = "(spectrum only)" if not carriers else "+".join(carriers)
        print(f"    sig,{lbl:34s} -> {nsp}", flush=True)
        prev = (carriers, nsp)

    # explicit witnesses: cospectral + equal simple-cycle vector, different H
    simple = ["c6", "c7", "c8"] if n == 8 else ["c6", "c7", "c8", "c9"]
    buckets = defaultdict(list)
    for sig, rs in classes.items():
        for r in rs:
            buckets[(sig,) + tuple(r[c] for c in simple)].append(r)
    witnesses = []
    for key, rs in buckets.items():
        Hs = {r["H"] for r in rs}
        if len(Hs) >= 2:
            extra = sorted({(r["H"], r["Q44"], r.get("T333", 0), r["TF"]) for r in rs})
            witnesses.append((key[1:], extra))
    print(f"\n  WITNESSES: cospectral + equal ({'+'.join(simple)}) but DIFFERENT H: "
          f"{len(witnesses)} buckets", flush=True)
    for sv, extra in witnesses[:8]:
        print(f"    {dict(zip(simple, sv))}: (H,Q44,T333,TF)={extra}", flush=True)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "form9"
    NS = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    rng = random.Random(20260615)
    if which == "lin8":
        run_lin8(NS, rng)
    elif which == "form9":
        run_form9(NS, rng)
    elif which == "dim9":
        run_dim9(NS, rng)
    elif which == "chain8":
        run_chain(8, NS, rng)
    elif which == "chain9":
        run_chain(9, NS, rng)


if __name__ == "__main__":
    main()
