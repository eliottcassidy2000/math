#!/usr/bin/env python3
"""
paley_R31_prediction_kpc1.py
kind-pasteur-2026-06-10-S1 -- Thread E (HYP-2371):
A FALSIFIABLE H(T_31) PREDICTION FROM THE PROVEN CLUSTER EXPANSION (HYP-2307/THM-438).

Background (proved / verified in monad-explorer-2026-06-07):
  R(p) = H(T_p) * 2^(p-1) / p!  =  (1/p!) sum_sigma prod_k (1 + chi(d_k))
       = sum over slot-subsets S of T(S)/p!,  and S decomposes into maximal runs.
  EXACT identity:  T(S) = A_joint(lengths) * (p - n)!   with n = sum(L_i + 1),
  where A_joint is the joint distinct-vertex character integral of the runs.
  Known exactly: A_2 = p(p-1); A_L = 0 for odd L (negation, chi(-1)=-1);
  A_{2k} = C_k p^{k+1} + O(p^k) (THM-438).  R(p) -> e.

THIS SCRIPT (all number theory in EXACT integer arithmetic):
 (1) Exact engine: distinct-tuple character sums via Moebius inversion over set
     partitions + integer tensor contraction (np.einsum, int64; bounds checked).
     Cross-validated against brute force and against monad-stored A_4/A_6/A_8.
 (2) Finite-p truncated predictors R_pred(p) built from:
       singles A_2,A_4,A_6,A_8 (exact), pair-collision corrections
       delta(2,2), delta(2,4), delta(3,3) (exact joints), cherry-triple
       correction delta3(2,2,2) (exact), with an EXACT placement DP
       (integer weights, exact Fractions; up to 3 disjoint marked pairs).
 (3) OUT-OF-SAMPLE VALIDATION: residuals vs the 5 known exact R(p)
     (p = 3,7,11,19,23), incl. a strict holdout test at p=23.
 (4) Prediction R(31) with an error bar derived from the observed residual
     trend (honest: bar widened if holdout fails / order is unstable).
 (5) Integrality sharpening: H(T_31) = R(31)*31!/2^30; |Aff(QR)| = 465 divides H
     (freeness, THM-F3 / Thread B), H/465 is ODD (orbit parity HYP-1713,
     re-verified here on all 5 known values)  =>  H = 465 * odd, i.e.
     H == 465 (mod 930).  Count admissible integers in the predicted interval.
 (6) Cross-check: direct R(p) = e(1 - C/p - D/p^2) fits (the ansatz THM-438
     ADDENDUM (4) justifies) -- reported, not used as the primary predictor.

Floats are used ONLY for residual analysis, extrapolation and display.
Run:  python3 04-computation/paley_R31_prediction_kpc1.py
"""
import sys, time
import numpy as np
from fractions import Fraction
from math import factorial, e as MATH_E
from itertools import permutations, combinations
from collections import defaultdict

H_KNOWN = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349}
VAL_PRIMES = [3, 7, 11, 19, 23]
P_TARGET = 31
ALL_PRIMES = VAL_PRIMES + [P_TARGET]

T0 = time.time()
def log(msg=""):
    print(msg, flush=True)

# ----------------------------------------------------------------------------
# Legendre symbol tables (exact)
# ----------------------------------------------------------------------------
def legendre(p):
    chi = [0] * p
    qr = set((x * x) % p for x in range(1, p))
    for d in range(1, p):
        chi[d] = 1 if d in qr else -1
    return chi

def chi_matrix(p, chi):
    # C[a,b] = chi[(b-a) mod p], int64
    a = np.arange(p)
    return np.array(chi, dtype=np.int64)[(a[None, :] - a[:, None]) % p]

# ----------------------------------------------------------------------------
# Set partitions (restricted growth strings)
# ----------------------------------------------------------------------------
def rgs_partitions(N):
    a = [0] * N
    b = [1] * N  # b[i] = max(a[:i]) + 1 (allowed ceiling for a[i]); b[0] unused
    while True:
        yield tuple(a)
        j = N - 1
        while j > 0 and a[j] >= b[j]:
            j -= 1
        if j == 0:
            return
        a[j] += 1
        nb = max(b[j], a[j] + 1)
        for i in range(j + 1, N):
            a[i] = 0
            b[i] = nb

BELL = {1:1,2:2,3:5,4:15,5:52,6:203,7:877,8:4140,9:21147,10:115975}

# ----------------------------------------------------------------------------
# Exact distinct-tuple character-sum engine
#   distinct_sum = sum over INJECTIVE (x_0..x_{N-1}) in F_p^N of
#                  prod_{(u,v) in edges} chi(x_v - x_u)
# Moebius over the partition lattice: D(0) = sum_sigma mu(0,sigma) * M_sigma,
# M_sigma = FREE sum with values constant on blocks  (mu = prod (-1)^{s-1}(s-1)!)
# Pruning (THM-438 'vanishing engine'): M_sigma = 0 if any edge is a self-loop
# (chi(0)=0) or any block has multigraph degree 1 (leaf: sum chi = 0).
# ----------------------------------------------------------------------------
LETTERS = 'abcdefghijklmnopqrstuvwxyz'

def precompute_structure(N, edges):
    """Return list of (mu, subscript) for surviving partitions (p-independent)."""
    surv = []
    npart = 0
    for part in rgs_partitions(N):
        npart += 1
        V = max(part) + 1
        deg = [0] * V
        ok = True
        for (u, v) in edges:
            if part[u] == part[v]:
                ok = False
                break
            deg[part[u]] += 1
            deg[part[v]] += 1
        if not ok or min(deg) <= 1:
            continue
        sizes = [0] * V
        for b in part:
            sizes[b] += 1
        mu = 1
        for s in sizes:
            mu *= (-1) ** (s - 1) * factorial(s - 1)
        sub = ','.join(LETTERS[part[u]] + LETTERS[part[v]] for (u, v) in edges) + '->'
        surv.append((mu, sub, V))
    assert npart == BELL[N], (N, npart)
    return surv

def eval_structure(p, C, N, surv):
    """Exact integer distinct-tuple sum at prime p (0 if N > p)."""
    if N > p:
        return 0
    total = 0
    ops_cache = {}
    for (mu, sub, V) in surv:
        nops = sub.count(',') + 1
        M = int(np.einsum(sub, *([C] * nops), optimize='greedy'))
        # overflow guard: |M| <= p^V << 2^63 for p<=31, V<=9
        assert abs(M) <= p ** V, "einsum overflow guard tripped"
        total += mu * M
    return total

def brute_distinct(p, N, edges, chi):
    tot = 0
    for tup in permutations(range(p), N):
        w = 1
        for (u, v) in edges:
            w *= chi[(tup[v] - tup[u]) % p]
            if w == 0:
                break
        tot += w
    return tot

def path_edges(L, off=0):
    return [(off + i, off + i + 1) for i in range(L)]

# structures: name -> (N, edges)
STRUCTS = {
    'A2':  (3, path_edges(2)),
    'A4':  (5, path_edges(4)),
    'A6':  (7, path_edges(6)),
    'A8':  (9, path_edges(8)),
    'J22': (6, path_edges(2) + path_edges(2, 3)),
    'J24': (8, path_edges(2) + path_edges(4, 3)),
    'J33': (8, path_edges(3) + path_edges(3, 4)),
    'J222': (9, path_edges(2) + path_edges(2, 3) + path_edges(2, 6)),
    # zero-checks (negation symmetry: total edge count odd => exact 0)
    'A5':  (6, path_edges(5)),
    'J23': (7, path_edges(2) + path_edges(3, 3)),
}

log("=" * 78)
log("PART 0: exact engine -- build partition tables and self-test")
log("=" * 78)
SURV = {}
for name, (N, edges) in STRUCTS.items():
    t = time.time()
    SURV[name] = precompute_structure(N, edges)
    log(f"  {name}: N={N} Bell={BELL[N]:>6} surviving-patterns={len(SURV[name]):>5}"
        f"  ({time.time()-t:.1f}s)")

# self-tests vs brute force at p=7 (and known monad values below)
chi7 = legendre(7)
C7 = chi_matrix(7, chi7)
for name in ['A2', 'A4', 'J22', 'A5']:
    N, edges = STRUCTS[name]
    bf = brute_distinct(7, N, edges, chi7)
    en = eval_structure(7, C7, N, SURV[name])
    assert bf == en, (name, bf, en)
    log(f"  self-test p=7 {name}: brute={bf} engine={en}  OK")
assert eval_structure(7, C7, 3, SURV['A2']) == 7 * 6

log("")
log("=" * 78)
log("PART 1: exact cluster integrals at p = 3,7,11,19,23,31 (exact integers)")
log("=" * 78)
INT = {}  # INT[p][name] = exact integer
for p in ALL_PRIMES:
    chi = legendre(p)
    C = chi_matrix(p, chi)
    INT[p] = {}
    for name, (N, edges) in STRUCTS.items():
        t = time.time()
        INT[p][name] = eval_structure(p, C, N, SURV[name])
        dt = time.time() - t
        if dt > 5:
            log(f"    [{name} @ p={p}: {dt:.1f}s]")
    # exact-law checks
    assert INT[p]['A2'] == p * (p - 1) if p >= 3 else True
    assert INT[p]['A5'] == 0 and INT[p]['J23'] == 0, "negation-symmetry zero failed"
log(f"  (elapsed {time.time()-T0:.1f}s)")

# cross-check vs monad-stored exact values
MONAD_A4 = {7: 336, 11: 1760, 19: 10944, 23: 20240, 31: 52080}
MONAD_A6 = {7: 1008, 11: 22880, 19: 361152, 23: 870320}
MONAD_A8 = {7: 0, 11: 154880}
for p, v in MONAD_A4.items():
    assert INT[p]['A4'] == v, ('A4', p, INT[p]['A4'], v)
for p, v in MONAD_A6.items():
    assert INT[p]['A6'] == v, ('A6', p, INT[p]['A6'], v)
for p, v in MONAD_A8.items():
    assert INT[p]['A8'] == v, ('A8', p, INT[p]['A8'], v)
log("  cross-checks vs monad-stored A_4 (p<=31), A_6 (p<=23), A_8 (p<=11): ALL OK")
log("  negation-symmetry exact zeros A_5 = J(2,3) = 0: OK at all primes")
log("")
hdr = f"{'p':>3} | {'A_2':>6} {'A_4':>8} {'A_6':>10} {'A_8':>13} | {'d22':>12} {'d24':>14} {'d33':>12} {'d3(222)':>14}"
log(hdr)
log("-" * len(hdr))
DELTA = {}
for p in ALL_PRIMES:
    I = INT[p]
    d22 = I['J22'] - I['A2'] * I['A2']
    d24 = I['J24'] - I['A2'] * I['A4']
    d33 = I['J33'] - 0                     # A_3 = 0 exactly
    d3 = I['J222'] - I['A2'] ** 3 - 3 * d22 * I['A2']
    DELTA[p] = {'d22': d22, 'd24': d24, 'd33': d33, 'd3': d3}
    log(f"{p:>3} | {I['A2']:>6} {I['A4']:>8} {I['A6']:>10} {I['A8']:>13} | "
        f"{d22:>12} {d24:>14} {d33:>12} {d3:>14}")
log("")
log("scaled (orders):   d22/p^3      d24/p^4      d33/p^4      d3/p^4   (expected O(1))")
for p in ALL_PRIMES[1:]:
    D = DELTA[p]
    log(f"{p:>3} | {D['d22']/p**3:>12.5f} {D['d24']/p**4:>12.5f} "
        f"{D['d33']/p**4:>12.5f} {D['d3']/p**4:>12.5f}")

# ----------------------------------------------------------------------------
# PART 2: exact placement DP
#   R_pred(p) = sum over run-configs (disjoint maximal runs in slots 1..p-1)
#     [ prod singles  +  sum_{matchings of <=3 disjoint marked pairs} prod delta
#       * prod singles(rest)  +  sum_{marked cherry triples} delta3 * prod
#       singles(rest) ] * (p-n)!/p!        -- all exact integers / Fractions.
# ----------------------------------------------------------------------------
def placement_pass(p, singles, pair_delta=None, triple=False, d3=0, max_pairs=3):
    """One DP pass. Returns dict n -> exact integer weight.
       pair_delta: dict (l1<=l2) -> delta  (pass collects only >=1 completed pair)
       triple: collect only configs with exactly one marked cherry-triple."""
    pair_delta = pair_delta or {}
    open_ls = sorted(set(l for k in pair_delta for l in k))
    lengths = sorted(set(list(singles.keys()) + open_ls + ([2] if triple else [])))
    # state: (n, pend tuple sorted, donepairs, tcnt) at slot pointer j
    layers = defaultdict(lambda: defaultdict(int))
    layers[1][(0, (), 0, 0)] = 1
    out = defaultdict(int)
    for j in range(1, p + 1):
        cur = layers.pop(j, None)
        if cur is None:
            continue
        if j == p:
            for (n, pend, dp_, tc), w in cur.items():
                if pend:
                    continue
                if pair_delta and dp_ == 0:
                    continue
                if triple and tc != 3:
                    continue
                if not pair_delta and not triple and (dp_ or tc):
                    continue
                out[n] += w
            continue
        for (n, pend, dp_, tc), w in cur.items():
            # skip slot j
            layers[j + 1][(n, pend, dp_, tc)] += w
            # place a run of length l starting at slot j
            for l in lengths:
                if j + l - 1 > p - 1:
                    continue
                jn = min(j + l + 1, p)  # forced gap slot j+l
                nn = n + l + 1
                if nn > p:
                    continue
                # plain run
                if l in singles and singles[l] != 0:
                    layers[jn][(nn, pend, dp_, tc)] += w * singles[l]
                # pair-opener
                if pair_delta and len(pend) < 2 and l in open_ls and dp_ < max_pairs:
                    layers[jn][(nn, tuple(sorted(pend + (l,))), dp_, tc)] += w
                # pair-closer (choose which pending partner; multiplicity)
                if pair_delta and pend:
                    for l1 in set(pend):
                        key = (min(l1, l), max(l1, l))
                        if key in pair_delta and dp_ < max_pairs:
                            mult = pend.count(l1)
                            np_ = list(pend)
                            np_.remove(l1)
                            layers[jn][(nn, tuple(np_), dp_ + 1, tc)] += \
                                w * mult * pair_delta[key]
                # triple-member (cherries only)
                if triple and l == 2 and tc < 3:
                    layers[jn][(nn, pend, dp_, tc + 1)] += w
    if triple:
        return {n: w * d3 for n, w in out.items()}
    return dict(out)

def R_predict(p, singles, pair_delta=None, triple_d3=None):
    W = defaultdict(int)
    for n, w in placement_pass(p, singles).items():
        W[n] += w
    if pair_delta:
        for n, w in placement_pass(p, singles, pair_delta=pair_delta).items():
            W[n] += w
    if triple_d3 is not None:
        for n, w in placement_pass(p, singles, triple=True, d3=triple_d3).items():
            W[n] += w
    R = Fraction(0)
    fp = factorial(p)
    for n, w in W.items():
        R += Fraction(w * factorial(p - n), fp)
    return R

# --- DP unit test vs brute force over all slot-subsets (fake weights) -------
def brute_placement(p, singles, pair_delta, d3):
    tot = defaultdict(int)
    for S in range(1 << (p - 1)):
        runs = []
        l = 0
        for b in range(p - 1):
            if (S >> b) & 1:
                l += 1
            else:
                if l:
                    runs.append(l)
                l = 0
        if l:
            runs.append(l)
        n = sum(r + 1 for r in runs)
        k = len(runs)
        # plain
        w = 1
        for r in runs:
            w *= singles.get(r, 0)
        tot[n] += w
        # matchings of size 1..3
        def matchings(idx):
            if len(idx) < 2:
                return [[]]
            res = [[]]
            first = idx[0]
            for jx in range(1, len(idx)):
                pairm = (first, idx[jx])
                rest = [x for x in idx if x != first and x != idx[jx]]
                for sub in matchings(rest):
                    res.append([pairm] + sub)
            # also matchings not using first
            for sub in matchings(idx[1:]):
                if sub:
                    res.append(sub)
            return res
        seen = set()
        for m in matchings(list(range(k))):
            if not (1 <= len(m) <= 3):
                continue
            key = tuple(sorted(m))
            if key in seen:
                continue
            seen.add(key)
            used = set()
            wm = 1
            ok = True
            for (i1, i2) in m:
                kk = (min(runs[i1], runs[i2]), max(runs[i1], runs[i2]))
                if kk not in pair_delta:
                    ok = False
                    break
                wm *= pair_delta[kk]
                used |= {i1, i2}
            if not ok:
                continue
            for i in range(k):
                if i not in used:
                    wm *= singles.get(runs[i], 0)
            tot[n] += wm
        # cherry triples
        ch = [i for i in range(k) if runs[i] == 2]
        for tri in combinations(ch, 3):
            wt = d3
            for i in range(k):
                if i not in tri:
                    wt *= singles.get(runs[i], 0)
            tot[n] += wt
    return {n: w for n, w in tot.items() if w}

log("")
log("=" * 78)
log("PART 2: placement DP -- unit test vs exhaustive slot-subset enumeration")
log("=" * 78)
for ptest in [7, 11, 13]:
    fake_singles = {2: 5, 4: 3, 6: 2, 8: 7}
    fake_pairs = {(2, 2): -11, (2, 4): 13, (3, 3): -17}
    fake_d3 = 23
    Wdp = defaultdict(int)
    for n, w in placement_pass(ptest, fake_singles).items():
        Wdp[n] += w
    for n, w in placement_pass(ptest, fake_singles, pair_delta=fake_pairs).items():
        Wdp[n] += w
    for n, w in placement_pass(ptest, fake_singles, triple=True, d3=fake_d3).items():
        Wdp[n] += w
    Wbf = brute_placement(ptest, fake_singles, fake_pairs, fake_d3)
    Wdp = {n: w for n, w in Wdp.items() if w}
    assert Wdp == Wbf, (ptest, Wdp, Wbf)
    log(f"  p={ptest}: DP == brute-force over all 2^{ptest-1} slot subsets "
        f"(incl. pair matchings <=3 and cherry triples)  OK")

# ----------------------------------------------------------------------------
# PART 3: predictor ladder + out-of-sample validation
# ----------------------------------------------------------------------------
log("")
log("=" * 78)
log("PART 3: predictor ladder vs exact R(p) (OUT-OF-SAMPLE validation)")
log("=" * 78)
VARIANTS = ['L2', 'L4', 'L6', 'L8', 'L8+pp', 'L8+P', 'L8+P+T']

def predict(p, variant):
    I = INT[p]
    D = DELTA[p]
    singles_all = {2: I['A2'], 4: I['A4'], 6: I['A6'], 8: I['A8']}
    if variant == 'L2':
        return R_predict(p, {2: I['A2']})
    if variant == 'L4':
        return R_predict(p, {2: I['A2'], 4: I['A4']})
    if variant == 'L6':
        return R_predict(p, {2: I['A2'], 4: I['A4'], 6: I['A6']})
    if variant == 'L8':
        return R_predict(p, singles_all)
    if variant == 'L8+pp':
        return R_predict(p, singles_all, pair_delta={(2, 2): D['d22']})
    pairs = {(2, 2): D['d22'], (2, 4): D['d24'], (3, 3): D['d33']}
    if variant == 'L8+P':
        return R_predict(p, singles_all, pair_delta=pairs)
    if variant == 'L8+P+T':
        return R_predict(p, singles_all, pair_delta=pairs, triple_d3=D['d3'])
    raise ValueError(variant)

R_EXACT = {p: Fraction(H_KNOWN[p] * 2 ** (p - 1), factorial(p)) for p in VAL_PRIMES}
PRED = {}
for p in ALL_PRIMES:
    PRED[p] = {v: predict(p, v) for v in VARIANTS}

log(f"{'p':>3} | {'R exact':>10} | " + " | ".join(f"{v:>10}" for v in VARIANTS))
for p in ALL_PRIMES:
    re_ = f"{float(R_EXACT[p]):>10.6f}" if p in R_EXACT else f"{'?':>10}"
    log(f"{p:>3} | {re_} | " + " | ".join(f"{float(PRED[p][v]):>10.6f}" for v in VARIANTS))
log("")
log("residuals eps = R_exact - R_pred:")
log(f"{'p':>3} | " + " | ".join(f"{v:>11}" for v in VARIANTS))
EPS = {}
for p in VAL_PRIMES:
    EPS[p] = {v: float(R_EXACT[p] - PRED[p][v]) for v in VARIANTS}
    log(f"{p:>3} | " + " | ".join(f"{EPS[p][v]:>11.3e}" for v in VARIANTS))
log("")
log("residual scalings for the FINAL variant (L8+P+T):  eps*p, eps*p^2, eps*p^3")
log(f"{'p':>3} | {'eps':>11} | {'eps*p':>9} | {'eps*p^2':>9} | {'eps*p^3':>10}")
for p in VAL_PRIMES:
    e_ = EPS[p]['L8+P+T']
    log(f"{p:>3} | {e_:>11.3e} | {e_*p:>9.4f} | {e_*p**2:>9.4f} | {e_*p**3:>10.3f}")
log("")
log("local empirical order gamma-hat = log(eps(p1)/eps(p2)) / log(p2/p1):")
ps = [p for p in VAL_PRIMES if p >= 7]
for (p1, p2) in zip(ps, ps[1:]):
    e1, e2 = EPS[p1]['L8+P+T'], EPS[p2]['L8+P+T']
    if e1 * e2 > 0:
        g = np.log(abs(e1) / abs(e2)) / np.log(p2 / p1)
        log(f"  ({p1:>2},{p2:>2}): gamma-hat = {g:.2f}")
    else:
        log(f"  ({p1:>2},{p2:>2}): SIGN CHANGE (eps {e1:.2e} -> {e2:.2e})")

# ----------------------------------------------------------------------------
# PART 4: HONEST VERDICT on the truncation ladder + the proven-form predictor
# ----------------------------------------------------------------------------
log("")
log("=" * 78)
log("PART 4: verdict on truncation; proven-form asymptotic fit; R(31) prediction")
log("=" * 78)
log("""
VERDICT ON THE FINITE-p TRUNCATED CLUSTER SUMS (the out-of-sample answer):
the depth-L / collision-order ladder does NOT converge at p <= 31. The exact
collision corrections have rapidly growing coefficients
   d22/p^3 -> ~-10,   d24/p^4 -> ~-28,   d3(222)/p^4 -> ~+262,
i.e. the multi-run inclusion-exclusion series in 1/p is ASYMPTOTIC with
factorially growing coefficients -- the same A088368 ~ e*k! growth that
THM-438 ADDENDUM-6 proved makes the cluster generating function Gevrey-1 /
resurgent, NOT algebraic. Successive truncations oscillate (L8: eps ~ -4;
L8+pp: eps ~ +-0.2-0.4; L8+P, L8+P+T: divergent oscillation). The expansion
PROVES the limit e but cannot deliver 1e-3 accuracy at p=31 by termwise
truncation. This is recorded as the honest negative half of HYP-2371.
""")
# cluster-only coarse bracket from the best-behaved truncation (L8+pp)
BEST_TRUNC = 'L8+pp'
res_bt = [abs(EPS[p][BEST_TRUNC]) for p in [11, 19, 23]]
bar_bt = max(res_bt)
R31_bt = float(PRED[31][BEST_TRUNC])
log(f"cluster-only coarse bracket (best truncation {BEST_TRUNC}, bar = max |residual|")
log(f"  over p=11,19,23 = {bar_bt:.3f}):   R(31) in [{R31_bt-bar_bt:.3f}, {R31_bt+bar_bt:.3f}]")
log("")
log("PRIMARY PREDICTOR: the PROVEN functional form. THM-438 ADDENDUM (4) proves")
log("A_{2k} = C_k p^{k+1} + O(p^k), hence R(p) - e is RELATIVE O(1/p); the")
log("asymptotic ansatz R(p) = e(1 - C/p - D/p^2 - E/p^3 - ...) is justified.")
log("We fit the first m coefficients on the m largest exact points (exact data,")
log("float solve), with an ORDER-LADDERED HOLDOUT at p=23 to calibrate the bar.")
log("")

def fit_predict(points, p_star):
    """Fit 1 - R/e = sum_j c_j / p^j (j=1..m) through the m exact points; predict R(p_star)."""
    m = len(points)
    M = np.array([[1.0 / p ** j for j in range(1, m + 1)] for p in points])
    rhs = np.array([1 - float(R_EXACT[p]) / MATH_E for p in points])
    c = np.linalg.solve(M, rhs)
    pred = MATH_E * (1 - sum(c[j] / p_star ** (j + 1) for j in range(m)))
    return pred, c

# holdout at p=23 (each model order, anchored on the largest points below 23)
log("HOLDOUT at p=23 (fit on points < 23, predict R(23) = 2.556980):")
hold_err = {}
for m, pts in [(1, [19]), (2, [11, 19]), (3, [7, 11, 19])]:
    pr, c = fit_predict(pts, 23)
    hold_err[m] = float(R_EXACT[23]) - pr
    log(f"  m={m} fit on {pts}: R(23)_pred = {pr:.6f}  err = {hold_err[m]:+.5f}"
        f"  coeffs = {np.array2string(c, precision=3)}")
log("  (coefficient growth |C|~1.6, |D|~5, |E|~30+ confirms the asymptotic/")
log("   resurgent character; adding orders does not monotonically help)")
log("")
# predictions at p=31 (anchored on the largest available points)
log("PREDICTIONS at p=31:")
preds31 = {}
for tag, pts in [('m=1 (23)', [23]), ('m=2 (19,23)', [19, 23]),
                 ('m=2 (11,23)', [11, 23]), ('m=3 (11,19,23)', [11, 19, 23]),
                 ('m=4 (7,11,19,23)', [7, 11, 19, 23])]:
    pr, c = fit_predict(pts, 31)
    preds31[tag] = pr
    log(f"  {tag:>18}: R(31) = {pr:.6f}   coeffs = {np.array2string(c, precision=3)}")
log("")
# central value and bar:
#   central = mean of the two highest-order top-anchored fits (m=2 (19,23), m=3)
#   bar     = (model spread)/2 + max over m of |holdout err| * (23/31)^(m+1),
#             then x1.5 safety for the asymptotic (non-convergent) tail.
core = [preds31['m=2 (19,23)'], preds31['m=3 (11,19,23)']]
spread_models = [preds31[t] for t in preds31]
spread = max(spread_models) - min(spread_models)
proj = max(abs(hold_err[m]) * (23 / 31) ** (m + 1) for m in hold_err)
R31_central = float(np.mean(core))
half31 = (spread / 2 + proj) * 1.5
R31_lo, R31_hi = R31_central - half31, R31_central + half31
log(f"model spread (all 5 fits)         = {spread:.5f}")
log(f"max holdout error projected to 31 = {proj:.5f}")
log(f"==> R(31) PREDICTED = {R31_central:.5f} +/- {half31:.5f}"
    f"   i.e.  R(31) in [{R31_lo:.5f}, {R31_hi:.5f}]")
log(f"    (bar = (spread/2 + projected holdout) x 1.5 safety; e - R(31) ~ "
    f"{MATH_E-R31_central:.4f}; (e-R)*31 = {(MATH_E-R31_central)*31:.2f}, cf. 3.71 at p=23,")
log(f"     monotone-increase consistency: R(23)=2.55698 < R(31) < e REQUIRED: "
    f"{'OK' if 2.556980 < R31_lo and R31_hi < MATH_E else 'VIOLATED'})")
in_bt = (R31_bt - bar_bt) <= R31_central <= (R31_bt + bar_bt)
log(f"    cluster-only coarse bracket contains the fit prediction: "
    f"{'YES' if in_bt else 'NO -- TENSION, widen!'}")

# ----------------------------------------------------------------------------
# PART 5: integrality sharpening and the falsifiable H(T_31) interval
# ----------------------------------------------------------------------------
log("")
log("=" * 78)
log("PART 5: integrality -- H(T_31) = 465 * ODD, and the predicted interval")
log("=" * 78)
log("orbit-parity re-verification (HYP-1713 PROVED: H odd by Redei; |Aff(QR)| =")
log("p(p-1)/2 odd since v2(p-1)=1 at p=3 mod 4; hence H/|Aff| is an ODD integer):")
for p in VAL_PRIMES:
    aff = p * (p - 1) // 2
    H = H_KNOWN[p]
    q, r = divmod(H, aff)
    log(f"  p={p:>2}: |Aff|={aff:>4}  H/|Aff| = {q}  remainder={r}  odd={q % 2 == 1}")
    assert r == 0 and q % 2 == 1 and H % 2 == 1
log("  all 5 known values: H/|Aff| is an odd integer -- CONFIRMED")
log("")
aff31 = 31 * 30 // 2
fac = Fraction(factorial(31), 2 ** 30)   # H = R * 31!/2^30
log(f"31!/2^30 = {float(fac):.6e} (exact: {fac.numerator}/{fac.denominator})")
# interval endpoints exactly via Fractions (floats only entered via the bar)
from math import floor as _floor, ceil as _ceil
R31_lo_fr = Fraction(int(round(R31_lo * 10 ** 12)), 10 ** 12)
R31_hi_fr = Fraction(int(round(R31_hi * 10 ** 12)), 10 ** 12)
H_lo = _floor(R31_lo_fr * fac)
H_hi = _ceil(R31_hi_fr * fac)
H_central = R31_central * float(fac)
width = H_hi - H_lo
# admissible: H == 465 mod 930  (H = 465*odd)
first = H_lo + ((465 - H_lo) % 930)
count_adm = 0 if first > H_hi else (H_hi - first) // 930 + 1
log(f"|Aff(QR_31)| = {aff31}  (= 465 = 3*5*31; divides H by freeness/THM-F3)")
log(f"H(T_31) = 465 * ODD  <=>  H == 465 (mod 930)")
log("")
log(f"PREDICTED:  H(T_31) ~= {H_central:.4e}")
log(f"  interval  [{H_lo}, {H_hi}]")
log(f"  width     {width}  ({float(width):.3e})")
log(f"  admissible integers (465*odd) inside: {count_adm}  ({float(count_adm):.3e})")
log(f"  predicted H/465 ('1729-analog') ~= {H_central/465:.4e}  (must be ODD)")
rel = half31 / R31_central
log(f"  relative half-width: {rel:.2e}  => ~{int(np.floor(-np.log10(2*rel)))} "
    f"significant digits honestly supported")
log("")
log(f"[total elapsed {time.time()-T0:.1f}s]")
