#!/usr/bin/env python3
"""
paley_ladder_fit_kpc1_verify.py
ADVERSARIAL VERIFIER, kind-pasteur-2026-06-10-S1, thread E-R31-prediction.
Verifies claims E4, E5, E6, E7.

INDEPENDENT METHOD (vs the worker's slot-pointer DP):
  R_variant(p) is assembled COMPOSITION-LEVEL: enumerate ordered sequences of
  run lengths (l_1..l_k), parts in the variant's allowed set, with
  sum(l)+k-1 <= p-1; the number of slot-placements of that ordered sequence is
  the closed form C((p-1) - sum(l) + 1, k) (stars-and-bars with mandatory
  gaps).  Weights: product of singles + sum over matchings (<=3 marked pairs)
  of delta-products + cherry-triple term.  All exact integers / Fractions.

  TWO pair-matching semantics are computed:
    'spec' = ALL matchings of <=3 disjoint marked pairs (the worker's stated
             definition and their brute-force reference);
    'dp'   = only matchings realizable left-to-right with at most 2 pending
             open runs (what the worker's placement_pass actually implements:
             its pair-opener branch requires len(pend) < 2).
  These coincide for p <= 13 (six matched runs do not fit) -- exactly the
  range of the worker's unit test -- but DIFFER for p >= 19.

  VALIDATION:
    (i)  identity sum over compositions of placements == 2^(p-1) - 1
         (every nonempty slot subset has a unique (sequence, placement));
    (ii) exhaustive enumeration over ALL 2^(p-1) slot subsets at p = 7..19,
         computing both semantics directly per subset.

Then: residual tables vs exact R(p) (E4), the L8+pp bracket (E5), the
asymptotic-form fit with exact-Fraction linear solves (E6), and the
integrality sharpening with exact endpoints (E7).
"""
import json, os, sys, time
from fractions import Fraction
from math import factorial, comb, e as MATH_E, log
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
T0 = time.time()

H_KNOWN = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349}
VAL_PRIMES = [3, 7, 11, 19, 23]
ALL_PRIMES = VAL_PRIMES + [31]
R_EXACT = {p: Fraction(H_KNOWN[p] * 2 ** (p - 1), factorial(p)) for p in VAL_PRIMES}

with open(os.path.join(ROOT, "05-knowledge", "results",
                       "paley_cluster_sums_kpc1_verify.json")) as f:
    INT = {int(k): v for k, v in json.load(f).items()}

# ---------------------------------------------------------------------------
# weight of one ordered run-length sequence under a variant
# ---------------------------------------------------------------------------
def matchings_weight(seq, singles, pairs, mode):
    """Sum over matchings of 1..3 disjoint marked pairs (keys in `pairs`) of
       prod(delta) * prod(singles) over unmatched runs.
       mode='dp': only matchings realizable with <=2 pending opens."""
    k = len(seq)
    idx = [i for i in range(k) if any(seq[i] in key for key in pairs)]
    total = 0
    def feasible_dp(m):
        opens = {}
        for (i, j) in m:
            opens[i] = 'o'
            opens[j] = 'c'
        pend = 0
        for t in sorted(opens):
            if opens[t] == 'o':
                if pend >= 2:
                    return False
                pend += 1
            else:
                pend -= 1
        return True
    def emit(m):
        nonlocal total
        if 1 <= len(m) <= 3 and (mode == 'spec' or feasible_dp(m)):
            w = 1
            used = set(x for pr in m for x in pr)
            for (i, j) in m:
                w *= pairs[(min(seq[i], seq[j]), max(seq[i], seq[j]))]
            for i in range(k):
                if i not in used:
                    w *= singles.get(seq[i], 0)
            total += w
    def rec(avail, m):
        # each matching reached by exactly ONE decision path: at the first
        # available position, either leave it unmatched forever or pair it.
        if not avail:
            emit(m)
            return
        first = avail[0]
        rec(avail[1:], m)                      # first stays unmatched
        if len(m) < 3:
            for jx in range(1, len(avail)):
                j = avail[jx]
                key = (min(seq[first], seq[j]), max(seq[first], seq[j]))
                if key in pairs:
                    rec([x for x in avail[1:] if x != j], m + [(first, j)])
    rec(idx, [])
    return total

def seq_weight(seq, singles, pairs, d3, mode):
    w = 1
    for l in seq:
        w *= singles.get(l, 0)
    if pairs:
        w += matchings_weight(seq, singles, pairs, mode)
    if d3 is not None:
        ch = sum(1 for l in seq if l == 2)
        if ch >= 3:
            rest = 1
            for l in seq:
                rest *= singles.get(l, 0) if l != 2 else 1
            w += comb(ch, 3) * d3 * rest * (singles.get(2, 0) ** (ch - 3))
    return w

# ---------------------------------------------------------------------------
# composition-level assembly (closed-form placements)
# ---------------------------------------------------------------------------
def R_variant(p, singles, pairs, d3, mode):
    parts = sorted(set(list(singles) + [l for key in (pairs or {}) for l in key]
                       + ([2] if d3 is not None else [])))
    N = p - 1
    total = Fraction(1)          # empty configuration
    fp = factorial(p)
    cache = {}
    def rec(seq, suml):
        nonlocal total
        if seq:
            k = len(seq)
            n = suml + k
            placements = comb(N - suml + 1, k)
            if placements:
                key = (tuple(seq), mode)
                if key not in cache:
                    cache[key] = seq_weight(seq, singles, pairs or {}, d3, mode)
                w = cache[key]
                if w:
                    total += Fraction(placements * w * factorial(p - n), fp)
        for l in parts:
            k = len(seq)
            span = suml + (k - 1 if k else 0)
            if span + l + (1 if k else 0) <= N and suml + l + k + 1 <= p:
                seq.append(l)
                rec(seq, suml + l)
                seq.pop()
    rec([], 0)
    return total

# ---------------------------------------------------------------------------
# exhaustive subset reference (validation)
# ---------------------------------------------------------------------------
def R_exhaustive(p, singles, pairs, d3, mode):
    N = p - 1
    fp = factorial(p)
    total = Fraction(0)
    Wn = {}
    cache = {}
    for S in range(1 << N):
        seq, l = [], 0
        for b in range(N):
            if (S >> b) & 1:
                l += 1
            else:
                if l:
                    seq.append(l)
                l = 0
        if l:
            seq.append(l)
        key = tuple(seq)
        if key not in cache:
            cache[key] = seq_weight(seq, singles, pairs or {}, d3, mode) if seq else 1
        n = sum(seq) + len(seq)
        Wn[n] = Wn.get(n, 0) + cache[key]
    for n, w in Wn.items():
        total += Fraction(w * factorial(p - n), fp)
    return total

VARIANTS = ['L2', 'L4', 'L6', 'L8', 'L8+pp', 'L8+P', 'L8+P+T']

def variant_args(p, v):
    I = INT[p]
    singles_all = {2: I['A2'], 4: I['A4'], 6: I['A6'], 8: I['A8']}
    if v == 'L2':
        return {2: I['A2']}, None, None
    if v == 'L4':
        return {2: I['A2'], 4: I['A4']}, None, None
    if v == 'L6':
        return {2: I['A2'], 4: I['A4'], 6: I['A6']}, None, None
    if v == 'L8':
        return singles_all, None, None
    if v == 'L8+pp':
        return singles_all, {(2, 2): I['d22']}, None
    pairs = {(2, 2): I['d22'], (2, 4): I['d24'], (3, 3): I['d33']}
    if v == 'L8+P':
        return singles_all, pairs, None
    if v == 'L8+P+T':
        return singles_all, pairs, I['d3']
    raise ValueError(v)

print("=" * 78)
print("VALIDATION (i): sum of closed-form placements over compositions = 2^N - 1")
print("=" * 78)
for N in [4, 7, 10, 13, 16]:
    tot = 0
    def reca(seq_k, suml):
        global tot
        k = seq_k
        if k:
            tot += comb(N - suml + 1, k)
        for l in range(1, N + 1):
            span = suml + (k - 1 if k else 0)
            if span + l + (1 if k else 0) <= N:
                reca(k + 1, suml + l)
    # iterate over part counts only (placements depend on (k, sum l) only via comb)
    tot = 0
    def recb(k, suml):
        global tot
        if k:
            c = comb(N - suml + 1, k)
            if c == 0:
                return
            tot += c
        for l in range(1, N - suml - (k if k else 0) + 1 + 1):
            if suml + l + (k - 1 if k else 0) + (1 if k else 0) <= N:
                recb(k + 1, suml + l)
    recb(0, 0)
    print(f"  N={N:>2}: sum placements = {tot}  ?= 2^N-1 = {2**N - 1}  "
          f"{'OK' if tot == 2**N - 1 else 'FAIL'}")
    assert tot == 2 ** N - 1

print()
print("=" * 78)
print("VALIDATION (ii): composition assembly vs exhaustive 2^(p-1) subsets")
print("=" * 78)
FAKE = ({2: 5, 4: -3, 6: 2, 8: 7},
        {(2, 2): -7, (2, 4): 11, (3, 3): -13}, 17)
for p in [7, 11, 13, 17, 19]:
    if p in INT:
        cases = [variant_args(p, v) for v in ['L8', 'L8+pp', 'L8+P+T']]
        kind = "real weights"
    else:
        cases = [FAKE]          # fake weights incl. nonzero (3,3) pair delta
        kind = "FAKE weights (exercises (3,3) pairs, d33 nonzero)"
    for (s_, pr_, d3_) in cases:
        for mode in ['spec', 'dp']:
            a = R_variant(p, s_, pr_, d3_, mode)
            b = R_exhaustive(p, s_, pr_, d3_, mode)
            assert a == b, (p, s_, mode, a, b)
    print(f"  p={p:>2}: composition == exhaustive ({kind})  OK"
          f"   ({time.time()-T0:.0f}s)")
print("  NOTE: at p=19 six cherry runs fit (span 17 <= 18), so the two")
print("  semantics genuinely differ there and BOTH were validated exhaustively.")

print()
print("=" * 78)
print("E4/E5: predictor ladder, exact Fractions, both pair-matching semantics")
print("=" * 78)
PRED = {}
for mode in ['spec', 'dp']:
    PRED[mode] = {}
    for p in ALL_PRIMES:
        PRED[mode][p] = {}
        for v in VARIANTS:
            s_, pr_, d3_ = variant_args(p, v)
            PRED[mode][p][v] = R_variant(p, s_, pr_, d3_, mode)
for mode in ['spec', 'dp']:
    print(f"\n--- semantics: {mode} "
          f"({'all matchings <=3 (stated spec)' if mode=='spec' else 'pending<=2 replica of worker DP'}) ---")
    print(f"{'p':>3} | {'R exact':>10} | " + " | ".join(f"{v:>10}" for v in VARIANTS))
    for p in ALL_PRIMES:
        re_ = f"{float(R_EXACT[p]):>10.6f}" if p in R_EXACT else f"{'?':>10}"
        print(f"{p:>3} | {re_} | " +
              " | ".join(f"{float(PRED[mode][p][v]):>10.6f}" for v in VARIANTS))
    print("residuals eps = R_exact - R_pred:")
    print(f"{'p':>3} | " + " | ".join(f"{v:>11}" for v in VARIANTS))
    for p in VAL_PRIMES:
        print(f"{p:>3} | " + " | ".join(
            f"{float(R_EXACT[p] - PRED[mode][p][v]):>11.3e}" for v in VARIANTS))

print()
print("difference (dp - spec) at the pair variants (implementation-semantics gap):")
for p in [19, 23, 31]:
    print(f"  p={p:>2}: " + "  ".join(
        f"{v}: {float(PRED['dp'][p][v] - PRED['spec'][p][v]):+.5f}"
        for v in ['L8+pp', 'L8+P', 'L8+P+T']))

print()
print("E4 checks (worker numbers come from THEIR dp semantics):")
for mode in ['spec', 'dp']:
    e19 = float(R_EXACT[19] - PRED[mode][19]['L8+P+T'])
    e23 = float(R_EXACT[23] - PRED[mode][23]['L8+P+T'])
    if e19 * e23 > 0:
        g = log(abs(e19) / abs(e23)) / log(23 / 19)
        print(f"  [{mode}] L8+P+T eps(19)={e19:+.4f} eps(23)={e23:+.4f} "
              f"gamma-hat(19,23) = {g:+.2f}   (claimed ~ -0.9)")
    else:
        print(f"  [{mode}] L8+P+T eps(19)={e19:+.4f} eps(23)={e23:+.4f}  SIGN CHANGE")
    allres = [abs(float(R_EXACT[p] - PRED[mode][p][v]))
              for p in [7, 11, 19, 23] for v in VARIANTS]
    print(f"  [{mode}] |residual| range over p=7..23, all variants: "
          f"{min(allres):.3f} .. {max(allres):.3f}   (claimed oscillation O(0.2-4.7))")

print()
print("E5 checks: L8+pp bracket")
for mode in ['spec', 'dp']:
    res = [abs(float(R_EXACT[p] - PRED[mode][p]['L8+pp'])) for p in [11, 19, 23]]
    bar = max(res)
    R31bt = float(PRED[mode][31]['L8+pp'])
    print(f"  [{mode}] residuals(11,19,23) = {[f'{r:.3f}' for r in res]}, "
          f"bar = {bar:.3f}, R31_bt = {R31bt:.3f}, "
          f"bracket = [{R31bt-bar:.3f}, {R31bt+bar:.3f}]")
print("  claimed: bar 0.369, R31_bt 2.619, bracket [2.250, 2.988]")

# ---------------------------------------------------------------------------
# E6: asymptotic-form fit, exact-Fraction linear algebra
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("E6: fit 1 - R/e = sum_{j=1..m} c_j p^-j  (exact-Fraction solve)")
print("=" * 78)
EF = Fraction(MATH_E)          # exact binary64 value of e (worker used float e)

def solve_frac(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(n + 1)]
    return [M[i][n] / M[i][i] for i in range(n)]

def fit_predict(points, p_star):
    m = len(points)
    A = [[Fraction(1, q ** j) for j in range(1, m + 1)] for q in points]
    b = [1 - R_EXACT[q] / EF for q in points]
    c = solve_frac(A, b)
    pred = EF * (1 - sum(c[j] / p_star ** (j + 1) for j in range(m)))
    return pred, c

print("holdout at p=23 (fit on smaller points, predict 23):")
hold = {}
for m, pts in [(1, [19]), (2, [11, 19]), (3, [7, 11, 19])]:
    pr, c = fit_predict(pts, 23)
    hold[m] = R_EXACT[23] - pr
    print(f"  m={m} on {pts}: R23_pred = {float(pr):.6f}  err = {float(hold[m]):+.5f}"
          f"  coeffs = {[f'{float(x):.3f}' for x in c]}")
print("  claimed: m=1 -0.0034, m=2 +0.0025, m=3 +0.0022, max |err| <= 3.4e-3")

print("\npredictions at p=31 (the worker's five anchor sets):")
preds = {}
for tag, pts in [('m=1 (23)', [23]), ('m=2 (19,23)', [19, 23]),
                 ('m=2 (11,23)', [11, 23]), ('m=3 (11,19,23)', [11, 19, 23]),
                 ('m=4 (7,11,19,23)', [7, 11, 19, 23])]:
    pr, c = fit_predict(pts, 31)
    preds[tag] = pr
    print(f"  {tag:>18}: R31 = {float(pr):.6f}  coeffs = {[f'{float(x):.3f}' for x in c]}")
print("  claimed: predictions across all 5 fits in 2.5937 - 2.5986;")
print("  claimed coefficients (m=4): C~1.4-1.6, then -5, -48, +427")

core = [preds['m=2 (19,23)'], preds['m=3 (11,19,23)']]
spread = max(preds.values()) - min(preds.values())
proj = max(abs(hold[m]) * Fraction(23, 31) ** (m + 1) for m in hold)
central = (core[0] + core[1]) / 2
half = (spread / 2 + proj) * Fraction(3, 2)
lo, hi = central - half, central + half
print(f"\nmodel spread (5 fits)      = {float(spread):.5f}   (claimed 0.0049)")
print(f"projected holdout to 31    = {float(proj):.5f}")
print(f"central (mean of 2 core)   = {float(central):.5f}   (claimed 2.59599)")
print(f"half-bar (spread/2+proj)*1.5 = {float(half):.5f}  (claimed 0.00650)")
print(f"interval                   = [{float(lo):.5f}, {float(hi):.5f}]"
      f"   (claimed [2.58949, 2.60249])")
print(f"consistency: 2.55698 < lo and hi < e: "
      f"{'OK' if float(lo) > 2.556980 and float(hi) < MATH_E else 'VIOLATED'}")
for q in [19, 23]:
    print(f"  (e - R({q})) * {q} = {(MATH_E - float(R_EXACT[q])) * q:.4f}")
print(f"  (e - R31_central) * 31 = {(MATH_E - float(central)) * 31:.2f}"
      f"   (claimed 3.79, continuing 3.63 / 3.71)")

# ---------------------------------------------------------------------------
# E7: integrality sharpening
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("E7: orbit parity + exact H(T_31) interval")
print("=" * 78)
print("parity on the 5 known values (Redei: H odd; |Aff| = p(p-1)/2 odd at p=3 mod 4):")
for p in VAL_PRIMES:
    aff = p * (p - 1) // 2
    q, r = divmod(H_KNOWN[p], aff)
    print(f"  p={p:>2}: |Aff|={aff:>4}  H/|Aff|={q}  rem={r}  odd={q % 2 == 1}")
    assert r == 0 and q % 2 == 1
print("  claimed quotients 1, 9, 1729, 6857869865, 62293308207033 all odd: "
      + str([H_KNOWN[p] // (p * (p - 1) // 2) for p in VAL_PRIMES]))

fac = Fraction(factorial(31), 2 ** 30)
print(f"\n31!/2^30 = {float(fac):.6e}")
# worker recipe: round interval ends to 12 decimal digits, then floor/ceil * fac
lo_f, hi_f = float(lo), float(hi)
lo_fr = Fraction(round(lo_f * 10 ** 12), 10 ** 12)
hi_fr = Fraction(round(hi_f * 10 ** 12), 10 ** 12)
import math as _m
H_lo = _m.floor(lo_fr * fac)
H_hi = _m.ceil(hi_fr * fac)
width = H_hi - H_lo
first = H_lo + ((465 - H_lo) % 930)
count = 0 if first > H_hi else (H_hi - first) // 930 + 1
print(f"MY interval : [{H_lo}, {H_hi}]")
print(f"   width    : {width}  ({float(width):.3e})")
print(f"   admissible (465*odd): {count}  ({float(count):.4e})")
print(f"   H central ~ {float(central * fac):.4e};  H/465 ~ {float(central*fac/465):.4e}")
rel = half / central
print(f"   relative half-width {float(rel):.2e}")

CL_LO = 19830629617139608462365775
CL_HI = 19930130881568868002912737
CL_CNT = 106990606913182301663
print(f"\nWORKER claimed interval [{CL_LO}, {CL_HI}], count {CL_CNT}")
print(f"  my endpoints equal claimed: {H_lo == CL_LO and H_hi == CL_HI}")
if H_lo != CL_LO or H_hi != CL_HI:
    print(f"  diff lo = {H_lo - CL_LO}, diff hi = {H_hi - CL_HI} "
          f"(float-pipeline noise tolerance: |diff| < ~1e16 is immaterial at 1e25)")
fc = CL_LO + ((465 - CL_LO) % 930)
cc = (CL_HI - fc) // 930 + 1
print(f"  internal consistency of claimed numbers: count from claimed endpoints "
      f"= {cc}  ?= {CL_CNT}  {'OK' if cc == CL_CNT else 'FAIL'}")
print(f"  claimed-implied R bounds: lo/fac = {float(Fraction(CL_LO)/fac):.10f}, "
      f"hi/fac = {float(Fraction(CL_HI)/fac):.10f}")
print(f"\n[elapsed {time.time()-T0:.1f}s]")
