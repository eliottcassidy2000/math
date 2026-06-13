#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPLORATORY (hypothesis-generating, NOT proof).

The code-side mirror of THM-485's pentagonal sign-rigidity.

BACKGROUND (THM-486):
  The deterministic EXTREMAL Type II weight enumerator at length n = 24m is
      W_n = sum_{j=0}^{m} c_j * W8^{3(m-j)} * P24^{j},   c_0 = 1,
  with W8 = W_e8hat = x^8 + 14 x^4 y^4 + y^8  and
       P24 = x^4 y^4 (x^4 - y^4)^4  (= 16 * eta(q^2)^24 = 16 Delta, the
       24th power of the pentagonal product).
  The c_j are FIXED by extremality (A_4 = ... = A_{4m} = 0). The leading
  correction is c_1(m) = -42m (PROVED, THM-486 part B); the deeper c_j carry
  an alternating "secular" stream whose signs maximally DELAY the first negative
  weight-enumerator coefficient A_w to the famous Zhang threshold n = 3696
  (residue-0 branch).

THIS SCRIPT — the random-sign analog (parallel to THM-485):
  Keep the MAGNITUDES |c_j| from the extremal solution but RANDOMIZE the signs
  (iid +-1, seeded). c_0 is kept +1 (it is the normalization that sets A_0 = 1;
  flipping it trivially makes A_0 = -1, which is not interesting). Then reassemble
      W_n^{(rand)} = sum_j  (s_j * |c_j|) * W8^{3(m-j)} * P24^j,  s_0 = +1, s_j in {+-1}.
  and find the FIRST length 24m at which some A_w < 0.

  We measure the DISTRIBUTION of the first-negativity length over many seeds, for
  m up to MAX_M (~60), and compare to the deterministic 3696. If random signs
  make negativity appear MUCH earlier than 3696, that is exploratory EVIDENCE
  that the deterministic extremal alternation is "maximally delaying" negativity
  -- the code-side avatar of THM-485's Euler-signs-are-maximally-subexponential
  rigidity.

ARITHMETIC: signs are random, but the magnitudes |c_j| and the assembled A_w are
  computed in EXACT big-int / exact-rational arithmetic (Fraction). No floats in
  the enumerator coefficients.

CAVEAT on the sign-randomization model:
  Within a single fixed m, the random-sign enumerator has |c_j| fixed and only the
  signs s_j (j>=1) flipped. The "first-negativity length" for a given seed is the
  smallest m (over the sweep) such that SOME A_w(m) < 0 under the seed's sign
  pattern at that m. Two natural conventions for "the same seed across m":
    MODE 'fresh'   : at each m draw a fresh sign pattern for c_1..c_m (length m).
    MODE 'prefix'  : draw one long sign stream once per seed; at length m use its
                     first m signs (so growing m only APPENDS signs -- a coherent
                     single random correction stream, closest to "randomize the
                     Burmann correction stream").
  We report BOTH. 'prefix' is the honest analog of one random correction stream.

Output: tee to 05-knowledge/results/random_extremal_enumerator_kps3_0611.out
"""

import sys
import time
from fractions import Fraction
from math import comb

# ----------------------------------------------------------------------------
# Polynomial machinery (homogeneous in (x,y), tracked as {y-degree: coeff})
# x-degree is implicit (= total degree - y-degree). Type II A_w lives at key w.
# ----------------------------------------------------------------------------

def poly_mul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            out[a + b] = out.get(a + b, 0) + ca * cb
    return out


def poly_pow(p, e):
    r = {0: 1}
    base = dict(p)
    while e:
        if e & 1:
            r = poly_mul(r, base)
        base = poly_mul(base, base)
        e >>= 1
    return r


W8 = {0: 1, 4: 14, 8: 1}  # W_e8hat, degree 8


def build_P24():
    base = {}
    for k in range(5):
        base[4 * k] = comb(4, k) * (-1) ** k
    return {d + 4: c for d, c in base.items()}  # * x^4 y^4 -> shift y-deg by 4


P24 = build_P24()


# ----------------------------------------------------------------------------
# Deterministic extremal solve: returns the EXACT c_j (Fraction) and the basis.
# Cache basis polynomials by m to avoid recompute (they are reused per seed).
# ----------------------------------------------------------------------------

_basis_cache = {}


def basis_for(m):
    if m in _basis_cache:
        return _basis_cache[m]
    basis = []
    for j in range(m + 1):
        bj = poly_mul(poly_pow(W8, 3 * (m - j)), poly_pow(P24, j))
        basis.append(bj)
    _basis_cache[m] = basis
    return basis


_cj_cache = {}


def extremal_cj(m):
    """Exact extremal coefficients c_0..c_m (Fraction). c_0 = 1."""
    if m in _cj_cache:
        return _cj_cache[m]
    basis = basis_for(m)
    c = [Fraction(0)] * (m + 1)
    c[0] = Fraction(1)
    for i in range(1, m + 1):
        s = Fraction(0)
        for j in range(i):
            s += c[j] * basis[j].get(4 * i, 0)
        bii = basis[i].get(4 * i, 0)
        c[i] = -s / bii
    _cj_cache[m] = c
    return c


def assemble(m, signs):
    """Assemble W_n with coefficient vector signs[j]*|c_j| (Fraction in, but
    extremal c_j are integers so result is integer). Returns {y-deg: int A_w}.
    signs has length m+1, signs[0] must be +1."""
    c = extremal_cj(m)
    basis = basis_for(m)
    W = {}
    for j in range(m + 1):
        cj = abs(c[j]) * signs[j]
        bj = basis[j]
        for d, coeff in bj.items():
            W[d] = W.get(d, 0) + cj * coeff
    Wi = {}
    for d, coeff in W.items():
        assert coeff.denominator == 1, f"non-integer A_{d}: {coeff}"
        if coeff != 0:
            Wi[d] = int(coeff)
    return Wi


def first_negative_weight(W):
    """Smallest weight w>0 with A_w<0, or None. (A_0 = +1 always here.)"""
    negs = [w for w, v in W.items() if w > 0 and v < 0]
    return min(negs) if negs else None


def has_negative(W):
    return any(v < 0 for v in W.values())


# ----------------------------------------------------------------------------
# RNG: a small seeded LCG so results are reproducible without numpy-version drift.
# ----------------------------------------------------------------------------

class SignStream:
    """Deterministic +-1 stream from a seed (splitmix64-style)."""
    def __init__(self, seed):
        self.state = (seed * 0x9E3779B97F4A7C15 + 0x1234567) & ((1 << 64) - 1)

    def _next64(self):
        self.state = (self.state + 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & ((1 << 64) - 1)
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & ((1 << 64) - 1)
        z = z ^ (z >> 31)
        return z

    def sign(self):
        return 1 if (self._next64() & 1) else -1


# ----------------------------------------------------------------------------
# Verification: reproduce the deterministic landmarks.
# ----------------------------------------------------------------------------

def det_signs(m):
    c = extremal_cj(m)
    return [1 if cj >= 0 else -1 for cj in c]


def verify_deterministic(window_m=60, full_neg_search=False, max_m_neg=170):
    """Reproduce THM-486 landmarks and confirm the deterministic extremal
    enumerator stays NONNEGATIVE through the random experiment's window (m<=window_m).
    The full first-negativity n=3696 (m=154) is PROVED/verified in THM-486 part C;
    re-deriving it costs a degree-3696 polynomial sweep, so we only optionally redo
    it (full_neg_search=True). Returns (det_n, det_m) for the first negativity:
    the cited 3696 if not re-derived, or the recomputed value if full_neg_search."""
    print("=== (V) deterministic sanity (mirrors THM-486 build) ===", flush=True)
    W24 = assemble(1, det_signs(1))
    golay = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    print(f"   m=1 Golay match: {W24 == golay}", flush=True)
    W72 = assemble(3, det_signs(3))
    print(f"   m=3 A_16 = {W72.get(16)} (Sloane 249849): {W72.get(16) == 249849}", flush=True)
    print(f"   m=3 all coeffs positive: {all(v > 0 for v in W72.values())}", flush=True)
    ok = all(extremal_cj(m)[1] == -42 * m for m in range(1, 12))
    print(f"   c_1(m) = -42m for m=1..11: {ok}", flush=True)

    # Confirm: deterministic extremal enumerator has NO negativity through window_m.
    det_neg_in_window = None
    for m in range(1, window_m + 1):
        W = assemble(m, det_signs(m))
        if has_negative(W):
            det_neg_in_window = m
            break
    if det_neg_in_window is None:
        print(f"   DETERMINISTIC extremal enumerator: ALL coeffs >= 0 through "
              f"m={window_m} (n={24*window_m}) -- negativity only at the cited "
              f"Zhang threshold n=3696 (m=154), THM-486 part C", flush=True)
    else:
        print(f"   (!) deterministic negativity at m={det_neg_in_window} "
              f"(unexpected vs THM-486 n=3696)", flush=True)

    if full_neg_search:
        first = None
        for m in range(1, max_m_neg):
            W = assemble(m, det_signs(m))
            if has_negative(W):
                first = (24 * m, m, first_negative_weight(W))
                break
        if first:
            n, m, w = first
            print(f"   DETERMINISTIC first negativity (recomputed): n={n} (m={m}), "
                  f"weight {w} -- Zhang n=3696: {n == 3696}", flush=True)
            return n, m
        print(f"   no deterministic negativity in m<{max_m_neg}", flush=True)
        return None, None
    # cite THM-486
    return 3696, 154


# ----------------------------------------------------------------------------
# Random-sign experiment.
# ----------------------------------------------------------------------------

def first_neg_length_for_seed(seed, max_m, mode):
    """Return the smallest m in 1..max_m s.t. the random-sign enumerator at length
    24m has a negative coefficient, or None. mode in {'prefix','fresh'}.
    'prefix': one sign stream per seed, length-m uses first m signs.
    'fresh' : fresh stream per (seed,m)."""
    if mode == 'prefix':
        stream = SignStream(seed)
        long_signs = [1]  # c_0 fixed +1
        for m in range(1, max_m + 1):
            # extend signs to length m+1 (append one sign for the new c_m)
            long_signs.append(stream.sign())
            W = assemble(m, long_signs[:m + 1])
            if has_negative(W):
                return m
        return None
    else:  # fresh
        for m in range(1, max_m + 1):
            stream = SignStream(seed * 1000003 + m)
            signs = [1] + [stream.sign() for _ in range(m)]
            W = assemble(m, signs)
            if has_negative(W):
                return m
        return None


def run_experiment(n_seeds, max_m, mode):
    print(f"\n=== random-sign experiment: mode='{mode}', "
          f"{n_seeds} seeds, max_m={max_m} (max length {24*max_m}) ===", flush=True)
    results = []
    none_count = 0
    t0 = time.time()
    for seed in range(n_seeds):
        m = first_neg_length_for_seed(seed, max_m, mode)
        if m is None:
            none_count += 1
        else:
            results.append(m)
        if (seed + 1) % max(1, n_seeds // 10) == 0:
            print(f"   ... {seed+1}/{n_seeds} seeds, {time.time()-t0:.1f}s", flush=True)
    results.sort()
    n_neg = len(results)
    print(f"   seeds with negativity by m={max_m}: {n_neg}/{n_seeds} "
          f"(never-negative: {none_count})", flush=True)
    if results:
        import statistics
        lengths = [24 * m for m in results]
        med_m = statistics.median(results)
        med_len = 24 * med_m
        mn, mx = min(results), max(results)
        mean_m = statistics.mean(results)
        print(f"   first-negativity m: min={mn}, median={med_m}, "
              f"mean={mean_m:.2f}, max={mx}", flush=True)
        print(f"   first-negativity LENGTH 24m: min={24*mn}, median={med_len}, "
              f"max={24*mx}", flush=True)
        # quantiles
        def q(p):
            idx = min(len(results) - 1, int(p * len(results)))
            return results[idx]
        print(f"   m-quantiles: 1%={q(0.01)}, 10%={q(0.10)}, 25%={q(0.25)}, "
              f"50%={q(0.50)}, 75%={q(0.75)}, 90%={q(0.90)}, 99%={q(0.99)}", flush=True)
        # histogram of m
        from collections import Counter
        hist = Counter(results)
        print("   histogram (m: count):", flush=True)
        for mm in sorted(hist):
            bar = "#" * min(60, hist[mm])
            print(f"      m={mm:3d} (n={24*mm:4d}): {hist[mm]:4d} {bar}", flush=True)
    print(f"   experiment time {time.time()-t0:.1f}s", flush=True)
    return results, none_count


def main():
    t0 = time.time()
    print("############################################################", flush=True)
    print("# random_extremal_enumerator_kps3_0611  (EXPLORATORY)", flush=True)
    print("# code-side mirror of THM-485 pentagonal sign rigidity", flush=True)
    print("# deterministic extremal first-negativity: n=3696 (m=154)", flush=True)
    print("############################################################", flush=True)

    # Parameters: prompt asks m up to ~60, many seeds.
    MAX_M = 60
    N_SEEDS = 400

    det_n, det_m = verify_deterministic(window_m=MAX_M, full_neg_search=False)

    res_prefix, none_p = run_experiment(N_SEEDS, MAX_M, mode='prefix')
    res_fresh, none_f = run_experiment(N_SEEDS, MAX_M, mode='fresh')

    print("\n=== COMPARISON / READOUT (exploratory) ===", flush=True)
    print(f"   deterministic extremal first negativity: n = {det_n} "
          f"(m = {det_n//24 if det_n else '?'})", flush=True)
    for label, res, none_c in (('prefix', res_prefix, none_p),
                               ('fresh', res_fresh, none_f)):
        if res:
            import statistics
            med = statistics.median(res)
            print(f"   [{label}] median random first-negativity: n = {24*med} "
                  f"(m = {med}); fraction reaching negativity by n={24*MAX_M}: "
                  f"{len(res)}/{N_SEEDS}", flush=True)
            if det_n:
                ratio = det_n / (24 * med)
                print(f"   [{label}] deterministic/median ratio = {ratio:.1f}x "
                      f"=> deterministic delays negativity ~{ratio:.0f}x longer "
                      f"than a typical random sign pattern", flush=True)
        else:
            print(f"   [{label}] NO seed reached negativity by m={MAX_M} "
                  f"(deterministic also doesn't until m=154; randomization did "
                  f"NOT accelerate within this window)", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)
    print("NOTE: exploratory evidence, not proof. c_0 kept +1 (normalization). "
          "Magnitudes |c_j| exact big-int; only signs randomized.", flush=True)


if __name__ == "__main__":
    main()
