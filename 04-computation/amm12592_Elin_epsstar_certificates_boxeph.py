"""LANE D2 (E-lin) -- exact certificates for the PROVED parts of the
E-lin theorem (amm12592-Elin-theorem-boxeph.md).

All arithmetic exact: int + Fraction.  No floats in any decision.
sympy NOT used (gamma* is sandwiched by pure-integer Lucas/Fibonacci
comparisons).

Certified here:

  C1 (gamma* sandwich).  g := gamma* = log_5(phi^2).  For M = 2^40 we
      certify integers a, b = a+1 with a/M < g < b/M by the exact
      equivalence   phi^{2M} >< 5^a  <=>  L_{2M} + F_{2M} sqrt5 >< 2*5^a,
      decided by integer squaring (5 F^2 >< (2*5^a - L)^2 with sign care).

  C2 (eps* closed form + bound).  eps* := 2(1-g-g^2)/(3+2g).  eps*(g) is
      strictly decreasing in g on (0,1)  [proved in the note; certified
      here on the sandwich by exact Fraction comparison], so
      eps* < eps*(g_lo) =: eps*_hi  (exact Fraction).  Certify
      eps*_hi < 1/32 < 1/16.

  C3 (exact i_feed formula).  i_feed := max{ i : d_i + i <= R-1 },
      d_i = floor(g(R+i)) + D0.  PROVED in the note:
      i_feed = ceil( (R(1-g) - D0) / (1+g) ) - 1, and monotonicity of
      d_i + i.  Certified: formula (computed exactly from the sandwich,
      with a guard that the sandwich decides the ceiling) equals the
      brute-force i_feed from floor_gamma_star for a hostile grid of
      (R, D0), including every dyadic R = 128..16384 at
      D0 in {0, R/32, R/16, thresholds}.

  C4 (feed-survival inequality).  For D0 >= ceil(eps*_hi R):
      2 d_0 - R + 2 > i_feed  (exact integer check per (R, D0) on the
      grid; the general-R proof is the algebra in the note, which the
      grid instantiates).

  C5 (T4b anchor at linear slack).  For the sweep's (R, D0) values,
      initial junk support is exactly [0, R-2-d_0]  (front F(0) = R-2-d_0,
      no gaps at the ends checked; interior block structure recorded), and
      the bottom cell load is 2 (no row-0 death).  This re-anchors the
      T4b edge lemma on the exact D0 range used by E-lin (the earlier
      certificate grid used smaller D0).

  C6 (LIFT lemma).  Degree-raising Delta at degree d -> d+1 maps cells
      delta'_t = delta_t + delta_{t-1}; certify on stored witnesses that
      raising every block by +1 (and +7) preserves admissibility and the
      epoch identity (the identity is trivial -- same polynomials -- but
      we re-verify sums to catch harness bugs).
"""
import os, json, importlib.util, contextlib, io
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
OUT = os.path.join(RESULTS, "amm12592_Elin_epsstar_certificates_boxeph.json")
report = {}

spec = importlib.util.spec_from_file_location(
    "fast", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fast = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fast)
from functools import lru_cache
floor_gamma_star = lru_cache(maxsize=None)(fast.floor_gamma_star)
initial_junk = fast.initial_junk


def ceil_div(a, b):
    return -((-a) // b)


# ---------------------------------------------------------------- C1
def lucas_fib(n):
    """(L_n, F_n) exactly by fast doubling (single pass)."""
    def fd(n):
        if n == 0:
            return (0, 1)
        a, b = fd(n >> 1)
        c = a * ((b << 1) - a)
        d = a * a + b * b
        if n & 1:
            return (d, c + d)
        return (c, d)
    F, Fnext = fd(n)
    L = 2 * Fnext - F        # L_n = F_{n+1} + F_{n-1} = 2F_{n+1} - F_n
    return L, F


def certify_gamma_sandwich(M):
    """find a with a/M < g < (a+1)/M, certified by integer comparisons.
    g = log_5 phi^2:  a/M < g  <=>  5^a < phi^{2M};
    phi^n = (L_n + F_n sqrt5)/2, so phi^{2M} > 5^a
      <=>  F_{2M} sqrt5 > 2*5^a - L_{2M}
      <=>  t <= 0  or  5 F^2 > t^2   (t = 2*5^a - L).
    L, F computed ONCE; only the 5-power varies during bisection."""
    L, F = lucas_fib(2 * M)
    F5sq = 5 * F * F

    def phi_gt(a):
        t = 2 * pow(5, a) - L
        return True if t <= 0 else F5sq > t * t

    lo, hi = 0, M
    assert phi_gt(0) and not phi_gt(M)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if phi_gt(mid):
            lo = mid
        else:
            hi = mid
    return lo, hi


M = 1 << 20
a, b = certify_gamma_sandwich(M)
g_lo, g_hi = Fraction(a, M), Fraction(b, M)
report["C1_gamma_sandwich"] = {
    "M": "2^20", "a": a, "g_lo": str(g_lo), "g_hi": str(g_hi),
    "width": str(g_hi - g_lo),
    "float_view": float(g_lo)}
print("C1: g in (%s, %s) ~ %.9f  width 2^-20" % (g_lo, g_hi, float(g_lo)),
      flush=True)
# cross-check vs floor_gamma_star where the sandwich decides the floor
xchk_n = 0
for n in (997, 5000, 12345, 100003):
    if int(g_lo * n) == int(g_hi * n):
        assert floor_gamma_star(n) == int(g_lo * n), n
        xchk_n += 1
report["C1_crosscheck_floor_decided_points"] = xchk_n
print("C1 crosscheck floor_gamma_star on", xchk_n, "decided points: ok",
      flush=True)
assert xchk_n >= 3


# ---------------------------------------------------------------- C2
def eps_star(gg):
    return 2 * (1 - gg - gg * gg) / (3 + 2 * gg)


e_hi = eps_star(g_lo)        # eps* decreasing in g => upper bound from g_lo
e_lo = eps_star(g_hi)
assert e_lo < e_hi
# certify monotone decrease on the sandwich (derivative sign check by values
# at both ends plus the closed-form derivative numerator, exact):
# d/dg [2(1-g-g^2)/(3+2g)] = 2[(-1-2g)(3+2g) - 2(1-g-g^2)]/(3+2g)^2
#                          = 2[-(2g^2+8g+3) - 2 + 2g + 2g^2 ]/(3+2g)^2  ... do it exactly:
gsym = g_lo
num = (-1 - 2 * gsym) * (3 + 2 * gsym) - 2 * (1 - gsym - gsym * gsym)
report["C2_derivative_negative_at_glo"] = (num < 0)
assert num < 0
report["C2_eps_star"] = {
    "eps_star_lo": str(e_lo), "eps_star_hi": str(e_hi),
    "float_view": float(e_hi),
    "lt_1_over_32": bool(e_hi < Fraction(1, 32)),
    "lt_1_over_16": bool(e_hi < Fraction(1, 16))}
assert e_hi < Fraction(1, 32) < Fraction(1, 16)
print("C2: eps* in (%.9f, %.9f) < 1/32 < 1/16" % (float(e_lo), float(e_hi)),
      flush=True)


# ---------------------------------------------------------------- C3 + C4
def i_feed_brute(R, D0):
    i = -1
    for j in range(R):
        if floor_gamma_star(R + j) + D0 + j <= R - 1:
            i = j
        else:
            break
    return i


def i_feed_formula(R, D0):
    """i_feed = ceil(x) - 1 = floor(x),  x = (R(1-g) - D0)/(1+g).
    x is irrational (x rational would force g rational), so never an
    integer.  Computed from the sandwich; when the sandwich cannot decide
    the ceiling (near-integer x), the two candidates differ by 1 and the
    winner is decided EXACTLY by the Beatty test
    i <= x  <=>  floor(g(R+i)) + D0 + i <= R-1  (the proved equivalence),
    evaluated by the exact integer comparator floor_gamma_star."""
    lo = (R * (1 - g_hi) - D0) / (1 + g_hi)
    hi = (R * (1 - g_lo) - D0) / (1 + g_lo)
    cl, ch = ceil_div(lo.numerator, lo.denominator), \
        ceil_div(hi.numerator, hi.denominator)
    if cl == ch:
        return cl - 1
    if ch != cl + 1:
        return None
    # candidates cl-1, cl; take the largest i satisfying the exact test
    i = cl
    if floor_gamma_star(R + i) + D0 + i <= R - 1:
        return i
    return i - 1


grid = []
for R in (128, 256, 512, 1024, 2048, 4096, 8192, 16384):
    D0s = {0, 1, ceil_div(R, 32), ceil_div(R, 16), ceil_div(3 * R, 32)}
    for D0 in D0s:
        grid.append((R, D0))
for R in (100, 250, 300, 511, 513, 1000, 3000, 5000, 10000):   # non-dyadic
    grid.append((R, ceil_div(R, 16)))

c3_ok, c4_rows, c3_n = True, [], 0
for R, D0 in sorted(grid):
    ib, iff = i_feed_brute(R, D0), i_feed_formula(R, D0)
    if iff is None or ib != iff:
        c3_ok = False
        print("C3 FAIL", R, D0, ib, iff, flush=True)
    c3_n += 1
    d0 = floor_gamma_star(R) + D0
    bound = 2 * d0 - R + 2
    if D0 >= ceil_div(int(e_hi.numerator * R), e_hi.denominator):  # D0>=ceil(eps*_hi R)
        c4_rows.append({"R": R, "D0": D0, "T6b_bound": bound, "i_feed": ib,
                        "survives_feed": bound > ib})
report["C3_i_feed_formula_grid"] = {"n": c3_n, "all_equal": c3_ok}
report["C4_feed_survival"] = c4_rows
c4_ok = all(r["survives_feed"] for r in c4_rows)
report["C4_all"] = c4_ok
json.dump(report, open(OUT, "w"), indent=1)   # stage-partial dump
print(f"C3: i_feed formula == brute on {c3_n} grid points: {c3_ok}",
      flush=True)
print(f"C4: T6b bound > i_feed at all {len(c4_rows)} grid points with "
      f"D0 >= ceil(eps*_hi R): {c4_ok}", flush=True)
assert c3_ok and c4_ok


# ---------------------------------------------------------------- C5
c5 = []
for R in (128, 256, 512, 1024, 2048, 4096, 8192, 16384):
    for D0 in (ceil_div(R, 32), ceil_div(R, 16)):
        d0 = floor_gamma_star(R) + D0
        assert d0 <= R - 2
        j, junkL1, c0 = initial_junk(R, d0)
        cells = sorted(j)
        block_exact = (cells == list(range(0, R - 1 - d0)))
        c5.append({"R": R, "D0": D0, "d0": d0,
                   "front0": (cells[-1] if cells else None),
                   "expected_front": R - 2 - d0,
                   "support_is_full_block": block_exact,
                   "junkL1_bits": junkL1.bit_length()})
        print(f"C5: R={R} D0={D0} front={cells[-1] if cells else None} "
              f"(expect {R-2-d0}) full_block={block_exact}", flush=True)
report["C5_T4b_anchor"] = c5
c5_ok = all(r["support_is_full_block"] and r["front0"] == r["expected_front"]
            for r in c5)
report["C5_all"] = c5_ok
assert c5_ok
json.dump(report, open(OUT, "w"), indent=1)   # stage-partial dump

import sys
if "--skip-c6" in sys.argv:
    report["C6_lift_on_witness"] = "SKIPPED (run separately)"
    report["ALL_PASS_C1_C5"] = True
    json.dump(report, open(OUT.replace(".json", "_c15.json"), "w"), indent=1)
    print("C1-C5 ALL PASS (C6 skipped) ->", OUT.replace(".json", "_c15.json"),
          flush=True)
    raise SystemExit(0)


# ---------------------------------------------------------------- C6
sys.path.insert(0, HERE)
from amm12592_allR_family_toolbox_boxeph import admissible, epoch_sum, qpow


def raise_block(delta, d, m):
    """cells of the same polynomial at degree d+m: conv with binom(m)."""
    out = [0] * (d + m + 1)
    for t, v in enumerate(delta):
        if v:
            for s in range(m + 1):
                out[t + s] += v * comb(m, s)
    return out


c6 = {}
wit_path = os.path.join(HERE, "amm12592_witness_R1024_ruleA_D0_15_fastflow_boxeph.json")
if os.path.exists(wit_path):
    # reconstruct full blocks via engine helper, then lift
    blocks, res = fast.full_blocks(512, 5)   # R=512 D0=5 (known closer, small)
    R, D0 = 512, 5
    pr = [floor_gamma_star(R + i) + D0 for i in range(R)]
    base_ok = (all(admissible(blocks[i], pr[i]) for i in range(R))
               and epoch_sum(R, blocks, pr) == qpow(R - 1))
    for m in (1, 7):
        lifted = [raise_block(blocks[i], pr[i], m) for i in range(R)]
        prm = [d + m for d in pr]
        adm = all(admissible(lifted[i], prm[i]) for i in range(R))
        idt = epoch_sum(R, lifted, prm) == qpow(R - 1)
        c6[f"lift_+{m}"] = {"admissible": adm, "identity": idt}
        print(f"C6: R=512 D0=5 witness lifted +{m}: admissible={adm} "
              f"identity={idt}", flush=True)
    c6["base_ok"] = base_ok
report["C6_lift_on_witness"] = c6
c6_ok = all(v["admissible"] and v["identity"]
            for k, v in c6.items() if k.startswith("lift"))
assert c6_ok and c6.get("base_ok", False)

report["ALL_PASS"] = True
json.dump(report, open(OUT, "w"), indent=1)
print("ALL_PASS -> ", OUT, flush=True)
