"""ANGLE B2 -- Lemma T4: closed-form Bernstein-cell decode of the initial error
e_{-1} = 2 G_R = q^{R-1} - E_{R-1}, and its exact margin profile vs the
asymmetric even box  c_t in [-2 binom(d-1,t), +2 binom(d-1,t-1)].

CLAIM (T4).  For 0 <= t <= d <= R-2:
    w_t := decode_d( trunc_d(2 G_R) )_t
         = (-1)^{d-t} binom(R-2-t, d-t) - binom(d+1, t+1) + 2 binom(d, t).
Proof route (generating functions, exact):
    [x^j] 2G_R = 2 (j=0);  (-1)^j binom(R-1,j) - 1  (1 <= j <= R-1).
    decode_d(P)_t = sum_j P_j binom(d-j, t)   [toolbox poly_to_bern].
    sum_{j>=0} (-1)^j binom(R-1,j) binom(d-j,t) = [x^{d-t}](1-x)^{R-2-t}
        = (-1)^{d-t} binom(R-2-t, d-t)                       (Vandermonde)
    sum_{j=0}^{d} binom(d-j,t) = binom(d+1, t+1)             (hockey stick)
    j=0 correction: a_0 = 2 = ((-1)^0 binom(R-1,0) - 1) + 2  -> +2 binom(d,t).

This script (all exact ints):
  1. verifies T4 against the toolbox decode for many (R, d) incl. non-dyadic;
  2. tabulates the exact in-box window and overflow bit-profile of row 0
     for R = 2^k, k = 3..13, at d = floor(gamma*(R)) + D0, D0 in {0,1,2,4,8,16};
  3. reports the saturation depth fraction t_in/d (exact Fractions).
"""
import sys, os, json, io, contextlib, importlib.util
from math import comb
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

from amm12592_allR_family_toolbox_boxeph import psub, qpow, poly_to_bern

spec = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(iref)

def Em(m): return [-1] + [1]*m if m >= 0 else []
def two_G(R): return psub(qpow(R - 1), Em(R - 1))

def w_closed(R, d, t):
    s = d - t
    return ((-1) ** s) * comb(R - 2 - t, s) - comb(d + 1, t + 1) + 2 * comb(d, t)

def verify_T4():
    bad = []
    for R in (6, 8, 12, 16, 25, 32, 64, 100, 128, 256):
        e = two_G(R)
        d0 = iref.floor_gamma_star(R)
        for d in {2, 3, d0 - 2, d0, d0 + 3, R - 2} - {0, 1}:
            if d < 2 or d > R - 2:
                continue
            direct = poly_to_bern(e[:d + 1], d)
            closed = [w_closed(R, d, t) for t in range(d + 1)]
            if direct != closed:
                bad.append((R, d))
    return bad

def profile(R, D0):
    """Exact row-0 margin profile at d = floor(gamma* R) + D0."""
    d = iref.floor_gamma_star(R) + D0
    t_in_lo = None          # least t with cells t..d ALL in box
    worst_bits = 0          # max over cells of overflow bit length
    n_out = 0
    sample = {}
    for t in range(d, -1, -1):
        w = w_closed(R, d, t)
        lo, hi = -2 * comb(d - 1, t), 2 * (comb(d - 1, t - 1) if t else 0)
        inbox = lo <= w <= hi
        if not inbox:
            n_out += 1
            if t_in_lo is None:
                t_in_lo = t + 1
            over = max(w - hi, lo - w)
            worst_bits = max(worst_bits, over.bit_length())
    if t_in_lo is None:
        t_in_lo = 0
    # a few sampled overflow bit lengths at relative depths
    for num, den in ((0, 1), (1, 4), (1, 2), (3, 4)):
        t = (num * d) // den
        w = w_closed(R, d, t)
        lo, hi = -2 * comb(d - 1, t), 2 * (comb(d - 1, t - 1) if t else 0)
        over = max(w - hi, lo - w, 0)
        cap_bits = max(-lo, hi).bit_length()
        sample[f"t={num}/{den}d"] = {"over_bits": over.bit_length(),
                                     "cap_bits": cap_bits}
    return {"R": R, "D0": D0, "d": d, "t_in_lo": t_in_lo,
            "in_frac": str(Fraction(d - t_in_lo + 1, d + 1)),
            "in_frac_milli": (1000 * (d - t_in_lo + 1)) // (d + 1),
            "n_out": n_out, "worst_over_bits": worst_bits,
            "sample": sample}

if __name__ == "__main__":
    bad = verify_T4()
    print(f"T4 closed form vs toolbox decode: "
          f"{'PASS (all R,d probed)' if not bad else f'FAIL {bad}'}", flush=True)
    tab = []
    for k in range(3, 14):
        R = 2 ** k
        for D0 in (0, 1, 2, 4, 8, 16):
            if iref.floor_gamma_star(R) + D0 > R - 2:
                continue        # T4 needs d <= R-2
            p = profile(R, D0)
            tab.append(p)
            print(f"R={R:5d} D0={D0:2d} d={p['d']:5d}: in-window t>={p['t_in_lo']:5d} "
                  f"(top {p['in_frac_milli']/10:.1f}% of cells in-box), "
                  f"worst overflow bits {p['worst_over_bits']}", flush=True)
    json.dump({"T4_verified": not bad, "profiles": tab},
              open(os.path.join(RESULTS,
                   "amm12592_transient_initial_decode_boxeph.json"), "w"),
              indent=1)
    print("saved profiles json", flush=True)
