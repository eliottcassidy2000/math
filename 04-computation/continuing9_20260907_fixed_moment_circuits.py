"""Exact fixed-first-two-moment Newton-circuit controls.

The all-degree conclusions are proved in the companion report. This finite
certificate reconstructs the q-polynomials, explicit perturbation radii,
all ternary words in the stated small banks, cubic image, and strict-Newton
hostiles. No root approximation is accepted as a sign certificate.
"""
from __future__ import annotations
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import sympy as sp

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
STEM = Path(__file__).stem
X, q, a, v, p = sp.symbols("X q a v p")
gates = 0


def need(test, label):
    global gates
    gates += 1
    if not bool(test):
        raise RuntimeError(label)


def center(d, qq):
    return [sp.binomial(d, k) * qq ** math.comb(k, 2) for k in range(d + 1)]


def poly(aa):
    return sp.Poly(sum(c * X**k for k, c in enumerate(aa)), X)


def eval_alt(aa, xx):
    return sum(c * (-xx)**k for k, c in enumerate(aa))


def ratios(aa):
    d = len(aa) - 1
    hh = [c / sp.binomial(d, k) for k, c in enumerate(aa)]
    rr = [sp.cancel(hh[k]**2 / (hh[k - 1] * hh[k + 1])) for k in range(1, d)]
    cc = [sp.cancel(rr[k] / rr[k - 1]) for k in range(1, d - 1)]
    return rr, cc


def chart(d, qq, cc):
    hh = []
    for k in range(d + 1):
        h = qq ** math.comb(k, 2)
        for j in range(2, k):
            h /= cc[j - 2] ** math.comb(k - j + 1, 2)
        hh.append(sp.cancel(h))
    return [sp.binomial(d, k) * hh[k] for k in range(d + 1)]


def encode(value):
    if isinstance(value, dict):
        return {str(k): encode(vv) for k, vv in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(vv) for vv in value]
    if isinstance(value, (sp.Basic,)):
        return str(value)
    return value


certificate = {"scope": "fixed S,Q circuit box and exact cubic image", "cases": []}

# Polynomial identities are checked coefficientwise, not by floating samples.
for d in range(1, 11):
    pp = poly(center(d, q)).as_expr()
    nxt = poly(center(d + 1, q)).as_expr()
    need(sp.expand(nxt - pp - X * pp.subs(X, q * X)) == 0, "q recurrence")
    need(sp.expand(sp.diff(pp, X) - d * poly(center(d - 1, q)).as_expr().subs(X, q * X)) == 0,
         "q derivative")
    rr, cc = ratios(center(d, q))
    need(all(sp.cancel(r - 1/q) == 0 for r in rr), "center R")
    need(all(c == 1 for c in cc), "center C")
    for k in range(d + 1):
        need(sum(math.comb(k - j + 1, 2) for j in range(2, k)) == math.comb(k, 3),
             "chart total exponent")

# A separate recursive reconstruction checks the closed coefficient chart.
for d in range(3, 9):
    cc = [sp.Rational(j + 2, j + 1) for j in range(d - 2)]
    hh = [sp.S.One, sp.S.One]
    r = sp.Rational(7, 5)
    for k in range(1, d):
        if k >= 2:
            r *= cc[k - 2]
        hh.append(sp.cancel(hh[-1]**2 / (r * hh[-2])))
    aa = chart(d, sp.Rational(5, 7), cc)
    need(aa == [sp.binomial(d, k)*hh[k] for k in range(d + 1)], "chart via recurrence")
    need(ratios(aa)[1] == cc, "chart exact C")

# Explicit center root isolators and whole-box sign samples. Rational intervals
# are a finite control of the separately proved q-separation induction.
banks = [(3, sp.Rational(3, 4)), (4, sp.Rational(7, 8)),
         (5, sp.Rational(275, 338)), (6, sp.Rational(1, 2)),
         (3, sp.Rational(1, 4)), (4, sp.Rational(9, 10))]
word_count = 0
for d, qq in banks:
    aa0 = center(d, qq)
    ff = sp.Poly(poly(aa0).as_expr().subs(X, -X), X)
    raw_intervals = ff.intervals(eps=sp.Rational(1, 10**10))
    need(len(raw_intervals) == d, "all center roots isolated")
    intervals = []
    for (lo, hi), multiplicity in raw_intervals:
        need(multiplicity == 1 and lo > 0 and hi >= lo, "positive simple center root")
        intervals.append((lo, hi))
    need(all(intervals[i+1][0] > intervals[i][1] / qq for i in range(d - 1)),
         "strict q separation")
    samples = [(intervals[i-1][1] + intervals[i][0])/2 for i in range(1, d)]
    samples.append(intervals[-1][1] + 1)
    values = [eval_alt(aa0, xx) for xx in samples]
    weights = [sum(math.comb(k, 3)*aa0[k]*xx**k for k in range(3, d + 1))
               for xx in samples]
    for i, val in enumerate(values, 1):
        need((-1)**i * val > 0, "alternating center samples")
    raw_bound = min([sp.Rational(1, 2), sp.Rational(1, 2*math.comb(d, 3))]
                    + [abs(val)/(4*ww) for val, ww in zip(values, weights)])
    delta = sp.Rational(1, 2)
    while delta > raw_bound:
        delta /= 2
    need(delta > 0 and delta <= raw_bound, "explicit dyadic radius")
    targets = []
    for word in itertools.product((-1, 0, 1), repeat=d - 2):
        word_count += 1
        cc = [1 + delta*w for w in word]
        aa = chart(d, qq, cc)
        need(aa[:3] == aa0[:3], "first two moments fixed")
        need(ratios(aa)[1] == cc, "target circuit exact")
        for k in range(3, d + 1):
            need(abs(aa[k]/aa0[k] - 1) <= 2*math.comb(k, 3)*delta,
                 "relative coefficient budget")
        vals = [eval_alt(aa, xx) for xx in samples]
        need(all((-1)**i * val > 0 for i, val in enumerate(vals, 1)),
             "every ternary word has all roots simple negative")
        need(all(abs(val-base) <= abs(base)/2 for val, base in zip(vals, values)),
             "sample sign margin")
        targets.append({"word": word, "C": cc, "coefficients": aa})
    case = {"d": d, "q": qq, "center_coefficients": aa0,
            "positive_roots_of_P_minus_X": intervals, "samples": samples,
            "sample_values": values, "weights": weights, "raw_bound": raw_bound,
            "delta": delta, "targets": targets}
    certificate["cases"].append(case)
    print(f"box d={d} q={qq}: delta={delta}; all {3**(d-2)} ternary words PASS")

# Exact cubic circuit image in the normalization S=3.
cubic = X**3 - 3*X**2 + 3*(1-v)*X - p
disc = sp.discriminant(cubic, X)
need(sp.expand(disc - 27*(4*v**3 - (p-(1-3*v))**2)) == 0, "cubic discriminant")
lo = (1+a)**2 * (1-2*a)
hi = (1-a)**2 * (1+2*a)
need(sp.expand(lo - (1-3*a*a-2*a**3)) == 0, "lower product endpoint")
need(sp.expand(hi - (1-3*a*a+2*a**3)) == 0, "upper product endpoint")
cmin = sp.cancel((1-a*a)**3/hi)
cmax = sp.cancel((1-a*a)**3/lo)
need(sp.cancel(cmin-1 + a**3*(2+a)/(1+2*a)) == 0, "lower circuit below one")
need(sp.cancel(cmax-1 - a**3*(2-a)/(1-2*a)) == 0, "upper circuit above one")
need(sp.cancel(cmin/(1-a*a) - (1+a)**2/(1+2*a)) == 0, "strict improvement over Newton")
cubic_controls = []
for aval in [sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(3, 4)]:
    vv, qq = aval**2, 1-aval**2
    lval, uval = lo.subs(a, aval), hi.subs(a, aval)
    pmid = (max(0, lval)+uval)/2
    pp = sp.Poly(cubic.subs({v: vv, p: pmid}), X)
    need(pp.discriminant() > 0 and pp.count_roots(0, sp.oo) == 3, "cubic interior iff")
    endpoint = sp.Poly(cubic.subs({v: vv, p: uval}), X)
    need(endpoint.discriminant() == 0 and endpoint.count_roots(0, sp.oo) == 2,
         "cubic positive repeated endpoint")
    need(sp.expand(endpoint.as_expr() - (X-(1-aval))**2*(X-(1+2*aval))) == 0,
         "upper product endpoint roots")
    need(sp.expand(cubic.subs({v: vv, p: lval})
                   - (X-(1+aval))**2*(X-(1-2*aval))) == 0,
         "lower product endpoint roots")
    cubic_controls.append({"a": aval, "q": qq, "product_lower": lval,
                           "product_upper": uval, "interior_product": pmid})

# Strict Newton inequalities do not characterize fixed-moment feasibility.
h3 = [sp.S.One, sp.Integer(3), sp.Rational(9, 4), sp.Rational(135, 256)]
h3N = sp.Poly(sum(ee*X**(3-k) for k, ee in enumerate(h3)), X)
r3, c3 = ratios(h3)
need(r3 == [sp.Rational(4, 3), sp.Rational(16, 15)], "cubic hostile all R")
need(c3 == [sp.Rational(4, 5)], "cubic hostile C")
need(h3N.discriminant() == -sp.Rational(25515, 65536), "cubic hostile not real rooted")
need(sp.Rational(4, 5) < sp.Rational(27, 32), "cubic hostile outside exact image")
h4N = sp.Poly(h3N.as_expr()*(X+1), X)
h4 = h4N.all_coeffs()
r4, c4 = ratios(h4)
need(h4 == [1, 4, sp.Rational(21, 4), sp.Rational(711, 256), sp.Rational(135, 256)],
     "quartic hostile factor")
need(all(rr > 1 for rr in r4), "quartic hostile all strict Newton")
need(h4[1]**2-2*h4[2] == sp.Rational(11, 2), "quartic hostile moments")
need(h4N.count_roots(-sp.oo, 0) == 2, "quartic hostile exactly two negative real roots")
certificate["hostiles"] = {"cubic_e": h3, "cubic_R": r3, "cubic_C": c3,
                            "cubic_discriminant": h3N.discriminant(),
                            "quartic_e": h4, "quartic_R": r4, "quartic_C": c4}
certificate["cubic"] = {"discriminant": disc, "product_lower": lo, "product_upper": hi,
                        "circuit_lower": cmin, "circuit_upper_if_a_lt_half": cmax,
                        "controls": cubic_controls}
print(f"strict-Newton cubic hostile: R={r3}, C={c3}, discriminant={h3N.discriminant()}")
print(f"strict-Newton quartic hostile: R={r4}, C={c4}; negative real roots=2")

# Actual degree-five first-two-moment anchor only: no C/D interlacer claim.
anchor = [ee*sp.Rational(13, 5)**k for k, ee in enumerate(center(5, sp.Rational(275, 338)))]
need(anchor[1] == 13 and anchor[1]**2-2*anchor[2] == 59, "actual S13 Q59")
need(ratios(anchor)[0] == [sp.Rational(338, 275)]*4, "actual anchor R")
need(ratios(anchor)[1] == [1]*3, "actual anchor ties")
certificate["anchor_S13_Q59"] = anchor
print("anchor S=13,Q=59: q=275/338, every R=338/275, every C=1")
print(f"finite universe: {len(banks)} boxes; {word_count} complete ternary words")
certificate["gates"] = gates + 1
payload = (json.dumps(encode(certificate), sort_keys=True, indent=2) + "\n").encode("utf-8")
outpath = HERE / (STEM + "_certificate.json")
outpath.write_bytes(payload)
need(b"\r" not in payload, "JSON LF")
print(f"certificate SHA256 {hashlib.sha256(payload).hexdigest()}")
print(f"PASS {gates} always-active exact gates")
