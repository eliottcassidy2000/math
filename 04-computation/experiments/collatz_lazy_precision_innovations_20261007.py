"""Exact joint dyadic atoms and lazy-reveal rules for paired odd 2-adic orbits.

U(x)=(3x+1)/2^v2(3x+1); ROOT is not made absorbing in this Haar model.
No orbit-independence, recurrence, or positive-integer convergence assertion.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


def integer(n, minimum=0):
    need(type(n) is int and n >= minimum, "exact integer outside domain")
    return n


def valuation(n):
    need(type(n) is int and n != 0, "nonzero exact integer")
    n = abs(n)
    return (n & -n).bit_length() - 1


def data(word):
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         "tuple of positive exact valuations")
    p, q, carry = 1, 1, 0
    for a in word:
        p, q, carry = 3 * p, q * (1 << a), 3 * carry + q
    residue = (q - carry) * pow(p, -1, 2 * q) % (2 * q)
    return p, q, carry, residue


def joint_atom(wy, wx, alpha=3, v=1, beta=1):
    """The exact Haar source cylinder for both prefixes, or None."""
    need(type(alpha) is int and alpha % 2 == 1, "odd exact affine unit")
    integer(v)
    need(type(beta) is int and beta % 2 == (1 if v else 0), "odd-output offset")
    py, qy, cy, ry = data(wy)
    px, qx, cx, rx = data(wx)
    a, b = qy.bit_length() - 1, qx.bit_length() - 1
    m = max(a, b - v)
    modulus = 1 << (m + 1)
    residue = ry
    if b + 1 > v:
        difference = rx - beta
        if difference % (1 << v):
            return None
        xmod = 1 << (b + 1 - v)
        pull = (difference // (1 << v)) * pow(alpha, -1, xmod) % xmod
        if b - v > a:
            residue = pull
    residue %= modulus
    if (residue - ry) % (2 * qy) or (alpha * (1 << v) * residue + beta - rx) % (2 * qx):
        return None
    ycur = (py * residue + cy) // qy
    xcur = (px * (alpha * (1 << v) * residue + beta) + cx) // qx
    return dict(residue=residue, precision=m, mass=F(1, 1 << m),
                a=a, b=b, y_residue=ycur % (1 << (m - a + 1)),
                x_residue=xcur % (1 << (m + v - b + 1)),
                y_bits=m - a, x_bits=m + v - b)


def query_rule(wy, wx, side, alpha=3, v=1, beta=1):
    """Recompute guards; no caller-supplied atom is trusted."""
    need(type(side) is str and side in ("y", "x"), "query side y or x")
    atom = joint_atom(wy, wx, alpha, v, beta)
    if atom is None:
        return None
    bits, residue = atom[side + "_bits"], atom[side + "_residue"]
    z = (3 * residue + 1) % (1 << (bits + 1))
    if z:
        return dict(kind="determined", valuation=valuation(z), bits=bits,
                    precision=atom["precision"], source_residue=atom["residue"])
    return dict(kind="fresh", shift=bits, bits=bits,
                precision=atom["precision"], source_residue=atom["residue"])


def extend(wy, wx, side, observed, alpha=3, v=1, beta=1):
    integer(observed, 1)
    rule = query_rule(wy, wx, side, alpha, v, beta)
    if rule is None:
        return None
    if rule["kind"] == "determined" and observed != rule["valuation"]:
        return None
    if rule["kind"] == "fresh" and observed <= rule["shift"]:
        return None
    return joint_atom(wy + (observed,) if side == "y" else wy,
                      wx + (observed,) if side == "x" else wx, alpha, v, beta)


def nb_above(count, threshold):
    """P(sum of count geometric(1/2) variables > threshold), exact."""
    integer(count)
    integer(threshold)
    if count == 0:
        return F(0)
    if threshold < count:
        return F(1)
    return F(sum(comb(threshold, k) for k in range(count)), 1 << threshold)


def nb_below(count, threshold):
    """Strict P(sum < threshold), including the count0 point mass."""
    integer(count)
    integer(threshold)
    if threshold == 0:
        return F(0)
    return 1 - nb_above(count, threshold - 1)


def overlap_bound(s, t, v, threshold):
    """Upper bound for P(A_s > B_(t-1)-v), no independence premise."""
    integer(s, 1)
    integer(t, 1)
    integer(v)
    integer(threshold)
    return min(F(1), nb_above(s, threshold)
               + (1 << v) * nb_below(t - 1, v + threshold))


def covariance_squared_bound(s, t, v, threshold):
    integer(v)
    need(type(t) is int and t >= v + 1, "X marginal must have flushed initial coset")
    return min(F(4), 38 * overlap_bound(s, t, v, threshold))


def literal(n, word):
    need(type(n) is int and n % 2, "signed odd exact integer")
    data(word)
    for a in word:
        z = 3 * n + 1
        need(valuation(z) == a, "literal prefix mismatch")
        n = z // (1 << a)
    return n


def main():
    checks = 0

    def check(condition, message):
        nonlocal checks
        checks += 1
        if not condition:
            raise RuntimeError(message)

    words = [()]
    for length in (1, 2):
        words += list(product(range(1, 4), repeat=length))
    couplings = [(3, 1, 1), (3, 2, 1), (-1, 0, 0), (5, 0, 2), (1, 3, 5)]
    atoms = fresh = determined = 0
    for alpha, v, beta in couplings:
        for wy in words:
            for wx in words:
                atom = joint_atom(wy, wx, alpha, v, beta)
                if atom is None:
                    continue
                atoms += 1
                m, c = atom["precision"], atom["residue"]
                check(min(atom["y_bits"], atom["x_bits"]) == 0, "one frontier exists")
                for side in ("y", "x"):
                    rule = query_rule(wy, wx, side, alpha, v, beta)
                    if rule["kind"] == "determined":
                        determined += 1
                        child = extend(wy, wx, side, rule["valuation"], alpha, v, beta)
                        check(child["precision"] == m and child["mass"] == atom["mass"],
                              "old-bit query has no information cost")
                        check(extend(wy, wx, side, rule["valuation"] + 1, alpha, v, beta) is None,
                              "wrong determined value rejected")
                    else:
                        fresh += 1
                        check(extend(wy, wx, side, max(1, rule["shift"]), alpha, v, beta) is None
                              if rule["shift"] else True, "fresh lower boundary")
                        for g in range(1, 6):
                            child = extend(wy, wx, side, rule["shift"] + g, alpha, v, beta)
                            check(child is not None and child["precision"] == m + g,
                                  "fresh precision increment")
                            check(child["mass"] / atom["mass"] == F(1, 1 << g),
                                  "conditional geometric law")
                    counts = {}
                    for j in range(32):
                        y = c + (1 << (m + 1)) * j
                        x = alpha * (1 << v) * y + beta
                        yend, xend = literal(y, wy), literal(x, wx)
                        n = yend if side == "y" else xend
                        a = valuation(3 * n + 1)
                        counts[a] = counts.get(a, 0) + 1
                        check(n % (1 << (rule["bits"] + 1)) == atom[side + "_residue"],
                              "independent signed literal coset")
                    if rule["kind"] == "determined":
                        check(counts == {rule["valuation"]: 32}, "literal determined query")
                    else:
                        for g in range(1, 6):
                            check(counts.get(rule["shift"] + g, 0) == 1 << (5 - g),
                                  "literal fresh geometric frequencies")

    # Same consumed-bit data, different next-query type.
    r1 = query_rule((1, 2), (1,), "x")
    r2 = query_rule((2, 1), (1,), "x")
    check(r1["precision"] == r2["precision"] == 3 and r1["bits"] == r2["bits"] == 3,
          "same information budget")
    check(r1["kind"] == "fresh" and r1["shift"] == 3
          and r2["kind"] == "determined" and r2["valuation"] == 1,
          "residue sidecar is essential")

    # v1: A1 versus B2. Exactly one is1, with geometric marginals.
    for a, b in product(range(1, 7), repeat=2):
        atom = joint_atom((a,), (1, b))
        expected = F(1, 1 << max(a, b)) if (a == 1) != (b == 1) else F(0)
        check((atom["mass"] if atom else F(0)) == expected, "near-diagonal dependence")

    # Negative-binomial formulas checked by an independent sum of point masses.
    for count in range(7):
        for threshold in range(21):
            direct = F(int(count == 0 and threshold > 0)) if count == 0 else sum(
                (F(comb(total - 1, count - 1), 1 << total)
                 for total in range(count, threshold)), F(0))
            check(nb_below(count, threshold) == direct, "strict lower-tail boundary")
    check(nb_below(0, 0) == 0 and nb_below(0, 1) == 1, "zero-length boundary")

    hostiles = [lambda: joint_atom([], ()), lambda: joint_atom((True,), ()),
                lambda: joint_atom((), (), 3.0, 1, 1), lambda: joint_atom((), (), 3, True, 1),
                lambda: joint_atom((), (), 3, 0, 1), lambda: query_rule((), (), True),
                lambda: covariance_squared_bound(1, 1, 1, 2), lambda: nb_below(2, -1)]
    for hostile in hostiles:
        try:
            hostile()
        except ValueError:
            check(True, "typed hostile")
        else:
            check(False, "typed hostile unexpectedly accepted")

    bounds = []
    for s, t, v in ((1, 20, 1), (10, 31, 1), (20, 61, 1), (40, 121, 1),
                    (20, 21, 1), (40, 41, 1)):
        h = s + t - 1
        bounds.append(dict(s=s, t=t, v=v, threshold=h,
                           overlap=str(overlap_bound(s, t, v, h)),
                           covariance_squared=str(covariance_squared_bound(s, t, v, h))))
    print("PROVED lazy precision innovations; no synchronous recurrence conclusion")
    print("universe=5 affine couplings; both prefix lengths0..2; valuations1..3")
    print("compatible_atoms=" + str(atoms) + "; fresh_queries=" + str(fresh)
          + "; determined_queries=" + str(determined))
    print("fresh_law=P(G=j)=2^-j; precision_increment=G; stopped-stream_scope=retained")
    print("same_budget_hostile=" + json.dumps([r1, r2], sort_keys=True))
    print("v1_A1_B2: exactly_one_equals1; E_product=3; covariance=-1")
    print("off_diagonal_bounds=" + json.dumps(bounds, sort_keys=True))
    print("typed_hostiles=" + str(len(hostiles)))
    print("checks=" + str(checks))


if __name__ == "__main__":
    main()
