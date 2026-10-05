"""Exact controls for the scoped -5 shadow obstruction; standard library only.

The proof is in the matching results note. Checks use explicit exceptions so
normal and optimized Python execute the same controls. No home oracle is used.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def natural(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(name + " must be an exact nonnegative integer")


def v2(value):
    if type(value) is not int or value == 0:
        raise ValueError("v2 needs a nonzero integer")
    value = abs(value)
    return (value & -value).bit_length() - 1


def step(value):
    """Actual odd Collatz step, including negative rational odd denominators."""
    value = F(value)
    if value.numerator % 2 != 1 or value.denominator % 2 != 1:
        raise ValueError("an odd rational is required")
    raw = 3 * value + 1
    if raw == 0:
        raise ValueError("zero has no finite valuation")
    exponent = v2(raw.numerator) - v2(raw.denominator)
    need(exponent >= 1, "odd input gives a positive valuation")
    return raw / (2**exponent), exponent


def stats(word):
    p, q, carry = 1, 1, 0
    for a in word:
        if type(a) is not int or a < 1:
            raise ValueError("valuation letters must be exact positive integers")
        p, carry, q = 3 * p, 3 * carry + q, q * 2**a
    return p, q, carry


def replay(value, word):
    for expected in word:
        value, actual = step(value)
        if actual != expected:
            return None
    return F(value)


@lru_cache(None)
def compositions(total, length):
    if length == 0:
        return ((),) if total == 0 else ()
    if total < length:
        return ()
    return tuple((a,) + tail for a in range(1, total - length + 2)
                 for tail in compositions(total - a, length - 1))


def cell(r_bound, s_bound):
    natural(r_bound, "R")
    natural(s_bound, "S")
    bits = 3 * r_bound // 2 + 1
    q, ternary = 2**bits, 3**s_bound
    residue = ternary * ((-5 * pow(ternary, -1, q)) % q)
    return bits, residue, q * ternary


def exponent_phase(bits):
    """Unique a=3 mod8 phase with 3**a=-5 mod2**bits, for bits>=5."""
    if type(bits) is not int or bits < 5:
        raise ValueError("bits must be an exact integer at least five")
    a, period = 3, 8
    for precision in range(6, bits + 1):
        modulus = 2**precision
        if pow(3, a, modulus) != (-5) % modulus:
            a += period
        period *= 2
        need(pow(3, a, modulus) == (-5) % modulus, "binary lift")
    return a, period


def contraction_word_controls():
    """Every child word under the necessary strict-slope cost bound, r,s<=5."""
    counts = Counter()
    for r in range(6):
        w = tuple(1 + i % 2 for i in range(r))
        p, q, bw = stats(w)
        anchor = F(-5 * p + bw, q)
        need(anchor in (F(-5), F(-7)), "signed anchor phase")
        for s in range(6):
            regime = "r<s" if r < s else "r=s" if r == s else "r>s"
            # lambda < 1 iff 2**D * 3**r < 2**A * 3**s.
            # This is a proved bound on every possible useful child cost,
            # not an arbitrary cutoff on individual valuation exponents.
            dmax = -1
            while 2**(dmax + 1) * p < q * 3**s:
                dmax += 1
            for cost in range(max(s, 0), dmax + 1):
                for v in compositions(cost, s):
                    pv, qv, bv = stats(v)
                    slope = F(p * qv, q * pv)
                    intercept = F(bw * qv - bv * q, q * pv)
                    y = F(anchor * qv - bv, pv)
                    need(slope < 1 and y == -5 * slope + intercept,
                         "independent inverse endpoint")
                    counts[regime] += 1
                    if replay(y, v) != anchor:
                        continue
                    counts["signed_actual"] += 1
                    if intercept.denominator & (intercept.denominator - 1):
                        continue
                    counts["dyadic_intercept"] += 1
                    need(intercept <= -1, "forbidden contracting shadow row")
    return counts


def literal_integer_controls():
    """An independent actual-graph path: no affine-word candidate generation."""
    sources = []
    for rb in range(5):
        for sb in range(5):
            _, residue, modulus = cell(rb, sb)
            for t in range(3):
                n = residue + t * modulus
                if n > 1:
                    sources.append((rb, sb, n))
    maximum = max(n for _, _, n in sources)
    inverse_index = defaultdict(list)
    paths = {}
    for n in range(1, maximum + 1, 2):
        x, word = F(n), ()
        path = [(x, word)]
        inverse_index[0, x].append((n, word))
        for depth in range(1, 5):
            x, a = step(x)
            word += (a,)
            path.append((x, word))
            inverse_index[depth, x].append((n, word))
        paths[n] = path
    joins = 0
    for rb, sb, n in sources:
        for r in range(rb + 1):
            endpoint, w = paths[n][r]
            p, q, bw = stats(w)
            for s in range(sb + 1):
                for h, v in inverse_index[s, endpoint]:
                    if h >= n:
                        continue
                    pv, qv, bv = stats(v)
                    slope = F(p * qv, q * pv)
                    intercept = F(bw * qv - bv * q, q * pv)
                    need(h == slope * n + intercept, "literal join identity")
                    need(intercept <= -1, "literal graph violates shadow theorem")
                    joins += 1
    return len(sources), maximum, len(paths), joins


def phase_and_fuel_controls():
    phases = []
    controls = 0
    for bits in range(5, 19):
        a, period = exponent_phase(bits)
        need(a % 8 == 3 and period == 2**(bits - 2), "phase type")
        need(pow(3, period, 2**bits) == 1, "period")
        need(pow(3, period // 2, 2**bits) != 1, "exact period")
        for j in range(4):
            need(pow(3, a + j * period, 2**bits) == (-5) % 2**bits,
                 "whole exponent progression")
            controls += 1
        if bits in (5, 7, 10, 12, 16, 18):
            phases.append((bits, a, period))
    # Direct complete residue search, independent of the digit-lift algorithm.
    for bits in range(5, 13):
        a, period = exponent_phase(bits)
        found = [j for j in range(period)
                 if pow(3, j, 2**bits) == (-5) % 2**bits]
        need(found == [a], "independent phase uniqueness")
    for valuation in range(1, 21):
        n = (1 if 2**valuation > 5 else 3) * 2**valuation - 5
        x, blocks = F(n), 0
        while True:
            y, a = step(x)
            if a != 1:
                break
            z, b = step(y)
            if b != 2:
                break
            x, blocks = z, blocks + 1
        need(blocks == (valuation - 1) // 3, "exact (12)-block fuel")
        need(x + 5 == F(9, 8)**blocks * (n + 5), "anchor distance")
    return phases, controls


def boundary_controls():
    need(step(F(-11, 3)) == (F(-5), 1), "rational-basin hostile")
    need(step(-3)[0] == -1 and step(-1)[0] == -1, "small signed basin")
    need(replay(3, (1,)) == 5, "outside ternary guard positive control")
    need(F(2, 3) * 5 - F(1, 3) == 3, "predecessor row")
    need(replay(9, (2,)) == 7, "outside binary guard positive control")
    need(F(3, 4) * 9 + F(1, 4) == 7, "direct descent row")
    need(F(9 * 3 + 5, 8) == 4 and replay(3, (1, 2)) is None,
         "integral formal endpoint is not an exact odd-word guard")
    for length in range(1, 9):
        for b in (1, 3):
            x, cost = F(b), 0
            for _ in range(length):
                x, a = step(x)
                cost += a
            expected = 2 * length if b == 1 else 1 if length == 1 else 2 * length + 1
            need(cost == expected and F(2**(cost + 1), 3**length) > 1,
                 "small-positive-intercept contradiction")
    # The successful later-checkpoint family remains present on a different phase.
    for t in range(20):
        n, h = 155 + 4096 * t, 111 + 2916 * t
        w = (1, 2, 1, 1, 1, 2, 3)
        need(replay(n, w) == replay(h, (1,)) and 0 < h < n,
             "positive depth-growing portfolio control")
        need(h == F(729 * n + 669, 1024), "positive family intercept")
    need(pow(3, 483, 4096) == 155, "positive family exponent phase")
    need(pow(3, 27107, 131072) == 155, "original discovery phase")
    invalid = 0
    for thunk in (lambda: cell(True, 0), lambda: cell(0, -1),
                  lambda: exponent_phase(4), lambda: exponent_phase(5.0),
                  lambda: stats((False,)), lambda: v2(0)):
        try:
            thunk()
        except ValueError:
            invalid += 1
        else:
            raise ArithmeticError("invalid exact domain accepted")
    return invalid


def main():
    word_counts = contraction_word_controls()
    integer_counts = literal_integer_controls()
    phases, phase_count = phase_and_fuel_controls()
    invalid = boundary_controls()
    print("PROVED statement in note: bounded depths R,S; affine intercept b>-1.")
    print("FINITE-EXACT word universe: all r,s=0..5 and every positive child word")
    print("under the necessary strict-slope total-cost bound; no individual-letter cap.")
    for key in ("r<s", "r=s", "r>s", "signed_actual", "dyadic_intercept"):
        print(f"word controls {key}: {word_counts[key]}")
    count, maximum, path_count, joins = integer_counts
    print(f"Independent literal graph: R,S=0..4; first three cell members >1.")
    print(f"Source instances {count}; largest source {maximum}; odd paths {path_count}.")
    print(f"Literal smaller joins found {joins}; every one has intercept <=-1.")
    print("Diagonal cells (R=S, bits, residue, modulus):")
    for bound in range(6):
        bits, residue, modulus = cell(bound, bound)
        print(f"  {bound}: {bits}, {residue}, {modulus}")
    print("Shadow power phases (bits, exponent residue, exponent period):")
    for row in phases:
        print(" ", row)
    print(f"Phase progression checks {phase_count}; independent complete phase searches 8.")
    print("Exact (12)-block fuel controls 20; small-positive-intercept controls 16.")
    print("Positive checkpoint-family controls 20; its exponent483 mod1024 is distinct.")
    print("Hostiles: rational -11/3 enters -5; n3 has integral but illegal formal12 endpoint.")
    print("Both missing-guard positive controls pass; empty-word cases are included.")
    print(f"Invalid exact-domain controls rejected {invalid}.")
    print("No arbitrary-intercept no-go, universal coverage, or new home proof claimed.")


if __name__ == "__main__":
    main()
