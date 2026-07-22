#!/usr/bin/env python3
"""lrc14_disproof_search_boxeph_S206.py -- corrected S206 audit

For a finite positive speed set S, write

    M(S) = max_t min_{v in S} ||v t||.

LRC(14) asks whether M(S) >= 1/14 for every thirteen-speed set.  This
program computes the certified rational-witness lower bound

    L_Q(S) = max_{2 <= q <= Q, 1 <= a < q} min_v ||a v/q|| <= M(S).

By THM-1002 (pair-sum-denominator file), section 1, the bound is exact when
Q >= 2 max(S): at a maximum of the piecewise-linear
lower envelope, an active rising branch of ||v t|| meets an active falling
branch of ||w t||.  Hence (v+w)t is integral.  The case v=w is the cusp of a
single triangular wave.  Thus a reduced maximizing denominator divides a
pair sum v+w <= 2 max(S).

Consequently L_Q(S) >= 1/14 is already a rigorous safe-phase certificate,
whether or not the full maximum has been found.  By contrast, L_Q(S) < 1/14
is a disproof only when exactness has separately been certified.

After displaying one strict bounded-scan underestimate, the finite bank is
replayed exactly on every pair-sum ruler.  These rows test AP, near-AP,
dilated, Fibonacci, and random covering shapes.  They do not prove a global
covering minimum, nor do they make anti-golden or higher-order-autocorrelation
heuristics necessary.

Tournament Analysis is intentionally not imposed here.  Runners, residue
classes, gaps, wall crossings, and witness denominators were considered as
vertices, but no canonical binary relation among them preserves the max-min
predicate; forcing a tournament would discard the phase height that decides
the certificate.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd
import random

THRESH = F(1, 14)
DEEPWELL_M = F(14, 183)

def rational_witness_lower_bound(S, Qmax=1500):
    """Return (L_Q(S), argmax (a,q), cutoff_complete)."""
    best = F(0)
    arg = None
    for q in range(2, Qmax + 1):
        for a in range(1, q):
            m = q
            for v in S:
                r = (v * a) % q
                d = min(r, q - r)
                if d < m:
                    m = d
                if m == 0:
                    break
            if m * best.denominator > best.numerator * q:
                best = F(m, q)
                arg = (a, q)
    return best, arg, Qmax >= 2 * max(S)


def exact_pair_sum_maximum(S):
    """Compute M(S) using every numerator over every pair-sum denominator."""
    best_num = 0
    best_den = 1
    arg = None
    denominators = sorted({v + w for i, v in enumerate(S) for w in S[i:]})
    for q in denominators:
        # f(1-t)=f(t), so one representative from each reflected pair suffices.
        for a in range(1, q // 2 + 1):
            m = q
            for v in S:
                r = (v * a) % q
                d = min(r, q - r)
                if d < m:
                    m = d
                if m == 0:
                    break
            if m * best_den > best_num * q:
                best_num = m
                best_den = q
                t = F(a, q)
                arg = (t.numerator, t.denominator)
    return F(best_num, best_den), arg


def is_speed_set(S):
    return len(S) == 13 and len(set(S)) == 13 and all(v > 0 for v in S)


def primitive_normalize(S):
    common = reduce(gcd, S)
    return sorted(v // common for v in S), common


def is_covering(S):
    """Necessary obstruction to the elementary witnesses t=1/q, 2<=q<=14."""
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def first_uncovered_denominator(S):
    return next((q for q in range(2, 15) if not any(v % q == 0 for v in S)), None)


# ---------- exact control: complete breakpoint range ----------
DW = list(range(1, 13)) + [182]
deep_lb, deep_arg, deep_complete = rational_witness_lower_bound(DW)
assert deep_complete and deep_lb == DEEPWELL_M
print("EXACT CONTROL (Q >= 2 max(S) contains every lower-envelope breakpoint):")
print(
    "  deep well {1..12,182}: covering=%s  M=%s (~%.5f) at t*=%s  cutoff_complete=%s"
    % (is_covering(DW), deep_lb, float(deep_lb), deep_arg, deep_complete)
)
print(
    "  M-1/14 = %s (~%.5f), hence this row is certified safe.\n"
    % (deep_lb - THRESH, float(deep_lb - THRESH))
)

# This row proves that an incomplete cutoff can be strictly below the true M.
strict_control = list(range(1, 13)) + [5460]
strict_lower, strict_lower_arg, strict_complete = rational_witness_lower_bound(strict_control, 1200)
strict_M, strict_M_arg = exact_pair_sum_maximum(strict_control)
assert not strict_complete and strict_lower == F(92, 1197)
assert strict_M == F(420, 5461) and strict_lower < strict_M
print("STRICT LOWER-BOUND CONTROL (pair-sum theorem):")
print(
    "  AP12 + far=5460: L_1200=%s at %s < M=%s at %s  cutoff_complete=%s\n"
    % (strict_lower, strict_lower_arg, strict_M, strict_M_arg, strict_complete)
)

# ---------- finite candidate families ----------
print("FINITE CANDIDATE AUDIT (exact pair-sum replay for every row):")
cands = {}
for far in [182, 364, 5040, 5460]:
    S = list(range(1, 13)) + [far]
    cands["AP12 + far=%d" % far] = S
cands["skip12: {1..11,13,182}"] = list(range(1, 12)) + [13, 182]
cands["skip7: {1..6,8..13,182}"] = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 182]
cands["2*AP: {2,4..24,182}"] = list(range(2, 25, 2)) + [182]
fib = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233]
cands["Fibonacci12 + 720720"] = fib + [720720]
cands["Fib-ish covering"] = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 182]

# Replay the original seeded completion rows, but now reject repeats so every
# audited row really has thirteen distinct speeds.
random.seed(2026)
found_cov = 0
for _ in range(400):
    base = sorted(random.sample(range(1, 60), 12))
    completion = 1
    for q in range(2, 15):
        if not any(v % q == 0 for v in base):
            completion = completion * q // gcd(completion, q)
    extra = completion * random.choice([1, 2, 3])
    S13 = sorted(base + [extra])
    if is_speed_set(S13) and is_covering(S13):
        found_cov += 1
        if found_cov <= 6:
            label = "random#%d %s" % (found_cov, S13[:4] + ["..", S13[-1]])
            cands[label] = S13

ranked = []
certified_safe = 0
disproofs = 0
for name, raw_S in cands.items():
    if not is_speed_set(raw_S):
        print("  %-42s INVALID (need 13 distinct positive speeds)" % name)
        continue

    S, dilation = primitive_normalize(raw_S)
    uncovered = first_uncovered_denominator(S)
    exact_M, arg = exact_pair_sum_maximum(S)
    if exact_M >= THRESH:
        verdict = "CERTIFIED SAFE"
        certified_safe += 1
    else:
        verdict = "DISPROOF"
        disproofs += 1

    covering_note = "covering" if uncovered is None else "misses q=%d" % uncovered
    dilation_note = "" if dilation == 1 else " normalized/gcd=%d" % dilation
    print(
        "  %-42s M=%-10s (~%.5f) %-14s witness=%-12s [%s%s]"
        % (
            name,
            str(exact_M),
            float(exact_M),
            verdict,
            str(arg),
            covering_note,
            dilation_note,
        )
    )
    ranked.append((exact_M, name, verdict, arg))

ranked.sort()
print("\n  smallest exact minima in the finite bank:")
for exact_M, name, verdict, arg in ranked[:5]:
    print("    %-42s M=%-10s witness=%s" % (name, exact_M, arg))
print(
    "\n  RESULT: all %d tested valid rows were computed exactly and are safe; disproofs=%d."
    % (certified_safe, disproofs)
)
print("  The finite audit does NOT prove that the deep well is the global covering minimizer.")
print("  Primitivity is used only as a WLOG normalization (M(d*S)=M(S)).")
print("  Anti-golden, far-blocker, and higher-order-autocorrelation conditions remain search heuristics.")
print("  ASSUMPTION CHALLENGE: a tournament quotient was rejected because it forgets resolved phase height.")
