#!/usr/bin/env python3
"""
codex-S109: the zeta -1/12 clue as a prove/disprove dialectic.

The user prompt asks to hold together:

  M({1..11,13}) = 1/12
  1 + 2 + 3 + ...  ->  -1/12

Incoming KPS S31p frames this as a resonance-killing game.  This script makes
the analogy exact enough to be useful for the first missing-killer core, then
tries to turn it into a counterexample.  The counterexample attempt fails for a
discrete reason:
the only one-speed extension of C={1..11,13} that kills both q-witnesses is
w=84m, and that infinite family has the exact formula

  M(C union {84m}) = 7m/(84m+5) > 1/14.

So the formal negative constant points to a real boundary resonance, but finite
positive integer speeds convert it into divisibility constraints and binding
pair margins.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial


F = Fraction
N = 14
CORE = tuple(range(1, 12)) + (13,)


def norm_frac(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def envelope_candidates(speeds: tuple[int, ...]) -> set[F]:
    speeds = tuple(sorted(set(speeds)))
    out: set[F] = {F(0), F(1, 2)}
    for v in speeds:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2):
                break
            out.add(t)
            k += 1
    for i, a in enumerate(speeds):
        for b in speeds[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                for k in range(1, d // 2 + 1):
                    out.add(F(k, d))
    return out


def lonely_constant(speeds: tuple[int, ...]) -> tuple[F, F, tuple[int, ...]]:
    best = F(-1)
    arg = F(0)
    active: tuple[int, ...] = ()
    for t in envelope_candidates(speeds):
        vals = [(norm_frac(F(v) * t), v) for v in speeds]
        val = min(x for x, _ in vals)
        if val > best:
            best = val
            arg = t
            active = tuple(v for x, v in vals if x == val)
    return best, arg, active


def covering_deficits(speeds: tuple[int, ...], n: int = N) -> tuple[int, ...]:
    return tuple(q for q in range(2, n + 1) if all(v % q for v in speeds))


def series_mul(a: list[F], b: list[F], deg: int) -> list[F]:
    out = [F(0) for _ in range(deg + 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j <= deg:
                out[i + j] += x * y
    return out


def series_inv(a: list[F], deg: int) -> list[F]:
    out = [F(0) for _ in range(deg + 1)]
    out[0] = 1 / a[0]
    for k in range(1, deg + 1):
        out[k] = -sum(a[i] * out[k - i] for i in range(1, k + 1) if i < len(a)) / a[0]
    return out


def exp_series(deg: int) -> list[F]:
    return [F(1, factorial(k)) for k in range(deg + 1)]


def zeta_minus_one_expansion(deg: int = 8) -> list[F]:
    """Return coefficients of e^x/(e^x-1)^2 = x^-2 * P(x)."""

    # A(x)=x/(e^x-1)=1/((e^x-1)/x).
    em1_over_x = [F(1, factorial(j + 1)) for j in range(deg + 1)]
    A = series_inv(em1_over_x, deg)
    return series_mul(exp_series(deg), series_mul(A, A, deg), deg)


def triangular(k: int) -> int:
    return k * (k + 1) // 2


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def print_zeta_anchor() -> None:
    print("=" * 96)
    print("A. The formal -1/12 anchor")
    print("=" * 96)
    coeffs = zeta_minus_one_expansion(8)
    print("Exact formal expansion:")
    print("  sum_{n>=1} n*exp(-epsilon*n) = epsilon^-2 + c0 + c2*epsilon^2 + ...")
    print(f"  c0 = {coeffs[2]}")
    print(f"  c2 = {coeffs[4]}")
    print(f"  c4 = {coeffs[6]}")
    print("Readout: the zeta-regularized negative mass is a Bernoulli boundary")
    print("constant.  It is not a finite negative runner; a finite LRC row has to")
    print("realize any cancellation through divisibility and binding pairs.")
    print()


def print_core() -> None:
    print("=" * 96)
    print("B. The positive 1/12 core")
    print("=" * 96)
    m, t, active = lonely_constant(CORE)
    print(f"C = {CORE}")
    print(f"M(C) = {fmt(m)}, attained at t={t}, active={active}")
    print("Proof sketch for the upper bound: speeds 1..11 make the 12 points")
    print("{0,t,2t,...,11t}; if all distances were >1/12, twelve circle gaps")
    print("would sum to >1.  The witness t=5/12 gives the lower bound, and speed")
    print("13 also lands at distance 5/12.")
    print()


def print_one_tail_classification() -> None:
    print("=" * 96)
    print("C. Prove/disprove dialectic for one added integer speed")
    print("=" * 96)
    print("Disproof pressure: add one speed w to kill the 1/12 core witness.")
    print("Proof response: the q-witness split is exhaustive.")
    print()
    print("  if 12 does not divide w:  q=12 witness remains, M>=1/12.")
    print("  if 14 does not divide w:  q=14 witness remains, M>=1/14.")
    print("  only if 84 divides w:     the row is genuinely q-covering.")
    print()
    print("For the covering tail w=84m there is an exact binding-pair formula:")
    print("  t_m = (35m+2)/(84m+5)")
    print("  M(C union {84m}) = 7m/(84m+5), active pair (5,84m)")
    print()
    print("Finite affine residue proof table at t_m.  Each numerator distance is")
    print("measured modulo D=84m+5 and is >=7m; equality occurs for speed 5.")
    table = {
        1: "35m+2",
        2: "14m+1",
        3: "21m+1",
        4: "28m+2",
        5: "7m",
        6: "42m+2",
        7: "7m+1",
        8: "28m+1",
        9: "21m+2",
        10: "14m",
        11: "35m+3",
        13: "35m+1",
        "84m": "7m",
    }
    for v, expr in table.items():
        print(f"  speed {str(v):>3}: distance numerator {expr}")
    print("Thus the infinite covering one-tail family is strictly lonely:")
    print("  7m/(84m+5) > 1/14  <=>  98m > 84m+5  <=>  14m>5.")
    print()


def print_exact_scans() -> None:
    print("=" * 96)
    print("D. Finite/discrete families: counterexample search versus proof split")
    print("=" * 96)
    target = F(1, 14)
    counts = {"q12": 0, "q14": 0, "covering_84m": 0}
    examples: dict[str, list[int]] = {key: [] for key in counts}
    for w in range(1, 337):
        if w in CORE:
            continue
        if w % 12:
            key = "q12"
        elif w % 14:
            key = "q14"
        else:
            key = "covering_84m"
        counts[key] += 1
        if len(examples[key]) < 8:
            examples[key].append(w)
    print("Classification of all one-tail rows C union {w}, 1<=w<=336:")
    print(f"  q=12 witness survives: {counts['q12']:3d} rows, examples {examples['q12']}")
    print(f"  q=14 witness survives: {counts['q14']:3d} rows, examples {examples['q14']}")
    print(
        f"  covering tail w=84m:   {counts['covering_84m']:3d} rows, "
        f"examples {examples['covering_84m']}"
    )
    print("This is a proof split, not a heuristic scan: every one-tail integer is")
    print("in exactly one of the three rows above, and each row has a certified")
    print("lower bound >=1/14.")
    print()

    audit_ws = (12, 24, 36, 84, 168, 252, 336)
    print("Exact envelope audit on representative disproof/proof rows:")
    for w in audit_ws:
        speeds = tuple(sorted(CORE + (w,)))
        m, t, active = lonely_constant(speeds)
        deficits = covering_deficits(speeds)
        tag = "TIGHT" if m == target else ("COUNTER?" if m < target else "loose")
        formula = F(7 * (w // 84), w + 5) if w % 84 == 0 else None
        extra = f", formula={formula}" if formula is not None else ""
        print(
            f"  w={w:3d}: M={str(m):>7s}, t={str(t):>7s}, "
            f"active={active!s:<10s}, deficits={deficits!s:<10s}{extra} {tag}"
        )
    print()

    print("Discrete zeta-flavored tails.  Triangular/Faulhaber choices do not")
    print("create negative mass; only their divisibility class matters.")
    families = [
        ("T_k", lambda k: triangular(k)),
        ("12*T_k", lambda k: 12 * triangular(k)),
        ("14*T_k", lambda k: 14 * triangular(k)),
        ("84*T_k", lambda k: 84 * triangular(k)),
    ]
    for name, fn in families:
        best: tuple[F, int, int, str] | None = None
        for k in range(1, 18):
            w = fn(k)
            if w in CORE:
                continue
            if w % 12:
                bound = F(1, 12)
                reason = "q=12"
            elif w % 14:
                bound = F(1, 14)
                reason = "q=14"
            else:
                m = w // 84
                bound = F(7 * m, 84 * m + 5)
                reason = "84m formula"
            if best is None or bound < best[0]:
                best = (bound, k, w, reason)
        assert best is not None
        m, k, w, reason = best
        print(f"  {name:<8s}: best k={k:2d}, w={w:5d}, certified M>={m} via {reason}")
    print()


def print_tournament_analysis() -> None:
    print("=" * 96)
    print("Tournament Analysis")
    print("=" * 96)
    vertices = [
        "formal_zeta_negative_constant",
        "finite_divisibility_realization",
        "q12_witness",
        "q14_witness",
        "covering_tail_84m",
        "binding_pair_formula",
        "multi_large_resonant_middle",
    ]
    print("Vertices are proof/disproof mechanisms, not runners.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("Observable: whether a formal cancellation survives as a finite integer")
    print("counterexample.  The switch is q-covering: non-covering rows are killed")
    print("by q-witnesses; the one-tail covering quotient has a positive formula.")
    print("Destroyed information: arbitrary multi-tail resonance geometry.  That is")
    print("the remaining OPEN-Q-108 middle, not this one-tail zeta boundary.")


def main() -> None:
    print("codex-S109 zeta -1/12 / LRC 1/12 core dialectic")
    print_zeta_anchor()
    print_core()
    print_one_tail_classification()
    print_exact_scans()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
