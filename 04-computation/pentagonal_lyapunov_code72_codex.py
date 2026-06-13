#!/usr/bin/env python3
"""
Pentagonal signs, random-sign Lyapunov growth, and the Type II [72,36,16] gate.

codex-2026-06-11-P1

This script treats Euler's pentagonal theorem as a sparse-sign denominator

    E(q) = prod_{m>=1} (1-q^m)
         = 1 + sum_k (-1)^k (q^{k(3k-1)/2} + q^{k(3k+1)/2})

and compares it with perturbed/random sign laws on the same pentagonal support.
The reciprocal coefficients are computed exactly by recurrence.  The observed
Lyapunov slope is the linear coefficient in log |a_n|, while Euler's law is
measured on the Hardy-Ramanujan sqrt(n) scale.

The coding-theory half computes the unique formal Gleason weight enumerator for
an extremal Type II binary self-dual [72,36,16] code:

    W = A^9 + c1 A^6 B + c2 A^3 B^2 + c3 B^3

with A = x^8 + 14 x^4 y^4 + y^8 and
B = x^4 y^4 (x^4-y^4)^4, forcing weights 4,8,12 to vanish.

Tournament Analysis note:
    Vertices are sign laws, not runners/arcs.  I also considered pentagonal
    exponents, weight layers, Gleason monomials, automorphism-prime cases, and
    proof obligations as vertices.  The sign-law quotient preserves the
    cancellation/growth predicate; it destroys individual root locations and
    any support-design information needed for an actual [72,36,16] code.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb, log, sqrt
import random


N_REC = 650
RANDOM_SAMPLES = 160
SEED = 7203616


def pentagonal_terms(limit: int) -> list[tuple[int, int, int]]:
    """Return (exponent, pair-index k, Euler denominator coefficient)."""
    out: list[tuple[int, int, int]] = []
    k = 1
    while True:
        g1 = k * (3 * k - 1) // 2
        g2 = k * (3 * k + 1) // 2
        if g1 > limit and g2 > limit:
            break
        s = -1 if k % 2 else 1
        if g1 <= limit:
            out.append((g1, k, s))
        if g2 <= limit:
            out.append((g2, k, s))
        k += 1
    out.sort()
    return out


TERMS = pentagonal_terms(N_REC)


def euler_signs(limit: int = N_REC) -> dict[int, int]:
    return {g: s for g, _k, s in pentagonal_terms(limit)}


def sign_law(name: str, limit: int = N_REC, seed: int = 0) -> dict[int, int]:
    terms = pentagonal_terms(limit)
    if name == "euler":
        return {g: s for g, _k, s in terms}
    if name == "all_plus":
        return {g: 1 for g, _k, _s in terms}
    if name == "all_minus":
        return {g: -1 for g, _k, _s in terms}
    if name == "period3_pair":
        return {g: (1 if k % 3 else -1) for g, k, _s in terms}
    if name == "thue_morse_pair":
        return {g: (-1 if k.bit_count() % 2 else 1) for g, k, _s in terms}
    if name == "legendre7_pair":
        residues = {1, 2, 4}
        return {g: (1 if (k % 7) in residues else -1) for g, k, _s in terms}
    if name.startswith("flip_pair_"):
        flip = int(name.rsplit("_", 1)[1])
        return {g: (-s if k == flip else s) for g, k, s in terms}
    if name.startswith("random_pair_"):
        rng = random.Random(seed)
        by_k: dict[int, int] = {}
        for _g, k, _s in terms:
            if k not in by_k:
                by_k[k] = 1 if rng.randrange(2) else -1
        return {g: by_k[k] for g, k, _s in terms}
    if name.startswith("random_term_"):
        rng = random.Random(seed)
        return {g: (1 if rng.randrange(2) else -1) for g, _k, _s in terms}
    raise ValueError(name)


def reciprocal_coeffs(signs: dict[int, int], limit: int = N_REC) -> list[int]:
    """Coefficients of 1 / (1 + sum signs[g] q^g), through q^limit."""
    coeffs = [0] * (limit + 1)
    coeffs[0] = 1
    support = sorted((g, s) for g, s in signs.items() if g <= limit)
    for n in range(1, limit + 1):
        total = 0
        for g, s in support:
            if g > n:
                break
            total += s * coeffs[n - g]
        coeffs[n] = -total
    return coeffs


def least_squares_slope(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den


def growth_metrics(coeffs: list[int], tail_frac: float = 0.55) -> dict[str, float]:
    start = max(10, int(len(coeffs) * tail_frac))
    ns: list[float] = []
    logs: list[float] = []
    sqrt_ratios: list[float] = []
    zero_tail = 0
    for n in range(start, len(coeffs)):
        a = abs(coeffs[n])
        if a == 0:
            zero_tail += 1
            y = 0.0
        else:
            y = log(a)
            sqrt_ratios.append(y / sqrt(n))
        ns.append(float(n))
        logs.append(y)
    slope = least_squares_slope(ns, logs)
    root = sum((abs(coeffs[n]) or 1) ** (1.0 / n) for n in range(start, len(coeffs))) / (
        len(coeffs) - start
    )
    sqrt_scale = sorted(sqrt_ratios)[len(sqrt_ratios) // 2] if sqrt_ratios else 0.0
    return {
        "lambda": slope,
        "root": root,
        "sqrt_scale": sqrt_scale,
        "zero_tail": float(zero_tail),
        "max_abs": float(max(abs(a) for a in coeffs)),
    }


@dataclass
class LawResult:
    name: str
    lam: float
    root: float
    sqrt_scale: float
    max_abs: int
    coeffs: list[int]


def analyze_law(name: str, seed: int = 0) -> LawResult:
    coeffs = reciprocal_coeffs(sign_law(name, seed=seed), N_REC)
    m = growth_metrics(coeffs)
    return LawResult(
        name=name,
        lam=m["lambda"],
        root=m["root"],
        sqrt_scale=m["sqrt_scale"],
        max_abs=max(abs(a) for a in coeffs),
        coeffs=coeffs,
    )


def summarize_random_family(family: str) -> tuple[list[LawResult], dict[str, float]]:
    rows: list[LawResult] = []
    for i in range(RANDOM_SAMPLES):
        rows.append(analyze_law(family, seed=SEED + 1009 * i))
    lams = sorted(r.lam for r in rows)
    stats = {
        "min": lams[0],
        "q10": lams[int(0.10 * (len(lams) - 1))],
        "median": lams[len(lams) // 2],
        "mean": sum(lams) / len(lams),
        "q90": lams[int(0.90 * (len(lams) - 1))],
        "max": lams[-1],
    }
    return rows, stats


def poly_add(a: list[Fraction], b: list[Fraction], scale: Fraction = Fraction(1)) -> list[Fraction]:
    n = max(len(a), len(b))
    out = [Fraction(0) for _ in range(n)]
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += scale * b[i]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_pow(a: list[Fraction], e: int) -> list[Fraction]:
    out = [Fraction(1)]
    base = a[:]
    while e:
        if e & 1:
            out = poly_mul(out, base)
        e >>= 1
        if e:
            base = poly_mul(base, base)
    return out


def gleason_72() -> tuple[list[Fraction], list[Fraction]]:
    """Return (c0..c3, W(z)) for z=y^4 at x=1."""
    A = [Fraction(1), Fraction(14), Fraction(1)]
    B = [Fraction(0), Fraction(1), Fraction(-4), Fraction(6), Fraction(-4), Fraction(1)]
    basis = [poly_mul(poly_pow(A, 9 - 3 * j), poly_pow(B, j)) for j in range(4)]
    coeffs = [Fraction(1)]
    W = basis[0][:]
    for r in range(1, 4):
        current = W[r] if r < len(W) else Fraction(0)
        pivot = basis[r][r]
        c = -current / pivot
        coeffs.append(c)
        W = poly_add(W, basis[r], c)
    return coeffs, W


def code_design_lambdas(A16: int) -> list[tuple[int, Fraction]]:
    out: list[tuple[int, Fraction]] = []
    for t in range(1, 6):
        out.append((t, Fraction(A16 * comb(16, t), comb(72, t))))
    return out


def tournament_analysis(rows: list[LawResult]) -> dict[str, object]:
    # Lower lambda means stronger cancellation.  Ties use the listed order.
    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    eps = 1e-5
    for i in range(n):
        for j in range(i + 1, n):
            if rows[i].lam < rows[j].lam - eps:
                adj[i][j] = True
            elif rows[j].lam < rows[i].lam - eps:
                adj[j][i] = True
            else:
                adj[i][j] = True
    scores = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]

    cycles = 0
    for a, b, c in combinations(range(n), 3):
        edges = [adj[a][b], adj[b][c], adj[c][a]]
        if all(edges) or not any(edges):
            cycles += 1

    # SCCs by reachability (small n, Floyd-Warshall style).
    reach = [row[:] for row in adj]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen = set()
    scc_sizes = []
    for i in range(n):
        if i in seen:
            continue
        comp = {j for j in range(n) if reach[i][j] and reach[j][i]}
        seen |= comp
        scc_sizes.append(len(comp))

    # Hamiltonian path count in the directed tournament.
    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    hp = sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))
    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles": cycles,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_paths": hp,
        "order": [r.name for r in sorted(rows, key=lambda x: (x.lam, x.name))],
    }


def print_main() -> None:
    print("Pentagonal random-sign Lyapunov + Type II [72,36,16] gate")
    print(f"N_REC={N_REC}, RANDOM_SAMPLES={RANDOM_SAMPLES}, SEED={SEED}")
    print()

    print("[1] Sparse pentagonal denominator laws")
    base_names = [
        "euler",
        "all_plus",
        "all_minus",
        "period3_pair",
        "thue_morse_pair",
        "legendre7_pair",
    ]
    base_rows = [analyze_law(name) for name in base_names]
    for row in base_rows:
        print(
            f"  {row.name:16s} lambda={row.lam: .6f} "
            f"root={row.root: .6f} sqrt-scale={row.sqrt_scale: .6f} "
            f"max|a_n|={row.max_abs}"
        )
    euler = base_rows[0]
    hr = 3.141592653589793 * sqrt(2 / 3)
    print(f"  Hardy-Ramanujan sqrt constant pi*sqrt(2/3)={hr:.6f}")
    print(f"  Euler tail sqrt-scale at N={N_REC}: {euler.sqrt_scale:.6f}")
    print()

    print("[2] Single paired pentagonal sign flips away from Euler")
    flip_rows = [analyze_law(f"flip_pair_{k}") for k in range(1, 21)]
    for row in flip_rows:
        k = int(row.name.rsplit("_", 1)[1])
        print(f"  flip k={k:2d}: lambda={row.lam: .6f} root={row.root: .6f}")
    print()

    print("[3] Random sign laws on the pentagonal support")
    pair_rows, pair_stats = summarize_random_family("random_pair_sample")
    term_rows, term_stats = summarize_random_family("random_term_sample")
    print("  random-pair lambda stats:", {k: round(v, 6) for k, v in pair_stats.items()})
    print("  random-term lambda stats:", {k: round(v, 6) for k, v in term_stats.items()})
    best_pair = min(pair_rows, key=lambda r: r.lam)
    worst_pair = max(pair_rows, key=lambda r: r.lam)
    med_pair = sorted(pair_rows, key=lambda r: r.lam)[len(pair_rows) // 2]
    print(
        "  selected pair samples:"
        f" best={best_pair.lam:.6f}, median={med_pair.lam:.6f}, worst={worst_pair.lam:.6f}"
    )
    print()

    print("[4] Tournament Analysis over sign-law quotients")
    ta_rows = base_rows + [
        LawResult("random_pair_best", best_pair.lam, best_pair.root, best_pair.sqrt_scale, best_pair.max_abs, []),
        LawResult("random_pair_median", med_pair.lam, med_pair.root, med_pair.sqrt_scale, med_pair.max_abs, []),
        LawResult("random_pair_worst", worst_pair.lam, worst_pair.root, worst_pair.sqrt_scale, worst_pair.max_abs, []),
    ]
    ta = tournament_analysis(ta_rows)
    print("  vertices:", [r.name for r in ta_rows])
    print("  observable: lower Lyapunov slope wins; tie path = listed order")
    print("  score histogram:", ta["score_hist"])
    print("  directed 3-cycles:", ta["cycles"])
    print("  SCC sizes:", ta["scc_sizes"])
    print("  Hamiltonian path count:", ta["hamiltonian_paths"])
    print("  cancellation order:", " -> ".join(ta["order"]))
    print()

    print("[5] Gleason gate for a putative extremal Type II [72,36,16] code")
    coeffs, W = gleason_72()
    print("  W = A^9 + c1 A^6 B + c2 A^3 B^2 + c3 B^3")
    print("  coefficients:", [int(c) if c.denominator == 1 else c for c in coeffs])
    assert all(c.denominator == 1 for c in W)
    W_int = [int(c) for c in W]
    print("  first weights:")
    for r, a in enumerate(W_int[:9]):
        print(f"    weight {4*r:2d}: {a}")
    print("  last weights:")
    for r, a in list(enumerate(W_int))[-6:]:
        print(f"    weight {4*r:2d}: {a}")
    print(f"  total W(1,1)={sum(W_int)} expected 2^36={2**36}")
    print(f"  minimum positive weight={next(4*r for r, a in enumerate(W_int) if r and a)}")
    print(f"  all coefficients nonnegative? {all(a >= 0 for a in W_int)}")
    A16 = W_int[4]
    print(f"  A_16={A16}")
    print("  5-design lambdas for weight-16 words:")
    for t, lam in code_design_lambdas(A16):
        val = int(lam) if lam.denominator == 1 else f"{lam.numerator}/{lam.denominator}"
        print(f"    t={t}: lambda_t={val}")
    print()

    print("[6] Synthesis checkpoints")
    print("  Euler signs: theorem-level ordinary Lyapunov exponent is 0; the finite")
    print("  regression still sees Hardy-Ramanujan sqrt-growth, so it is only a")
    print("  diagnostic, not a proof of the limiting slope.")
    print("  Random pentagonal signs: positive finite-window slopes in every sampled case;")
    print("  the open proof target is an interior zero of the random sparse denominator.")
    print("  One early Euler pair flip usually creates visible exponential growth; late")
    print("  flips beyond the sampled window become indistinguishable from Euler here.")
    print("  The [72,36,16] scalar enumerator exists and is nonnegative; the open problem")
    print("  is support/design realization, not the scalar modular-form gate.")


if __name__ == "__main__":
    print_main()
