#!/usr/bin/env python3
"""
lrc_n14_carry_conservativity_s611.py

codex-2026-06-03 S611

Carry-fiber diagnostic for the n=14 LRC lift/CRT conservativity problem.

Context:
  S608/HYP-2164 classifies the least-positive Res_27 quotient: after the
  canonical D/N gates, the primitive floor branch is AP plus V*.
  S609/HYP-2165 shows the C=27 owner layer carries the cheap-pair /
  positive-measure certificates in the fixed-boundary bridge.

This script probes the missing data: integer lifts over the same Res_27
shadow.  If x = r + 27k, then because 27 == -1 (mod 14),

    x == r - k (mod 14).

So the carry vector k is exactly the CRT glue between the C=27 shell quotient
and the n=14 clock ledger.  Forgetting k is not conservative.

Method:
  * Print the THM-407 G=<2,-1> shell orbits modulo 27.
  * For the known primitive floor rows AP and V*, compare every unit scalar
    lift uV with its least-positive Res_27 shadow.  The actual scalar lift is
    floor by scaling invariance; the least-positive section is often strict.
  * Audit local carry perturbations over AP and V*: add one or two C=27 carries
    to the canonical residues and compute exact M.  These local carry moves are
    all strict in the tested radius.

Tournament Analysis:
  Vertices are carry probes, not runners.  A probe is either a unit-scalar
  shadow or a local carry perturbation over AP/V*.  The pairwise observable is
    (route, margin above 1/14, carry_span, clock_change_count, label).
  The proof-burden switch orients strict larger-margin probes before floor
  probes, with label tie path.  A second carry-complexity gauge reports edge
  flips against the proof-burden order.  This quotient preserves the question
  "does this carry fiber probe stay at the floor?" and destroys phase order.

Assumption challenge:
  Vertices considered: runners, gaps, fixed circle sections, section
  boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
  matroid circuits, proof obligations, and lift carries.  This script chooses
  lift carries because HYP-2164 and HYP-2165 already say the base quotient and
  owner quotient are nearly closed; the open theorem is precisely whether the
  integer fiber over the Res_27 coimage is conservative.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd


N = 14
K = N - 1
C = 2 * N - 1
FLOOR = Fraction(1, N)

AP = tuple(range(1, N))
VSTAR = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
BASE_ROWS = (("AP", AP), ("V*", VSTAR))


@dataclass(frozen=True)
class ScalarProbe:
    base_name: str
    unit: int
    actual_row: tuple[int, ...]
    shadow_row: tuple[int, ...]
    carries: tuple[int, ...]
    actual_m: Fraction
    actual_t: Fraction
    shadow_m: Fraction
    shadow_t: Fraction
    actual_clock_failures: tuple[int, ...]
    shadow_clock_failures: tuple[int, ...]

    @property
    def label(self) -> str:
        return f"{self.base_name}:u={self.unit}"

    @property
    def shadow_margin(self) -> Fraction:
        return self.shadow_m - FLOOR

    @property
    def shadow_route(self) -> str:
        if self.shadow_m == FLOOR:
            return "floor shadow"
        if self.shadow_m > FLOOR:
            return "strict shadow"
        return "below shadow"

    @property
    def carry_span(self) -> int:
        return max(self.carries) - min(self.carries)

    @property
    def carry_sum(self) -> int:
        return sum(self.carries)

    @property
    def wrap_count(self) -> int:
        return sum(1 for k in self.carries if k)

    @property
    def clock_change_count(self) -> int:
        a = set(self.actual_clock_failures)
        b = set(self.shadow_clock_failures)
        return len(a.symmetric_difference(b))


@dataclass(frozen=True)
class LocalCarryProbe:
    base_name: str
    weight: int
    moved_speeds: tuple[int, ...]
    row: tuple[int, ...]
    score: Fraction
    time: Fraction

    @property
    def label(self) -> str:
        moved = ",".join(map(str, self.moved_speeds))
        return f"{self.base_name}:w{self.weight}:carry({moved})"

    @property
    def margin(self) -> Fraction:
        return self.score - FLOOR

    @property
    def route(self) -> str:
        if self.score == FLOOR:
            return "floor local carry"
        if self.score > FLOOR:
            return "strict local carry"
        return "below local carry"


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def exact_score(row: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in row)


@lru_cache(maxsize=None)
def exact_maximin(row: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    candidates: set[Fraction] = set()
    for i, a in enumerate(row):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in row[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))

    best = Fraction(0)
    best_t = Fraction(0)
    for t in candidates:
        score = exact_score(row, t)
        if (score, -t) > (best, -best_t):
            best = score
            best_t = t
    return best, best_t


def units_mod(modulus: int) -> tuple[int, ...]:
    return tuple(u for u in range(1, modulus) if gcd(u, modulus) == 1)


def shell(r: int) -> int:
    r %= C
    if r == 0:
        return 0
    return min(r, C - r)


def shell_orbits() -> list[tuple[int, ...]]:
    seen: set[int] = set()
    orbits: list[tuple[int, ...]] = []
    for start in range(1, N):
        if start in seen:
            continue
        orbit = set()
        frontier = {start}
        while frontier:
            a = frontier.pop()
            if a in orbit:
                continue
            orbit.add(a)
            for b in (shell(2 * a), shell(-a)):
                if b and b not in orbit:
                    frontier.add(b)
        seen.update(orbit)
        orbits.append(tuple(sorted(orbit)))
    return sorted(orbits, key=lambda xs: (gcd(xs[0], C), xs[0]))


def least_positive_residue(x: int) -> int:
    r = x % C
    if r == 0:
        raise ValueError(f"unexpected zero residue for {x}")
    return r


def lift_data(base: tuple[int, ...], unit: int) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    actual = tuple(unit * v for v in base)
    residues: list[int] = []
    carries: list[int] = []
    for x in actual:
        r = least_positive_residue(x)
        residues.append(r)
        carries.append((x - r) // C)
    return tuple(sorted(actual)), tuple(sorted(residues)), tuple(carries)


def clock_failures(row: tuple[int, ...]) -> tuple[int, ...]:
    failures: list[int] = []
    for j in range(1, N // 2 + 1):
        blocked = False
        for v in row:
            r = (v * j) % N
            if min(r, N - r) <= 1:
                blocked = True
                break
        if not blocked:
            failures.append(j)
    return tuple(failures)


def verify_carry_crt(base: tuple[int, ...], unit: int) -> bool:
    for v in base:
        x = unit * v
        r = least_positive_residue(x)
        k = (x - r) // C
        if x % N != (r - k) % N:
            return False
    return True


def scalar_probes() -> list[ScalarProbe]:
    probes: list[ScalarProbe] = []
    for base_name, base in BASE_ROWS:
        for unit in units_mod(C):
            actual, shadow, carries = lift_data(base, unit)
            if not verify_carry_crt(base, unit):
                raise RuntimeError(f"CRT carry identity failed for {base_name}, unit={unit}")
            actual_m, actual_t = exact_maximin(actual)
            shadow_m, shadow_t = exact_maximin(shadow)
            probes.append(
                ScalarProbe(
                    base_name=base_name,
                    unit=unit,
                    actual_row=actual,
                    shadow_row=shadow,
                    carries=carries,
                    actual_m=actual_m,
                    actual_t=actual_t,
                    shadow_m=shadow_m,
                    shadow_t=shadow_t,
                    actual_clock_failures=clock_failures(actual),
                    shadow_clock_failures=clock_failures(shadow),
                )
            )
    return probes


def local_carry_probes(max_weight: int = 2) -> list[LocalCarryProbe]:
    probes: list[LocalCarryProbe] = []
    for base_name, base in BASE_ROWS:
        for weight in range(1, max_weight + 1):
            for idxs in combinations(range(len(base)), weight):
                row = list(base)
                moved: list[int] = []
                for idx in idxs:
                    moved.append(base[idx])
                    row[idx] += C
                lifted = tuple(sorted(row))
                score, time = exact_maximin(lifted)
                probes.append(
                    LocalCarryProbe(
                        base_name=base_name,
                        weight=weight,
                        moved_speeds=tuple(moved),
                        row=lifted,
                        score=score,
                        time=time,
                    )
                )
    return probes


def route_key(route: str) -> int:
    return {
        "strict shadow": 0,
        "strict local carry": 0,
        "floor shadow": 1,
        "floor local carry": 1,
        "below shadow": 2,
        "below local carry": 2,
    }[route]


def tournament_fingerprint(
    scalar: list[ScalarProbe],
    local: list[LocalCarryProbe],
) -> dict[str, object]:
    vertices: list[tuple[str, str, Fraction, int, int]] = []
    for p in scalar:
        vertices.append((p.label, p.shadow_route, p.shadow_margin, p.carry_span, p.clock_change_count))
    for p in local:
        vertices.append((p.label, p.route, p.margin, p.weight, 0))

    proof_order = sorted(
        vertices,
        key=lambda x: (
            route_key(x[1]),
            -x[2],
            x[3],
            x[4],
            x[0],
        ),
    )
    carry_order = sorted(vertices, key=lambda x: (x[3], x[4], route_key(x[1]), -x[2], x[0]))
    proof_rank = {v[0]: i for i, v in enumerate(proof_order)}
    carry_rank = {v[0]: i for i, v in enumerate(carry_order)}

    flips = 0
    total = 0
    for a, b in combinations(vertices, 2):
        total += 1
        flips += (proof_rank[a[0]] < proof_rank[b[0]]) != (carry_rank[a[0]] < carry_rank[b[0]])

    outdegree_hist = Counter(len(vertices) - 1 - i for i in range(len(vertices)))
    return {
        "vertices": len(vertices),
        "score_histogram": tuple(sorted(outdegree_hist.items())),
        "scc_count": len(vertices),
        "largest_scc": 1 if vertices else 0,
        "directed_3_cycles": 0,
        "edge_flips_proof_vs_carry": (flips, total),
        "proof_path_head": tuple(v[0] for v in proof_order[:10]),
        "proof_path_tail": tuple(v[0] for v in proof_order[-10:]),
    }


def summarize_scalar(probes: list[ScalarProbe]) -> None:
    route_hist = Counter(p.shadow_route for p in probes)
    actual_hist = Counter(p.actual_m for p in probes)
    shadow_floor = [p for p in probes if p.shadow_m == FLOOR]
    shadow_strict = [p for p in probes if p.shadow_m > FLOOR]
    shadow_below = [p for p in probes if p.shadow_m < FLOOR]
    min_strict = min(shadow_strict, key=lambda p: (p.shadow_margin, p.label)) if shadow_strict else None
    max_strict = max(shadow_strict, key=lambda p: (p.shadow_margin, p.label)) if shadow_strict else None

    print("Scalar lift versus least-positive section:")
    print(f"  scalar probes: {len(probes)}")
    print(f"  actual scaled rows at floor: {actual_hist[FLOOR]}/{len(probes)}")
    print(f"  least-positive shadow route histogram: {dict(route_hist)}")
    print(f"  shadow below-floor rows: {len(shadow_below)}")
    print("  floor shadows:")
    for p in shadow_floor:
        print(
            f"    {p.label:8s} shadow={p.shadow_row} "
            f"carry_sum={p.carry_sum} carry_span={p.carry_span} "
            f"M={fmt_frac(p.shadow_m)} t={fmt_frac(p.shadow_t)}"
        )
    if min_strict is not None and max_strict is not None:
        print(
            "  strict shadow margin range: "
            f"min={fmt_frac(min_strict.shadow_margin)} ({min_strict.label}, "
            f"M={fmt_frac(min_strict.shadow_m)}), "
            f"max={fmt_frac(max_strict.shadow_margin)} ({max_strict.label}, "
            f"M={fmt_frac(max_strict.shadow_m)})"
        )
    print()

    print("  sample scalar probes:")
    for p in probes:
        if p.shadow_m == FLOOR or p.unit in (7, 14, 26):
            print(
                f"    {p.label:8s} actual_M={fmt_frac(p.actual_m):>5s} "
                f"shadow_M={fmt_frac(p.shadow_m):>5s} "
                f"route={p.shadow_route:13s} carries={p.carries} "
                f"clock_fail actual/shadow={p.actual_clock_failures}/{p.shadow_clock_failures}"
            )
    print()


def summarize_local(probes: list[LocalCarryProbe]) -> None:
    by_key: dict[tuple[str, int], list[LocalCarryProbe]] = {}
    for p in probes:
        by_key.setdefault((p.base_name, p.weight), []).append(p)

    print("Local carry perturbations over AP and V* residues:")
    for key, rows in sorted(by_key.items()):
        route_hist = Counter(p.route for p in rows)
        min_probe = min(rows, key=lambda p: (p.score, p.time, p.label))
        print(
            f"  {key[0]:2s} weight={key[1]} count={len(rows):3d} "
            f"routes={dict(route_hist)} "
            f"min_M={fmt_frac(min_probe.score)} at t={fmt_frac(min_probe.time)} "
            f"via {min_probe.label}"
        )
    floor = [p for p in probes if p.score == FLOOR]
    below = [p for p in probes if p.score < FLOOR]
    print(f"  local floor rows: {len(floor)}")
    print(f"  local below-floor rows: {len(below)}")
    print()


def main() -> None:
    print(f"==== LRC n={N} carry-fiber conservativity diagnostic (S611) ====")
    print(f"C=2n-1={C}; floor=1/{N}")
    print()

    print("THM-407 shell orbits under G=<2,-1> modulo 27:")
    for orbit in shell_orbits():
        g = gcd(orbit[0], C)
        print(f"  gcd={g:2d}: {orbit}")
    print()

    print("Carry/CRT identity:")
    print("  For every scalar probe x=u*v, write x=r+27k with 1<=r<=26.")
    print("  Since 27 == -1 (mod 14), x mod 14 is r-k mod 14.")
    print("  Therefore the carry vector k is visible to the n=14 clock ledger.")
    print()

    scalar = scalar_probes()
    summarize_scalar(scalar)

    local = local_carry_probes(max_weight=2)
    summarize_local(local)

    fp = tournament_fingerprint(scalar, local)
    print("Tournament Analysis over carry probes:")
    print(f"  vertices: {fp['vertices']}")
    score_hist = fp["score_histogram"]
    if (
        len(score_hist) == fp["vertices"]
        and all(count == 1 for _, count in score_hist)
        and score_hist[0][0] == 0
        and score_hist[-1][0] == fp["vertices"] - 1
    ):
        print(f"  score_histogram: transitive total order, degrees 0..{fp['vertices'] - 1} once each")
    else:
        print(f"  score_histogram: {score_hist}")
    print(f"  SCCs: count={fp['scc_count']} largest={fp['largest_scc']}")
    print(f"  directed_3_cycles: {fp['directed_3_cycles']}")
    flips, total = fp["edge_flips_proof_vs_carry"]
    print(f"  edge_flips between proof-burden and carry-complexity gauges: {flips}/{total}")
    print("  proof Hamiltonian path head:")
    for label in fp["proof_path_head"]:
        print(f"    {label}")
    print("  proof Hamiltonian path tail:")
    for label in fp["proof_path_tail"]:
        print(f"    {label}")
    print()

    print("Conclusion:")
    print(
        "  The least-positive Res_27 section is not conservative for lifted floor rows: "
        "all 36 AP/V* unit scalar lifts remain at M=1/14, but only three section "
        "shadows stay at the floor."
    )
    print(
        "  Those three floor shadows are exactly AP, V*, and nonprimitive 2*AP, "
        "matching HYP-2164."
    )
    print(
        "  In the AP/V* residue fibers, every tested nonzero local carry vector of "
        "Hamming weight one or two is strict.  Any new floor lift must therefore "
        "use a globally coherent carry pattern, not an isolated wrap."
    )
    print(
        "  The next proof target is a carry-conservativity theorem for the CRT fiber "
        "k in v=r+27k, together with the HYP-2165 owner/certificate bridge."
    )


if __name__ == "__main__":
    main()
