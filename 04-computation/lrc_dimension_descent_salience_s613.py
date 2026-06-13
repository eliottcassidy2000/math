#!/usr/bin/env python3
"""
lrc_dimension_descent_salience_s613.py

codex-2026-06-03 S613

Dimension-descent salience audit for the LRC proof program.

Question:
  If the n=14 residue quotient is hard, can the proof be reduced in size from
  n=14 to n=13 to n=12 to n=11, etc.?  What actually shrinks?

Method:
  For each LRC floor parameter n, set C=2n-1 and K=n-1 runners.  Work in the
  THM-401 / THM-407 residue quotient modulo C.  Rows are K-subsets of the
  nonzero residues that hit every unit antipodal shell {a,C-a}; this is the
  canonical unit-shell quotient used by the n=14 work.

  Each row is passed through the same cheap necessary gates used in the n=14
  quotient tower:
    D gate: every divisor-sized clock q=2..n-1 is hit by some v == 0 mod q.
    N gate: every n-clock frequency j<=floor(n/2) has a runner within one
            residue of the origin.
    Pinch gate: some pair-sum pinch t=m/(a+b) gives min_i ||v_i t|| >= 1/n.

  The audit is not a new proof of all n.  It is a burden comparison: which
  proof obligations remain after the same quotient pipeline is applied as n
  descends?

Tournament Analysis:
  Vertices are n-level quotient summaries, not runners.  The pairwise
  observable is the proof-burden tuple

    (below, survivors, unit_rows, primitive_floor, floor,
     nonunit_shells, shell_orbits, n).

  The switch orients lower burden before higher burden, with n as the
  Hamiltonian tie path.  A second "dimension descent" gauge orients smaller n
  before larger n; edge flips between the two gauges mark places where naive
  size descent is not the real proof descent.

Assumption challenge:
  Candidate tournament vertices considered: runners, gaps, fixed circle
  sections, section boundaries, wall-crossing events, residues, cover arcs,
  Fourier modes, matroid circuits, proof obligations, shell orbits, and
  n-level summaries.  This script chooses n-level summaries because the user
  question is explicitly about reducing the proof size from 14 to 13 to 12 to
  11.  The quotient preserves the LRC floor predicate seen by the D/N/pinch
  gates and destroys the actual integer carry fiber, so S611 remains the
  separate lift theorem.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import reduce
from itertools import combinations
from math import comb, gcd


MIN_N = 3
MAX_N = 14


@dataclass(frozen=True)
class PinchRoute:
    route: str
    witness_m: int | None
    witness_q: int | None


@dataclass(frozen=True)
class NReport:
    n: int
    C: int
    factors: tuple[tuple[int, int], ...]
    raw_subsets: int
    shell_orbits: tuple[tuple[int, ...], ...]
    orbit_gcds: tuple[int, ...]
    unit_shells: int
    nonunit_shells: int
    unit_rows: int
    d_pass: int
    dn_pass: int
    survivors: int
    strict: int
    floor: int
    primitive_floor: int
    below: int
    floor_examples: tuple[str, ...]
    survivor_type_hist: tuple[tuple[str, int], ...]

    @property
    def max_prime_power_depth(self) -> int:
        return max(exp for _, exp in self.factors)

    @property
    def burden(self) -> tuple[int, int, int, int, int, int, int, int]:
        return (
            self.below,
            self.survivors,
            self.unit_rows,
            self.primitive_floor,
            self.floor,
            self.nonunit_shells,
            len(self.shell_orbits),
            self.n,
        )


def factorize(x: int) -> tuple[tuple[int, int], ...]:
    out: list[tuple[int, int]] = []
    d = 2
    while d * d <= x:
        if x % d == 0:
            e = 0
            while x % d == 0:
                x //= d
                e += 1
            out.append((d, e))
        d += 1 if d == 2 else 2
    if x > 1:
        out.append((x, 1))
    return tuple(out)


def factor_str(factors: tuple[tuple[int, int], ...]) -> str:
    parts = []
    for p, e in factors:
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(parts)


def shell_label(r: int, C: int) -> int:
    r %= C
    if r == 0:
        return 0
    return min(r, C - r)


def shell_orbits(n: int) -> tuple[tuple[int, ...], ...]:
    C = 2 * n - 1
    seen: set[int] = set()
    orbits: list[tuple[int, ...]] = []
    for start in range(1, n):
        if start in seen:
            continue
        orbit: set[int] = set()
        frontier = {start}
        while frontier:
            a = frontier.pop()
            if a in orbit:
                continue
            orbit.add(a)
            b = shell_label(2 * a, C)
            if b and b not in orbit:
                frontier.add(b)
        seen.update(orbit)
        orbits.append(tuple(sorted(orbit)))
    return tuple(sorted(orbits, key=lambda xs: (gcd(xs[0], C), xs[0])))


def iter_unit_shell_rows(n: int):
    C = 2 * n - 1
    K = n - 1
    shells = tuple((a, C - a, gcd(a, C) == 1) for a in range(1, n))
    suffix_min_unit = [0] * (len(shells) + 1)
    for i in range(len(shells) - 1, -1, -1):
        suffix_min_unit[i] = suffix_min_unit[i + 1] + (1 if shells[i][2] else 0)

    def rec(i: int, chosen: tuple[int, ...]):
        chosen_len = len(chosen)
        if chosen_len > K:
            return
        remaining = len(shells) - i
        if chosen_len + suffix_min_unit[i] > K:
            return
        if chosen_len + 2 * remaining < K:
            return
        if i == len(shells):
            if chosen_len == K:
                yield tuple(sorted(chosen))
            return

        a, b, is_unit = shells[i]
        options = ((a,), (b,), (a, b)) if is_unit else ((), (a,), (b,), (a, b))
        for opt in options:
            yield from rec(i + 1, chosen + opt)

    yield from rec(0, ())


def d_failures(row: tuple[int, ...], n: int) -> tuple[int, ...]:
    return tuple(q for q in range(2, n) if all(v % q != 0 for v in row))


def n_clock_failures(row: tuple[int, ...], n: int) -> tuple[int, ...]:
    failures: list[int] = []
    for j in range(1, n // 2 + 1):
        hit = False
        for v in row:
            r = (v * j) % n
            if min(r, n - r) <= 1:
                hit = True
                break
        if not hit:
            failures.append(j)
    return tuple(failures)


def pinch_route(row: tuple[int, ...], n: int) -> PinchRoute:
    denominators = sorted({a + b for a, b in combinations(row, 2)})
    floor_witness: tuple[int, int] | None = None

    for q in denominators:
        for m in range(1, q // 2 + 1):
            local = q
            for v in row:
                r = (v * m) % q
                d = min(r, q - r)
                if d < local:
                    local = d
                    if n * local < q:
                        break
            cmp = n * local - q
            if cmp > 0:
                return PinchRoute("strict", m, q)
            if cmp == 0 and floor_witness is None:
                floor_witness = (m, q)

    if floor_witness is not None:
        return PinchRoute("floor", floor_witness[0], floor_witness[1])
    return PinchRoute("below", None, None)


def row_gcd(row: tuple[int, ...]) -> int:
    return reduce(gcd, row)


def row_name(row: tuple[int, ...], n: int) -> str:
    C = 2 * n - 1
    ap = tuple(range(1, n))
    if row == ap:
        return "AP"
    two_ap = tuple(2 * a for a in range(1, n))
    if row == two_ap:
        return "2*AP"
    for a in range(1, n):
        b = 2 * a
        if b > C - 1 or b in ap:
            continue
        candidate = tuple(sorted(x for x in ap if x != a) + [b])
        if row == candidate:
            return f"AP[{a}->{b}]"
    return str(row)


def survivor_type(row: tuple[int, ...], n: int) -> str:
    C = 2 * n - 1
    doubled = 0
    missed_nonunit = 0
    full_nonunit = 0
    for a in range(1, n):
        pair = {a, C - a}
        hits = sum(1 for x in pair if x in row)
        if hits == 2:
            doubled += 1
        if gcd(a, C) > 1:
            if hits == 0:
                missed_nonunit += 1
            elif hits == 2:
                full_nonunit += 1
    primitive = "prim" if row_gcd(row) == 1 else "imprim"
    return (
        f"{primitive};double={doubled};"
        f"miss_nonunit={missed_nonunit};full_nonunit={full_nonunit}"
    )


def audit_n(n: int) -> NReport:
    C = 2 * n - 1
    factors = factorize(C)
    orbits = shell_orbits(n)
    orbit_gcds = tuple(gcd(orb[0], C) for orb in orbits)
    unit_shells = sum(1 for a in range(1, n) if gcd(a, C) == 1)
    nonunit_shells = (n - 1) - unit_shells

    unit_rows = 0
    d_pass = 0
    dn_pass = 0
    survivors = 0
    strict = 0
    floor = 0
    primitive_floor = 0
    below = 0
    floor_examples: list[str] = []
    type_hist: Counter[str] = Counter()

    for row in iter_unit_shell_rows(n):
        unit_rows += 1
        if d_failures(row, n):
            continue
        d_pass += 1
        if n_clock_failures(row, n):
            continue
        dn_pass += 1
        survivors += 1
        type_hist[survivor_type(row, n)] += 1
        route = pinch_route(row, n)
        if route.route == "strict":
            strict += 1
        elif route.route == "floor":
            floor += 1
            if row_gcd(row) == 1:
                primitive_floor += 1
            if len(floor_examples) < 12:
                floor_examples.append(
                    f"{row_name(row, n)} via {route.witness_m}/{route.witness_q}"
                )
        else:
            below += 1

    return NReport(
        n=n,
        C=C,
        factors=factors,
        raw_subsets=comb(C - 1, n - 1),
        shell_orbits=orbits,
        orbit_gcds=orbit_gcds,
        unit_shells=unit_shells,
        nonunit_shells=nonunit_shells,
        unit_rows=unit_rows,
        d_pass=d_pass,
        dn_pass=dn_pass,
        survivors=survivors,
        strict=strict,
        floor=floor,
        primitive_floor=primitive_floor,
        below=below,
        floor_examples=tuple(floor_examples),
        survivor_type_hist=tuple(type_hist.most_common(6)),
    )


def score_histogram(reports: tuple[NReport, ...]) -> Counter[int]:
    scores: Counter[int] = Counter()
    for a in reports:
        out = 0
        for b in reports:
            if a is b:
                continue
            if a.burden < b.burden:
                out += 1
        scores[out] += 1
    return scores


def directed_3_cycles(reports: tuple[NReport, ...]) -> int:
    cycles = 0
    for a, b, c in combinations(reports, 3):
        ab = a.burden < b.burden
        bc = b.burden < c.burden
        ca = c.burden < a.burden
        ba = not ab
        cb = not bc
        ac = not ca
        if (ab and bc and ca) or (ba and cb and ac):
            cycles += 1
    return cycles


def edge_flips_vs_dimension(reports: tuple[NReport, ...]) -> int:
    flips = 0
    total = 0
    for a, b in combinations(reports, 2):
        total += 1
        burden_says_a_easier = a.burden < b.burden
        dimension_says_a_easier = a.n < b.n
        if burden_says_a_easier != dimension_says_a_easier:
            flips += 1
    return flips, total


def adjacent_note(high: NReport, low: NReport) -> str:
    notes: list[str] = []
    if high.max_prime_power_depth > low.max_prime_power_depth:
        notes.append(
            f"prime-power depth drops {high.max_prime_power_depth}->{low.max_prime_power_depth}"
        )
    elif high.max_prime_power_depth < low.max_prime_power_depth:
        notes.append(
            f"prime-power depth rises {high.max_prime_power_depth}->{low.max_prime_power_depth}"
        )
    if high.nonunit_shells and not low.nonunit_shells:
        notes.append("nonunit shells vanish")
    elif not high.nonunit_shells and low.nonunit_shells:
        notes.append("nonunit shells reappear")
    if len(high.shell_orbits) != len(low.shell_orbits):
        notes.append(
            f"shell orbits {len(high.shell_orbits)}->{len(low.shell_orbits)}"
        )
    if high.survivors and low.survivors:
        if high.survivors > low.survivors:
            notes.append(f"survivors shrink {high.survivors}->{low.survivors}")
        elif high.survivors < low.survivors:
            notes.append(f"survivors grow {high.survivors}->{low.survivors}")
    if high.floor != low.floor:
        notes.append(f"floor atoms {high.floor}->{low.floor}")
    return "; ".join(notes) if notes else "mostly cardinality descent"


def print_report(reports: tuple[NReport, ...]) -> None:
    print("==== LRC dimension-descent salience audit (S613) ====")
    print("Rows are unit-shell quotient rows modulo C=2n-1.")
    print("Pinch route is after D/N gates; floor is M=1/n in the quotient.")
    print()
    print(
        " n   C factors raw_subsets shell_orb orbit_gcds unit/nonunit "
        "unit_rows D_pass DN_pass survivors strict floor prim_floor below"
    )
    print("-" * 128)
    for r in sorted(reports, key=lambda x: -x.n):
        print(
            f"{r.n:2d} {r.C:3d} {factor_str(r.factors):>7} "
            f"{r.raw_subsets:11d} {len(r.shell_orbits):9d} "
            f"{str(r.orbit_gcds):>14} {r.unit_shells:4d}/{r.nonunit_shells:<7d} "
            f"{r.unit_rows:9d} {r.d_pass:6d} {r.dn_pass:7d} "
            f"{r.survivors:9d} {r.strict:6d} {r.floor:5d} "
            f"{r.primitive_floor:10d} {r.below:5d}"
        )

    print()
    print("Floor atoms seen by the quotient:")
    for r in sorted(reports, key=lambda x: -x.n):
        examples = ", ".join(r.floor_examples) if r.floor_examples else "none"
        print(f"  n={r.n:2d}: {examples}")

    print()
    print("Shell orbits under <2,-1> on antipodal shells:")
    for r in sorted(reports, key=lambda x: -x.n):
        orbit_text = "; ".join(
            f"gcd {gcd(orb[0], r.C)}:{orb}" for orb in r.shell_orbits
        )
        print(f"  n={r.n:2d}, C={r.C:2d}: {orbit_text}")

    print()
    print("Adjacent descent ledger:")
    ordered = sorted(reports, key=lambda x: -x.n)
    for high, low in zip(ordered, ordered[1:]):
        print(f"  {high.n:2d}->{low.n:2d}: {adjacent_note(high, low)}")

    print()
    print("Top survivor obligation types:")
    for r in ordered:
        hist = ", ".join(f"{k}:{v}" for k, v in r.survivor_type_hist)
        print(f"  n={r.n:2d}: {hist}")

    print()
    print("Tournament Analysis over n-level quotient summaries:")
    hist = score_histogram(reports)
    flips, total = edge_flips_vs_dimension(reports)
    burden_path = sorted(reports, key=lambda x: x.burden)
    print(f"  vertices: {len(reports)}")
    print(f"  score_histogram: {dict(sorted(hist.items()))}")
    print(f"  directed_3_cycles: {directed_3_cycles(reports)}")
    print("  SCCs: transitive burden order, singleton components")
    print(f"  edge_flips_vs_naive_dimension: {flips}/{total}")
    print(
        "  burden Hamiltonian path easiest->hardest: "
        + " -> ".join(f"n={r.n}" for r in burden_path)
    )

    print()
    print("Salient conclusions:")
    print(
        "  * Descent is controlled by C=2n-1 arithmetic, not by n alone. "
        "The quotient gets easier when nonunit shells and prime-power depth "
        "drop, and harder when they reappear."
    )
    print(
        "  * n=12 is the clean reset in the 14->13->12->11 descent: "
        "C=23 is prime, all shells are unit shells, and the unit quotient is "
        "exactly one choice from each antipodal shell."
    )
    print(
        "  * n=14 is arithmetically special because C=27=3^3 has a "
        "three-level nonunit tower gcd 1/3/9; this is the shell-depth behind "
        "the owner and carry-fiber work."
    )
    print(
        "  * n=11 is smaller but not monotone-easier: C=21 is composite, so "
        "nonunit obligations reappear after the prime C=23 reset."
    )
    print(
        "  * The actual lift theorem cannot descend by projecting "
        "27->25->23: the modulus changes.  What descends is the obligation "
        "ledger, while S611's carry fiber remains a separate n=14 issue."
    )


def main() -> None:
    reports = tuple(audit_n(n) for n in range(MIN_N, MAX_N + 1))
    print_report(reports)


if __name__ == "__main__":
    main()
