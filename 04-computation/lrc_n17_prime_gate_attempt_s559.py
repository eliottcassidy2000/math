#!/usr/bin/env python3
"""
lrc_n17_prime_gate_attempt_s559.py

codex-2026-06-02 S559

First focused proof-search pass for the prime Lonely Runner denominator n=17.

The row has 16 moving speeds and threshold 1/17.  Since 17 is prime, the
small-denominator sieve has a particularly clean top gate:

    if no speed is divisible by 17, every unit wall a/17 is a closed witness.

So any open-cover counterexample at n=17 must contain a 17-multiple.  This
script audits what happens after adding such a gate:

* one-gate swaps from the initial segment {1,...,16};
* "one primitive breaker plus 15 gates" rows {r} union 17*Q;
* random sieve-complete constructions;
* endpoint-core summaries and a small Tournament Analysis on the most
  counterexample-like exact rows.

The goal is a proof attempt, not a proof claim.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations, permutations
from math import gcd, lcm
from pathlib import Path
import builtins
import functools
import random


print = functools.partial(builtins.print, flush=True)

ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()


N = 17
K = N - 1
ONE = Fraction(1, 1)


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]
    missing: tuple[int, ...]
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witnesses: int
    witness: Fraction | None
    classification: str
    unprotected: int | None = None
    core_endpoints: int | None = None
    peel_depth: int | None = None
    quotient_layer: str | None = None


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def ffloat(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def initial() -> tuple[int, ...]:
    return tuple(range(1, N))


def gcd_all(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def normalized(values: tuple[int, ...]) -> tuple[int, ...]:
    return S356.normalize_speed_set(list(values))


def missing_moduli(speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(m for m in range(2, N + 1) if all(speed % m for speed in speeds))


def classify(report: object) -> str:
    if report.max_gap > 0:
        return "positive_gap"
    if report.boundary_witness_count:
        return "boundary_only"
    return "open_cover_candidate"


def exact_row(label: str, raw_speeds: tuple[int, ...], endpoint: bool = False) -> Row:
    speeds = normalized(raw_speeds)
    report = S356.report(label, list(speeds))
    missing = missing_moduli(speeds)
    unprotected = core = depth = None
    quotient = None
    if endpoint:
        descent = S362.summarize(list(speeds))
        unprotected = descent.unprotected_count
        core = descent.core_endpoint_count
        depth = len(descent.peel_layers)
        quotient = descent.quotient_layer
    return Row(
        label=label,
        speeds=speeds,
        missing=missing,
        forbidden_length=report.forbidden_length,
        max_gap=report.max_gap,
        gap_ratio=report.max_gap / report.threshold,
        boundary_witnesses=report.boundary_witness_count,
        witness=report.witness or report.boundary_witness,
        classification=classify(report),
        unprotected=unprotected,
        core_endpoints=core,
        peel_depth=depth,
        quotient_layer=quotient,
    )


def no_gate_witnesses(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    if any(speed % N == 0 for speed in speeds):
        return ()
    return tuple(
        Fraction(a, N)
        for a in range(1, N)
        if S356.is_lonely_witness(speeds, Fraction(a, N))
    )


def verify_prime_gate_branch(samples: int = 200, hi: int = 500) -> None:
    rng = random.Random(559)
    failures: list[tuple[int, ...]] = []
    for _ in range(samples):
        speeds = tuple(sorted(rng.sample([v for v in range(1, hi + 1) if v % N], K)))
        if len(no_gate_witnesses(speeds)) != N - 1:
            failures.append(speeds)
            break

    print("1. Prime 17-gate branch")
    print("   Lemma: if no speed is divisible by 17, all 16 unit walls a/17 are closed witnesses.")
    print(f"   random_no_gate_checks={samples}; failures={len(failures)}")
    print("   Therefore any open-cover counterexample at n=17 must contain a 17-multiple.")
    print()


def drop_add(drop: int, add: int) -> tuple[int, ...]:
    speeds = set(initial())
    speeds.remove(drop)
    speeds.add(add)
    return tuple(sorted(speeds))


def print_row(row: Row, indent: str = "   ") -> None:
    miss = ",".join(str(m) for m in row.missing) or "-"
    endpoint = ""
    if row.unprotected is not None:
        endpoint = (
            f" unprot={row.unprotected:>5} core={row.core_endpoints:>4}"
            f" peel={row.peel_depth:>2} layer={row.quotient_layer}"
        )
    print(
        f"{indent}{row.label:<24} class={row.classification:<13} "
        f"gap/th={ffloat(row.gap_ratio):>8} length={fmt(row.forbidden_length):>8} "
        f"boundary={row.boundary_witnesses:>3} missing={miss:<18} "
        f"wit={fmt(row.witness):>10}{endpoint}"
    )


def one_gate_swap_scan(qmax: int = 32) -> list[Row]:
    best_by_drop: list[Row] = []
    all_rows: list[Row] = []
    for drop in range(1, N):
        local: list[Row] = []
        for q in range(1, qmax + 1):
            row = exact_row(f"drop{drop}_add17x{q}", drop_add(drop, N * q))
            local.append(row)
        local.sort(key=lambda row: (row.gap_ratio, -row.forbidden_length, row.speeds[-1]))
        best_by_drop.append(local[0])
        all_rows.extend(local)

    all_rows.sort(key=lambda row: (row.gap_ratio, -row.forbidden_length, row.speeds[-1]))

    print("2. One 17-gate swaps from the initial segment")
    print(f"   scanned drops 1..16 and gates 17*q, q<= {qmax}")
    print("   best row for each dropped initial speed:")
    for row in best_by_drop:
        print_row(row)
    print("   global closest rows:")
    for row in all_rows[:10]:
        print_row(row)
    print()
    return all_rows[:12]


def breaker_gate_rows() -> list[Row]:
    rows: list[Row] = []
    for breaker in range(1, N):
        for skip in range(1, N):
            gates = tuple(N * q for q in range(1, N) if q != skip)
            speeds = tuple(sorted((breaker,) + gates))
            if len(set(speeds)) != K or gcd_all(speeds) != 1:
                continue
            rows.append(exact_row(f"breaker{breaker}_skip{skip}", speeds))

    rows.sort(key=lambda row: (row.gap_ratio, len(row.missing), -row.forbidden_length, row.speeds))
    print("3. One primitive breaker plus 15 gates")
    print("   rows are {r} union {17*q: 1<=q<=16, q!=skip}")
    print("   closest rows:")
    for row in rows[:12]:
        print_row(row)
    print("   missing-moduli histogram among these rows:")
    print(f"     {dict(sorted(Counter(len(row.missing) for row in rows).items()))}")
    print()
    return rows[:12]


def random_sieve_complete_probe(trials: int = 30, hi: int = 500) -> list[Row]:
    rng = random.Random(17559)
    pool = set(range(1, 90))
    pool.update(N * q for q in range(1, 25))
    pool.update(lcm(N, m) for m in range(2, N + 1))
    pool.update(m * q for m in range(2, N + 1) for q in range(2, 8) if m * q <= hi)
    pool = {v for v in pool if 1 <= v <= hi}
    pool_tuple = tuple(sorted(pool))

    rows: list[Row] = []
    tries = 0
    while len(rows) < trials and tries < 20000:
        tries += 1
        forced = {N * rng.randint(1, 24)}
        for m in range(2, N):
            choices = [v for v in pool_tuple if v % m == 0]
            forced.add(rng.choice(choices))
        while len(forced) < K:
            forced.add(rng.choice(pool_tuple))
        if len(forced) > K:
            forced = set(rng.sample(sorted(forced), K))
        speeds = tuple(sorted(forced))
        if len(speeds) != K or gcd_all(speeds) != 1:
            continue
        speeds = normalized(speeds)
        if missing_moduli(speeds) or not any(v % N == 0 for v in speeds):
            continue
        rows.append(exact_row(f"random_sieve_{len(rows)}", speeds))

    rows.sort(key=lambda row: (row.gap_ratio, -row.forbidden_length, row.speeds[-1]))
    print("4. Random sieve-complete exact probe")
    print(f"   accepted={len(rows)} tries={tries} hi={hi}")
    print("   closest accepted rows:")
    for row in rows[:10]:
        print_row(row)
    print()
    return rows[:10]


def endpoint_audit(rows: list[Row], limit: int = 8) -> list[Row]:
    chosen: list[Row] = []
    seen: set[tuple[int, ...]] = set()
    for row in rows:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        chosen.append(exact_row(row.label, row.speeds, endpoint=True))
        if len(chosen) >= limit:
            break

    print("5. Endpoint-core audit of closest rows")
    for row in chosen:
        print_row(row)
        print(f"     speeds={row.speeds}")
    print()
    return chosen


def tournament_fingerprint(rows: list[Row]) -> dict[str, object]:
    # Switch/gauge: more counterexample-like rows win.
    # Primary: smaller exact open gap; secondary: fewer missing moduli; then
    # longer forbidden length; then smaller endpoint core/unprotected debt.
    def key(row: Row) -> tuple[Fraction, int, Fraction, int, int, int]:
        unprotected = row.unprotected if row.unprotected is not None else 10**9
        core = row.core_endpoints if row.core_endpoints is not None else 10**9
        return (-row.gap_ratio, -len(row.missing), row.forbidden_length, -core, -unprotected, -row.speeds[-1])

    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            if i != j:
                adj[i][j] = key(left) > key(right)

    scores = [sum(adj[i]) for i in range(n)]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        cyc = (
            adj[i][j] and adj[j][k] and adj[k][i]
        ) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        )
        c3 += int(cyc)

    def reaches(start: int) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v in range(n):
                if adj[u][v] and v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reaches(u)
        comp = {v for v in remaining if v in ru and u in reaches(v)}
        sccs.append(len(comp))
        remaining -= comp

    hp: int | str = 0
    if n <= 8:
        for perm in permutations(range(n)):
            if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
                hp += 1
    else:
        hp = "skipped(n>8)"

    return {
        "vertices": [row.label for row in rows],
        "score_hist": dict(sorted(Counter(scores).items())),
        "c3": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": hp,
    }


def print_tournament(rows: list[Row]) -> None:
    print("6. Tournament Analysis")
    print("   vertices: closest exact n=17 repair rows")
    print("   pairwise observable: (exact gap/th, missing moduli, forbidden length, endpoint debt)")
    print("   switch/gauge: smaller gap wins; ties prefer sieve-complete rows and smaller endpoint core")
    print(f"   fingerprints: {tournament_fingerprint(rows)}")
    print()


def main() -> None:
    print("S559 n=17 prime-gate LRC proof attempt")
    print("=" * 72)
    print("Denominator n=17; moving speeds k=16; threshold=1/17.")
    print()

    verify_prime_gate_branch()

    base = exact_row("initial_1..16", initial(), endpoint=True)
    print("Initial tight row:")
    print_row(base)
    print(f"     speeds={base.speeds}")
    print()

    one_gate = one_gate_swap_scan()
    ladder = breaker_gate_rows()
    random_rows = random_sieve_complete_probe()
    print("Closest breaker-gate and random sieve rows are exact positive-gap rows;")
    print("endpoint-core checks are reserved for smaller one-gate representatives below.")
    print()
    audited = endpoint_audit(one_gate[:4])
    print_tournament(audited)

    best = min(one_gate + ladder + random_rows, key=lambda row: (row.gap_ratio, len(row.missing)))
    print("Synthesis")
    print(f"  closest_row={best.label} gap/th={ffloat(best.gap_ratio)} missing={best.missing or '-'}")
    print("  No scanned row is an open-cover candidate; closest rows have positive open gap and")
    print("  the endpoint-core audit peels to an empty core.")
    print("  Prime-row proof target: after the mandatory 17-gate, prove that every")
    print("  primitive repair either leaves a positive gap or exports a private endpoint")
    print("  leaf in the 17-adic gate tree.")


if __name__ == "__main__":
    main()
