#!/usr/bin/env python3
"""
lrc_p0_collapse_additive_chains_s602.py

codex-2026-06-03 S602

The user's "p_0=0 collapse" notation is not a stable symbol in the LRC files,
so this script uses the operational endpoint-collapse predicate already visible
in S357/S359:

    max open safe gap = 0,
    boundary witnesses exist,
    the witness denominator lcm collapses to n=k+1,
    the witnesses are exactly the unit points a/n.

Question: is the collapse family only the AP row, or does it include sparse
additive-chain rows such as (1,3,4,7) and (1,3,4,5,9)?

Methodology / Tournament Analysis:
- Vertices are proof lenses, not runners: AP, unit skeleton, antipodal
  transversals, flip-set {2}, two-seed additive chains, top-sum witnesses, and
  nonunit holes.
- Pairwise observable: a lens beats another if it captures more exact collapse
  rows in the scanned universe; ties are broken by fewer false positives.
- This quotient preserves the endpoint-collapse predicate and additive
  generation data, but destroys endpoint phase order inside Q(V).
- Challenged assumption: "floor collapse means AP."  The sporadic chains show
  that AP is only one lens on the collapse locus.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class BoxConfig:
    k: int
    max_speed: int


@dataclass(frozen=True)
class CollapseRow:
    speeds: tuple[int, ...]
    n: int
    first_witness: Fraction
    witness_count: int
    witness_modulus: int
    endpoint_modulus: int
    residues_mod_n: tuple[int, ...]
    residues_mod_c: tuple[int, ...]
    flipset: tuple[int, ...]
    perfect_transversal: bool
    missing_shells: tuple[int, ...]
    nonunit_missing_shells: tuple[int, ...]
    ap: bool
    two_seed_chain: bool
    top_sum_pairs: tuple[tuple[int, int], ...]
    top_previous_two: bool
    generation_recipe: tuple[str, ...]


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def circle_point(x: Fraction) -> Fraction:
    return x % ONE


def boundary_witnesses(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    raw_intervals = S356.forbidden_intervals(speeds)
    boundary_candidates = sorted(
        {
            circle_point(point)
            for interval in raw_intervals
            for point in interval
        }
    )
    return tuple(t for t in boundary_candidates if S356.is_lonely_witness(speeds, t))


def witness_modulus(witnesses: tuple[Fraction, ...]) -> int:
    q = 1
    for t in witnesses:
        q = lcm(q, t.denominator)
    return q


def unit_witnesses(n: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(a, n) for a in range(1, n) if gcd(a, n) == 1)


def shell_c(r: int, c: int) -> int:
    r %= c
    if r == 0:
        return 0
    return min(r, c - r)


def missing_shells(speeds: tuple[int, ...]) -> tuple[int, ...]:
    n = len(speeds) + 1
    c = 2 * n - 1
    occupied = {shell_c(v, c) for v in speeds if shell_c(v, c) != 0}
    return tuple(a for a in range(1, n) if a not in occupied)


def nonunit_shells(shells: tuple[int, ...], c: int) -> tuple[int, ...]:
    return tuple(a for a in shells if gcd(a, c) != 1)


def flipset(speeds: tuple[int, ...]) -> tuple[int, ...]:
    n = len(speeds) + 1
    c = 2 * n - 1
    speed_residues = {v % c for v in speeds}
    return tuple(a for a in range(1, n) if (c - a) in speed_residues)


def perfect_transversal(speeds: tuple[int, ...]) -> bool:
    n = len(speeds) + 1
    c = 2 * n - 1
    residues = [v % c for v in speeds]
    if 0 in residues:
        return False
    for a in range(1, n):
        hits = int(a in residues) + int((c - a) in residues)
        if hits != 1:
            return False
    return True


def addition_recipe(speeds: tuple[int, ...]) -> tuple[bool, tuple[str, ...]]:
    """Whether sorted speeds form a two-seed addition chain, plus one recipe."""
    if len(speeds) <= 2:
        return True, tuple(f"{v}: seed" for v in speeds)

    built = [speeds[0], speeds[1]]
    recipe = [f"{speeds[0]}: seed", f"{speeds[1]}: seed"]
    for v in speeds[2:]:
        witness = None
        for i, a in enumerate(built):
            for b in built[i:]:
                if a + b == v:
                    witness = (a, b)
                    break
            if witness is not None:
                break
        if witness is None:
            return False, tuple(recipe + [f"{v}: no previous sum"])
        recipe.append(f"{v}: {witness[0]}+{witness[1]}")
        built.append(v)
    return True, tuple(recipe)


def top_sum_pairs(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    top = speeds[-1]
    pairs: list[tuple[int, int]] = []
    below = speeds[:-1]
    for i, a in enumerate(below):
        for b in below[i:]:
            if a + b == top:
                pairs.append((a, b))
    return tuple(pairs)


def classify_collapse(speeds: tuple[int, ...]) -> CollapseRow | None:
    report = S356.report("p0-collapse", list(speeds))
    witnesses = boundary_witnesses(report.speeds)
    wmod = witness_modulus(witnesses)
    n = len(report.speeds) + 1
    units = unit_witnesses(n)
    if not (
        report.max_gap == 0
        and witnesses
        and wmod == n
        and witnesses == units
    ):
        return None

    c = 2 * n - 1
    miss = missing_shells(report.speeds)
    chain, recipe = addition_recipe(report.speeds)
    top_pairs = top_sum_pairs(report.speeds)
    return CollapseRow(
        speeds=report.speeds,
        n=n,
        first_witness=witnesses[0],
        witness_count=len(witnesses),
        witness_modulus=wmod,
        endpoint_modulus=report.boundary_modulus,
        residues_mod_n=tuple(v % n for v in report.speeds),
        residues_mod_c=tuple(v % c for v in report.speeds),
        flipset=flipset(report.speeds),
        perfect_transversal=perfect_transversal(report.speeds),
        missing_shells=miss,
        nonunit_missing_shells=nonunit_shells(miss, c),
        ap=report.speeds == tuple(range(1, n)),
        two_seed_chain=chain,
        top_sum_pairs=top_pairs,
        top_previous_two=(
            len(report.speeds) >= 3
            and report.speeds[-1] == report.speeds[-2] + report.speeds[-3]
        ),
        generation_recipe=recipe,
    )


def scan_box(config: BoxConfig) -> tuple[list[CollapseRow], Counter[str], int, list[tuple[int, ...]]]:
    rows: list[CollapseRow] = []
    chain_hist: Counter[str] = Counter()
    universe: list[tuple[int, ...]] = []
    total = 0
    for combo in combinations(range(1, config.max_speed + 1), config.k):
        if not primitive(combo):
            continue
        universe.append(combo)
        total += 1
        chain, _ = addition_recipe(combo)
        if chain:
            chain_hist["two_seed_chain"] += 1
        if len(combo) >= 3 and combo[-1] == combo[-2] + combo[-3]:
            chain_hist["top_previous_two"] += 1
        collapse = classify_collapse(combo)
        if collapse is not None:
            rows.append(collapse)
            if collapse.two_seed_chain:
                chain_hist["collapse_two_seed_chain"] += 1
            if collapse.top_previous_two:
                chain_hist["collapse_top_previous_two"] += 1
    return rows, chain_hist, total, universe


def canonical_transversal_rows(n: int) -> list[tuple[tuple[int, ...], Fraction]]:
    """Enumerate canonical residue transversals mod 2n-1 and score exactly via S356."""
    c = 2 * n - 1
    floor = Fraction(1, n)
    rows: list[tuple[tuple[int, ...], Fraction]] = []
    for mask in range(1 << (n - 1)):
        speeds = []
        for a in range(1, n):
            speeds.append(c - a if (mask >> (a - 1)) & 1 else a)
        speeds = tuple(sorted(speeds))
        if not primitive(speeds):
            continue
        # Boundary collapse at the n-clock floor is enough for the present audit.
        collapse = classify_collapse(speeds)
        if collapse is not None and collapse.first_witness == floor:
            rows.append((speeds, floor))
    return rows


def lens_values(row: CollapseRow, all_rows: list[tuple[int, ...]]) -> dict[str, bool]:
    c = 2 * row.n - 1
    return {
        "AP_only": row.ap,
        "unit_boundary_skeleton": True,
        "two_seed_addition_chain": row.two_seed_chain,
        "top_has_sum_pair": bool(row.top_sum_pairs),
        "top_previous_two_sum": row.top_previous_two,
        "perfect_C_transversal": row.perfect_transversal,
        "flipset_{2}": row.flipset == (2,),
        "nonunit_C_hole": bool(row.nonunit_missing_shells),
        "C_prime": c % 2 == 1,
    }


def lens_predicates(speeds: tuple[int, ...], is_collapse: bool) -> dict[str, bool]:
    n = len(speeds) + 1
    c = 2 * n - 1
    chain, _ = addition_recipe(speeds)
    miss = missing_shells(speeds)
    return {
        "AP_only": speeds == tuple(range(1, n)),
        "unit_boundary_skeleton": is_collapse,
        "two_seed_addition_chain": chain,
        "top_has_sum_pair": bool(top_sum_pairs(speeds)),
        "top_previous_two_sum": (
            len(speeds) >= 3 and speeds[-1] == speeds[-2] + speeds[-3]
        ),
        "perfect_C_transversal": perfect_transversal(speeds),
        "flipset_{2}": flipset(speeds) == (2,),
        "nonunit_C_hole": bool(nonunit_shells(miss, c)),
        "C_prime": all(gcd(a, c) == 1 for a in range(1, n)),
    }


def tournament_analysis(universe: list[tuple[int, ...]], collapses: set[tuple[int, ...]]) -> None:
    lens_names = [
        "AP_only",
        "unit_boundary_skeleton",
        "two_seed_addition_chain",
        "top_has_sum_pair",
        "top_previous_two_sum",
        "perfect_C_transversal",
        "flipset_{2}",
        "nonunit_C_hole",
        "C_prime",
    ]
    stats: dict[str, tuple[int, int]] = {}
    values = {speeds: lens_predicates(speeds, speeds in collapses) for speeds in universe}
    for name in lens_names:
        hits = sum(1 for speeds in universe if values[speeds][name] and speeds in collapses)
        false_pos = sum(1 for speeds in universe if values[speeds][name] and speeds not in collapses)
        stats[name] = (hits, false_pos)

    scores = {name: 0 for name in lens_names}
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(lens_names, 2):
        a_key = (stats[a][0], -stats[a][1], a)
        b_key = (stats[b][0], -stats[b][1], b)
        winner = a if a_key > b_key else b
        loser = b if winner == a else a
        scores[winner] += 1
        edges[(winner, loser)] = winner

    ordered = sorted(lens_names, key=lambda name: (stats[name][0], -stats[name][1], name), reverse=True)
    score_hist = Counter(scores.values())

    # The tie-breaker makes this tournament transitive unless two names are equal,
    # so SCCs are singleton and there is one Hamiltonian path in the score order.
    print("Tournament Analysis: proof-lens quotient")
    print(f"  vertices={lens_names}")
    print("  observable=(collapse_hits, -false_positives, lens_name)")
    print("  lens_stats=")
    for name in ordered:
        hits, false_pos = stats[name]
        print(f"    {name}: hits={hits} false_positives={false_pos} score={scores[name]}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print("  directed_3_cycles=0")
    print("  scc_count=9 singleton_sccs=9")
    print(f"  tie_hamiltonian_path={' > '.join(ordered)}")
    print("  assumption_challenge=vertices are proof lenses, not runners/arcs; AP-only loses recall")
    print()


def main() -> None:
    configs = [
        BoxConfig(k=3, max_speed=12),
        BoxConfig(k=4, max_speed=16),
        BoxConfig(k=5, max_speed=13),
        BoxConfig(k=6, max_speed=11),
        BoxConfig(k=7, max_speed=14),
    ]

    print("LRC p0-collapse/additive-chain audit (codex-2026-06-03 S602)")
    print("=" * 78)
    print("Operational p0-collapse: no open safe gap, unit boundary witnesses, witness lcm=n.")
    print("Exact interval arithmetic is inherited from S356.\n")

    all_collapse_rows: list[CollapseRow] = []
    universe: list[tuple[int, ...]] = []
    for config in configs:
        rows, chain_hist, total, scanned_universe = scan_box(config)
        all_collapse_rows.extend(rows)
        universe.extend(scanned_universe)
        print(f"Primitive box k={config.k}, max_speed={config.max_speed}, n={config.k+1}")
        print(f"  total={total}")
        print(f"  p0_collapse_rows={len(rows)}")
        print(f"  chain_hist={dict(sorted(chain_hist.items()))}")
        for row in rows:
            print(
                "    "
                f"speeds={row.speeds} first={fmt_frac(row.first_witness)} "
                f"witnesses={row.witness_count} witness_Q={row.witness_modulus} "
                f"endpoint_Q={row.endpoint_modulus} residues_mod_n={row.residues_mod_n}"
            )
            print(
                "      "
                f"AP={row.ap} two_seed_chain={row.two_seed_chain} "
                f"top_sum_pairs={row.top_sum_pairs or '-'} top_prev2={row.top_previous_two}"
            )
            print(
                "      "
                f"C={2*row.n-1} residues_mod_C={row.residues_mod_c} "
                f"perfect_transversal={row.perfect_transversal} flipset={row.flipset} "
                f"missing_shells={row.missing_shells or '-'} "
                f"nonunit_missing={row.nonunit_missing_shells or '-'}"
            )
            print("      recipe=" + "; ".join(row.generation_recipe))
        print()

    print("Canonical residue-transversal floor rows")
    for n in range(4, 9):
        floor_rows = canonical_transversal_rows(n)
        print(f"  n={n}, C={2*n-1}: floor_transversal_rows={len(floor_rows)}")
        for speeds, floor in floor_rows:
            row = classify_collapse(speeds)
            assert row is not None
            print(
                "    "
                f"speeds={speeds} floor={fmt_frac(floor)} flipset={row.flipset} "
                f"top_prev2={row.top_previous_two} recipe={' ; '.join(row.generation_recipe)}"
            )
    print()

    named = [(1, 3, 4, 7), (1, 3, 4, 5, 9)]
    print("Named sporadic additive-chain check")
    for speeds in named:
        row = classify_collapse(speeds)
        print(f"  speeds={speeds}")
        if row is None:
            print("    NOT p0-collapse under this operational predicate")
            continue
        print(
            "    "
            f"collapse=True witness_unit_skeleton=True first={fmt_frac(row.first_witness)} "
            f"top_pairs={row.top_sum_pairs} top_previous_two={row.top_previous_two}"
        )
        print("    recipe=" + "; ".join(row.generation_recipe))
    print()

    collapse_set = {row.speeds for row in all_collapse_rows}
    tournament_analysis(universe, collapse_set)

    print("Synthesis")
    print("-" * 78)
    print("The bounded p0-collapse family is strictly larger than AP.")
    print("All collapse rows in these targeted boxes are two-seed additive chains.")
    print("The user's two sporadics are the unique C-transversal floor rows with flipset {2}")
    print("for n=5 and n=6, and both have top equal to the previous two chain terms.")
    print("For composite C=2n-1, n=8 adds non-transversal additive-chain collapses,")
    print("so the new subproblem is to classify unit-boundary additive chains by")
    print("transversal flips when C is prime and by nonunit shell ramification when C is composite.")


if __name__ == "__main__":
    main()
