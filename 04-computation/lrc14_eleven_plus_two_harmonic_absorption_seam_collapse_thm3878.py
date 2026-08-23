#!/usr/bin/env python3
"""Exact hostile audit of the THM-3818/THM-668 dispatch join.

The checker independently rebuilds the 5,855-ratio atlas twice, reconstructs
the 46,837 cyclic-gluing seam triples, applies the harmonic-pack absorption
criterion, and classifies the residual scale-two two-branch geometry by exact
rational wall cells.
"""

from __future__ import annotations

import ast
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md",
        "9e59db3780f0897e8b978f29a8c4a9a9f640c952a31705d54d988f0f30342dcf",
    ),
    (
        "01-canon/theorems/THM-668-detuned-harmonic-dispatch.md",
        "865a8dbcfdd3ecd58d5cc98dbc7cb592ef46a59ac4290c80136a78dd1225f9e2",
    ),
    (
        "04-computation/lrc14_two_component_cyclic_gluing_extension_thm3818.py",
        "57f4ba57204a8e987f48dce46f846247fb06a0d5c5d3eb2c9c2cd7664e78d0a9",
    ),
    (
        "05-knowledge/results/lrc14_two_component_cyclic_gluing_extension_thm3818.out",
        "576995ac92b8677dea53faf811c96a7614a9f06f1257c2f7d1f4a21f2fc72586",
    ),
)
EXPECTED_SEMANTIC_SHA256 = "58a950770c4984d4a1a3f4a4031a360042ca471ed406d529dc28aad7812c1de6"
GATES = 0


def require(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def factor_trial(n: int) -> tuple[tuple[int, int], ...]:
    result = []
    prime = 2
    while prime * prime <= n:
        if n % prime == 0:
            exponent = 0
            while n % prime == 0:
                exponent += 1
                n //= prime
            result.append((prime, exponent))
        prime += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def admissible_sum_trial(n: int) -> bool:
    factors = factor_trial(n)
    return bool(factors) and all(p % 3 == 2 and e <= 2 for p, e in factors)


def atlas_trial() -> tuple[tuple[int, int], ...]:
    return tuple(
        (p, q)
        for p in range(1, 356)
        for q in range(p + 1, 357 - p)
        if p + q <= 356
        and gcd(p, q) == 1
        and admissible_sum_trial(p + q)
    )


def smallest_prime_factors(limit: int) -> list[int]:
    spf = list(range(limit + 1))
    for p in range(2, int(limit ** 0.5) + 1):
        if spf[p] == p:
            for multiple in range(p * p, limit + 1, p):
                if spf[multiple] == multiple:
                    spf[multiple] = p
    return spf


def admissible_sum_sieve(n: int, spf: list[int]) -> bool:
    while n > 1:
        p = spf[n]
        exponent = 0
        while n % p == 0:
            exponent += 1
            n //= p
        if p % 3 != 2 or exponent > 2:
            return False
    return True


def atlas_sieve() -> tuple[tuple[int, int], ...]:
    spf = smallest_prime_factors(356)
    result = []
    for total in range(3, 357):
        if not admissible_sum_sieve(total, spf):
            continue
        for p in range(1, (total + 1) // 2):
            q = total - p
            if p < q and gcd(p, total) == 1:
                result.append((p, q))
    return tuple(sorted(result))


def divisors_at_least_two(n: int) -> tuple[int, ...]:
    return tuple(d for d in range(2, n + 1) if n % d == 0)


def seams_scan(atlas: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (p, q, s)
        for p, q in atlas
        for s in range(2, q + 1)
        if s == 2 or p % s == 0 or q % s == 0
    )


def seams_divisor_union(
    atlas: tuple[tuple[int, int], ...]
) -> tuple[tuple[int, int, int], ...]:
    result = []
    for p, q in atlas:
        scales = {2}
        scales.update(divisors_at_least_two(p))
        scales.update(divisors_at_least_two(q))
        result.extend((p, q, s) for s in sorted(scales))
    return tuple(result)


def distance(x: Fraction) -> Fraction:
    x %= 1
    return min(x, 1 - x)


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return (a, 1, 0)
    common, x1, y1 = extended_gcd(b, a % b)
    return (common, y1, x1 - (a // b) * y1)


def branch_clearance(p: int, q: int, z: Fraction, branch: int) -> Fraction:
    phase = z + Fraction(branch, 2)
    return min(distance(p * phase), distance(q * phase))


def both_branches_bad(p: int, q: int, z: Fraction) -> bool:
    return (
        branch_clearance(p, q, z, 0) < Fraction(1, 14)
        and branch_clearance(p, q, z, 1) < Fraction(1, 14)
    )


def exact_wall_witness(p: int, q: int) -> Fraction | None:
    """Return a bad open-cell midpoint, or None if every phase has a safe branch."""
    endpoints = {Fraction(0), Fraction(1)}
    for speed in (p, q):
        for branch in (0, 1):
            for integer in range(speed):
                for sign in (-1, 1):
                    endpoint = (
                        (Fraction(integer) + Fraction(sign, 14)) / speed
                        - Fraction(branch, 2)
                    ) % 1
                    endpoints.add(endpoint)
    walls = sorted(endpoints)
    for left, right in zip(walls, walls[1:]):
        if left < right:
            midpoint = (left + right) / 2
            if both_branches_bad(p, q, midpoint):
                return midpoint
    # The inserted 0 and 1 already split the circular wrap cell.
    return None


def wall_phase_witnesses(pairs: tuple[tuple[int, int], ...]):
    return tuple((p, q, exact_wall_witness(p, q)) for p, q in pairs)


def semantic_digest(value: object) -> str:
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    own_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(own_tree)), "no assert")
    for relative, expected in DEPENDENCIES:
        require(sha256(ROOT / relative) == expected, ("dependency hash", relative))

    atlas_a = atlas_trial()
    atlas_b = atlas_sieve()
    require(atlas_a == atlas_b, "independent atlas equality")
    require(len(atlas_a) == 5855, "atlas count")

    seams_a = seams_scan(atlas_a)
    seams_b = seams_divisor_union(atlas_b)
    require(seams_a == seams_b, "independent seam equality")
    require(len(seams_a) == 46837, "original positive-scale seam count")
    require(sum(s == 2 for _, _, s in seams_a) == 5855, "scale-two seam count")
    require(sum(s >= 3 for _, _, s in seams_a) == 40982, "large-scale seam count")

    dispatch_closed = []
    residual_scale_two = []
    for p, q, s in seams_a:
        if s >= 3:
            require(p % s == 0 or q % s == 0, ("tariff divisor", p, q, s))
        if p % s == 0:
            # The absorbed coordinate joins the eleven-speed component; the
            # other remains genuinely detuned for every unit t modulo s.
            require(gcd(q, s) == 1, ("other coordinate coprime", p, q, s))
            dispatch_closed.append((p, q, s, "p"))
        elif q % s == 0:
            require(gcd(p, s) == 1, ("other coordinate coprime", p, q, s))
            dispatch_closed.append((p, q, s, "q"))
        else:
            require(s == 2 and p % 2 == 1 and q % 2 == 1,
                    ("only odd scale-two residual", p, q, s))
            residual_scale_two.append((p, q))

    require(len(dispatch_closed) == 45186, "THM-668 dispatch closures")
    require(len(residual_scale_two) == 1651, "odd scale-two residual")
    require(len(set(residual_scale_two)) == 1651, "odd residual distinct")

    wall_witnesses = wall_phase_witnesses(tuple(residual_scale_two))
    witness_by_pair = {(p, q): witness for p, q, witness in wall_witnesses}
    universal_wall = tuple(
        (p, q) for p, q, witness in wall_witnesses if witness is None
    )
    universal_formula = tuple(
        (p, q) for p, q in residual_scale_two if p + q <= 7
    )
    require(universal_wall == universal_formula, "wall/formula classification")
    require(universal_wall == ((1, 3),), "unique universal atlas pair")

    # Formula audit: for coprime odd p,q, the closest p-danger and shifted
    # q-danger centres are 1/(2pq) apart; open radii sum to
    # (p+q)/(14pq).  The wall sweep independently checks the consequence.
    for p, q in residual_scale_two:
        require(gcd(2 * q, 2 * p) == 2, ("centre gcd", p, q))
        require((p + 1) % 2 == 0, ("Bezout target parity", p, q))
        common, x, y = extended_gcd(q, p)
        require(common == 1, ("primitive pair", p, q))
        multiplier = (p + 1) // 2
        a = x * multiplier
        b = -y * multiplier
        require(2 * q * a - p * (2 * b + 1) == 1,
                ("unit centre gap certificate", p, q))
        min_gap = Fraction(1, 2 * p * q)
        radius_sum = Fraction(p + q, 14 * p * q)
        predicted_obstruction = min_gap < radius_sum
        observed_obstruction = witness_by_pair[p, q] is not None
        require(predicted_obstruction == observed_obstruction,
                ("gap/wall obstruction", p, q))

    hostile_pair = (3, 7)
    hostile_phase = Fraction(13, 20)
    hostile_clearances = (
        branch_clearance(*hostile_pair, hostile_phase, 0),
        branch_clearance(*hostile_pair, hostile_phase, 1),
    )
    require(hostile_clearances == (Fraction(1, 20), Fraction(1, 20)),
            "minimal hostile clearances")
    require(both_branches_bad(*hostile_pair, hostile_phase), "minimal hostile")
    hostile_candidates = [
        (max(p, q), p + q, p, q)
        for p, q in residual_scale_two
        if witness_by_pair[p, q] is not None
    ]
    require(min(hostile_candidates) == (7, 10, 3, 7), "minimal hostile order")

    # Realize the method hostile at an actual eleven-speed good point.  This
    # shows only that the selected two lifts can both fail; the full row is
    # explicitly safe elsewhere and hence is not an LRC counterexample.
    hostile_core = tuple(range(1, 10)) + (11, 12)
    hostile_y = Fraction(3, 10)
    require(len(hostile_core) == 11, "hostile core size")
    require(all(distance(u * hostile_y) >= Fraction(1, 14) for u in hostile_core),
            "hostile core good point")
    hostile_row = tuple(2 * u for u in hostile_core) + hostile_pair
    hostile_lifts = (hostile_y / 2, (hostile_y + 1) / 2)
    hostile_lift_clearances = tuple(
        min(distance(n * time) for n in hostile_row) for time in hostile_lifts
    )
    require(hostile_lift_clearances == (Fraction(1, 20), Fraction(1, 20)),
            "actual-row failed selected lifts")
    safe_control_time = Fraction(1, 20)
    safe_control_clearance = min(
        distance(n * safe_control_time) for n in hostile_row
    )
    require(safe_control_clearance == Fraction(1, 10),
            "actual-row off-selector safe time")

    # Off-atlas positive control: geometry and cube-atlas membership are
    # deliberately separate predicates.
    require(exact_wall_witness(1, 5) is None, "off-atlas universal control")
    require((1, 5) not in set(atlas_a), "off-atlas factor-three exclusion")

    unresolved_scale_two = tuple(
        (p, q) for p, q in residual_scale_two if (p, q) not in universal_wall
    )
    require(len(unresolved_scale_two) == 1650, "final scale-two count")
    final_residual_count = len(atlas_a) + len(unresolved_scale_two)
    require(final_residual_count == 7505, "s1 plus s2 residual count")

    semantic = (
        atlas_a,
        seams_a,
        tuple(dispatch_closed),
        tuple(residual_scale_two),
        tuple(
            (p, q, None if witness is None else str(witness))
            for p, q, witness in wall_witnesses
        ),
        universal_wall,
        unresolved_scale_two,
        (hostile_pair, str(hostile_phase), tuple(map(str, hostile_clearances))),
        (
            hostile_core,
            str(hostile_y),
            hostile_row,
            tuple(map(str, hostile_lift_clearances)),
            str(safe_control_time),
            str(safe_control_clearance),
        ),
        final_residual_count,
    )
    digest = semantic_digest(semantic)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(digest == EXPECTED_SEMANTIC_SHA256, "frozen semantic digest")

    print("LRC14_THM3818_DISPATCH_JOIN_20260823")
    print("status=PROVED_BOUNDED_JOIN+FINITE_EXACT;LRC14=OPEN")
    print("source=THM3818_11+2_cyclic_seams;target=THM668_12+1_harmonic_dispatch")
    print("map=s_divides_pair_coordinate=>absorb_that_runner_into_s_times_12_pack")
    print("preserved=all_11_component_clearances+absorbed_pair_runner;lost=chosen_LRC13_time+branch_owner+arrival")
    print("sidecar=coprime_other_coordinate_and_complete_s_branch_orbit")
    print(f"dependencies={DEPENDENCIES}")
    print("atlas=(ratios=5855,positive_seams=46837,s_ge_3=40982,s_eq_2=5855)")
    print("dispatch=(closed=45186,s2_even_coordinate=4204,s2_odd_residual=1651)")
    print("scale_two_geometry=universal_iff_p_plus_q_le_7_for_coprime_odd_pair")
    print("atlas_universal=((1,3),);unresolved_s2=1650;s1=5855;combined_residual=7505")
    print("minimal_method_hostile=(p,q,s,z)=(3,7,2,13/20);branch_clearances=(1/20,1/20)")
    print("hostile_row=2*(1,2,3,4,5,6,7,8,9,11,12)+(3,7);selected_lifts=(3/20,13/20)_give_1/20;off_selector_t=1/20_gives_1/10")
    print("hostile_scope=two_branch_pair_selector_only;not_an_LRC_counterexample")
    print("factorial_connection=real_selector_analogy_only;local_divisor_filter_fails_until_adaptive_sidecar;here_divisibility_changes_the_object_by_absorption_and_is_proof_forcing")
    print(f"semantic_sha256={digest}")
    print(f"gates={GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
