#!/usr/bin/env python3
"""Independent hostile audit of the THM-3818 / detuned-dispatch join.

This implementation deliberately does not import the primary checker.  It
builds the atlas by admissible totals and Euler's totient identity, counts
seams with divisor-count formulas, constructs every obstructed odd-pair phase
from a modular inverse, and computes the full hostile-row maximum on pair-sum
events.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
ROOT = Path(__file__).resolve().parents[1]
GATES = 0


def require(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def factor(n: int) -> tuple[tuple[int, int], ...]:
    answer: list[tuple[int, int]] = []
    d = 2
    while d * d <= n:
        exponent = 0
        while n % d == 0:
            n //= d
            exponent += 1
        if exponent:
            answer.append((d, exponent))
        d = 3 if d == 2 else d + 2
    if n > 1:
        answer.append((n, 1))
    return tuple(answer)


def admissible_total(n: int) -> bool:
    return all(prime % 3 == 2 and exponent <= 2 for prime, exponent in factor(n))


def phi_from_factorization(n: int) -> int:
    answer = n
    for prime, _ in factor(n):
        answer = answer // prime * (prime - 1)
    return answer


def divisor_count(n: int) -> int:
    answer = 1
    for _, exponent in factor(n):
        answer *= exponent + 1
    return answer


def proper_nonunit_divisors(n: int) -> set[int]:
    return {d for d in range(2, n + 1) if n % d == 0}


def build_atlas() -> tuple[tuple[int, int], ...]:
    atlas: list[tuple[int, int]] = []
    for total in range(3, 357):
        if not admissible_total(total):
            continue
        block = tuple(
            (p, total - p)
            for p in range(1, (total + 1) // 2)
            if gcd(p, total) == 1
        )
        require(len(block) == phi_from_factorization(total) // 2,
                ("totient block", total))
        atlas.extend(block)
    return tuple(atlas)


def dist(x: Fraction) -> Fraction:
    x %= 1
    return min(x, 1 - x)


def pair_branch_clearance(p: int, q: int, z: Fraction, branch: int) -> Fraction:
    phase = z + Fraction(branch, 2)
    return min(dist(p * phase), dist(q * phase))


def obstructing_phase(p: int, q: int) -> Fraction:
    """Construct a point in D_p intersect (D_q-1/2) when p+q>7.

    For p,q odd and coprime, choose centres c_p=a/p and
    c_q=(2b+1)/(2q) with c_p-c_q=1/(2pq).  Weighting the two centres by the
    opposite radii puts z strictly inside both open intervals exactly when
    their radii sum exceeds that gap.
    """

    target = (p + 1) // 2
    if p == 1:
        a = 0
    else:
        a = (target * pow(q, -1, p)) % p
    b_numerator = q * a - target
    require(b_numerator % p == 0, ("modular inverse certificate", p, q))
    b = b_numerator // p
    c_p = Fraction(a, p)
    c_q = Fraction(2 * b + 1, 2 * q)
    require(c_p - c_q == Fraction(1, 2 * p * q),
            ("nearest cross centres", p, q))
    radius_p = Fraction(1, 14 * p)
    radius_q = Fraction(1, 14 * q)
    require(c_p - c_q < radius_p + radius_q,
            ("strict cross overlap", p, q))
    z = (radius_q * c_p + radius_p * c_q) / (radius_p + radius_q)
    require(abs(z - c_p) < radius_p, ("inside p danger", p, q))
    require(abs(z - c_q) < radius_q, ("inside shifted q danger", p, q))
    z %= 1
    require(pair_branch_clearance(p, q, z, 0) < Fraction(1, 14),
            ("branch zero killed", p, q))
    require(pair_branch_clearance(p, q, z, 1) < Fraction(1, 14),
            ("branch one killed", p, q))
    return z


def event_maximum(row: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    best_clearance = Fraction(0)
    best_time = Fraction(0)
    for i, left in enumerate(row):
        for right in row[i:]:
            modulus = left + right
            for numerator in range(1, modulus):
                time = Fraction(numerator, modulus)
                clearance = min(dist(speed * time) for speed in row)
                if clearance > best_clearance:
                    best_clearance = clearance
                    best_time = time
    return best_clearance, best_time


def decoder_components(row: tuple[int, ...], atlas: set[tuple[int, int]]) -> int:
    parent = list(range(len(row)))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for i, left in enumerate(row):
        for j in range(i + 1, len(row)):
            right = row[j]
            common = gcd(left, right)
            ratio = tuple(sorted((left // common, right // common)))
            if ratio in atlas:
                a, b = root(i), root(j)
                parent[a] = b
    return len({root(vertex) for vertex in range(len(row))})


def main() -> None:
    own_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(own_tree)),
            "assertion-free")

    atlas = build_atlas()
    require(len(atlas) == 5855, "atlas count")
    require(len(set(atlas)) == len(atlas), "atlas uniqueness")

    seams: list[tuple[int, int, int]] = []
    formula_count = 0
    for p, q in atlas:
        scales = {2} | proper_nonunit_divisors(p) | proper_nonunit_divisors(q)
        seams.extend((p, q, scale) for scale in sorted(scales))
        local_formula = divisor_count(p) + divisor_count(q) - 2
        if p % 2 and q % 2:
            local_formula += 1
        require(len(scales) == local_formula, ("divisor formula", p, q))
        formula_count += local_formula
    require(formula_count == len(seams) == 46837, "seam count")

    large = tuple(seam for seam in seams if seam[2] >= 3)
    require(len(large) == 40982, "large-scale count")
    require(all(p % scale == 0 or q % scale == 0 for p, q, scale in large),
            "large scale is a divisor seam")

    dispatch = tuple(
        (p, q, scale)
        for p, q, scale in seams
        if p % scale == 0 or q % scale == 0
    )
    odd_scale_two = tuple(
        (p, q)
        for p, q, scale in seams
        if scale == 2 and p % 2 and q % 2
    )
    require(len(dispatch) == 45186, "dispatch count")
    require(len(odd_scale_two) == 1651, "odd scale-two count")
    require(len(dispatch) + len(odd_scale_two) == len(seams), "seam partition")

    # Audit the divisibility direction needed by THM-668.  If scale divides
    # one primitive coordinate, it is coprime to the other.  Any t coprime to
    # scale therefore keeps the other coordinate genuinely detuned and even
    # gives a complete branch orbit.
    for p, q, scale in dispatch:
        if p % scale == 0:
            require(gcd(scale, q) == 1, ("p absorbed, q detuned", p, q, scale))
        else:
            require(q % scale == 0, ("q divisor", p, q, scale))
            require(gcd(scale, p) == 1, ("q absorbed, p detuned", p, q, scale))

    universal = tuple((p, q) for p, q in odd_scale_two if p + q <= 7)
    obstructed = tuple((p, q) for p, q in odd_scale_two if p + q > 7)
    require(universal == ((1, 3),), "unique universal atlas pair")
    require(len(obstructed) == 1650, "obstructed pair count")

    # The two same-speed terms are always disjoint: their centre gap is
    # 1/(2w), versus total open radius 1/(7w).  Cross terms are disjoint iff
    # 1/(2pq) >= (p+q)/(14pq), i.e. p+q<=7.
    for p, q in odd_scale_two:
        for speed in (p, q):
            require(Fraction(1, 2 * speed) > Fraction(1, 7 * speed),
                    ("same-speed half-shift", speed))
        cross_gap = Fraction(1, 2 * p * q)
        radius_sum = Fraction(p + q, 14 * p * q)
        require((cross_gap >= radius_sum) == (p + q <= 7),
                ("open-boundary criterion", p, q))
        if p + q > 7:
            obstructing_phase(p, q)

    # Positive/off-atlas and first-hostile controls straddle the strict wall.
    require(Fraction(1, 2 * 1 * 5) >= Fraction(1 + 5, 14 * 1 * 5),
            "(1,5) universally safe geometry")
    require((1, 5) not in set(atlas), "(1,5) excluded only by atlas")
    require(Fraction(1, 2 * 3 * 5) < Fraction(3 + 5, 14 * 3 * 5),
            "sum-eight strict overlap control")

    hostile_order = min((q, p + q, p, q) for p, q in obstructed)
    require(hostile_order == (7, 10, 3, 7), "minimal hostile order")
    z = Fraction(13, 20)
    require(tuple(pair_branch_clearance(3, 7, z, branch) for branch in (0, 1))
            == (Fraction(1, 20), Fraction(1, 20)), "named hostile phase")

    core = tuple(range(1, 10)) + (11, 12)
    y = Fraction(3, 10)
    row = tuple(2 * speed for speed in core) + (3, 7)
    require(all(dist(speed * y) >= Fraction(1, 14) for speed in core),
            "eleven-speed good point")
    selected = (y / 2, (y + 1) / 2)
    require(tuple(min(dist(speed * time) for speed in row) for time in selected)
            == (Fraction(1, 20), Fraction(1, 20)), "both selected lifts fail")
    exact_maximum, exact_time = event_maximum(row)
    require((exact_maximum, exact_time) == (Fraction(1, 10), Fraction(1, 20)),
            "full row exact positive control")

    # Crucial scope check: this concrete positive control is not an example
    # in the rank-eleven two-component branch; its packet decoder is connected.
    components = decoder_components(row, set(atlas))
    require(components == 1, "hostile control decoder is connected")

    semantic = {
        "atlas": len(atlas),
        "seams": len(seams),
        "large": len(large),
        "dispatch": len(dispatch),
        "odd_scale_two": len(odd_scale_two),
        "universal": universal,
        "obstructed": len(obstructed),
        "combined_residual": len(atlas) + len(obstructed),
        "minimal_hostile": hostile_order,
        "hostile_row_maximum": str(exact_maximum),
        "hostile_row_time": str(exact_time),
        "hostile_row_decoder_components": components,
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    print("LRC14_THM3818_DISPATCH_JOIN_INDEPENDENT_AUDIT_20260823")
    print("status=PASS;scope=bounded_rank11_11+2_join;LRC14=OPEN")
    print("counts=5855_atlas;46837_positive_seams;45186_dispatch_closed;1651_odd_s2;1_universal;7505_final_residual")
    print("logic=s_divides_one_coordinate+primitive_pair+gcd(s,t)=1=>other_runner_is_detuned_with_complete_s_orbit")
    print("geometry=self_terms_disjoint;cross_gap=1/(2pq);open_radius_sum=(p+q)/(14pq);universal_iff_p+q<=7")
    print("hostile=(3,7,2),z=13/20;both_selected_clearances=1/20")
    print("actual_row_control=M=1/10_at_1/20;decoder_components=1;not_a_rank11_two_component_example_and_not_an_LRC_hostile")
    print(f"semantic_sha256={digest}")
    print(f"gates={GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
