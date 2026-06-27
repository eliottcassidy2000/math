#!/usr/bin/env python3
"""Sixth-power equal-sum collision scout for LRC14 proof carriers.

This is not a Diophantine proof and not a proof of LRC14.  It extends the
S244/HYP-3076 sixth-power collision sidecar with a bounded
certificate-ledger test for the two equations named in the prompt:

    a^6 + b^6 = d^6 + e^6
    a^6 + b^6 + c^6 = d^6 + e^6 + f^6

The LRC use is controlled forgetting.  S244 treats a 3-vs-3 equality as native
support-six relation data and a 2-vs-2 equality as a rank-lowered square-cube
shadow that needs padding if lifted to six slots.  S248 adds the certificate
ledger around that split: a two-lane equality retains unordered pair identity,
primitive gcd, and residue payload rather than trusting the scalar sum, while a
three-lane equality carries a collision certificate before a route triple is
medianized.

Tournament Analysis vertices are proof obligations / sidecar carriers, not
runners, not arcs, and not the numbers a,b,c,d,e,f.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement
from math import gcd
from typing import Dict, Iterable, List, Sequence, Tuple


PAIR_BOUND = 250
TRIPLE_BOUND = 80
SAMPLE_COLLISIONS = 12


Pair = Tuple[int, int]
Triple = Tuple[int, int, int]


def tuple_gcd(values: Iterable[int]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def residue_word(values: Sequence[int], modulus: int) -> Tuple[int, ...]:
    return tuple(sorted(value % modulus for value in values))


def power_residue_word(values: Sequence[int], modulus: int) -> Tuple[int, ...]:
    return tuple(sorted(pow(value, 6, modulus) for value in values))


def sixth_power_residue_alphabet(modulus: int) -> Tuple[int, ...]:
    return tuple(sorted({pow(value, 6, modulus) for value in range(modulus)}))


def find_pair_collisions(bound: int) -> Dict[int, List[Pair]]:
    sums: Dict[int, List[Pair]] = defaultdict(list)
    powers = {i: i**6 for i in range(1, bound + 1)}
    for a, b in combinations_with_replacement(range(1, bound + 1), 2):
        sums[powers[a] + powers[b]].append((a, b))
    return {total: pairs for total, pairs in sums.items() if len(pairs) > 1}


def find_triple_collisions(bound: int) -> Dict[int, List[Triple]]:
    sums: Dict[int, List[Triple]] = defaultdict(list)
    powers = {i: i**6 for i in range(1, bound + 1)}
    for triple in combinations_with_replacement(range(1, bound + 1), 3):
        sums[sum(powers[x] for x in triple)].append(triple)
    return {total: triples for total, triples in sums.items() if len(triples) > 1}


def collision_records(
    triple_collisions: Dict[int, List[Triple]], limit: int
) -> List[Dict[str, object]]:
    records: List[Dict[str, object]] = []
    for total in sorted(triple_collisions):
        triples = triple_collisions[total]
        for left, right in combinations(triples, 2):
            if left == right:
                continue
            shared = sorted((Counter(left) & Counter(right)).elements())
            record = {
                "sum": total,
                "left": left,
                "right": right,
                "all_gcd": tuple_gcd(left + right),
                "left_gcd": tuple_gcd(left),
                "right_gcd": tuple_gcd(right),
                "shared_terms": tuple(shared),
                "shared_term_count": len(shared),
                "primitive_all_terms": tuple_gcd(left + right) == 1,
                "mod14_base": (residue_word(left, 14), residue_word(right, 14)),
                "mod14_pow6": (
                    power_residue_word(left, 14),
                    power_residue_word(right, 14),
                ),
                "mod27_pow6": (
                    power_residue_word(left, 27),
                    power_residue_word(right, 27),
                ),
                "mod41_pow6": (
                    power_residue_word(left, 41),
                    power_residue_word(right, 41),
                ),
            }
            records.append(record)
            if len(records) >= limit:
                return records
    return records


@dataclass(frozen=True)
class Carrier:
    name: str
    payload: frozenset[str]
    note: str


AXES = (
    "lrc_predicate",
    "exact_integer_relation",
    "lane_rank",
    "primitive_gcd",
    "shared_term_filter",
    "collision_graph",
    "mod14_packet",
    "mod27_lift",
    "mod41_pressure",
    "height_lattice",
    "route_triple_median",
    "sidecar_discharge",
    "raw_scalar",
)


CARRIERS = [
    Carrier(
        "labelled_lrc_packet_sheaf",
        frozenset(AXES) - {"raw_scalar"},
        "keeps the proof predicate above every quotient",
    ),
    Carrier(
        "sixth_power_collision_certificate",
        frozenset(
            {
                "lrc_predicate",
                "exact_integer_relation",
                "lane_rank",
                "primitive_gcd",
                "shared_term_filter",
                "collision_graph",
                "mod14_packet",
                "mod27_lift",
                "mod41_pressure",
                "sidecar_discharge",
            }
        ),
        "names the equality, lane rank, gcd, shared terms, and CRT payload",
    ),
    Carrier(
        "two_lane_rigidity_gate",
        frozenset(
            {
                "lrc_predicate",
                "exact_integer_relation",
                "lane_rank",
                "primitive_gcd",
                "shared_term_filter",
                "mod14_packet",
                "sidecar_discharge",
            }
        ),
        "uses the bounded absence of two-lane collisions as a no-hidden-pair gate",
    ),
    Carrier(
        "three_lane_resonance_graph",
        frozenset(
            {
                "exact_integer_relation",
                "lane_rank",
                "collision_graph",
                "primitive_gcd",
                "shared_term_filter",
                "route_triple_median",
                "sidecar_discharge",
            }
        ),
        "turns triple equal sums into route-triple sidecars",
    ),
    Carrier(
        "route_state_median_sidecar",
        frozenset(
            {
                "lrc_predicate",
                "route_triple_median",
                "sidecar_discharge",
                "collision_graph",
                "lane_rank",
            }
        ),
        "tests whether the collision certificate has a legal center after closure",
    ),
    Carrier(
        "ramanujan_crt_residue_word",
        frozenset(
            {
                "exact_integer_relation",
                "mod14_packet",
                "mod27_lift",
                "mod41_pressure",
                "sidecar_discharge",
            }
        ),
        "keeps sixth-power residue words as exact-period packet payload",
    ),
    Carrier(
        "roth_minkowski_height_fence",
        frozenset(
            {
                "exact_integer_relation",
                "primitive_gcd",
                "height_lattice",
                "sidecar_discharge",
            }
        ),
        "records low-height walls and lattice height before using finiteness pressure",
    ),
    Carrier(
        "fermat_catalan_hyperbolic_pressure",
        frozenset(
            {
                "lane_rank",
                "primitive_gcd",
                "height_lattice",
                "mod41_pressure",
            }
        ),
        "treats sixth powers as hyperbolic pressure, not as a theorem shortcut",
    ),
    Carrier(
        "partial_cube_simplex_lane",
        frozenset(
            {
                "lane_rank",
                "route_triple_median",
                "collision_graph",
                "sidecar_discharge",
            }
        ),
        "interfaces the three-lane equation with the Moser/fibbinary/simplex route",
    ),
    Carrier(
        "raw_equal_sum_scalar_shadow",
        frozenset({"raw_scalar"}),
        "remembers only the value of the sum and is intentionally unsafe",
    ),
]


RIGIDITY_WEIGHTS = {
    "lrc_predicate": 5,
    "exact_integer_relation": 4,
    "lane_rank": 4,
    "primitive_gcd": 3,
    "shared_term_filter": 4,
    "mod14_packet": 2,
    "mod27_lift": 1,
    "mod41_pressure": 1,
    "collision_graph": 1,
    "route_triple_median": 1,
    "sidecar_discharge": 3,
    "height_lattice": 1,
    "raw_scalar": -2,
}


RESONANCE_WEIGHTS = {
    "lrc_predicate": 4,
    "exact_integer_relation": 3,
    "lane_rank": 2,
    "primitive_gcd": 2,
    "shared_term_filter": 2,
    "mod14_packet": 2,
    "mod27_lift": 2,
    "mod41_pressure": 2,
    "collision_graph": 5,
    "route_triple_median": 5,
    "sidecar_discharge": 4,
    "height_lattice": 2,
    "raw_scalar": -2,
}


def carrier_score(carrier: Carrier, weights: Dict[str, int]) -> int:
    return sum(weights.get(axis, 0) for axis in carrier.payload)


def build_tournament(weights: Dict[str, int]) -> Dict[Tuple[str, str], str]:
    order = {carrier.name: idx for idx, carrier in enumerate(CARRIERS)}
    scores = {carrier.name: carrier_score(carrier, weights) for carrier in CARRIERS}
    edges: Dict[Tuple[str, str], str] = {}
    for left, right in combinations(CARRIERS, 2):
        ls = scores[left.name]
        rs = scores[right.name]
        if ls > rs or (ls == rs and order[left.name] < order[right.name]):
            winner, loser = left.name, right.name
        else:
            winner, loser = right.name, left.name
        edges[(winner, loser)] = winner
    return edges


def score_histogram(edges: Dict[Tuple[str, str], str]) -> Dict[int, int]:
    wins = Counter({carrier.name: 0 for carrier in CARRIERS})
    for winner in edges.values():
        wins[winner] += 1
    return dict(sorted(Counter(wins.values()).items()))


def directed_3_cycles(edges: Dict[Tuple[str, str], str]) -> int:
    names = [carrier.name for carrier in CARRIERS]

    def beats(a: str, b: str) -> bool:
        return edges.get((a, b)) == a or edges.get((b, a)) == a

    count = 0
    for a, b, c in combinations(names, 3):
        if (
            beats(a, b)
            and beats(b, c)
            and beats(c, a)
            or beats(a, c)
            and beats(c, b)
            and beats(b, a)
        ):
            count += 1
    return count


def scc_sizes(edges: Dict[Tuple[str, str], str]) -> List[int]:
    names = [carrier.name for carrier in CARRIERS]
    adj = {name: [] for name in names}
    rev = {name: [] for name in names}
    for winner, loser in edges:
        adj[winner].append(loser)
        rev[loser].append(winner)

    seen = set()
    finish: List[str] = []

    def dfs(node: str) -> None:
        seen.add(node)
        for nxt in adj[node]:
            if nxt not in seen:
                dfs(nxt)
        finish.append(node)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: List[int] = []

    def rdfs(node: str) -> int:
        seen.add(node)
        size = 1
        for nxt in rev[node]:
            if nxt not in seen:
                size += rdfs(nxt)
        return size

    for name in reversed(finish):
        if name not in seen:
            sizes.append(rdfs(name))
    return sizes


def hamiltonian_path_count(edges: Dict[Tuple[str, str], str]) -> int:
    names = [carrier.name for carrier in CARRIERS]

    def beats(a: str, b: str) -> bool:
        return edges.get((a, b)) == a or edges.get((b, a)) == a

    count = 0
    stack = deque([(name, frozenset({name})) for name in names])
    while stack:
        last, used = stack.pop()
        if len(used) == len(names):
            count += 1
            continue
        for nxt in names:
            if nxt not in used and beats(last, nxt):
                stack.append((nxt, used | {nxt}))
    return count


def edge_flip_count(
    left_edges: Dict[Tuple[str, str], str], right_edges: Dict[Tuple[str, str], str]
) -> int:
    flips = 0
    names = [carrier.name for carrier in CARRIERS]
    for a, b in combinations(names, 2):
        left_winner = left_edges.get((a, b)) or left_edges.get((b, a))
        right_winner = right_edges.get((a, b)) or right_edges.get((b, a))
        if left_winner != right_winner:
            flips += 1
    return flips


def tournament_path(edges: Dict[Tuple[str, str], str]) -> List[str]:
    names = [carrier.name for carrier in CARRIERS]
    wins = Counter({name: 0 for name in names})
    for winner in edges.values():
        wins[winner] += 1
    return [name for name, _ in sorted(wins.items(), key=lambda item: (-item[1], item[0]))]


def print_collision_record(record: Dict[str, object]) -> None:
    print(
        "  sum={sum} left={left} right={right} gcd_all={all_gcd} "
        "shared={shared_terms} mod14_pow6={mod14_pow6}".format(**record)
    )
    print(
        "    mod27_pow6={} mod41_pow6={} primitive={}".format(
            record["mod27_pow6"],
            record["mod41_pow6"],
            record["primitive_all_terms"],
        )
    )


def main() -> None:
    pair_collisions = find_pair_collisions(PAIR_BOUND)
    triple_collisions = find_triple_collisions(TRIPLE_BOUND)
    records = collision_records(triple_collisions, SAMPLE_COLLISIONS)

    pair_count = PAIR_BOUND * (PAIR_BOUND + 1) // 2
    triple_count = TRIPLE_BOUND * (TRIPLE_BOUND + 1) * (TRIPLE_BOUND + 2) // 6
    nontrivial_triple_pairs = sum(
        len(list(combinations(triples, 2))) for triples in triple_collisions.values()
    )
    primitive_records = sum(1 for record in records if record["primitive_all_terms"])
    shared_records = sum(1 for record in records if record["shared_term_count"])

    print("LRC14 S248 SIXTH-POWER CERTIFICATE EXTENSION SCOUT")
    print("================================================")
    print()
    print("Bounded Diophantine carrier scan")
    print(f"  pair equation bound: positive 1..{PAIR_BOUND}")
    print(f"  pair unordered inputs checked: {pair_count}")
    print(f"  pair sums with more than one unordered pair: {len(pair_collisions)}")
    if pair_collisions:
        first_pair_sum = min(pair_collisions)
        print(f"  first pair collision: sum={first_pair_sum} pairs={pair_collisions[first_pair_sum]}")
    else:
        print("  first pair collision: none found in this positive bounded scan")
    print()
    print(f"  triple equation bound: positive 1..{TRIPLE_BOUND}")
    print(f"  triple unordered inputs checked: {triple_count}")
    print(f"  triple sums with more than one unordered triple: {len(triple_collisions)}")
    print(f"  nontrivial triple-pair collision certificates: {nontrivial_triple_pairs}")
    print(f"  sample primitive certificates: {primitive_records}/{len(records)}")
    print(f"  sample certificates with shared terms: {shared_records}/{len(records)}")
    print()
    print("Sixth-power residue alphabets")
    for modulus in (14, 27, 41):
        print(f"  mod {modulus}: {sixth_power_residue_alphabet(modulus)}")
    print()
    print("First triple collision certificates")
    for record in records:
        print_collision_record(record)
    print()
    print("LRC carrier interpretation")
    print("  two_lane_rigidity_gate:")
    print(
        "    The bounded pair scan found no positive non-permutation collision, so a "
        "two-lane equality should be treated as a no-hidden-pair gate with gcd, "
        "residue, and shared-term sidecars retained."
    )
    print("  three_lane_resonance_graph:")
    print(
        "    The triple scan has bounded collisions, so a route triple can balance "
        "sixth-power mass without identifying lanes.  A legal proof state must "
        "carry a collision certificate before medianization."
    )
    print("  controlled_forgetting_rule:")
    print(
        "    Forgetting from the certificate to the raw scalar sum destroys lane "
        "rank, primitive gcd, shared-term cancellation, and CRT residue words. "
        "The loss must be reconstructed, dual-annihilated, descended, stopped at "
        "AP/GW boundary, or routed to named F7/THM-572 debt."
    )
    print()
    print("Tournament Analysis")
    print("  vertices: proof obligations / sidecar carriers, not runners or numbers")
    print("  pairwise observable: weighted retained payload axes")
    print("  tie path: listed carrier order in the script")
    rigidity = build_tournament(RIGIDITY_WEIGHTS)
    resonance = build_tournament(RESONANCE_WEIGHTS)
    for name, edges in (("rank2_rigidity_gauge", rigidity), ("rank3_resonance_gauge", resonance)):
        print(f"  {name}:")
        print(f"    score_histogram={score_histogram(edges)}")
        print(f"    directed_3_cycles={directed_3_cycles(edges)}")
        print(f"    scc_sizes={scc_sizes(edges)}")
        print(f"    hamiltonian_path_count={hamiltonian_path_count(edges)}")
        print(f"    path={' > '.join(tournament_path(edges))}")
    print(f"  edge_flips_between_gauges={edge_flip_count(rigidity, resonance)}")
    print()
    print("Assumption challenge")
    print(
        "  Do not map the carrier directly to runners or powers.  Useful vertices "
        "are the two-lane rigidity gate, the three-lane resonance graph, primitive "
        "gcd filters, shared-term cancellation filters, mod-14/mod-27/mod-41 "
        "residue packets, height/lattice fences, route-state median sidecars, and "
        "raw scalar shadows.  The LRC predicate preserved is legal discharge of "
        "packet/route/certificate/sidecar states; the destroyed information is "
        "the lane identity and residue/gcd payload behind equal sixth-power sums."
    )


if __name__ == "__main__":
    main()
