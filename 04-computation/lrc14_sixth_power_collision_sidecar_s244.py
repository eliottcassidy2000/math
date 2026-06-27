#!/usr/bin/env python3
"""Sixth-power collision sidecar scout for LRC14.

The prompt equations are

    a^6 + b^6 + c^6 = d^6 + e^6 + f^6
    a^6 + b^6       = d^6 + e^6

This script treats them as proof-interface carriers, not as imported number
theory theorems.  The 3-vs-3 equation is a native six-term signed relation and
therefore fits the THM-538 support-six / coimage-tail layer.  The 2-vs-2
equation is a rank-lowered square-cube shadow because x^6=(x^3)^2; to enter
the support-six layer it must be padded by a cancelling pair, which should be
recorded as a degenerate sidecar rather than a genuine six-term wall.

Tournament Analysis vertices are proof carriers and sidecar fields, not
runners or raw bases.  The pairwise observable is retained LRC payload:
native support-six status, owner visibility, residue phase, low-height wall
use, legal route-state closure, and degeneracy protection.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement
from math import gcd
from typing import Dict, Iterable, List, Sequence, Tuple


MODULI = (7, 9, 13, 14, 21, 27, 28, 41)


def tuple_gcd(values: Iterable[int]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def sixth_residue_set(modulus: int, units_only: bool = False) -> Tuple[int, ...]:
    vals = []
    for x in range(modulus):
        if units_only and gcd(x, modulus) != 1:
            continue
        vals.append(pow(x, 6, modulus))
    return tuple(sorted(set(vals)))


def sum_bucket_profile(modulus: int, width: int) -> Dict[str, object]:
    residues = sixth_residue_set(modulus)
    buckets: Dict[int, List[Tuple[int, ...]]] = defaultdict(list)
    for combo in combinations_with_replacement(residues, width):
        buckets[sum(combo) % modulus].append(combo)
    sizes = sorted((len(v) for v in buckets.values()), reverse=True)
    return {
        "modulus": modulus,
        "width": width,
        "sixth_residues": residues,
        "bucket_count": len(buckets),
        "collision_buckets": sum(1 for v in buckets.values() if len(v) > 1),
        "max_bucket": sizes[0],
        "top_bucket_sizes": sizes[:5],
    }


def pair_collisions(bound: int) -> List[Dict[str, object]]:
    seen: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    hits: List[Dict[str, object]] = []
    for a in range(1, bound + 1):
        a6 = a**6
        for b in range(a, bound + 1):
            pair = (a, b)
            total = a6 + b**6
            for prev in seen[total]:
                if set(prev) != set(pair):
                    vals = prev + pair
                    hits.append(
                        {
                            "left": prev,
                            "right": pair,
                            "sum": total,
                            "gcd": tuple_gcd(vals),
                        }
                    )
            seen[total].append(pair)
    hits.sort(key=lambda h: (h["sum"], h["left"], h["right"]))
    return hits


def triple_collisions(bound: int) -> List[Dict[str, object]]:
    seen: Dict[int, List[Tuple[int, int, int]]] = defaultdict(list)
    hits: List[Dict[str, object]] = []
    for a in range(1, bound + 1):
        a6 = a**6
        for b in range(a, bound + 1):
            ab = a6 + b**6
            for c in range(b, bound + 1):
                triple = (a, b, c)
                total = ab + c**6
                for prev in seen[total]:
                    if set(prev) != set(triple):
                        vals = prev + triple
                        hits.append(
                            {
                                "left": prev,
                                "right": triple,
                                "sum": total,
                                "gcd": tuple_gcd(vals),
                                "primitive": tuple_gcd(vals) == 1,
                            }
                        )
                seen[total].append(triple)
    hits.sort(key=lambda h: (h["sum"], h["left"], h["right"]))
    return hits


def residue_signature(left: Sequence[int], right: Sequence[int]) -> Dict[int, Dict[str, object]]:
    out: Dict[int, Dict[str, object]] = {}
    for modulus in MODULI:
        left6 = tuple(pow(x, 6, modulus) for x in left)
        right6 = tuple(pow(x, 6, modulus) for x in right)
        signed = left6 + tuple((-x) % modulus for x in right6)
        out[modulus] = {
            "left_base": tuple(x % modulus for x in left),
            "right_base": tuple(x % modulus for x in right),
            "left_sixth": left6,
            "right_sixth": right6,
            "signed_sixth": signed,
            "left_sum": sum(left6) % modulus,
            "right_sum": sum(right6) % modulus,
            "signed_sum": sum(signed) % modulus,
            "left_mask": tuple(sorted(Counter(left6).items())),
            "right_mask": tuple(sorted(Counter(right6).items())),
        }
    return out


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: Tuple[int, int, int, int, int, int, int]


SCORE_LABELS = (
    "preserves_lrc_predicate",
    "native_support6",
    "owner_visibility",
    "residue_phase_visible",
    "low_height_wall_use",
    "route_state_legal",
    "guards_degeneracy",
)


CARRIERS = (
    Carrier("labelled_packet_sheaf", (5, 5, 5, 5, 5, 5, 5)),
    Carrier("native_three_vs_three_support6_collision", (5, 5, 5, 5, 4, 4, 5)),
    Carrier("sixth_power_residue_phase_mask", (4, 4, 4, 5, 4, 4, 5)),
    Carrier("owner_gcd_common_factor_gate", (4, 3, 5, 3, 4, 4, 5)),
    Carrier("low_height_wall_ledger", (4, 4, 4, 4, 5, 4, 4)),
    Carrier("route_state_closure_sidecar", (5, 3, 4, 4, 4, 5, 4)),
    Carrier("rank_lowered_two_vs_two_square_cube_shadow", (3, 1, 3, 3, 2, 3, 5)),
    Carrier("padded_support6_canceling_pair", (3, 2, 3, 3, 2, 3, 4)),
    Carrier("raw_equal_sixth_power_scalar", (1, 0, 0, 1, 0, 1, 0)),
)


def carrier_total(carrier: Carrier) -> int:
    return sum(carrier.scores)


def orient(a: Carrier, b: Carrier) -> bool:
    if carrier_total(a) != carrier_total(b):
        return carrier_total(a) > carrier_total(b)
    if a.scores != b.scores:
        return a.scores > b.scores
    return CARRIERS.index(a) < CARRIERS.index(b)


def tournament_fingerprint() -> Dict[str, object]:
    n = len(CARRIERS)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(CARRIERS):
        for j, b in enumerate(CARRIERS):
            if i != j:
                adj[i][j] = orient(a, b)
    score_hist = Counter(sum(row) for row in adj)
    cycles3 = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles3 += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cycles3 += 1

    # The score/lexicographic orientation is transitive in this scout, so the
    # Hamiltonian path is the total carrier order.  Count it by tiny DP anyway.
    dp: Dict[Tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + ways
    full = (1 << n) - 1
    hpaths = sum(dp.get((full, i), 0) for i in range(n))
    path = sorted(CARRIERS, key=lambda c: (carrier_total(c), c.scores), reverse=True)
    return {
        "score_labels": SCORE_LABELS,
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles3,
        "scc_sizes": [1] * n,
        "hamiltonian_path_count": hpaths,
        "hamiltonian_path": [c.name for c in path],
    }


def print_residue_profiles() -> None:
    print("SIXTH-POWER RESIDUE BUCKET PROFILES")
    for modulus in MODULI:
        residues = sixth_residue_set(modulus)
        unit_residues = sixth_residue_set(modulus, units_only=True)
        prof2 = sum_bucket_profile(modulus, 2)
        prof3 = sum_bucket_profile(modulus, 3)
        print(
            f"mod {modulus}: residues={residues} units={unit_residues} "
            f"2-sum buckets={prof2['bucket_count']} colliding={prof2['collision_buckets']} "
            f"max={prof2['max_bucket']} | "
            f"3-sum buckets={prof3['bucket_count']} colliding={prof3['collision_buckets']} "
            f"max={prof3['max_bucket']}"
        )


def print_collision(title: str, hit: Dict[str, object]) -> None:
    left = hit["left"]
    right = hit["right"]
    print(title)
    print(f"  left={left} right={right} sum={hit['sum']} gcd={hit['gcd']}")
    if "primitive" in hit:
        print(f"  primitive={hit['primitive']}")
    sig = residue_signature(left, right)
    for modulus in MODULI:
        row = sig[modulus]
        print(
            f"  mod {modulus}: base={row['left_base']}|{row['right_base']} "
            f"sixth={row['left_sixth']}|{row['right_sixth']} "
            f"signed_sum={row['signed_sum']} masks={row['left_mask']}|{row['right_mask']}"
        )


def main() -> None:
    print("LRC14 SIXTH-POWER COLLISION SIDECAR S244")
    print("Equation A: a^6+b^6+c^6=d^6+e^6+f^6")
    print("Equation B: a^6+b^6=d^6+e^6")
    print()

    print_residue_profiles()
    print()

    pair_bound = 220
    pair_hits = pair_collisions(pair_bound)
    print(f"2-vs-2 bounded search: bound={pair_bound} nontrivial_hits={len(pair_hits)}")
    if pair_hits:
        print_collision("first 2-vs-2 hit", pair_hits[0])
    else:
        print(
            "  no nontrivial 2-vs-2 hit found; treat as a rank-lowered "
            "square-cube shadow, not a native support-six wall"
        )
    print()

    triple_bound = 80
    triple_hits = triple_collisions(triple_bound)
    primitive_hits = [h for h in triple_hits if h["primitive"]]
    print(
        f"3-vs-3 bounded search: bound={triple_bound} "
        f"collision_pairs={len(triple_hits)} primitive_pairs={len(primitive_hits)}"
    )
    for idx, hit in enumerate(triple_hits[:5], start=1):
        print(
            f"  hit {idx}: {hit['left']} = {hit['right']} "
            f"sum={hit['sum']} gcd={hit['gcd']} primitive={hit['primitive']}"
        )
    print()
    if primitive_hits:
        print_collision("first primitive 3-vs-3 support-six wall", primitive_hits[0])
    print()

    print("SIDECAR INTERPRETATION")
    print(
        "  3-vs-3: native six-term relation; use as support6_power_collision "
        "with residue masks, owner gcd, mod-7 unit collapse, mod-9 valuation, "
        "mod-13 character phase, and route-state closure fields."
    )
    print(
        "  2-vs-2: four-term equality of squares of cubes; if lifted to six "
        "slots it carries a canceling pair and must be marked "
        "degenerate_padded_support6 before entering THM-538/coimage ledgers."
    )
    print(
        "  proof target: a remaining LRC14 relation-lattice residual may use "
        "a sixth-power sidecar only after exact M, endpoint owner, topology, "
        "cycle image, and discharge route remain attached."
    )
    print()

    fp = tournament_fingerprint()
    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof carriers / sidecar fields, not runners or bases")
    print("  score_labels=", fp["score_labels"])
    print("  score_hist=", fp["score_hist"])
    print("  directed_3cycles=", fp["directed_3cycles"])
    print("  scc_sizes=", fp["scc_sizes"])
    print("  hamiltonian_path_count=", fp["hamiltonian_path_count"])
    print("  hamiltonian_path=", " > ".join(fp["hamiltonian_path"]))


if __name__ == "__main__":
    main()
