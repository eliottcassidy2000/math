#!/usr/bin/env python3
"""Exact wreath-centralizer and LCA-packet audit for THM-3539.

The script checks only finite group/combinatorial claims.  It does not
compute the Keller decomposition group, and it does not test any all-level
coordinate resultant.  Its roles are:

* verify the order formula identifying the centralizer of a supported bottom
  transposition with the stabilizer of its third sibling;
* generate the full predecessor-leaf stabilizer and independently recover its
  point and unordered-pair orbits through predecessor depth five;
* check every symbolic LCA orbit-size formula and the exact ``n^2`` packet
  count; and
* retain the trivial-image hostile showing that the proved inertia subgroup
  alone gives no exponential reduction.
"""

from __future__ import annotations

from collections import defaultdict, deque
from hashlib import sha256
from itertools import combinations, product
import json


EXPECTED_SEMANTIC_SHA256 = "96b5ba9f5c0c66e29c5b40d20f77b11242de706a7c8b98c66647cddc283efbea"


Word = tuple[int, ...]
Permutation = tuple[int, ...]
Pair = tuple[int, int]


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def words(depth: int) -> list[Word]:
    return list(product(range(3), repeat=depth))


def lca_depth(left: Word, right: Word) -> int:
    depth = 0
    while depth < len(left) and left[depth] == right[depth]:
        depth += 1
    return depth


def wreath_order(depth: int) -> int:
    return 6 ** ((3**depth - 1) // 2)


def centralizer_row(n: int) -> dict[str, int | str]:
    """Return the symbolic order ledger for a bottom transposition in W_n."""

    require(n >= 1, ("bad depth", n))
    predecessor_depth = n - 1
    block_count = 3**predecessor_depth
    total_order = wreath_order(n)
    predecessor_order = wreath_order(predecessor_depth)
    predecessor_leaf_stabilizer_order = predecessor_order // block_count
    centralizer_order = (
        2 * 6 ** (block_count - 1) * predecessor_leaf_stabilizer_order
    )
    leaf_stabilizer_order = total_order // 3**n

    require(
        centralizer_order == leaf_stabilizer_order,
        ("centralizer/leaf-stabilizer mismatch", n),
    )
    require(
        total_order // centralizer_order == 3**n,
        ("bottom-transposition conjugacy census", n),
    )
    require(
        3 * block_count == 3**n,
        ("three transpositions per bottom block", n),
    )
    centralizer_six_exponent = (3**n - 3) // 2
    return {
        "n": n,
        "predecessor_depth": predecessor_depth,
        "predecessor_blocks": block_count,
        "wreath_order": f"6^{(3**n - 1) // 2}",
        "centralizer_order": (
            f"2^{centralizer_six_exponent + 1}*"
            f"3^{centralizer_six_exponent - (n - 1)}"
        ),
        "centralizer_index": 3**n,
        "internal_orbits_full_gate": n,
        "pair_orbits_full_gate": n * (n - 1),
        "total_packets_full_gate": n * n,
        "raw_factors_trivial_image": block_count * (block_count + 1) // 2,
    }


def predecessor_stabilizer_generators(depth: int) -> tuple[list[Word], list[Permutation]]:
    """Generate the full W_depth stabilizer of the marked leaf 0^depth.

    At a vertex on the marked ray the local S_3 is reduced to the S_2 fixing
    the next zero digit.  Every off-ray vertex retains its full local S_3.
    """

    leaves = words(depth)
    positions = {leaf: index for index, leaf in enumerate(leaves)}
    generators: list[Permutation] = []
    for level in range(depth):
        marked_prefix = (0,) * level
        for prefix in product(range(3), repeat=level):
            if prefix == marked_prefix:
                local_generators = ((0, 2, 1),)
            else:
                local_generators = ((1, 0, 2), (1, 2, 0))
            for local in local_generators:
                permutation: list[int] = []
                for leaf in leaves:
                    image = list(leaf)
                    if leaf[:level] == prefix:
                        image[level] = local[leaf[level]]
                    permutation.append(positions[tuple(image)])
                require(
                    sorted(permutation) == list(range(3**depth)),
                    ("bad stabilizer generator", depth, level, prefix),
                )
                require(permutation[0] == 0, ("marked leaf moved", depth, prefix))
                generators.append(tuple(permutation))
    require(
        len(generators) == 3**depth - 1 - depth,
        ("stabilizer generator census", depth, len(generators)),
    )
    return leaves, generators


def generated_orbits(
    items: list[int] | list[Pair],
    generators: list[Permutation],
    action,
) -> list[frozenset[int] | frozenset[Pair]]:
    unseen = set(items)
    answer: list[frozenset[int] | frozenset[Pair]] = []
    while unseen:
        seed = min(unseen)
        seen = {seed}
        queue = deque([seed])
        while queue:
            item = queue.popleft()
            for generator in generators:
                image = action(generator, item)
                if image not in seen:
                    seen.add(image)
                    queue.append(image)
        unseen.difference_update(seen)
        answer.append(frozenset(seen))
    return answer


def block_signature(leaves: list[Word], position: int) -> tuple[str, int]:
    if position == 0:
        return ("marked", len(leaves[0]))
    return ("shell", lca_depth(leaves[0], leaves[position]))


def pair_signature(
    leaves: list[Word], pair: Pair
) -> tuple[str, int] | tuple[str, int, int] | tuple[str, int, int, int]:
    left, right = pair
    if left == 0:
        return ("marked_pair", lca_depth(leaves[0], leaves[right]))
    if right == 0:
        return ("marked_pair", lca_depth(leaves[0], leaves[left]))
    first = lca_depth(leaves[0], leaves[left])
    second = lca_depth(leaves[0], leaves[right])
    if first > second:
        first, second = second, first
    mutual = lca_depth(leaves[left], leaves[right])
    if first < second:
        require(mutual == first, ("ultrametric asymmetric failure", pair))
        return ("asymmetric", first, second)
    require(mutual >= first, ("ultrametric equal-shell failure", pair))
    return ("equal_shell", first, mutual)


def expected_block_orbit_size(depth: int, signature: tuple[str, int]) -> int:
    kind, value = signature
    if kind == "marked":
        require(value == depth, ("marked signature", depth, signature))
        return 1
    require(kind == "shell" and 0 <= value < depth, ("block signature", signature))
    return 2 * 3 ** (depth - value - 1)


def expected_pair_orbit_size(depth: int, signature: tuple) -> int:
    kind = signature[0]
    if kind == "marked_pair":
        level = signature[1]
        return 2 * 3 ** (depth - level - 1)
    if kind == "asymmetric":
        first, second = signature[1:]
        require(0 <= first < second < depth, ("asymmetric signature", signature))
        return 4 * 3 ** (2 * depth - first - second - 2)
    require(kind == "equal_shell", ("pair signature kind", signature))
    shell, mutual = signature[1:]
    require(0 <= shell <= mutual < depth, ("equal-shell signature", signature))
    if shell == mutual:
        return 3 ** (2 * (depth - shell - 1))
    return 2 * 3 ** (2 * depth - mutual - shell - 2)


def serialize_signature(signature: tuple) -> str:
    return ":".join(map(str, signature))


def orbit_audit(depth: int) -> dict[str, object]:
    leaves, generators = predecessor_stabilizer_generators(depth)
    degree = len(leaves)
    block_orbits = generated_orbits(
        list(range(degree)), generators, lambda generator, item: generator[item]
    )
    pair_items = list(combinations(range(degree), 2))
    pair_orbits = generated_orbits(
        pair_items,
        generators,
        lambda generator, item: tuple(
            sorted((generator[item[0]], generator[item[1]]))
        ),
    )

    block_ledger: dict[tuple, int] = {}
    for orbit in block_orbits:
        signatures = {block_signature(leaves, item) for item in orbit}
        require(len(signatures) == 1, ("mixed block orbit", depth, signatures))
        signature = next(iter(signatures))
        require(signature not in block_ledger, ("split block signature", depth, signature))
        expected = expected_block_orbit_size(depth, signature)
        require(len(orbit) == expected, ("block orbit size", depth, signature, len(orbit)))
        block_ledger[signature] = len(orbit)

    pair_ledger: dict[tuple, int] = {}
    for orbit in pair_orbits:
        signatures = {pair_signature(leaves, item) for item in orbit}
        require(len(signatures) == 1, ("mixed pair orbit", depth, signatures))
        signature = next(iter(signatures))
        require(signature not in pair_ledger, ("split pair signature", depth, signature))
        expected = expected_pair_orbit_size(depth, signature)
        require(len(orbit) == expected, ("pair orbit size", depth, signature, len(orbit)))
        pair_ledger[signature] = len(orbit)

    require(len(block_ledger) == depth + 1, ("block orbit count", depth))
    require(len(pair_ledger) == depth * (depth + 1), ("pair orbit count", depth))
    require(sum(block_ledger.values()) == degree, ("block partition", depth))
    require(
        sum(pair_ledger.values()) == degree * (degree - 1) // 2,
        ("pair partition", depth),
    )
    require(
        len(block_ledger) + len(pair_ledger) == (depth + 1) ** 2,
        ("quadratic packet census", depth),
    )
    return {
        "predecessor_depth": depth,
        "degree": degree,
        "generator_count": len(generators),
        "block_orbits": {
            serialize_signature(key): block_ledger[key]
            for key in sorted(block_ledger, key=serialize_signature)
        },
        "pair_orbits": {
            serialize_signature(key): pair_ledger[key]
            for key in sorted(pair_ledger, key=serialize_signature)
        },
        "total_packets": (depth + 1) ** 2,
        "raw_factors": degree * (degree + 1) // 2,
    }


def symbolic_census(max_n: int = 12) -> list[dict[str, int]]:
    rows: list[dict[str, int]] = []
    for n in range(1, max_n + 1):
        depth = n - 1
        blocks = 3**depth
        internal = n
        pairs = n * (n - 1)
        packets = n * n
        raw = blocks * (blocks + 1) // 2
        require(internal == depth + 1, ("internal formula", n))
        require(pairs == depth * (depth + 1), ("pair formula", n))
        require(packets == internal + pairs, ("packet formula", n))
        require(packets <= raw, ("packet floor exceeds raw census", n))
        rows.append(
            {
                "n": n,
                "blocks": blocks,
                "internal_orbits": internal,
                "pair_orbits": pairs,
                "full_gate_packets": packets,
                "trivial_image_packets": raw,
            }
        )
    return rows


centralizer_ledger = [centralizer_row(n) for n in range(1, 8)]
orbit_ledgers = [orbit_audit(depth) for depth in range(0, 6)]
census = symbolic_census()

semantic_payload = {
    "theorem": "THM-3539",
    "group": "W_n=S3 wr ... wr S3 on ternary words",
    "supported_inertia": "I=<t> with t one bottom-block transposition",
    "proved_decomposition_bound": "I <= D <= N_W(I)=C_W(t)=Stab_W(third sibling)",
    "unproved_gate": "image(D->W_(n-1)) has full point-and-pair stabilizer orbits",
    "full_gate_counts": {
        "internal": "n",
        "cross_block": "n(n-1)",
        "total": "n^2",
    },
    "unconditional_orbit_floor": "at least n^2 D-orbits; at most 3^(n-1)(3^(n-1)+1)/2",
    "centralizer_ledger": centralizer_ledger,
    "orbit_ledgers": orbit_ledgers,
    "symbolic_census": census,
}
semantic_text = json.dumps(semantic_payload, sort_keys=True, separators=(",", ":"))
semantic_hash = sha256(semantic_text.encode("utf-8")).hexdigest()
require(
    semantic_hash == EXPECTED_SEMANTIC_SHA256,
    ("semantic hash", semantic_hash, EXPECTED_SEMANTIC_SHA256),
)

print("status=VERIFIED_EXACT_GROUP_THEORY_CONDITIONAL_KELLER_PACKET_REDUCTION")
print("scope=centralizer+normalizer+LCA_orbits; decomposition_saturation_and_all_level_index_zero_OPEN")
print("decomposition_bound=I<=D<=N_W(I)=C_W(t)=Stab_W(third_sibling)")
print("full_gate_packets=internal:n cross:n*(n-1) total:n^2")
print("unconditional_floor=n^2; raw_upper=3^(n-1)*(3^(n-1)+1)/2")
print("centralizer_ledger=" + json.dumps(centralizer_ledger, sort_keys=True, separators=(",", ":")))
print("orbit_ledgers=" + json.dumps(orbit_ledgers, sort_keys=True, separators=(",", ":")))
print("symbolic_census=" + json.dumps(census, sort_keys=True, separators=(",", ":")))
print("semantic_sha256=" + semantic_hash)
