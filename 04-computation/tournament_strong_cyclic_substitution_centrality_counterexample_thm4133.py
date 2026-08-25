#!/usr/bin/env python3
"""Exact structured counterexample audit for THM-4133.

The load-bearing response evaluator is the independently audited THM-4131
primary implementation.  This script pins that source, builds a transparent
substitution family, and prints the exact first failure inside that family.
"""

from hashlib import sha256
import importlib.util
import os


BASE_PATH = os.path.join(
    os.path.dirname(__file__),
    "tournament_strong_centrality_through_order_eight_thm4131.py",
)
EXPECTED_BASE_SHA256 = "6b195b0379d1ae3e5d215aa1c495f7180daeecae189df86269d07ef855867881"
EXPECTED_SEMANTIC = "2d41d1a1bb6f8a936c6f8104d149cba898ae51d17371981dbf2d1332643d2873"


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def file_digest(path):
    hasher = sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


require(file_digest(BASE_PATH) == EXPECTED_BASE_SHA256, "THM-4131 evaluator hash")
SPEC = importlib.util.spec_from_file_location("thm4131", BASE_PATH)
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


# Quotient blocks are (a,B,c,d,e).  The ten displayed arcs specify Q exactly.
QUOTIENT_OUT = (
    frozenset((1, 2, 4)),       # a -> B,c,e
    frozenset((2, 3, 4)),       # B -> c,d,e
    frozenset((3,)),            # c -> d
    frozenset((0, 4)),          # d -> a,e
    frozenset((2,)),            # e -> c
)


def regular_cyclic_minus_zero(order):
    require(order >= 3 and order % 2 == 1, "odd cyclic order")
    old_vertices = tuple(range(1, order))
    new_index = {vertex: index for index, vertex in enumerate(old_vertices)}
    adjacency = [0] * (order - 1)
    for vertex in old_vertices:
        for step in range(1, (order + 1) // 2):
            target = (vertex + step) % order
            if target in new_index:
                adjacency[new_index[vertex]] |= 1 << new_index[target]
    return tuple(adjacency)


def substitute_block(block):
    sizes = (1, len(block), 1, 1, 1)
    starts = []
    cursor = 0
    for size in sizes:
        starts.append(cursor)
        cursor += size
    adjacency = [0] * cursor
    for left, row in enumerate(block):
        for right in range(len(block)):
            if row & (1 << right):
                adjacency[1 + left] |= 1 << (1 + right)
    for source_block, targets in enumerate(QUOTIENT_OUT):
        for target_block in targets:
            for source in range(starts[source_block], starts[source_block] + sizes[source_block]):
                for target in range(starts[target_block], starts[target_block] + sizes[target_block]):
                    adjacency[source] |= 1 << target
    return tuple(adjacency)


def compact_family_row(cyclic_order):
    adjacency = substitute_block(regular_cyclic_minus_zero(cyclic_order))
    require(BASE.is_strong(adjacency), f"strong substitution at R_{cyclic_order}")
    record = BASE.analyze(adjacency)
    central = {0} if len(adjacency) % 2 == 0 else {-1, 1}
    best_central = max(
        layer["coset_floor"] for layer in record["layers"]
        if layer["t"] in central
    )
    best_outer = max(
        layer["coset_floor"] for layer in record["layers"]
        if layer["t"] not in central
    )
    return {
        "cyclic_order": cyclic_order,
        "order": len(adjacency),
        "adjacency": adjacency,
        "packet": (record["H"], BASE.pair(record["W"]), BASE.pair(record["D4"]),
                   BASE.pair(record["C_hd"])),
        "rho": BASE.pair(record["normalized_tilt"]),
        "rational_t": record["rational_t"],
        "coset_t": record["coset_t"],
        "actual_t": record["actual_t"],
        "coset_central_minus_outer": best_central - best_outer,
        "layers": tuple(
            (
                layer["size"], layer["t"], BASE.pair(layer["floor"]),
                layer["coset_floor"], layer["maximum"], layer["lattice"],
            )
            for layer in record["layers"]
        ),
    }


def main():
    family = tuple(compact_family_row(order) for order in (3, 5, 7, 9))
    expected_rhos = (
        (13, 557),
        (14809, 29733),
        (97613103, 107670631),
        (53092739331, 40435524866),
    )
    require(tuple(row["rho"] for row in family) == expected_rhos, "family rho row")
    require(all(row["rational_t"] == (0,) and row["coset_t"] == (0,)
                for row in family[:-1]), "central proper-prefix family")

    witness = family[-1]
    require(witness["order"] == 12, "witness order")
    require(witness["adjacency"] ==
            (3070, 3644, 3704, 3824, 4064, 4032, 3970, 3846, 3598, 1024, 2049, 512),
            "frozen witness adjacency")
    require(witness["packet"] ==
            (27759, (506085, 1), (80871049732, 1), (-23596773036, 1)),
            "witness exposure packet")
    require(witness["rational_t"] == (-2,), "noncentral rational optimizer")
    require(witness["coset_t"] == (-2,), "noncentral coset optimizer")
    require(witness["actual_t"] == (-6,), "actual optimizer sidecar")
    require(witness["coset_central_minus_outer"] == -2224,
            "strict noncentral coset margin")
    require(witness["rho"][0] > witness["rho"][1], "rho exceeds one")

    ledger = {
        "theorem": "THM-4133",
        "quotient_out": tuple(tuple(sorted(row)) for row in QUOTIENT_OUT),
        "family": family,
        "scope": (
            "explicit strong order-12 counterexample to all-order rational and exact-coset "
            "Johnson centrality; first failure only within this cyclic-deletion substitution family"
        ),
    }
    semantic = BASE.digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")
    print("status=PASS")
    print("theorem=THM-4133 structured strong centrality counterexample")
    print(f"quotient_out={ledger['quotient_out']}")
    print("family=R_(3,5,7,9)-zero substituted into block B")
    print(f"family_rows={tuple((row['order'], row['rho'], row['rational_t'], row['coset_t']) for row in family)}")
    print(f"witness_adjacency={witness['adjacency']}")
    print(f"witness_packet={witness['packet']}")
    print(f"witness_rho={witness['rho']}")
    print(f"witness_optimizers=rational:{witness['rational_t']};coset:{witness['coset_t']};actual:{witness['actual_t']}")
    print(f"witness_coset_central_minus_outer={witness['coset_central_minus_outer']}")
    print(f"witness_layers={witness['layers']}")
    print("scope=explicit n=12 refutation; no claim of globally minimal order or family classification")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
