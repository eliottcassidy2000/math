#!/usr/bin/env python3
"""Exact companion for THM-3491.

The script constructs the physical height-21 staircase orbit, its length-20
overlap graph, and the deterministic integer Collatz certificate

    4 A v <= 7 v.

It contains no assertion statements, so ``python3`` and ``python3 -O`` run
the same explicit gates.  Sparse matrix-vector products use signed 64-bit
integer arithmetic; all intermediate values are independently bounded far
below that range by the committed certificate statistics.
"""

from __future__ import annotations

import hashlib
from itertools import product

try:
    import numpy as np
    from scipy.sparse import coo_matrix
except ImportError as exc:  # pragma: no cover - an explicit environment gate
    raise RuntimeError("THM-3491 requires numpy and scipy") from exc


A_STATE = 0
B_STATE = 1
C_STATE = 2

GRAPH_HEIGHT = 21
EXPECTED_COUNTS_H14_21 = (
    18619,
    32885,
    57741,
    100901,
    175680,
    304714,
    526563,
    906525,
)
EXPECTED_ITERATIONS = 5789
EXPECTED_CERT_MIN = 11
EXPECTED_CERT_MAX = 2303
EXPECTED_CERT_SUM = 253229011
EXPECTED_CERT_SHA256 = "cf7adb12a00c294013135ebd19e606c0ff553eb2f27d65f300fbfbae45dc284a"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exclusive_prefix_parity(e_mask: int, height: int) -> int:
    """Return bit i = parity of e_mask bits strictly below i."""
    mask = (1 << height) - 1
    inclusive = e_mask
    shift = 1
    while shift < height:
        inclusive ^= (inclusive << shift) & mask
        shift <<= 1
    return (inclusive << 1) & mask


def mask_step(code: int, height: int, incoming: int) -> int:
    """One exact transition on the (E,C) two-mask staircase encoding."""
    require(incoming in (0, 1), "nonbinary staircase input")
    mask = (1 << height) - 1
    e_mask = code & mask
    c_mask = code >> height
    carrier_mask = exclusive_prefix_parity(e_mask, height)
    if incoming:
        carrier_mask ^= mask
    new_c = carrier_mask
    new_e = carrier_mask | c_mask
    return new_e | (new_c << height)


def decode_state(code: int, height: int) -> tuple[int, ...]:
    mask = (1 << height) - 1
    e_mask = code & mask
    c_mask = code >> height
    return tuple(
        C_STATE if (c_mask >> i) & 1 else B_STATE if (e_mask >> i) & 1 else A_STATE
        for i in range(height)
    )


def encode_state(state: tuple[int, ...]) -> int:
    height = len(state)
    e_mask = 0
    c_mask = 0
    for i, q in enumerate(state):
        require(q in (A_STATE, B_STATE, C_STATE), "bad ternary state")
        if q != A_STATE:
            e_mask |= 1 << i
        if q == C_STATE:
            c_mask |= 1 << i
    return e_mask | (c_mask << height)


def tuple_step(state: tuple[int, ...], incoming: int) -> tuple[int, ...]:
    """Independent serial Mealy implementation of the same transition."""
    carrier = incoming
    new_state = []
    for q in state:
        emitted = carrier ^ int(q != A_STATE)
        if carrier:
            q_new = C_STATE
        elif q == C_STATE:
            q_new = B_STATE
        else:
            q_new = A_STATE
        new_state.append(q_new)
        carrier = emitted
    return tuple(new_state)


def reachable_codes(height: int) -> set[int]:
    """Exact orbit of A^height under the two external input transitions."""
    seen = {0}
    queue = [0]
    cursor = 0
    while cursor < len(queue):
        code = queue[cursor]
        cursor += 1
        for incoming in (0, 1):
            successor = mask_step(code, height, incoming)
            if successor not in seen:
                seen.add(successor)
                queue.append(successor)
    return seen


def prefix_code(code: int, source_height: int, target_height: int) -> int:
    require(0 <= target_height <= source_height, "bad prefix height")
    source_mask = (1 << source_height) - 1
    target_mask = (1 << target_height) - 1
    e_mask = code & source_mask
    c_mask = code >> source_height
    return (e_mask & target_mask) | ((c_mask & target_mask) << target_height)


def suffix_code(code: int, source_height: int) -> int:
    """Drop coordinate zero from a source_height word."""
    require(source_height >= 1, "empty suffix request")
    source_mask = (1 << source_height) - 1
    e_mask = code & source_mask
    c_mask = code >> source_height
    target_height = source_height - 1
    return (e_mask >> 1) | ((c_mask >> 1) << target_height)


def word_code(word: str) -> int:
    state = tuple("ABC".index(letter) for letter in word)
    return encode_state(state)


def ceil_log2(value: int) -> int:
    require(value >= 1, "ceil_log2 requires a positive integer")
    return (value - 1).bit_length()


def main() -> None:
    # Small independent implementation gate.
    small_checks = 0
    for height in range(1, 9):
        states = reachable_codes(height)
        for code in states:
            state = decode_state(code, height)
            require(encode_state(state) == code, "two-mask round trip failed")
            for incoming in (0, 1):
                require(
                    mask_step(code, height, incoming)
                    == encode_state(tuple_step(state, incoming)),
                    "tuple/mask staircase transition mismatch",
                )
                small_checks += 1

    # One physical orbit supplies R_21 and, by prefix projection, every R_h
    # below it.  Prefix surjectivity follows because the same driving word acts
    # on every taller cascade.
    edges = reachable_codes(GRAPH_HEIGHT)
    counts = []
    projected: dict[int, set[int]] = {}
    for height in range(14, GRAPH_HEIGHT):
        states = {
            prefix_code(code, GRAPH_HEIGHT, height)
            for code in edges
        }
        projected[height] = states
        counts.append(len(states))
    counts.append(len(edges))
    require(tuple(counts) == EXPECTED_COUNTS_H14_21, "reachable count mismatch")

    vertex_height = GRAPH_HEIGHT - 1
    vertices = sorted(projected[vertex_height])
    require(len(vertices) == 526563, "G_21 vertex count mismatch")
    index = {code: i for i, code in enumerate(vertices)}

    rows = np.empty(len(edges), dtype=np.int32)
    columns = np.empty(len(edges), dtype=np.int32)
    for edge_index, code in enumerate(edges):
        prefix = prefix_code(code, GRAPH_HEIGHT, vertex_height)
        suffix = suffix_code(code, GRAPH_HEIGHT)
        require(suffix in index, "R_21 suffix is not in R_20")
        rows[edge_index] = index[prefix]
        columns[edge_index] = index[suffix]

    adjacency = coo_matrix(
        (np.ones(len(edges), dtype=np.int64), (rows, columns)),
        shape=(len(vertices), len(vertices)),
        dtype=np.int64,
    ).tocsr()
    require(adjacency.nnz == len(edges), "duplicate G_21 overlap edge")
    require(
        int(adjacency.data.min()) == 1 and int(adjacency.data.max()) == 1,
        "G_21 adjacency is not zero-one",
    )

    # Canonical monotone construction of the positive Collatz vector.
    certificate = np.ones(len(vertices), dtype=np.int64)
    iterations = 0
    while True:
        iterations += 1
        successor_sums = adjacency.dot(certificate)
        required = (4 * successor_sums + 6) // 7
        next_certificate = np.maximum(certificate, required)
        if np.array_equal(next_certificate, certificate):
            break
        certificate = next_certificate
        require(iterations < 10000, "Collatz certificate did not stabilize")
        require(int(certificate.max()) < 1 << 40, "certificate overflow guard failed")

    slack = 4 * adjacency.dot(certificate) - 7 * certificate
    cert_min = int(certificate.min())
    cert_max = int(certificate.max())
    cert_sum = int(certificate.sum())
    require(iterations == EXPECTED_ITERATIONS, "certificate iteration mismatch")
    require(cert_min == EXPECTED_CERT_MIN, "certificate minimum mismatch")
    require(cert_max == EXPECTED_CERT_MAX, "certificate maximum mismatch")
    require(cert_sum == EXPECTED_CERT_SUM, "certificate sum mismatch")
    require(int(slack.max()) == 0, "4 A v <= 7 v certificate failed")
    require(int(slack.min()) == -6, "certificate minimum slack mismatch")

    # Canonical payload: ascending two-mask code followed by its weight, both
    # unsigned 64-bit little-endian integers.
    payload = np.empty((len(vertices), 2), dtype="<u8")
    payload[:, 0] = np.asarray(vertices, dtype=np.uint64)
    payload[:, 1] = certificate.astype(np.uint64)
    certificate_sha = hashlib.sha256(payload.tobytes(order="C")).hexdigest()
    require(certificate_sha == EXPECTED_CERT_SHA256, "certificate digest mismatch")

    # Exact path prefactor and entropy coarsening.
    prefactor_left = cert_sum * 4**20
    prefactor_right = 318 * cert_min * 7**20
    require(prefactor_left < prefactor_right, "318*(7/4)^h prefactor failed")
    require(7**16 < 2**45, "log2(7/4) < 13/16 gate failed")
    require(3**16 < 2**26, "ternary marker exponent gate failed")

    # Charged word-budget samples.  The displayed inequalities prove the
    # asymptotic choice; samples catch floor, padding, and small-regime errors.
    tariff_samples = []
    for w in (128, 256, 512, 1024, 2048, 4096):
        d = ceil_log2(w)
        ell = w - 8 * d - 8
        m = (ell - 9) // 26
        h = 16 * m
        b = (ell - 9) // 2
        require(m >= 1 and b >= 1, "tariff sample below large-word regime")
        require(13 * m + 9 + b <= ell, "state/chunk address exceeds word budget")
        require(26 * m <= ell - 9, "ternary marker budget failed")
        tariff_samples.append((w, ell, h, b, 13 * m + 9 + b))

    # Finite-only lower-bound signal.  This is deliberately not promoted to an
    # all-height inclusion or a lower entropy bound.
    block_zero = "CABC"
    block_one = "BCBC"
    lower_checks = 0
    for block_count in range(1, 6):
        height = 4 * block_count
        states = projected[height] if height >= 14 else {
            prefix_code(code, GRAPH_HEIGHT, height) for code in edges
        }
        for bits in product((0, 1), repeat=block_count):
            word = "".join(block_one if bit else block_zero for bit in bits)
            require(word_code(word) in states, "finite block-code signal failed")
            lower_checks += 1
    require(lower_checks == 62, "finite block-code check count mismatch")

    print("THM-3491 Rule 30 seven-four staircase entropy exact companion")
    print(f"tuple_mask_transition_checks={small_checks}")
    print("reachable_counts_h14_21=" + repr(tuple(counts)))
    print(f"G21_vertices={len(vertices)} G21_edges={len(edges)} adjacency=zero_one")
    print(
        "certificate=4*A*v<=7*v "
        f"iterations={iterations} min={cert_min} max={cert_max} sum={cert_sum} "
        f"slack_range=({int(slack.min())},{int(slack.max())})"
    )
    print("certificate_payload=ascending_(two_mask_code_u64le,weight_u64le)_pairs")
    print(f"certificate_sha256={certificate_sha}")
    print(
        f"prefactor_check={prefactor_left}<{prefactor_right}; "
        "N_h<=318*(7/4)^h_for_h>=20"
    )
    print(f"entropy_check=7^16={7**16}<2^45={2**45}; eta<13/16")
    print(f"marker_check=3^16={3**16}<2^26={2**26}")
    print("tariff_samples_(w,ell,h,b,address_bits)=" + repr(tuple(tariff_samples)))
    print("compiler_query_tariff=(13/4+o(1))*n^2/w^2")
    print(
        "finite_lower_signal=CABC/BCBC_full_block_code_through_m=5 "
        f"checks={lower_checks} universal_claim=0"
    )
    print("scope=upper_bound_and_compiler_only prizes_solved=0 literature_novelty_claimed=0")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
