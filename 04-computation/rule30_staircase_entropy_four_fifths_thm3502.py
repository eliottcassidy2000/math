#!/usr/bin/env python3
"""Exact companion for provisional THM-3502.

The script reconstructs the physical height-23 Rule 30 staircase orbit, its
length-22 overlap graph, and the deterministic integer Collatz certificate

    50 A_23 v <= 87 v.

It also exhausts the height-24 orbit for the explicit spurious overlap-path
hostile and for the finite-only ``{CABC,BCBC}^6`` lower signal.  There are no
``assert`` statements, so ordinary and optimized Python execute the same
explicit gates.
"""

from __future__ import annotations

import gc
import hashlib
from itertools import product

try:
    import numpy as np
    from scipy.sparse import coo_matrix
except ImportError as exc:  # pragma: no cover - explicit environment gate
    raise RuntimeError("THM-3502 requires numpy and scipy") from exc


GRAPH_HEIGHT = 23
VERTEX_HEIGHT = 22
HOSTILE_HEIGHT = 24
NUMERATOR = 87
DENOMINATOR = 50

EXPECTED_COUNTS_H14_23 = (
    18619,
    32885,
    57741,
    100901,
    175680,
    304714,
    526563,
    906525,
    1555372,
    2660178,
)
EXPECTED_N24 = 4535965
EXPECTED_ITERATIONS = 3502
EXPECTED_CERT_MIN = 6
EXPECTED_CERT_MAX = 1313
EXPECTED_CERT_SUM = 389727182
EXPECTED_SLACK_MIN = -86
EXPECTED_CERT_SHA256 = (
    "a4ce932f067e0ff21e90ccd3b4216be4dbbba9440e657b4cebc7ff834b302e33"
)
SPURIOUS_WORD = "ABCBCBCBCBCAAAAAAAAAAAAA"
EXPECTED_SPURIOUS_CODE = 22884124670


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
    """One exact transition in THM-3480's (E,C) two-mask encoding."""
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


def serial_step(code: int, height: int, incoming: int) -> int:
    """Independent stage-by-stage Mealy transition without broadword parity."""
    require(incoming in (0, 1), "nonbinary serial input")
    mask = (1 << height) - 1
    e_mask = code & mask
    c_mask = code >> height
    new_e = 0
    new_c = 0
    carrier = incoming
    for index in range(height):
        old_c = (c_mask >> index) & 1
        if carrier:
            new_c |= 1 << index
        if carrier or old_c:
            new_e |= 1 << index
        carrier ^= (e_mask >> index) & 1
    return new_e | (new_c << height)


def reachable_codes(height: int, independent_check: bool = False) -> set[int]:
    """Exact orbit of A^height under the two external input transitions."""
    seen = {0}
    queue = [0]
    cursor = 0
    while cursor < len(queue):
        code = queue[cursor]
        cursor += 1
        for incoming in (0, 1):
            successor = mask_step(code, height, incoming)
            if independent_check:
                require(
                    successor == serial_step(code, height, incoming),
                    "serial/mask staircase transition mismatch",
                )
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
    height = len(word)
    e_mask = 0
    c_mask = 0
    for index, letter in enumerate(word):
        require(letter in "ABC", "bad staircase symbol")
        if letter != "A":
            e_mask |= 1 << index
        if letter == "C":
            c_mask |= 1 << index
    return e_mask | (c_mask << height)


def decode_word(code: int, height: int) -> str:
    mask = (1 << height) - 1
    e_mask = code & mask
    c_mask = code >> height
    return "".join(
        "C" if (c_mask >> index) & 1 else "B" if (e_mask >> index) & 1 else "A"
        for index in range(height)
    )


def ceil_log2(value: int) -> int:
    require(value >= 1, "ceil_log2 requires a positive integer")
    return (value - 1).bit_length()


def main() -> None:
    # The physical orbit and an exhaustive independent transition gate.
    edges = reachable_codes(GRAPH_HEIGHT, independent_check=True)
    projected_counts = []
    vertices_set: set[int] | None = None
    for height in range(14, GRAPH_HEIGHT):
        states = {
            prefix_code(code, GRAPH_HEIGHT, height)
            for code in edges
        }
        projected_counts.append(len(states))
        if height == VERTEX_HEIGHT:
            vertices_set = states
    projected_counts.append(len(edges))
    counts = tuple(projected_counts)
    require(counts == EXPECTED_COUNTS_H14_23, "reachable count mismatch")
    require(vertices_set is not None, "R_22 projection missing")

    # Build the zero-one G_23 overlap adjacency in ascending two-mask order.
    vertices = sorted(vertices_set)
    index = {code: position for position, code in enumerate(vertices)}
    rows = np.empty(len(edges), dtype=np.int32)
    columns = np.empty(len(edges), dtype=np.int32)
    for edge_index, code in enumerate(edges):
        prefix = prefix_code(code, GRAPH_HEIGHT, VERTEX_HEIGHT)
        suffix = suffix_code(code, GRAPH_HEIGHT)
        require(suffix in index, "R_23 suffix is not in R_22")
        rows[edge_index] = index[prefix]
        columns[edge_index] = index[suffix]
    adjacency = coo_matrix(
        (np.ones(len(edges), dtype=np.int64), (rows, columns)),
        shape=(len(vertices), len(vertices)),
        dtype=np.int64,
    ).tocsr()
    require(adjacency.nnz == len(edges), "duplicate G_23 overlap edge")
    require(
        int(adjacency.data.min()) == 1 and int(adjacency.data.max()) == 1,
        "G_23 adjacency is not zero-one",
    )

    # Canonical monotone construction of the positive Collatz vector.
    certificate = np.ones(len(vertices), dtype=np.int64)
    iterations = 0
    while True:
        iterations += 1
        successor_sums = adjacency.dot(certificate)
        required = (
            DENOMINATOR * successor_sums + NUMERATOR - 1
        ) // NUMERATOR
        next_certificate = np.maximum(certificate, required)
        if np.array_equal(next_certificate, certificate):
            break
        certificate = next_certificate
        require(iterations < 10000, "Collatz certificate did not stabilize")
        require(int(certificate.max()) < 1 << 40, "certificate overflow guard failed")

    slack = DENOMINATOR * adjacency.dot(certificate) - NUMERATOR * certificate
    cert_min = int(certificate.min())
    cert_max = int(certificate.max())
    cert_sum = int(certificate.sum())
    require(iterations == EXPECTED_ITERATIONS, "certificate iteration mismatch")
    require(cert_min == EXPECTED_CERT_MIN, "certificate minimum mismatch")
    require(cert_max == EXPECTED_CERT_MAX, "certificate maximum mismatch")
    require(cert_sum == EXPECTED_CERT_SUM, "certificate sum mismatch")
    require(int(slack.min()) == EXPECTED_SLACK_MIN, "minimum slack mismatch")
    require(int(slack.max()) == 0, "50 A v <= 87 v certificate failed")

    payload = np.empty((len(vertices), 2), dtype="<u8")
    payload[:, 0] = np.asarray(vertices, dtype=np.uint64)
    payload[:, 1] = certificate.astype(np.uint64)
    certificate_sha = hashlib.sha256(payload.tobytes(order="C")).hexdigest()
    require(certificate_sha == EXPECTED_CERT_SHA256, "certificate digest mismatch")

    # Exact path prefactor, entropy coarsening, and charged marker inequality.
    prefactor_left = cert_sum * DENOMINATOR**22
    prefactor_right = cert_min * 2**9 * NUMERATOR**22
    require(prefactor_left < prefactor_right, "2^9 path prefactor failed")
    entropy_left = NUMERATOR**5
    entropy_right = 2**4 * DENOMINATOR**5
    require(entropy_left < entropy_right, "log2(87/50) < 4/5 gate failed")
    require(3**5 < 2**8, "charged ternary marker inequality failed")

    # Floor-sensitive samples for h=5m, b=floor((L-9)/2).
    tariff_samples = []
    for word_size in (128, 256, 512, 1024, 2048, 4096):
        d = ceil_log2(word_size)
        ell = word_size - 8 * d - 8
        m = (ell - 9) // 8
        height = 5 * m
        chunk = (ell - 9) // 2
        require(m >= 1 and chunk >= 1, "tariff sample below large-word regime")
        require(4 * m + 9 + chunk <= ell, "state/chunk address exceeds word budget")
        require(8 * m <= ell - 9, "ternary marker budget failed")
        tariff_samples.append((word_size, ell, height, chunk, 4 * m + 9 + chunk))

    # Recheck the inherited finite lower signal through five macro symbols.
    lower_checks_h20 = 0
    for block_count in range(1, 6):
        height = 4 * block_count
        states = {
            prefix_code(code, GRAPH_HEIGHT, height)
            for code in edges
        }
        for bits in product((0, 1), repeat=block_count):
            word = "".join("BCBC" if bit else "CABC" for bit in bits)
            require(word_code(word) in states, "finite block-code signal failed")
            lower_checks_h20 += 1
    require(lower_checks_h20 == 62, "finite height-20 lower check count mismatch")

    # Record that both length-23 factors of the hostile are physical before
    # freeing G_23 and exhausting the height-24 physical orbit.
    spurious_code = word_code(SPURIOUS_WORD)
    require(len(SPURIOUS_WORD) == HOSTILE_HEIGHT, "spurious word has wrong height")
    require(spurious_code == EXPECTED_SPURIOUS_CODE, "spurious code mismatch")
    spurious_prefix = prefix_code(spurious_code, HOSTILE_HEIGHT, GRAPH_HEIGHT)
    spurious_suffix = suffix_code(spurious_code, HOSTILE_HEIGHT)
    require(spurious_prefix in edges, "spurious prefix is not in R_23")
    require(spurious_suffix in edges, "spurious suffix is not in R_23")

    graph_vertices = len(vertices)
    graph_edges = len(edges)
    slack_min = int(slack.min())
    slack_max = int(slack.max())
    del adjacency, certificate, columns, edges, index, payload, rows, slack
    del successor_sums, vertices, vertices_set
    gc.collect()

    states24 = reachable_codes(HOSTILE_HEIGHT)
    require(len(states24) == EXPECTED_N24, "N_24 mismatch")
    require(spurious_code not in states24, "local overlap hostile is physically reachable")

    lower_words_h24 = {
        word_code("".join("BCBC" if bit else "CABC" for bit in bits))
        for bits in product((0, 1), repeat=6)
    }
    require(len(lower_words_h24) == 64, "height-24 block-code collision")
    require(lower_words_h24 <= states24, "height-24 finite block-code signal failed")

    print("THM-3502 Rule 30 four-fifths staircase entropy exact companion")
    print("reachable_counts_h14_23=" + repr(counts))
    print(f"G23_vertices={graph_vertices} G23_edges={graph_edges} adjacency=zero_one")
    print(f"serial_mask_transition_checks={2 * graph_edges}")
    print(
        "certificate=50*A*v<=87*v "
        f"iterations={iterations} min={cert_min} max={cert_max} sum={cert_sum} "
        f"slack_range=({slack_min},{slack_max})"
    )
    print("certificate_payload=ascending_(two_mask_code_u64le,weight_u64le)_pairs")
    print(f"certificate_sha256={certificate_sha}")
    print(
        f"prefactor_check={prefactor_left}<{prefactor_right}; "
        "N_h<2^9*(87/50)^h_for_h>=22"
    )
    print(
        f"entropy_check=87^5={entropy_left}<16*50^5={entropy_right}; "
        "eta<4/5"
    )
    print(f"marker_check=3^5={3**5}<2^8={2**8}")
    print("tariff_samples_(w,ell,h,b,address_bits)=" + repr(tuple(tariff_samples)))
    print("compiler_query_tariff=(16/5+o(1))*n^2/w^2")
    print(
        f"spurious_path=N24:{len(states24)} code:{spurious_code} "
        f"word:{decode_word(spurious_code, HOSTILE_HEIGHT)} "
        f"prefix:{decode_word(spurious_prefix, GRAPH_HEIGHT)} "
        f"suffix:{decode_word(spurious_suffix, GRAPH_HEIGHT)} physical_word=0"
    )
    print(
        "finite_lower_signal={CABC,BCBC}^m_subset_S_4m_through_m=6 "
        f"checks={lower_checks_h20 + len(lower_words_h24)} universal_claim=0"
    )
    print("scope=upper_bound_and_compiler_only prizes_solved=0 literature_novelty_claimed=0")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
