#!/usr/bin/env python3
"""Exact companion for THM-3480.

The script uses no assertion statements, so ``python3`` and ``python3 -O``
execute the same gates.  It constructs the physical zero-boundary staircase
transducer, its reachable languages R_h, the exact k=13 SFT certificate, and
finite Myhill--Nerode and chunk-permutivity hostiles.
"""

from __future__ import annotations

import hashlib
from itertools import product


A = 0  # elementary right context 00
B = 1  # elementary right context 01
C = 2  # elementary right context 1*

EXPECTED_COUNTS = (3, 7, 16, 35, 71, 141, 272, 517, 971, 1792, 3263, 5873, 10483)
EXPECTED_R7_SHA256 = "97539561bb1476d314fc5f9db9ea179c5333e96d04d29122b93cf297854395c2"
EXPECTED_R12_SHA256 = "6d0de0cfed05a12fb0f278fbfe292fe75433d854a3b1d8ee13fa41392b8fe8e7"
EXPECTED_R13_SHA256 = "f10c8e76404828e265bcb66f34b06edaa0cb8cba7a615bf6a2b38018a0bad770"
EXPECTED_CERT_SHA256 = "52d2b23844f52fb6af9666de6034440887aafeb67f07ec2ffc9a95202a576a94"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def elementary_step(q: int, a: int) -> tuple[int, int]:
    """Consume a bit on the left; return (new state, emitted right output)."""
    require(q in (A, B, C) and a in (0, 1), "bad elementary input")
    emitted = a ^ (q != A)
    if a == 1:
        new_q = C
    elif q == C:
        new_q = B
    else:
        new_q = A
    return new_q, int(emitted)


def cascade_step(state: tuple[int, ...], a: int) -> tuple[tuple[int, ...], int]:
    """One external transition through a serial cascade of Rule-30 stages."""
    carrier = a
    new_state = []
    for q in state:
        q_new, carrier = elementary_step(q, carrier)
        new_state.append(q_new)
    return tuple(new_state), carrier


def mask_step(state: tuple[int, ...], a: int) -> tuple[tuple[int, ...], int]:
    """Independent check of the (E,C)-mask staircase recurrence."""
    e_mask = 0
    c_mask = 0
    for i, q in enumerate(state):
        if q != A:
            e_mask |= 1 << i
        if q == C:
            c_mask |= 1 << i

    carrier = a
    e_new = 0
    c_new = 0
    for i in range(len(state)):
        bit = 1 << i
        if carrier:
            e_new |= bit
            c_new |= bit
        elif c_mask & bit:
            e_new |= bit
        if e_mask & bit:
            carrier ^= 1

    decoded = []
    for i in range(len(state)):
        bit = 1 << i
        if c_new & bit:
            decoded.append(C)
        elif e_new & bit:
            decoded.append(B)
        else:
            decoded.append(A)
    return tuple(decoded), carrier


def direct_center(context: tuple[int, ...], h: int) -> int:
    require(len(context) == 2 * h + 1, "context has wrong width")
    row = list(context)
    for _ in range(h):
        row = [row[j] ^ (row[j + 1] | row[j + 2]) for j in range(len(row) - 2)]
    require(len(row) == 1, "direct cone did not close")
    return row[0]


def reachable_states(h: int) -> tuple[tuple[int, ...], ...]:
    """Orbit R_h of A^h under the two external bit transitions."""
    start = (A,) * h
    seen = {start}
    queue = [start]
    cursor = 0
    while cursor < len(queue):
        state = queue[cursor]
        cursor += 1
        for a in (0, 1):
            successor, _ = cascade_step(state, a)
            if successor not in seen:
                seen.add(successor)
                queue.append(successor)
    return tuple(sorted(seen))


def language_sha(states: tuple[tuple[int, ...], ...]) -> str:
    payload = "\n".join("".join(map(str, state)) for state in states) + "\n"
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def nerode_history(states: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    """Exact partition refinement for the reachable binary Mealy machine."""
    index = {state: i for i, state in enumerate(states)}
    transitions = []
    outputs = []
    for state in states:
        s0, z0 = cascade_step(state, 0)
        s1, z1 = cascade_step(state, 1)
        require(s0 in index and s1 in index, "reachable set is not transition closed")
        transitions.append((index[s0], index[s1]))
        outputs.append((z0, z1))

    signatures: dict[tuple[int, int], int] = {}
    classes = []
    for signature in outputs:
        if signature not in signatures:
            signatures[signature] = len(signatures)
        classes.append(signatures[signature])
    history = [len(signatures)]

    while True:
        signatures2: dict[tuple[tuple[int, int], int, int], int] = {}
        classes2 = []
        for i in range(len(states)):
            signature2 = (
                outputs[i],
                classes[transitions[i][0]],
                classes[transitions[i][1]],
            )
            if signature2 not in signatures2:
                signatures2[signature2] = len(signatures2)
            classes2.append(signatures2[signature2])
        history.append(len(signatures2))
        if len(signatures2) == history[-2]:
            return tuple(history)
        classes = classes2


def chunk_action(state: tuple[int, ...], word: int, b: int) -> tuple[tuple[int, ...], int]:
    """Consume b bits in scan order (least significant bit first)."""
    output_word = 0
    current = state
    for i in range(b):
        current, emitted = cascade_step(current, (word >> i) & 1)
        output_word |= emitted << i
    return current, output_word


def ceil_log2(value: int) -> int:
    require(value >= 1, "ceil_log2 requires a positive integer")
    return (value - 1).bit_length()


def main() -> None:
    expected_elementary = (
        ((A, 0), (C, 1)),
        ((A, 1), (C, 0)),
        ((B, 1), (C, 0)),
    )
    actual_elementary = tuple(
        tuple(elementary_step(q, a) for a in (0, 1)) for q in (A, B, C)
    )
    require(actual_elementary == expected_elementary, "elementary Mealy table mismatch")

    direct_contexts = 0
    for h in range(1, 7):
        width = 2 * h + 1
        for encoded in range(1 << width):
            context = tuple((encoded >> j) & 1 for j in range(width))
            state = (A,) * h
            emitted = 0
            for a in reversed(context):
                state, emitted = cascade_step(state, a)
            require(emitted == direct_center(context, h), "cascade/direct cone mismatch")
            direct_contexts += 1

    reachable: dict[int, tuple[tuple[int, ...], ...]] = {}
    counts = []
    mask_checks = 0
    for h in range(1, 14):
        states = reachable_states(h)
        reachable[h] = states
        counts.append(len(states))
        for state in states:
            for a in (0, 1):
                require(cascade_step(state, a) == mask_step(state, a), "mask recurrence mismatch")
                mask_checks += 1
    require(tuple(counts) == EXPECTED_COUNTS, "reachable-state count mismatch")
    require(language_sha(reachable[7]) == EXPECTED_R7_SHA256, "R7 digest mismatch")
    require(language_sha(reachable[12]) == EXPECTED_R12_SHA256, "R12 digest mismatch")
    require(language_sha(reachable[13]) == EXPECTED_R13_SHA256, "R13 digest mismatch")

    factor_checks = 0
    allowed_factors = {k: set(reachable[k]) for k in range(1, 14)}
    for h in range(2, 14):
        for state in reachable[h]:
            for k in range(1, h + 1):
                for start in range(h - k + 1):
                    require(
                        state[start : start + k] in allowed_factors[k],
                        "factorial-language gate failed",
                    )
                    factor_checks += 1

    nerode = []
    for h in range(1, 14):
        history = nerode_history(reachable[h])
        require(history[-1] == len(reachable[h]), "reachable states admit an unexpected quotient")
        nerode.append(history)

    vertices = reachable[12]
    vertex_index = {state: i for i, state in enumerate(vertices)}
    adjacency: list[list[int]] = [[] for _ in vertices]
    for word in reachable[13]:
        prefix = word[:-1]
        suffix = word[1:]
        require(prefix in vertex_index and suffix in vertex_index, "R13 edge leaves R12")
        adjacency[vertex_index[prefix]].append(vertex_index[suffix])
    require(sum(map(len, adjacency)) == 10483, "R13 edge count mismatch")

    # A stored, reproducible Collatz certificate: the deterministic monotone
    # construction is committed by EXPECTED_CERT_SHA256 and then checked with
    # exact integer arithmetic.  There is no floating-point eigenvalue gate.
    certificate = [1] * len(vertices)
    certificate_iterations = 0
    while True:
        certificate_iterations += 1
        next_certificate = []
        for i, successors in enumerate(adjacency):
            successor_sum = sum(certificate[j] for j in successors)
            required_weight = (64 * successor_sum + 116) // 117
            next_certificate.append(max(certificate[i], required_weight))
        if next_certificate == certificate:
            break
        certificate = next_certificate
        require(certificate_iterations < 5000, "certificate construction did not stabilize")

    certificate_payload = "\n".join(
        "".join(map(str, state)) + ":" + str(certificate[i])
        for i, state in enumerate(vertices)
    ) + "\n"
    certificate_sha = hashlib.sha256(certificate_payload.encode("ascii")).hexdigest()
    require(certificate_sha == EXPECTED_CERT_SHA256, "certificate digest mismatch")

    maximum_slack = None
    for i, successors in enumerate(adjacency):
        slack = 64 * sum(certificate[j] for j in successors) - 117 * certificate[i]
        maximum_slack = slack if maximum_slack is None else max(maximum_slack, slack)
        require(slack <= 0, "64 A v <= 117 v certificate failed")
    require(min(certificate) == 15, "certificate minimum mismatch")
    require(max(certificate) == 457, "certificate maximum mismatch")
    require(sum(certificate) == 854695, "certificate sum mismatch")
    require(maximum_slack == 0, "certificate boundary slack mismatch")

    # Exact prefactor and dyadic simplification used in the theorem.
    require(
        sum(certificate) * (64**12) <= 41 * min(certificate) * (117**12),
        "41*(117/64)^h prefix bound failed",
    )
    require(117**8 < 2**55, "117/64 is not below 2^(7/8)")
    for h in range(12, 14):
        require(len(reachable[h]) <= 41 * (117**h) // (64**h) + 1, "finite entropy bound failed")

    chunk_checks = 0
    for h in range(1, 7):
        for state in reachable[h]:
            for b in range(1, 6):
                outputs = set()
                for word in range(1 << b):
                    _, output_word = chunk_action(state, word, b)
                    outputs.add(output_word)
                    chunk_checks += 1
                require(len(outputs) == (1 << b), "fixed-state chunk map is not a permutation")

    tariff_samples = []
    for w in (128, 256, 512, 1024, 2048, 4096):
        d = ceil_log2(w)
        ell = w - 8 * d - 8
        m = (ell - 6) // 14
        b = (ell - 6) // 2
        h = 8 * m
        require(m >= 2 and b >= 1, "tariff sample below the large-word regime")
        require(7 * m + b + 6 <= ell, "chunk address exceeds its table budget")
        require(3 ** (8 * m) < 2 ** (13 * m), "ternary marker-array bound failed")
        tariff_samples.append((w, ell, h, b, 7 * m + b + 6))

    print("THM-3480 Rule 30 staircase-transducer exact companion")
    print("elementary_table=" + repr(actual_elementary))
    print(f"direct_cone_contexts_h1_6={direct_contexts}")
    print("reachable_counts_h1_13=" + repr(tuple(counts)))
    print(f"mask_transition_checks={mask_checks}")
    print(f"factor_language_checks={factor_checks}")
    print("nerode_histories_h1_13=" + repr(tuple(nerode)))
    print(
        "R12_vertices=5873 R13_edges=10483 "
        f"certificate_iterations={certificate_iterations} certificate_min=15 "
        "certificate_max=457 certificate_sum=854695 max_slack=0"
    )
    print(f"certificate_sha256={certificate_sha}")
    print("exact_entropy_certificate=64*A*v<=117*v; eta<=log2(117/64)<7/8")
    print(f"fixed_state_chunk_permutation_checks={chunk_checks}")
    print("tariff_samples_(w,ell,h,b,address_bits)=" + repr(tuple(tariff_samples)))
    print("compiler_query_tariff=(7/2+o(1))*n^2/w^2; prizes_solved=0")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
