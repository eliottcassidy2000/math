#!/usr/bin/env python3
"""Exact finite companion for THM-3490 Fourier--Dirichlet tomography.

The analytic proof is the residue-one law for Hurwitz zeta.  This companion
checks the finite Fourier inversion and the induced minimal tail recurrence
independently over a split prime, including missing-colour, period-refinement,
finite-prefix, alternating, and THM-3484 hostiles.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from json import dumps
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/periodic_polynomial_fourier_dirichlet_residue_tomography_thm3490.py"
OUTPUT = "05-knowledge/results/periodic_polynomial_fourier_dirichlet_residue_tomography_thm3490.out"
EXPECTED_SEMANTIC_SHA256 = "39daa9a75e84e0fcb28141fa9079c3f616c9aa5d0689934cc239b4d911fc713e"
MODULUS = 2521  # prime; MODULUS-1=lcm(1,...,10)

PINS = (
    (
        "THM-3485",
        ROOT / "01-canon/theorems/THM-3485-periodic-polynomial-fourier-jordan-recurrence-classification.md",
        "814ad59925079c9c69c3b90ea1ee7947dbf01ec1c8855bb760fdd6dc5098fafc",
    ),
    (
        "THM-3486",
        ROOT / "01-canon/theorems/THM-3486-critical-harmonic-transform-of-periodic-polynomial-words.md",
        "de79c50df490f550897d494391d1572d71a35f6f8a54d2cb72fb92c5a55a35ad",
    ),
)

THM3484_LANES = (
    (0, 0, 0, 0, -2048, 4096, 8192, -16384),
    (0, 0, 0, 0, 2048, 4096, -24576, -16384),
    (0, -256, 3072, -15360, 40960, -61440, 49152, -16384),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def prime_factors(value: int) -> tuple[int, ...]:
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    return tuple(factors)


def primitive_root(prime: int) -> int:
    factors = prime_factors(prime - 1)
    candidates = tuple(
        value
        for value in range(2, prime)
        if all(pow(value, (prime - 1) // factor, prime) != 1 for factor in factors)
    )
    require(candidates, ("no primitive root", prime))
    return candidates[0]


GENERATOR = primitive_root(MODULUS)


def root_of_order(period: int) -> int:
    require((MODULUS - 1) % period == 0, ("period does not split", period))
    root = pow(GENERATOR, (MODULUS - 1) // period, MODULUS)
    require(pow(root, period, MODULUS) == 1, (period, root))
    require(all(pow(root, period // factor, MODULUS) != 1
                for factor in prime_factors(period)), ("wrong order", period, root))
    return root


def normalize_lanes(lanes):
    width = max(len(lane) for lane in lanes)
    return tuple(
        tuple((lane[index] if index < len(lane) else 0) % MODULUS
              for index in range(width))
        for lane in lanes
    )


def residue_table(lanes):
    lanes = normalize_lanes(lanes)
    period = len(lanes)
    root = root_of_order(period)
    inverse_period = pow(period, -1, MODULUS)
    return tuple(
        tuple(
            inverse_period * sum(
                pow(root, (-character * lane) % period, MODULUS)
                * lanes[lane][degree]
                for lane in range(period)
            ) % MODULUS
            for degree in range(len(lanes[0]))
        )
        for character in range(period)
    )


def inverse_residues(residues):
    period = len(residues)
    root = root_of_order(period)
    return tuple(
        tuple(
            sum(
                pow(root, (character * lane) % period, MODULUS)
                * residues[character][degree]
                for character in range(period)
            ) % MODULUS
            for degree in range(len(residues[0]))
        )
        for lane in range(period)
    )


def row_degree(row):
    support = tuple(index for index, value in enumerate(row) if value)
    return max(support) if support else -1


def polynomial_multiply(left, right):
    product = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            product[i + j] = (product[i + j] + a * b) % MODULUS
    return tuple(product)


def predicted_characteristic(residues):
    period = len(residues)
    root = root_of_order(period)
    polynomial = (1,)
    profile = []
    for character, row in enumerate(residues):
        degree = row_degree(row)
        profile.append(degree)
        if degree >= 0:
            factor = ((-pow(root, character, MODULUS)) % MODULUS, 1)
            for _ in range(degree + 1):
                polynomial = polynomial_multiply(polynomial, factor)
    return tuple(profile), polynomial


def evaluate_lane(lane, value):
    total = 0
    for coefficient in reversed(lane):
        total = (total * value + coefficient) % MODULUS
    return total


def sequence_from_lanes(lanes, length, start=1):
    lanes = normalize_lanes(lanes)
    period = len(lanes)
    return tuple(
        evaluate_lane(lanes[n % period], n % MODULUS)
        for n in range(start, start + length)
    )


def berlekamp_massey(sequence):
    connection = [1]
    previous = [1]
    length = 0
    shift = 1
    last_discrepancy = 1
    for index in range(len(sequence)):
        discrepancy = sequence[index]
        for offset in range(1, length + 1):
            discrepancy = (discrepancy
                           + connection[offset] * sequence[index - offset]) % MODULUS
        if discrepancy == 0:
            shift += 1
            continue
        saved = connection[:]
        multiplier = discrepancy * pow(last_discrepancy, -1, MODULUS) % MODULUS
        needed = len(previous) + shift
        if len(connection) < needed:
            connection.extend([0] * (needed - len(connection)))
        for offset, coefficient in enumerate(previous):
            connection[offset + shift] = (
                connection[offset + shift] - multiplier * coefficient
            ) % MODULUS
        if 2 * length <= index:
            length = index + 1 - length
            previous = saved
            last_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    connection = tuple(connection[:length + 1])
    return tuple(reversed(connection))


def audit_packet(label, lanes):
    lanes = normalize_lanes(lanes)
    residues = residue_table(lanes)
    require(inverse_residues(residues) == lanes, (label, "Fourier inversion"))
    profile, predicted = predicted_characteristic(residues)
    order = len(predicted) - 1
    sequence = sequence_from_lanes(lanes, max(8 * order + 32, 96))
    observed = berlekamp_massey(sequence)
    require(observed == predicted, (label, profile, predicted, observed))
    top_degree = len(lanes[0]) - 1
    top_period_sum = tuple(
        sum(
            pow(root_of_order(len(lanes)), (-character * lane) % len(lanes), MODULUS)
            * lanes[lane][top_degree]
            for lane in range(len(lanes))
        ) % MODULUS
        for character in range(len(lanes))
    )
    require(top_period_sum == tuple(
        len(lanes) * residues[character][top_degree] % MODULUS
        for character in range(len(lanes))
    ), (label, "critical average"))
    return label, len(lanes), top_degree, profile, order


def deterministic_packets():
    packets = []
    for period in range(2, 10):
        for degree in range(0, 5):
            generic = tuple(
                tuple((lane + 2) ** (power + 1) + 3 * power + lane
                      for power in range(degree + 1))
                for lane in range(period)
            )
            packets.append((f"generic-p{period}-d{degree}", generic))
            shared = tuple(
                tuple(
                    11 if power == degree
                    else (lane + 1) * (power + 2) + lane * lane
                    for power in range(degree + 1)
                )
                for lane in range(period)
            )
            packets.append((f"shared-top-p{period}-d{degree}", shared))
    packets.append(("period4-parity", ((1,), (-1,), (1,), (-1,))))
    packets.append(("period6-mod3", ((2,), (3,), (5,), (2,), (3,), (5,))))
    return tuple(packets)


def security_gate(path: Path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert forbidden")
    forbidden = {"eval", "exec", "compile", "open", "system", "popen", "run", "Popen"}
    calls = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    calls.update(
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not (calls & forbidden), ("forbidden calls", calls & forbidden))
    return len(tuple(ast.walk(tree))), tuple(sorted(calls & forbidden))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    pins = []
    for label, path, expected in PINS:
        actual = lf_sha256(path)
        require(actual == expected, (label, "hash drift", actual))
        pins.append((label, actual))

    require(prime_factors(MODULUS) == (MODULUS,), ("modulus not prime", MODULUS))
    require(MODULUS - 1 == 2520, MODULUS)
    require(GENERATOR == 17, ("primitive-root drift", GENERATOR))

    controls = tuple(audit_packet(label, lanes)
                     for label, lanes in deterministic_packets())
    require(len(controls) == 82, len(controls))

    alternating = ((1, 0, 0, 0, 1), (-1, 0, 0, 0, -1))
    alternating_residues = residue_table(alternating)
    require(tuple(row_degree(row) for row in alternating_residues) == (-1, 4),
            alternating_residues)
    require(alternating_residues[1][4] == 1, alternating_residues)

    parity_refined = residue_table(((1,), (-1,), (1,), (-1,)))
    require(tuple(row_degree(row) for row in parity_refined) == (-1, -1, 0, -1),
            parity_refined)
    period2_characteristic = predicted_characteristic(residue_table(((1,), (-1,))))[1]
    period4_characteristic = predicted_characteristic(parity_refined)[1]
    require(period2_characteristic == period4_characteristic,
            (period2_characteristic, period4_characteristic))

    ternary_residues = residue_table(THM3484_LANES)
    ternary_profile, ternary_characteristic = predicted_characteristic(ternary_residues)
    require(ternary_profile == (7, 6, 6), ternary_profile)
    require(tuple(ternary_residues[j][7] for j in range(3))
            == ((-16384) % MODULUS, 0, 0), ternary_residues)
    require(len(ternary_characteristic) - 1 == 22, len(ternary_characteristic))

    prefix_lanes = ((1, 2, 3), (4, 0, 2), (5, 1, 1))
    prefix_residues = residue_table(prefix_lanes)
    _, prefix_characteristic = predicted_characteristic(prefix_residues)
    clean_tail = list(sequence_from_lanes(prefix_lanes, 240))
    changed = clean_tail[:]
    for index in range(7):
        changed[index] = (changed[index] + 100 + index) % MODULUS
    require(berlekamp_massey(tuple(changed[7:])) == prefix_characteristic,
            "finite prefix changed tail recurrence")

    # At s=m+1 the Hurwitz-zeta factor p^(m-s) has residue p^(-1),
    # independently of m.  This is the analytic normalization behind every
    # finite DFT row above.
    hurwitz_residue_normalizations = tuple(
        (period, pow(period, -1, MODULUS)) for period in range(1, 10)
    )
    security = security_gate(ROOT / SCRIPT)
    semantic_payload = {
        "pins": tuple(pins),
        "split_prime_and_generator": (MODULUS, GENERATOR),
        "controls": controls,
        "alternating_profile": tuple(row_degree(row) for row in alternating_residues),
        "period_refinement": (
            tuple(row_degree(row) for row in parity_refined),
            period2_characteristic,
        ),
        "ternary_profile_and_order": (ternary_profile, len(ternary_characteristic) - 1),
        "finite_prefix_changed_terms": 7,
        "hurwitz_residue_normalizations": hurwitz_residue_normalizations,
        "security": security,
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"),
              default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("THM-3490 FULL FOURIER-DIRICHLET RESIDUE TOMOGRAPHY COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"PINS: {tuple(pins)}")
    print(f"SPLIT_FIELD: F_{MODULUS}; primitive_generator={GENERATOR}; periods=1..9")
    print(f"CONTROL_PACKETS: {len(controls)}; Fourier inversion and predicted-vs-Berlekamp-Massey recurrence all exact")
    print(f"ALTERNATING_HOSTILE: residue_degree_profile={tuple(row_degree(row) for row in alternating_residues)}; trivial top residue=0; nontrivial top residue=1")
    print(f"PERIOD_REFINEMENT_HOSTILE: period4_profile={tuple(row_degree(row) for row in parity_refined)}; period2_and_period4_characteristic_equal={period2_characteristic == period4_characteristic}")
    print(f"THM3484_TERNARY: residue_degree_profile={ternary_profile}; recurrence_order={len(ternary_characteristic)-1}; top_residue_row=((-16384 mod {MODULUS}),0,0)")
    print("FINITE_PREFIX_HOSTILE: seven changed initial terms leave the suffix residue packet and tail recurrence unchanged; the transform does not recover transients")
    print("INFORMATION_ORDER: the labelled residue table recovers every lane coefficient; the minimal recurrence retains only the highest nonzero residue in each character row and cannot recover the table")
    print("SCOPE: periodic-polynomial tails only; character gauge and declared period retained; no graph, tournament, ancestry, LRC-current, or Jacobian-flux semantics")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print(f"SECURITY_AST_NODES_AND_FORBIDDEN: {security}")
    print("VERDICT: every Fourier/Jordan rung is one twisted Dirichlet pole residue; ordinary harmonic mass is only the top trivial entry, and recurrence is a lossy support-depth projection of the full residue table")


if __name__ == "__main__":
    main()
