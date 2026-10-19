#!/usr/bin/env python3
"""Exact radix-block transport, with source, carry and integer clocks retained.

No assertion here is a universal Collatz test. The all-parameter identities
are proved in the companion note; this file checks finite controls by literal
odd-map replay, independent affine composition, and modular readers.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def v2(n):
    if n <= 0:
        raise ValueError("positive valuation argument required")
    return (n & -n).bit_length() - 1


def route(n, k):
    if n <= 0 or not n & 1:
        raise ValueError("positive odd source required")
    letters, states = [], [n]
    for _ in range(k):
        z = 3 * n + 1
        a = v2(z)
        n = z >> a
        letters.append(a)
        states.append(n)
    return tuple(letters), states


def affine(word):
    # C is updated before the current letter enters the halving clock.
    C, A = 0, 0
    for a in word:
        C = 3 * C + (1 << A)
        A += a
    return 3 ** len(word), C, A


def order(base, modulus):
    if base % 3 == 0:
        raise ValueError("this unchanged-anchor compiler requires a 3-unit base")
    z, L = base % modulus, 1
    while z != 1:
        z = z * base % modulus
        L += 1
        if L > modulus:
            raise ValueError("order search exceeded unit-group size")
    return L


def digits(n, base):
    alphabet = "0123456789abcdefghijklmnopqrstuvwxyz"
    if not 2 <= base <= len(alphabet) or n < 0:
        raise ValueError("unsupported positional display")
    out = ""
    while n:
        n, d = divmod(n, base)
        out = alphabet[d] + out
    return out or "0"


def pump(n, C, p, B):
    require((C * (B - 1)) % p == 0, "nonintegral source translation")
    return B * n + C * (B - 1) // p


def modular_family(j, modulus):
    # Divide only after retaining the full denominator modulus.
    residue = (pow(2, 18 * j + 11, 27 * modulus) - 23) % (27 * modulus)
    require(residue % 27 == 0, "full-denominator residue lost divisibility")
    return residue // 27


def affine_modular_family(j, modulus):
    # Independent binary powering of x -> 262144*x + 223307.
    a, c, x = 262144 % modulus, 223307 % modulus, 75 % modulus
    while j:
        if j & 1:
            x = (a * x + c) % modulus
        a, c = a * a % modulus, c * (a + 1) % modulus
        j >>= 1
    return x


def run():
    report = {"status": "FINITE-EXACT controls for proved scoped identities"}
    bases = [2, 4, 5, 7, 8, 10, 11, 16]
    configurations, replays, shifted_examples = 0, 0, []
    block_cache = {}
    for n in range(1, 4096, 2):
        for k in range(1, 5):
            word, states = route(n, k)
            p, C, A = affine(word)
            require(p * n + C == (1 << A) * states[-1], "affine source identity")
            require(C % 3 != 0, "carry must be a ternary unit")
            q, r = divmod(C, p)
            for base in bases:
                key = (base, p)
                if key not in block_cache:
                    L = order(base, p)
                    block_cache[key] = (L, base ** L)
                L, B = block_cache[key]
                e, u = v2(B), B >> v2(B)
                block = r * (B - 1) // p
                require(0 < block < B, "shifted coordinate has an honest radix block")
                x = n
                for j in range(1, 4):
                    x = pump(x, C, p, B)
                    closed = (B ** j * (p * n + C) - C) // p
                    require(x == closed, "iterative/closed source mismatch")
                    got, got_states = route(x, k)
                    expected = word[:-1] + (word[-1] + j * e,)
                    require(got == expected, "changed an earlier valuation or lost clock")
                    require(got_states[-1] == u ** j * states[-1], "lost odd endpoint factor")
                    require(sum(got) == A + j * e, "integer halving clock mismatch")
                    if j == 1:
                        expected_digits = digits(n + q, base) + digits(block, base).zfill(L)
                        require(digits(x + q, base) == expected_digits, "radix carry/shift mismatch")
                    replays += 1
                if q and len(shifted_examples) < 4 and base == 2:
                    shifted_examples.append({"source": n, "word": word, "offset": q,
                                             "base": base, "block_length": L, "block": block})
                configurations += 1
    report["general_controls"] = {
        "odd_sources": "1..4095", "word_lengths": [1, 2, 3, 4],
        "bases": bases, "repeat_counts": [1, 2, 3],
        "configurations": configurations, "literal_odd_map_replays": replays,
        "shifted_coordinate_examples": shifted_examples,
    }

    binary = []
    for j in list(range(101)) + [1000, 10000, 1000000]:
        n = ((1 << (18 * j + 11)) - 23) // 27
        word, states = route(n, 3)
        require(word == (1, 2, 18 * j + 8), "binary-family word")
        require(states[-1] == 1 and all(s > 1 for s in states[:-1]), "first-hit ROOT")
        require(n.bit_length() == 18 * j + 7, "binary length")
        require(27 * n + 23 == 1 << (18 * j + 11), "binary-family source identity")
        if j < 20:
            require(format(n, "b") == "1001011" + "110110100001001011" * j,
                    "literal appended binary block")
        for modulus in [2, 3, 7, 9, 27, 64, 81, 125, 1024]:
            require(modular_family(j, modulus) == n % modulus, "modular/materialized source")
        if j < 4 or j >= 1000:
            row = {"j": j, "bits": n.bit_length(), "word": word,
                   "odd_endpoint": 1, "shortcut_clock": sum(word),
                   "ordinary_3n_plus_1_clock": sum(word) + 3}
            if j < 4:
                row["source"] = n
            else:
                row["source_little_endian_sha256"] = hashlib.sha256(
                    n.to_bytes((n.bit_length() + 7) // 8, "little")).hexdigest()
            binary.append(row)
    report["binary_family"] = binary
    huge_j = 10 ** 12
    huge_moduli = [2, 3, 7, 9, 27, 64, 81, 125, 1024, 19683, 1000000007]
    huge_residues = {}
    for modulus in huge_moduli:
        a, b = modular_family(huge_j, modulus), affine_modular_family(huge_j, modulus)
        require(a == b, "two independent giant-source readers disagree")
        huge_residues[str(modulus)] = a
    report["symbolic_trillion_repeats"] = {
        "status": "all-j theorem specialization; source not materialized",
        "j": huge_j, "bits": 18 * huge_j + 7,
        "odd_word": [1, 2, 18 * huge_j + 8], "endpoint": 1,
        "shortcut_clock": 18 * huge_j + 11,
        "ordinary_3n_plus_1_clock": 18 * huge_j + 14,
        "independently_compared_source_residues": huge_residues,
    }

    decimal = []
    for j in range(51):
        n = int("75" + "851" * j)
        word, states = route(n, 3)
        require(27 * n + 23 == 2048 * 1000 ** j, "decimal source identity")
        require(word == (1, 2, 8 + 3 * j), "decimal word")
        require(states[-1] == 125 ** j, "decimal odd-factor debt")
        if j <= 3:
            decimal.append({"j": j, "source": n, "word": word, "endpoint": states[-1]})
    require(decimal[1]["endpoint"] != 1, "hostile must reject unpaid decimal ROOT claim")
    report["decimal_family"] = {"tested_repeats": "0..50", "examples": decimal,
                                 "all_j_convergence": "NOT CLAIMED"}

    # Same-anchor pumps commute, but changing the anchor changes the source.
    for base1 in bases:
        B = base1 ** order(base1, 27)
        for base2 in bases:
            D = base2 ** order(base2, 27)
            left = pump(pump(75, 23, 27, D), 23, 27, B)
            right = pump(pump(75, 23, 27, B), 23, 27, D)
            require(left == right, "same-anchor commutativity")
            word, states = route(left, 3)
            require(word == (1, 2, 8 + v2(B) + v2(D)), "composed integer clocks")
            require(states[-1] == (B >> v2(B)) * (D >> v2(D)), "composed odd-factor cost")
    P = lambda n: 4 * n + 1
    Q = lambda n: 262144 * n + 223307
    pq, qp = P(Q(75)), Q(P(75))
    require(pq - qp == 407778, "changed-anchor translation defect")
    require(Fraction(pq - qp) == 3 * 262143 * (Fraction(23, 27) - Fraction(1, 3)),
            "general commutator identity")
    w_pq, s_pq = route(pq, 3)
    w_qp, s_qp = route(qp, 3)
    require(w_pq == (3, 2, 26) and s_pq[-1] == 1, "legal composition")
    require(P(75) % 16 != 11, "guard failure must be visible")
    require(w_qp == (1, 2, 16) and s_qp[-1] == 4075, "wrong-order terminal hostile")
    report["generator_composition"] = {
        "same_anchor_base_pairs": len(bases) ** 2,
        "changed_anchor": {"P_Q_75": pq, "Q_P_75": qp, "difference": pq - qp,
                            "first_word": w_pq, "first_endpoint": s_pq[-1],
                            "second_word": w_qp, "second_endpoint": s_qp[-1],
                            "failed_intermediate_guard_source": P(75)},
    }

    # The same long leading binary prefix permits growth and descent.
    for K in range(3, 129):
        n1, n2 = (1 << K) + 1, (1 << K) + 3
        w1, s1 = route(n1, 1)
        w2, s2 = route(n2, 1)
        require(w1 == (2,) and s1[-1] < n1, "leading-digit descending control")
        require(w2 == (1,) and s2[-1] > n2, "leading-digit growing control")
        require(Fraction(n2 - n1, 1 << (K + 1)) == Fraction(1, 1 << K), "mantissa distance")
    require(Fraction(1, 2) == Fraction(2, 4), "length-erasure control")
    # Correct interval identity is rho(Tx)=2*rho(x)-parity(x).
    endpoint_cases = []
    for n, rho, rho_next in [(-2, Fraction(1, 2), Fraction(1)),
                              (-1, Fraction(1), Fraction(1)),
                              (0, Fraction(0), Fraction(0)),
                              (1, Fraction(2, 3), Fraction(1, 3)),
                              (2, Fraction(1, 3), Fraction(2, 3))]:
        require(2 * rho - (n % 2) == rho_next, "marked interval shift")
        ordinary_doubling = 2 * rho - (2 * rho).__floor__()
        require((ordinary_doubling == rho_next) == (n not in [-2, -1]), "dyadic endpoint exception")
        require((rho_next - ordinary_doubling).denominator == 1, "circle semiconjugacy")
        endpoint_cases.append({"source": n, "rho": str(rho), "rho_next": str(rho_next),
                               "fractional_part_doubling": str(ordinary_doubling)})
    report["hostile_controls"] = {"leading_digit_pairs": 126,
                                    "parity_number_endpoint_cases": endpoint_cases}
    report["checks"] = CHECKS
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.verify:
        if args.verify.read_text() != encoded:
            raise ValueError("saved receipt differs from recomputation")
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({"checks": result["checks"],
                      "replays": result["general_controls"]["literal_odd_map_replays"],
                      "materialized_max_bits": result["binary_family"][-1]["bits"],
                      "symbolic_bits": result["symbolic_trillion_repeats"]["bits"],
                      "status": "PASS"}, sort_keys=True))


if __name__ == "__main__":
    main()
