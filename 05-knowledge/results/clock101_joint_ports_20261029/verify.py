"""Exact clock lifts, finite power-free windows, and joint-port controls.

Proofs are in README.md; this finite checker is not a Collatz solver.
Python 3.10+, SymPy. All gates survive -O. No floating-point mathematics.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from pathlib import Path
import argparse
import importlib.util
import json
import sympy as sp

GATES = Counter()


def need(ok, name):
    GATES[name] += 1
    if not ok:
        raise ArithmeticError(name)


def compile_word(a, word):
    h, c = 0, 0
    prefixes = []
    for j, v in enumerate(word):
        prefixes.append((a ** j, c, 1 << h))
        c = a * c + (1 << h)
        h += v
    return a ** len(word), c, 1 << h, prefixes


def affine_power(a, b, k, modulus, x=0):
    """Binary composition, not an orbit table."""
    while k:
        if k & 1:
            x = (a * x + b) % modulus
        a, b = a * a % modulus, b * (a + 1) % modulus
        k //= 2
    return x


def geometric_power(A, C, D, k, modulus):
    """Independent division of a geometric numerator at extra precision."""
    gap = D - A
    large = modulus * gap
    alpha = A * pow(D, -1, large) % large
    numerator = (1 - pow(alpha, k, large)) % large
    need(numerator % gap == 0, "geometric exact division")
    return C * (numerator // gap) % modulus


def lift_clock(A, C, D, root_num, root_den, p, e):
    i, q = 0, 1
    digits = []
    for level in range(e):
        modulus = q * p
        alpha = A * pow(D, -1, modulus) % modulus
        beta = C * pow(D, -1, modulus) % modulus
        target = root_num * pow(root_den, -1, modulus) % modulus
        x = affine_power(alpha, beta, i, modulus)
        need((target - x) % q == 0, "compatible Hensel prefix")
        digit = (target - x) // q * pow(beta, -1, p) % p
        i += digit * q
        digits.append(digit)
        need(affine_power(alpha, beta, i, modulus) == target, "lifted clock hits root")
        need(geometric_power(A, C, D, i, modulus) == target, "independent geometric hit")
        q = modulus
    return i, digits


def capacities(times, modulus, phase, k=3):
    times = sorted(times)
    return [(right - 1 - phase) // k - (left - phase) // k
            for left, right in zip(times, times[1:] + [modulus])]


def survival(caps, r):
    return sum(max(c - r, 0) for c in caps)


def gf2_gcd(a, b):
    while b:
        while a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def cyclic_three(times, modulus):
    outdegrees = [sum(0 < (b - a) % modulus <= modulus // 2 for b in times)
                  for a in times]
    return outdegrees == [1, 1, 1]


def brute_clock(a, word, p, e, expected):
    """Literal native one-step iteration; all states and all survival lengths."""
    N, M = p ** e, 3 * p ** e
    x = 0
    states = [set() for _ in range(3)]
    zero_times = []
    for t in range(M):
        j = t % 3
        need(x not in states[j], "joint clock does not repeat early")
        states[j].add(x)
        if x == 0:
            zero_times.append(t)
        x = (a * x + 1) * pow(1 << word[j], -1, N) % N
    need(x == 0 and all(len(s) == N for s in states), "joint clock closes with all states")
    need(zero_times == expected["joint_zero_times"], "native orbit vs lifted roots")
    # Distance to the next zero, including a zero at the current node.
    anchors = zero_times + [M]
    pointer = 0
    histogram = [Counter() for _ in range(3)]
    for t in range(M):
        while anchors[pointer] < t:
            pointer += 1
        distance = anchors[pointer] - t
        histogram[t % 3][(distance - 1) // 3] += 1
    for j in range(3):
        count = 0
        for r in range(N, -1, -1):
            count += histogram[j][r]
            need(count == survival(expected["phase_capacities"][j], r), "all-length survival vs literal native clock")
    # Mod-2 certificates imply at least these rational ranks; the proofs give equality.
    block_polynomial = sum(1 << i for i in expected["named_block_indices"])
    joint_polynomial = sum(1 << i for i in zero_times)
    need(gf2_gcd((1 << N) | 1, block_polynomial) == 1, "finite full block rank mod two")
    need(gf2_gcd((1 << M) | 1, joint_polynomial) == 7, "finite joint rank exactly M minus two mod two")


def audit_towers():
    output = {}
    definitions = [("115", 3, (1, 1, 5), 101),
                   ("133", 3, (1, 3, 3), 101),
                   ("7x1_135", 7, (1, 3, 5), 13)]
    for label, a, word, p in definitions:
        A, C, D, prefixes = compile_word(a, word)
        need((D - A) % p == 0 and gcd(C * D, p) == 1, "unit-carry resonance")
        phases = []
        for j in range(3):
            Aj, Cj, Dj, _ = compile_word(a, word[j:] + word[:j])
            phases.append(str(Fraction(Cj, Dj - Aj)))
        rows = []
        for e in range(1, 13):
            N, M = p ** e, 3 * p ** e
            indices, roots, digits = [], [], []
            for slope, intercept, denominator in prefixes:
                i, ds = lift_clock(A, C, D, -intercept, slope, p, e)
                indices.append(i)
                digits.append(ds)
                roots.append(-intercept * pow(slope, -1, N) % N)
            need(len(set(roots)) == 3, "three distinct roots")
            times = sorted(3 * indices[j] + j for j in range(3))
            differences = Counter((y - x) % M for x in times for y in times if x != y)
            need(len(differences) == 6 and set(differences.values()) == {1},
                 "all six joint-clock differences distinct")
            offset = 3 * indices[1] + 1
            need({t for t in times if (t + offset) % M in times} == {0},
                 "two timed ports isolate named zero state")
            caps = [capacities(times, M, j) for j in range(3)]
            # Independently recompile each rotated word from its own zero source.
            # This catches the phase-order changes first visible at exponent 3.
            for j in range(3):
                rotated = word[j:] + word[:j]
                Ar, Cr, Dr, pr = compile_word(a, rotated)
                rotated_indices = sorted(lift_clock(Ar, Cr, Dr, -c0, a0, p, e)[0]
                                         for a0, c0, d0 in pr)
                rotated_caps = [rotated_indices[1], rotated_indices[2] - rotated_indices[1],
                                N - rotated_indices[2] - 1]
                need(sorted(rotated_caps) == sorted(caps[j]), "independent rotated-word capacities")
            limits = [max(c) - 1 for c in caps]
            need(max(limits) - min(limits) <= 1, "native phase limits differ at most one")
            need(gcd(gcd(times[1], times[2]), p) == 1, "joint Fourier zeros only primitive thirds")
            for j in range(3):
                need(survival(caps[j], limits[j]) > 0 and survival(caps[j], limits[j] + 1) == 0,
                     "sharp local boundary")
            row = {"exponent": e, "modulus": N, "roots": roots, "named_block_indices": indices,
                   "clock_digits": digits, "joint_zero_times": times, "phase_capacities": caps,
                   "last_possible_blocks_by_phase": limits,
                   "two_port_delta_offset": offset,
                   "cyclic_root_triple": cyclic_three(times, M),
                   "block_Q_rank": N, "joint_Q_rank": M - 2, "joint_deterministic_states": M}
            if e <= 2:
                brute_clock(a, word, p, e, row)
            rows.append(row)
        output[label] = {"multiplier": a, "word": word, "prime": p, "A_C_D": [A, C, D],
                         "rational_phases": phases, "rows": rows}
    a, b = output["115"]["rows"], output["133"]["rows"]
    need(a[1]["last_possible_blocks_by_phase"][0] == 4985 and b[1]["last_possible_blocks_by_phase"][0] == 8600,
         "inherited squarefree limits")
    need(a[2]["last_possible_blocks_by_phase"] == [576469, 576470, 576469], "115 cube-free phase limits")
    need(b[2]["last_possible_blocks_by_phase"] == [979116, 979115, 979115], "133 cube-free phase limits")
    need(a[3]["last_possible_blocks_by_phase"][0] == 64332491 and b[3]["last_possible_blocks_by_phase"][0] == 45324462,
         "fourth-free ranking reversal")
    need(output["7x1_135"]["rows"][1]["last_possible_blocks_by_phase"][0] == 63, "recovered 79 over 169 bound")
    return output


def deck(S, N, order):
    return Counter(tuple((v - base) % N for v in rest)
                   for base in S for rest in product(S, repeat=order - 1))


def translate_normal_form(S, N):
    return min(tuple(sorted((x - base) % N for x in S)) for base in S)


def audit_ports(towers):
    rows = {}
    for label in ("115", "133"):
        for e in (1, 2, 3, 4):
            data = towers[label]["rows"][e - 1]
            N = data["modulus"]
            S = set(data["named_block_indices"])
            reflection = {(-x) % N for x in S}
            two, three = deck(S, N, 2), deck(S, N, 3)
            need(two == deck(reflection, N, 2), "pair deck reflection blindness")
            need(three != deck(reflection, N, 3), "third port separates reflection")
            witness = next(key for key in sorted(three) if three[key] != deck(reflection, N, 3)[key])
            support = {k[0] for k in two if k != (0,)}
            candidates = {translate_normal_form({0, u, v}, N) for u, v in combinations(support, 2)
                          if deck({0, u, v}, N, 2) == two}
            need(candidates == {translate_normal_form(S, N), translate_normal_form(reflection, N)},
                 "pair deck candidates exactly mirror pair in specified examples")
            need(deck({(x + 17) % N for x in S}, N, 3) == three, "third deck loses absolute address")
            rows[label + "/" + str(e)] = {"third_port_witness": witness,
                "values": [three[witness], deck(reflection, N, 3)[witness]],
                "pair_deck_candidate_orbits": len(candidates)}
    N = 101 ** 2
    need(deck(set(towers["115"]["rows"][1]["named_block_indices"]), N, 2) !=
         deck(set(towers["133"]["rows"][1]["named_block_indices"]), N, 2), "two ports already distinguish 115 and 133")
    # Small full census independently exercises rank and three-deck reconstruction.
    X = sp.Symbol("X")
    for p, e in ((5, 1), (7, 1), (3, 2), (5, 2)):
        N = p ** e
        records = {}
        for size in (1, 2, 3):
            if size % p == 0:
                continue
            for rest in combinations(range(1, N), size - 1):
                S = {0, *rest}
                polynomial = sp.Poly(sum(X ** t for t in S), X)
                need(sp.gcd(polynomial, sp.Poly(X ** N - 1, X)).degree() == 0,
                     "small prime-power observer rank census")
                key = tuple(sorted(deck(S, N, 3).items()))
                form = translate_normal_form(S, N)
                need(key not in records or records[key] == form, "small three-deck uniqueness census")
                records[key] = form
    hostile = {0, 35, 70, 1, 22, 43, 64, 85}
    polynomial = sp.Poly(sum(X ** t for t in hostile), X)
    common = sp.gcd(polynomial, sp.Poly(X ** 105 - 1, X))
    need(gcd(len(hostile), 105) == 1 and common.degree() == 56, "composite clock full-rank hostile")
    need(common == sp.Poly(sp.cyclotomic_poly(15, X) * sp.cyclotomic_poly(105, X), X),
         "composite obstruction exact cyclotomic factors")
    rows["composite_hostile"] = {"N": 105, "mask": sorted(hostile), "deterministic_states": 105,
                                 "Q_rank": 49, "missing_factors": [15, 105]}
    return rows


def audit_arithmetic_controls():
    need(int(sp.n_order(2, 79)) == 39 and int(sp.n_order(2, 19)) == 18,
         "historical power-tournament distinction")
    need(all(int(sp.n_order(p, 55)) == 10 and pow(p, 5, 55) == 54 for p in (19, 79)),
         "historical conjugation Frobenius clocks")
    Z = sp.Symbol("Z")
    need(sp.rem((3 - Z) ** 2 - (8 - 7 * Z), Z ** 2 + Z + 1, Z) == 0 and 512 - 343 == 13 ** 2,
         "historical Eisenstein square clock")
    need(119 ** 2 + 120 ** 2 == 169 ** 2, "historical Pell hypotenuse")
    ring_field = []
    for a, w, p, phase_prime in ((3, (1, 1, 5), 101, 13), (3, (1, 3, 3), 101, 7),
                                  (7, (1, 3, 5), 13, 19)):
        A, C, D, prefixes = compile_word(a, w)
        for slope, intercept, denominator in prefixes:
            value = Fraction(slope * C + intercept * (D - A), denominator * (D - A))
            need(max(sp.factorint(value.numerator).values(), default=0) <= 2,
                 "all nonresonant numerator valuations at most two")
        modulus_ring = p * p
        alpha = A * pow(D, -1, modulus_ring) % modulus_ring
        beta = C * pow(D, -1, modulus_ring) % modulus_ring
        c = (alpha - 1) // p
        beta0, beta1 = beta % p, beta // p
        for x in range(modulus_ring):
            low, high = x % p, x // p
            low_next = (low + beta0) % p
            high_next = (high + beta1 + c * low + (low + beta0) // p) % p
            need((alpha * x + beta) % modulus_ring == low_next + p * high_next,
                 "two-digit carry skew product")
        need(affine_power(alpha, beta, p, modulus_ring) != 0,
             "prime steps do not close ring clock")
        need((p * beta0) % p == 0 and beta0 != 0, "field translation order exactly prime")
        ring_field.append({"word": w, "p": p, "alpha": alpha, "beta": beta, "ring_return_period": p * p,
                           "field_return_period": p, "field_number_of_cycles": p})
        modulus = phase_prime ** 2
        f = C * pow(D - A, -1, modulus) % modulus
        x = (f + phase_prime) % modulus
        start = x
        steps = 0
        while True:
            for slope, intercept, denominator in prefixes:
                need((slope * x + intercept) * pow(denominator, -1, modulus) % modulus != 0,
                     "exceptional square-prime repair all phases")
            x = (A * x + C) * pow(D, -1, modulus) % modulus
            steps += 1
            if x == start:
                break
            need(steps < modulus, "exceptional repair terminates finite orbit")
    fractions = []
    for word in ((1, 1, 5, 1), (1, 3, 3, 1)):
        A, C, D, _ = compile_word(3, word)
        fractions.append(Fraction(C, D - A))
    need(fractions == [Fraction(37, 35), Fraction(221, 175)], "append 35 versus 175 retained")
    path = Path(__file__).resolve().parent.parent / "four_vertex_glue_fermat_20261029" / "verify.py"
    spec = importlib.util.spec_from_file_location("inherited_four_vertex", path)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    return {"append_fractions": list(map(str, fractions)), "hostile_81_27": old.hostile_control(),
            "ring_vs_field": ring_field}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    towers = audit_towers()
    output = {"towers": towers, "ports": audit_ports(towers), "controls": audit_arithmetic_controls()}
    output["gates"] = dict(sorted(GATES.items()))
    output["total_gates"] = sum(GATES.values())
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"total_gates": output["total_gates"], "tower_depth": 12,
                      "finite_orbit_universe": "all phase/residue states through exponent 2 for three words"}))


if __name__ == "__main__":
    main()
