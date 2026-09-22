"""Exact controls for affine, accelerated, odd-part, and half-step conjugacies.

Standard library only; every check survives -O. The finite affine search is
a control for the universal proofs, not a proof outside its stated box.
"""

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def v2(n):
    check(n != 0, "valuation of zero")
    n = abs(n)
    return (n & -n).bit_length() - 1


def raw(n, b):
    return 3 * n + b


def oddpart(n, b):
    value = raw(n, b)
    check(value != 0, "zero odd-part numerator")
    return value // 2**v2(value)


def accelerated(n, b):
    check(n % 2 and b % 2, "odd-accelerated domain")
    return oddpart(n, b)


def shortcut(n, b):
    check(b % 2, "shortcut needs odd parameter")
    return n // 2 if n % 2 == 0 else (3 * n + b) // 2


def swapped(n, d, e):
    check(d % 2 == 0 and e % 2, "swapped branch integrality")
    return (3 * n + d) // 2 if n % 2 == 0 else (n + e) // 2


def shell(n, b, r):
    check(v2(n) == r and v2(b) == r, "common shell domain")
    value = raw(n, b)
    return value // 2**(v2(value) - r)


def candidate_holds(a, c, b, d, inputs, operation):
    for n in inputs:
        if 3 * n + b == 0:
            continue
        image = a * n + c
        if 3 * image + d == 0:
            return False
        if operation(image, d) != a * operation(n, b) + c:
            return False
    return True


def main():
    # Raw maps: coefficient classification checked independently on a finite box.
    raw_checks = 0
    for a in range(-8, 9):
        if a == 0:
            continue
        for c in range(-5, 6):
            for b in range(-7, 8):
                d = a * b - 2 * c
                for n in range(-15, 16):
                    check(raw(a * n + c, d) == a * raw(n, b) + c, "raw affine identity")
                    raw_checks += 1

    raw_pairs = []
    for n in range(1, 101):
        h = 3 - 4 * n
        check(h < 0 and h % 4 == 3, "first raw image progression")
        check(raw(h, -2) == 3 - 4 * raw(n, -1), "first requested raw map")
        m = -n
        image = -3 - 4 * m
        check(image > 0 and image % 4 == 1, "second raw image progression")
        check(raw(image, 2) == -3 - 4 * raw(m, 1), "second requested raw map")
        if n <= 4:
            raw_pairs.append({"positive_source": n, "negative_target": h,
                              "negative_source": m, "positive_target": image})
    check(accelerated(1, -1) == 1 and oddpart(-1, -2) == -5, "first accelerated hostile")
    check(accelerated(-1, 1) == -1 and oddpart(1, 2) == 5, "second accelerated hostile")

    # Independent finite classifier: no formula is used to select survivors.
    odd_inputs = [s * n for n in range(1, 130, 2) for s in (1, -1)]
    accelerated_candidates = accelerated_survivors = 0
    for a in range(-7, 8):
        if a == 0:
            continue
        for c in range(-6, 7):
            if (a + c) % 2 != 1:
                continue
            for b in range(-7, 8, 2):
                for d in range(-7, 8, 2):
                    accelerated_candidates += 1
                    actual = candidate_holds(a, c, b, d, odd_inputs, accelerated)
                    predicted = c == 0 and a % 2 != 0 and d == a * b
                    check(actual == predicted, "accelerated affine finite classifier")
                    accelerated_survivors += actual

    full_candidates = full_survivors = 0
    full_inputs = list(range(-80, 81))
    for a in range(-5, 6):
        if a == 0:
            continue
        for c in range(-4, 5):
            for b in range(-6, 7):
                for d in range(-6, 7):
                    full_candidates += 1
                    actual = candidate_holds(a, c, b, d, full_inputs, oddpart)
                    predicted = c == 0 and a % 2 != 0 and d == a * b
                    check(actual == predicted, "full odd-part affine finite classifier")
                    full_survivors += actual

    # Full odd-part loses even scaling; retaining an exact shell restores it.
    shell_checks = decomposition_checks = 0
    for r in range(6):
        for b in range(-9, 10, 2):
            for n in range(-61, 62, 2):
                if 3 * n + b == 0:
                    continue
                check(shell(2**r * n, 2**r * b, r) == 2**r * accelerated(n, b),
                      "shell conjugacy")
                check(oddpart(2**r * n, 2**r * b) == accelerated(n, b),
                      "odd-part erases shell scale")
                shell_checks += 1
    for r in range(5):
        for s in range(5):
            for u in range(-9, 10, 2):
                for beta in range(-9, 10, 2):
                    n, b = 2**r * u, 2**s * beta
                    if 3 * n + b == 0:
                        continue
                    predicted = (3 * u + 2**(s-r) * beta if r < s else
                                 3 * 2**(r-s) * u + beta if r > s else accelerated(u, beta))
                    check(oddpart(n, b) == predicted, "valuation-stratum decomposition")
                    decomposition_checks += 1

    even_parameter_cycles = []
    for b in range(-12, 13, 2):
        fixed = -b // 2
        valid_fixed = fixed % 2 != 0
        for m in range(-51, 52, 2):
            x = m
            for t in range(1, 9):
                x = oddpart(x, b)
                check(x == 3**t * (m + b // 2) - b // 2, "even-parameter closed orbit")
                check((x == m) == (m == fixed), "even-parameter periodic-point classification")
        even_parameter_cycles.append({"b": b, "only_cycle": [fixed] if valid_fixed else None})
    for n in range(1, 102, 2):
        check(oddpart(-n, -2) < -n and oddpart(n, 2) > n, "requested sign sectors escape")
    check(oddpart(2, -2) == oddpart(1, -2) == 1, "even nonfixed start can enter fixed basin")
    check(oddpart(-2, 2) == oddpart(-1, 2) == -1, "reflected even preimage of fixed point")
    fixed_basin_checks = 0
    for beta in range(-9, 10, 2):
        for k in range(13):
            numerator = -beta * (2**k + 2)
            if numerator % 3 == 0:
                check(oddpart(numerator // 3, 2 * beta) == -beta, "even parameter fixed-basin formula")
                fixed_basin_checks += 1

    # Half-step conjugacies require transforming both branches.
    halfstep_checks = 0
    for a in range(-5, 6, 2):
        for b in range(-7, 8, 2):
            for e in range(-5, 6, 2):
                d = a * b - e
                for n in range(-50, 51):
                    check(swapped(a * n + e, d, e) == a * shortcut(n, b) + e,
                          "parity-swapped half-step affine conjugacy")
                    halfstep_checks += 1
    shortcut_candidates = shortcut_survivors = 0
    swapped_candidates = swapped_survivors = 0
    for a in range(-5, 6):
        if not a:
            continue
        for c in range(-4, 5):
            for b in range(-5, 6, 2):
                for d in range(-5, 6, 2):
                    actual = all(shortcut(a*n+c, d) == a*shortcut(n, b)+c for n in range(-8, 9))
                    predicted = a % 2 != 0 and c == 0 and d == a*b
                    check(actual == predicted, "standard shortcut affine finite classifier")
                    shortcut_candidates += 1
                    shortcut_survivors += actual
                for d in range(-4, 5, 2):
                    for e in range(-3, 4, 2):
                        actual = all(swapped(a*n+c, d, e) == a*shortcut(n, b)+c for n in range(-8, 9))
                        predicted = a % 2 != 0 and c == e and d == a*b-c
                        check(actual == predicted, "swapped shortcut affine finite classifier")
                        swapped_candidates += 1
                        swapped_survivors += actual
    check(swapped(2, -2, 1) == 2 and swapped(-2, 2, -1) == -2, "shifted fixed point controls")
    check(oddpart(1, 1) == 1 and oddpart(2, 1) == 7,
          "odd-part input quotient does not define a common forward value")

    # Commuting permutations of the 20 explicitly known tagged cycle states.
    core = [(j, side, i) for j in (1, 2, 7) for side in (0, 1) for i in range(j)]
    def advance(state):
        j, side, i = state
        return j, side, (i + 1) % j
    choices = [[(swap, r0, r1) for swap in (0, 1) for r0 in range(j) for r1 in range(j)]
               for j in (1, 2, 7)]
    centralizer = third_power_identity = 0
    for blocks in product(*choices):
        rules = dict(zip((1, 2, 7), blocks))
        def permute(state):
            j, side, i = state
            swap, r0, r1 = rules[j]
            return j, side ^ swap, (i + (r0 if side == 0 else r1)) % j
        images = {state: permute(state) for state in core}
        check(len(set(images.values())) == 20, "centralizer permutation is bijective")
        check(all(advance(images[state]) == images[advance(state)] for state in core),
              "centralizer commutes")
        centralizer += 1
        third_power_identity += all(images[images[images[state]]] == state for state in core)
    check(centralizer == 1568 and third_power_identity == 1, "known-core no order-three action")

    return {
        "status": "FINITE-EXACT controls for the proved affine and valuation classifications",
        "raw": {"checks": raw_checks, "universe": "nonzero -8<=a<=8; -5<=c<=5; -7<=b<=7; -15<=n<=15; d=ab-2c",
                "requested_sign_pair_examples": raw_pairs, "smallest_absolute_negative_raw_slope": 4},
        "accelerated_affine_search": {"candidates": accelerated_candidates, "survivors": accelerated_survivors,
                "universe": "nonzero -7<=a<=7; -6<=c<=6; a+c odd; odd b,d in [-7,7]; odd inputs +/-1..129",
                "survivor_rule": "c=0, a odd, d=ab"},
        "full_oddpart_affine_search": {"candidates": full_candidates, "survivors": full_survivors,
                "universe": "nonzero -5<=a<=5; -4<=c<=4; b,d in [-6,6]; inputs -80..80; zero source numerator omitted",
                "survivor_rule": "c=0, a odd, d=ab"},
        "shell_controls": {"checks": shell_checks, "universe": "0<=r<=5; odd -9<=b<=9; odd -61<=n<=61; numerator nonzero"},
        "valuation_decomposition": {"checks": decomposition_checks, "universe": "0<=r,s<=4; odd u,beta in [-9,9]; numerator nonzero"},
        "even_parameters": {"universe": "even -12<=b<=12; odd -51<=m<=51; iterations 1..8", "cycles": even_parameter_cycles,
                            "fixed_basin_checks": fixed_basin_checks,
                            "fixed_basin_universe": "odd beta in [-9,9]; 0<=k<=12; n=-beta(2^k+2)/3 when integral",
                            "nonfixed_even_preimage_hostiles": [[2, -2, 1], [-2, 2, -1]]},
        "halfsteps": {"direct_checks": halfstep_checks,
                "direct_universe": "odd a,e in [-5,5]; odd b in [-7,7]; -50<=n<=50; d=ab-e",
                "standard_candidates": shortcut_candidates, "standard_survivors": shortcut_survivors,
                "swapped_candidates": swapped_candidates, "swapped_survivors": swapped_survivors,
                "classifier_universe": "nonzero a in [-5,5]; c in [-4,4]; odd b in [-5,5]; n in [-8,8]; standard odd d in [-5,5]; swapped even d in [-4,4], odd e in [-3,3]"},
        "hostiles": {"raw_to_accelerated_first": {"source": 1, "H_source": -1, "H_U_source": -1, "O_target": -5},
                     "raw_to_accelerated_second": {"source": -1, "H_source": 1, "H_U_source": 1, "O_target": 5},
                     "oddpart_input_identification": {"inputs": [1, 2], "same_oddpart": 1, "O1_outputs": [1, 7]}},
        "known_cycle_core": {"tagged_states": 20, "cycle_type": "1^2 2^2 7^2",
                            "centralizer_size": centralizer, "permutations_with_cube_identity": third_power_identity,
                            "scope": "Only the three known cycles and their parameter-reflected partners"},
        "source_sha256_lf": sha256(Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(main(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8", newline="\n")
    else:
        print(rendered, end="")
