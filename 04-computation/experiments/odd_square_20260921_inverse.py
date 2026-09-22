"""Exact guarded inverse, edge-content, and primitive parameter-five controls.

The mixed completion universe is all 13 prefixes of length 0..2 with
exponents in {1,2,3}, six named targets, every residue modulo 9, and every
unit residue modulo 25. This is not a cycle census or a convergence test.
Standard library only; all checks survive -O.
"""

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from math import gcd
from pathlib import Path


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def valuation(n, p):
    check(n != 0, "zero valuation")
    n = abs(n)
    if p == 2:
        return (n & -n).bit_length() - 1
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def step(n, b):
    check(n % 2 and b % 2, "odd domain")
    numerator = 3 * n + b
    k = valuation(numerator, 2)
    return numerator // 2**k, k


def word_data(word):
    total, carry = 0, 0
    for k in word:
        check(k >= 1, "positive exponent")
        carry = 3 * carry + 2**total
        total += k
    return total, carry


def triangle(x, y):
    return abs(x * y), abs(x*x-y*y) // 2, (x*x+y*y) // 2


@lru_cache(maxsize=None)
def base_two_logs(prime, exponent):
    check(prime in (3, 5) and exponent >= 1, "log modulus")
    modulus = prime**exponent
    period = (prime - 1) * prime**(exponent - 1)
    logs = {}
    power = 1
    for k in range(period):
        check(power not in logs, "power-two order too short")
        logs[power] = k
        power = (2 * power) % modulus
    check(power == 1 and len(logs) == period, "full unit group not generated")
    return modulus, period, logs


def completion(prefix, target, s, residue3, t, residue5, family_index=0):
    """Construct prefix+(k,ell), k=1 or2, ending at target under U_-5."""
    check(target % 2 and gcd(target, 15) == 1, "primitive admissible target")
    check(s >= 0 and t >= 1 and gcd(residue5, 5) == 1, "source residue domain")
    total, carry = word_data(prefix)
    length = len(prefix) + 2
    modulus3, period3, logs3 = base_two_logs(3, length + s)
    modulus5, period5, logs5 = base_two_logs(5, t)
    possible = []
    for k in (1, 2):
        c = 2**(total+k) + 3 * 2**total + 9 * carry
        rhs3 = ((3**length * residue3 - 5*c) * pow(target, -1, modulus3)) % modulus3
        rhs5 = ((3**length * residue5 - 5*c) * pow(target, -1, modulus5)) % modulus5
        check(rhs3 in logs3 and rhs5 in logs5, "nonunit completion residue")
        e3, e5 = logs3[rhs3], logs5[rhs5]
        if (e3 - e5) % 2 == 0:
            possible.append((k, c, e3, e5))
    check(len(possible) == 1, "one inserted exponent parity must work")
    k, c, e3, e5 = possible[0]
    check(gcd(period3, period5) == 2, "exponent CRT overlap")
    j = (((e5-e3)//2) * pow(period3//2, -1, period5//2)) % (period5//2)
    period = period3 * (period5//2)
    exponent = (e3 + period3*j) % period
    if exponent <= total + k:
        exponent += ((total+k-exponent)//period + 1)*period
    exponent += family_index * period
    while True:
        ell = exponent - total - k
        full_word = tuple(prefix) + (k, ell)
        numerator = 2**exponent * target + 5*c
        check(numerator % 3**length == 0, "source integrality")
        source = numerator // 3**length
        # Independent backward reconstruction checks every intermediate /3 guard.
        backward = [target]
        for power in reversed(full_word):
            value = 2**power * backward[-1] + 5
            check(value % 3 == 0, "inverse guard")
            value //= 3
            check(value % 2, "odd inverse node")
            backward.append(value)
        check(backward[-1] == source, "inverse composition source")
        nodes = tuple(reversed(backward))
        if all(n * target > 0 for n in nodes):
            break
        exponent += period
    # Independent forward divisions verify the prescribed actual valuations.
    for i, power in enumerate(full_word):
        check(step(nodes[i], -5) == (nodes[i+1], power), "actual forward valuation")
        check(gcd(nodes[i], nodes[i+1]) == 1, "primitive completed edge")
    check(source % 3**s == residue3 % 3**s, "prescribed ternary residue")
    check(source % 5**t == residue5 % 5**t, "prescribed quinary residue")
    cylinder_modulus = 2**(total+1)
    cylinder = ((2**total + 5*carry) * pow(3**len(prefix), -1, cylinder_modulus)) % cylinder_modulus
    check(source % cylinder_modulus == cylinder, "preserved prefix cylinder")
    return {"source": source, "nodes": nodes, "word": full_word,
            "inserted_k": k, "terminal_ell": ell, "exponent_period": period}


def main():
    edge_checks = 0
    growth_hostiles = []
    for b in range(-21, 22, 2):
        for x in range(-201, 202, 2):
            if 3*x+b == 0:
                continue
            y, k = step(x, b)
            g = gcd(x, b)
            check(gcd(x, y) == g, "sharp edge content")
            check(gcd(y, b) == gcd(3*x, b), "next state content")
            check(valuation(gcd(y, b), 3) == min(valuation(x, 3)+1, valuation(b, 3)),
                  "three-primary content update")
            a, leg, hyp = triangle(x, y)
            check(a*a + leg*leg == hyp*hyp, "Pythagorean identity")
            check(gcd(gcd(a, leg), hyp) == g*g, "triangle content is edge content squared")
            primitive = triangle(x//g, y//g)
            check((a, leg, hyp) == tuple(g*g*z for z in primitive), "primitive triangle dilation")
            check(step(x//g, b//g) == (y//g, k), "edge content parameter reduction")
            if gcd(y, b) != g and len(growth_hostiles) < 4:
                growth_hostiles.append({"b": b, "x": x, "y": y, "edge_content": g,
                                        "next_edge_content": gcd(y, b)})
            edge_checks += 1
    check(step(1, 3) == (3, 1), "minimal three-primary growth witness")

    inverse_checks = 0
    for b in range(-15, 16, 2):
        for target in range(-31, 32, 2):
            for k in range(1, 9):
                admissible = (2**k * target - b) % 3 == 0
                predicted = ((target % 3 == 0) if b % 3 == 0 else
                             (target % 3 != 0 and pow(2, k, 3)*target % 3 == b % 3))
                check(admissible == predicted, "inverse target/parity guard")
                if admissible:
                    x = (2**k * target-b)//3
                    check(step(x, b) == (target, k), "inverse arrow exact valuation")
                    check(step(4*x+b, b) == (target, k+2), "inverse braid exact lift")
                inverse_checks += 1

    cycle_nodes = [
        (5,), (25,35), (85,125,185,275,205,305,455), (-1,), (-5,),
        (-19,-31,-49), (-23,-37,-29),
        (-187,-283,-427,-643,-967,-1453,-1091,-1639,-2461,-1847,-2773,-2081,-781,-587,-883,-1327,-1993),
        (-347,-523,-787,-1183,-1777,-667,-1003,-1507,-2263,-3397,-2549,-1913,-359,-541,-407,-613,-461),
    ]
    cycles = []
    for nodes in cycle_nodes:
        word = []
        content = gcd(nodes[0], 5)
        for i, n in enumerate(nodes):
            target, k = step(n, -5)
            check(target == nodes[(i+1) % len(nodes)], "inherited cycle replay")
            check(gcd(n, target) == content, "cycle edge content")
            word.append(k)
        total, carry = word_data(word)
        gap = 2**total-3**len(word)
        denominator = abs(gap)//gcd(carry, abs(gap))
        check(gap*nodes[0] == -5*carry and denominator == 5//content, "cycle carry/content gate")
        reduced = tuple(n//content for n in nodes)
        for i, n in enumerate(reduced):
            check(step(n, -5//content)[0] == reduced[(i+1) % len(reduced)], "reduced cycle")
        cycles.append({"nodes": nodes, "word": word, "L": len(word), "K": total,
                       "B": carry, "Delta": gap, "q": denominator, "content": content,
                       "reduced_parameter": -5//content, "reduced_nodes": reduced})

    prefixes = [word for length in range(3) for word in product((1,2,3), repeat=length)]
    targets = (1,-1,-19,-23,-187,-347)
    residues5 = [r for r in range(25) if r % 5]
    count = next_family_count = 0
    max_bits = 0
    inserted_counts = {1: 0, 2: 0}
    for prefix in prefixes:
        for target in targets:
            for r3 in range(9):
                for r5 in residues5:
                    first = completion(prefix, target, 2, r3, 2, r5)
                    second = completion(prefix, target, 2, r3, 2, r5, family_index=2)
                    check(abs(second["source"]) > abs(first["source"]), "distinct unbounded completion family")
                    check(first["word"][:-1] == second["word"][:-1], "only final exponent changes")
                    inserted_counts[first["inserted_k"]] += 1
                    max_bits = max(max_bits, abs(first["source"]).bit_length())
                    count += 1
                    next_family_count += 1
    check(count == 14040, "completion universe count")
    sharp = completion((), 1, 0, 0, 1, 4)
    check(sharp["source"] == 459 and sharp["word"] == (2,10), "two-tail sharp witness")
    check(sharp["nodes"] == (459,343,1), "sharp witness direct trajectory")
    one_step_residues = {(2**k+5)//3 % 5 for k in range(2, 26, 2)}
    check(one_step_residues == {2,3} and 4 not in one_step_residues | {1},
          "zero-or-one-step completion cannot cover all unit residues")
    check(step(5,-5) == (5,1) and step(1,-5) == (-1,1), "content and sign portal hostiles")
    check(step(11,1)[0] == 17 and 11 % 6 == 17 % 6 == 5, "no forced mod-six alternation")
    check((3*5+1)//2 == 8, "one halving does not guarantee an odd output")

    return {
        "status": "PROVED statements in note; FINITE-EXACT controls here; no all-cycle census",
        "edge_controls": {"count": edge_checks, "universe": "odd b in [-21,21], odd x in [-201,201], 3x+b nonzero",
                          "three_primary_examples": growth_hostiles,
                          "minimal_three_primary_hostile": {"b":3,"x":1,"y":3}},
        "inverse_controls": {"count": inverse_checks, "universe": "odd b in [-15,15], odd target in [-31,31], k=1..8"},
        "cycle_witnesses": {"scope": "Nine inherited b=-5 cycles, each directly replayed; no completeness assertion", "cycles": cycles},
        "primitive_mixed_completions": {"count": count, "additional_family_checks": next_family_count,
             "prefixes": prefixes, "targets": targets, "source_modulus3":9,"source_modulus5":25,
             "source_residues3":list(range(9)),"source_residues5":residues5,
             "maximum_first_source_bit_length":max_bits,"inserted_k_counts":inserted_counts},
        "two_tail_sharpness": sharp,
        "one_step_source_residues_mod5_for_target1": sorted(one_step_residues),
        "source_sha256_lf":sha256(Path(__file__).read_bytes().replace(b"\r\n",b"\n")).hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    result = main()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps({"status":"PASS", "output":str(args.output),
                      "mixed_completions":result["primitive_mixed_completions"]["count"],
                      "source_sha256_lf":result["source_sha256_lf"]}))
