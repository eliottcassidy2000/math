"""Exact constructive controls for finite-word completion in 3n+b graphs.

No trajectory census or claim that every starting integer reaches a cycle.
Checks use explicit exceptions and remain active under python -O.
"""
from itertools import product
from pathlib import Path
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def word_data(word):
    k_total, carry = 0, 0
    for k in word:
        require(k >= 1, "positive halving exponents required")
        carry = 3 * carry + 2 ** k_total
        k_total += k
    return k_total, carry


def log4_principal(value, exponent):
    """Unique t mod 3^(exponent-1) with 4^t=value mod 3^exponent."""
    modulus = 3 ** exponent
    require(exponent >= 1 and value % 3 == 1, "principal unit required")
    t, period = 0, 1
    for depth in range(2, exponent + 1):
        target_modulus = 3 ** depth
        hits = [t + a * period for a in range(3)
                if pow(4, t + a * period, target_modulus)
                == value % target_modulus]
        require(len(hits) == 1, "ternary lift not unique")
        t = hits[0]
        period *= 3
    require(pow(4, t, modulus) == value % modulus, "log verification")
    return t


def complete(word, b, target, residue=0, ternary_depth=0):
    """Complete word to target, prescribing source mod3^ternary_depth.

    Source has the same sign as target; b is odd and coprime to3.
    Returns compact exponent certificate, plus exact integer trajectory.
    """
    require(b % 2 and b % 3 and target % 2 and target % 3,
            "odd b,target, both coprime to3 required")
    length = len(word)
    total, carry = word_data(word)
    k0 = 1 if (2 * target - b) % 3 == 0 else 2
    precision = length + 1 + ternary_depth
    modulus = 3 ** precision
    # Numerator = 2^total * (2^k0 * 4^t * target - b) -3*b*carry.
    # Prescribe numerator /3^(length+1) modulo3^ternary_depth directly.
    rhs = (b * (2 ** total + 3 * carry)
           + 3 ** (length + 1) * residue)
    unit = rhs * pow(2 ** (total + k0) * target, -1, modulus) % modulus
    t = log4_principal(unit, precision)
    period = 3 ** (length + ternary_depth)
    while True:
        predecessor = (2 ** (k0 + 2 * t) * target - b) // 3
        numerator = 2 ** total * predecessor - b * carry
        require(numerator % 3 ** length == 0, "source integrality")
        source = numerator // 3 ** length
        if source * target > 0 and predecessor * target > 0:
            break
        t += period
    states = [source]
    exponents = []
    for _ in range(length + 1):
        image = 3 * states[-1] + b
        require(image != 0, "zero reached")
        k = (abs(image) & -abs(image)).bit_length() - 1
        exponents.append(k)
        states.append(image // 2 ** k)
    require(tuple(exponents[:-1]) == tuple(word), "incorrect prefix")
    require(states[-2] == predecessor and states[-1] == target,
            "incorrect endpoint")
    require(exponents[-1] == k0 + 2 * t, "incorrect final exponent")
    require(source % 3 ** ternary_depth == residue % 3 ** ternary_depth,
            "incorrect ternary residue")
    require(all(n % 2 for n in states), "nonodd node")
    return {"b": b, "target": target, "word": list(word), "k0": k0,
            "t": t, "K": total, "B": carry, "source": source,
            "states": states, "exponents": exponents}


def compatible_word(odd_residue, binary_depth, b):
    """Read a realizable prefix which determines at least binary_depth bits."""
    n = odd_residue
    word = []
    total = 0
    while total < binary_depth - 1:
        image = 3 * n + b
        require(image != 0, "zero encountered in prefix extraction")
        k = (abs(image) & -abs(image)).bit_length() - 1
        word.append(k)
        total += k
        n = image // 2 ** k
    return tuple(word)


def small_same_class_witnesses(modulus=72, limit=200000):
    """Independent direct iteration, not the completion construction."""
    roots = {1: 1, 5: 5, 7: 5,
             17: 17, 25: 17, 37: 17, 55: 17, 41: 17, 61: 17, 91: 17}
    found = {}
    for source in range(1, limit + 1, 2):
        if source % modulus != 1:
            continue
        n, seen = source, set()
        for step in range(10000):
            if n in roots:
                root = roots[n]
                if root not in found:
                    found[root] = {"source": source, "hit": n, "steps": step}
                break
            if n in seen:
                break
            seen.add(n)
            image = 3 * n - 1
            k = (image & -image).bit_length() - 1
            n = image // 2 ** k
        if len(found) == 3:
            break
    require(len(found) == 3, "three direct witnesses missing")
    return found


def run():
    checks, max_bits = 0, 0
    # Full rectangular universe; no accepted-word filter.
    for length in range(5):
        for word in product(range(1, 4), repeat=length):
            for b, target in [(1, 1), (-1, 1), (-1, 5), (-1, 17),
                              (-5, 5), (-5, -1), (5, 1)]:
                for residue in range(9):
                    row = complete(word, b, target, residue, 2)
                    checks += 1
                    max_bits = max(max_bits, abs(row["source"]).bit_length())
    cylinder_checks = 0
    for b, target in [(1, 1), (-1, 1), (-1, 5), (-1, 17)]:
        for depth in range(1, 6):
            for r2 in range(1, 2 ** depth, 2):
                word = compatible_word(r2, depth, b)
                for r3 in range(9):
                    row = complete(word, b, target, r3, 2)
                    require(row["source"] % 2 ** depth == r2,
                            "binary cylinder not preserved")
                    cylinder_checks += 1
    # Genuine exclusion: b=-5, target5 ancestors have gcd(n,5)=5.
    # Thus 2/3-adic density cannot be advertised as density modulo every M.
    hostile = complete((1, 2, 1), -5, 5, 7, 2)
    require(hostile["source"] % 5 == 0, "gcd stratum lost")
    examples = []
    for target in [1, 5, 17]:
        row = complete((1, 1), -1, target, 1, 1)
        examples.append(row)
    return {
        "status": "FINITE-EXACT; general statements proved in companion note",
        "word_universe": "L=0..4, exponents1..3, seven (b,target) pairs, all9 residues mod9",
        "word_completions": checks, "maximum_source_bit_length": max_bits,
        "mixed_cylinder_universe": "H=1..5,s=2, all odd residues, four (b,target) pairs",
        "mixed_cylinder_checks": cylinder_checks,
        "direct_same_class_mod72": small_same_class_witnesses(),
        "common_prefix_examples": examples,
        "hostile_control": "b=-5 ancestors of5 stay divisible by5",
        "limitations": ["No completeness theorem for any cycle census",
                        "No natural density or height-efficient cover proved",
                        "A different completion integer is used for each finite prefix"]}


if __name__ == "__main__":
    result = run()
    destination = Path(__file__).with_suffix(".json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
