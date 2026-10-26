"""Exact finite controls for the finite-center observer obstruction.

No imports from prior mathematical producers. Assertions are not used, so
normal and optimized Python execute the same gates. The universal statements
are proved in README.md; these bounded checks do not establish nonhalting.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def v2(x):
    x = F(x)
    if not x:
        return math.inf
    n, d = abs(x.numerator), x.denominator
    return (n & -n).bit_length() - (d & -d).bit_length()


def carry(word):
    out = 0
    for i, bit in enumerate(word):
        out = 3 * out + int(bit) * 2**i
    return out


def center(word):
    return -F(carry(word), 3**len(word) - 2**len(word))


def forward(word, z):
    return F(3**len(word) * z + carry(word), 2**len(word))


def backward(word, z):
    return F(2**len(word) * z - carry(word), 3**len(word))


def balanced(word):
    q = 0
    for bit in word:
        q += 3 * int(bit) - 1
        if q < 0:
            return False
    return bool(word) and q == 0


def residue(word):
    modulus = 2**len(word)
    return -carry(word) * pow(3**len(word), -1, modulus) % modulus


def crt_two_power_odd(a, power, b, odd):
    require(power >= 0 and odd >= 1 and odd % 2, "CRT domain")
    two = 2**power
    out = a % two
    out += two * ((b - out) * pow(two, -1, odd) % odd)
    return out, two * odd


def positive_lift(a, modulus, floor=3):
    if a < floor:
        a += ((floor - a + modulus - 1) // modulus) * modulus
    return a


def native_prefix(z, q, length):
    bits = []
    for _ in range(length):
        bit = z % 2
        require(not (q == 1 and bit == 0), "unexpected blocked macro")
        require(z >= 3 and q >= 1, "physical sector or credit")
        bits.append(str(bit))
        z = (3 * z + bit) // 2
        q += 3 * bit - 1
    return z, q, "".join(bits)


def halt_word(z, q, word):
    executed = 0
    for i, letter in enumerate(word):
        bit = z % 2
        require(bit == int(letter), "literal parity differs")
        require(z >= 3 and q >= 1, "positive physical state")
        if q == 1 and bit == 0:
            require(i == len(word) - 1, "blocked before final symbol")
            return executed
        require(i != len(word) - 1, "last symbol was executed")
        z = (3 * z + bit) // 2
        q += 3 * bit - 1
        executed += 1
    raise RuntimeError("no final blocked symbol")


@lru_cache(None)
def trees(internal):
    if internal == 0:
        return ("0",)
    return tuple(
        "1" + u + v + w
        for a in range(internal)
        for b in range(internal - a)
        for u in trees(a)
        for v in trees(b)
        for w in trees(internal - 1 - a - b)
    )


def observers(sensors, z):
    return tuple(v2(F(z) - s) for s in sensors)


def run():
    p, qword = "100", "101000"
    rp, rq = center(p), center(qword)
    require((rp, rq) == (-F(9, 19), -F(351, 665)), "base centers")
    rho = backward(p, rq)
    require(rho == -F(977, 1995), "missing pullback")
    require(rho.denominator % 3 == 0 and rho.numerator % 3 != 0,
            "pullback is not 3-integral")

    words = tuple(w[:-1] for e in range(1, 5) for w in trees(e))
    require(len(words) == 71 and all(map(balanced, words)), "tree bank")
    bank = tuple(sorted({F(0)} | {center(w) for w in words}))
    require(len(bank) == 67, "distinct sensor bank")
    require(all(s.denominator % 3 for s in bank), "periodic centers 3-integral")
    h = max(v2(rho - s) for s in bank)
    require(h == 11, "bank separation depth")

    # Test the transport identity, including zero/positive sensors and infinity.
    transport = 0
    for word in (p, qword, p + p, p + qword, qword + p):
        for j in range(1, 17):
            z = positive_lift(residue(word), 2**len(word)) + j * 2**len(word)
            y = forward(word, z)
            require(y.denominator == 1, "cylinder integrality")
            for sensor in (F(0), F(7), -F(5, 7), F(1, 3), rp, rq, y):
                left = v2(y - sensor)
                right = v2(F(z) - backward(word, sensor)) - len(word)
                require(left == right, "exact pullback identity")
                transport += 1

    # Generic necessity: arbitrary fixed M; compatible binary class is chosen
    # from rho, not claimed equal to an independently fixed integer source.
    generic_cases = 0
    moduli = (1, 2, 8, 15, 30, 2**18 * 105)
    sensor_sets = ((F(0),), (F(7),), (-F(5, 7),),
                   (rp, rq, F(0), F(1, 3), F(7)), bank)
    for sensors in sensor_sets:
        # Target sensor belongs to the bank; pullback is deliberately missing.
        target = next(s for s in sensors if backward(p, s) not in sensors)
        missing = backward(p, target)
        separation = max(v2(missing - s) for s in sensors)
        a, b = missing.numerator, missing.denominator
        require(b % 2, "odd pullback denominator")
        for modulus in moduli:
            e = int(v2(modulus))
            odd = modulus >> e
            base = max(int(separation) + 1, len(p), e)
            signatures = []
            for n in (base, base + 1, base + 4, base + 9):
                binary = (a + 2**n) * pow(b, -1, 2**(n + 1))
                z, period = crt_two_power_odd(binary, n + 1, 5, odd)
                z = positive_lift(z, period, 10**4)
                require(v2(F(z) - missing) == n, "exact missing depth")
                y = forward(p, z)
                require(y.denominator == 1, "missing cylinder")
                require(v2(y - target) == n - len(p), "unbounded target")
                signatures.append((z % modulus, observers(sensors, z)))
                for initial_q in (1, 5):
                    end, credit, got = native_prefix(z, initial_q, len(p))
                    require(got == p and end == y and credit == initial_q,
                            "balanced legality")
                generic_cases += 1
            require(len(set(signatures)) == 1, "observer fibre not constant")

    # Equal-center/powers is a genuine positive survivor, not a counterexample.
    powers_cases = 0
    for power in range(1, 6):
        word = p * power
        require(center(word) == rp and backward(word, rp) == rp,
                "powers have common center")
        z = positive_lift(residue(word), 2**len(word))
        require(v2(forward(word, z) - rp) == v2(F(z) - rp) - len(word),
                "common-center fuel")
        powers_cases += 1
    require(backward(p, F(0)) != 0, "zero singleton is not closed")
    require(backward(p, rq) != rq and backward(qword, rp) != rp,
            "different centers are not closed")
    require({backward(p, s) for s in ()} == set(), "empty bank is vacuous")
    require(forward(p, 21) != forward(p + p, 21),
            "equal center must not erase clock")

    # Only 86 actual macros are needed. This is not an extended orbit census.
    actual, credit, prefix = native_prefix(7, 1, 86)
    require((actual, credit) == (10952205219944621, 47), "actual marked state")
    actual_after = forward(p, actual)
    require(actual_after == 36963692617313097, "actual P output")
    require((v2(F(actual) - rp), v2(actual_after - rp),
             v2(actual_after - rq)) == (5, 2, 9), "actual refill")
    require(v2(F(actual) - rho) == 12, "actual missing-center depth")
    base_signature = observers(bank, actual)
    require(base_signature == observers(bank, rho), "actual bank signature")
    histogram = dict(sorted(Counter(base_signature).items()))
    require(histogram == {0: 1, 1: 30, 2: 22, 4: 6, 5: 3,
                          6: 1, 7: 1, 8: 1, 11: 2}, "sensor histogram")

    # Strong actual-prefix hostile: choose the future cylinder and an odd CRT
    # condition that makes the entire 86-bit pullback an ordinary integer.
    # It preserves exactly 98 initial binary bits, never the original integer7.
    odd_sidecar = 15015
    future_cases = []
    for repetitions in (2, 3, 4, 8, 16, 32):
        future = p + qword * repetitions + "0" * credit
        target_residue = residue(future)
        odd_modulus = 3**len(prefix) * odd_sidecar
        entry, period = crt_two_power_odd(target_residue, len(future),
                                          actual, odd_modulus)
        entry = positive_lift(entry, period, actual + 1)
        source = backward(prefix, entry)
        require(source.denominator == 1 and source >= 3 and source != 7,
                "positive distinct ordinary source")
        source = int(source)
        require((source - 7) % odd_sidecar == 0, "original odd sidecar")
        require(v2(source - 7) == 98, "exact initial binary agreement")
        require(v2(entry - actual) == 12, "actual entry binary limit")
        require(entry % (2**12 * odd_sidecar) == actual % (2**12 * odd_sidecar),
                "compatible actual observation class")
        require(entry % 2**13 != actual % 2**13, "higher precision hostile")
        require(observers(bank, entry) == base_signature, "full actual sensor bank")
        require(v2(F(entry) - rho) == 3 + 6 * repetitions, "exact future prefix")
        after, end_q, first = native_prefix(entry, credit, len(p))
        require(first == p and end_q == credit, "P native credit")
        require(v2(F(after) - rq) == 6 * repetitions, "refilled periodic fuel")
        got_entry, got_credit, got_prefix = native_prefix(source, 1, 86)
        require((got_entry, got_credit, got_prefix) == (entry, credit, prefix),
                "whole exact native history")
        count = halt_word(source, 1, prefix + future)
        require(count == 86 + 3 + 6 * repetitions + 46,
                "blocked symbol is not a macro")
        future_cases.append({"m": repetitions, "entry_precision": 3 + 6 * repetitions,
                             "Q_fuel_after_P": 6 * repetitions,
                             "initial_binary_agreement": 98,
                             "executed_macros": count})

    result = {
        "status": "PASS",
        "scope": "finite controls; no nonhalting or source7 alternative claim",
        "transport_checks": transport,
        "generic_constant_fibre_sources": generic_cases,
        "moduli": list(moduli),
        "equal_center_controls": powers_cases,
        "balanced_words": len(words),
        "distinct_sensors_including_zero": len(bank),
        "separation_H": h,
        "actual_native_macro_index": 86,
        "actual_entry": actual,
        "actual_live_credit": credit,
        "sensor_histogram": histogram,
        "actual_prefix_hostiles": future_cases,
    }
    encoded = json.dumps(result, sort_keys=True, separators=(",", ":"))
    result["payload_sha256"] = hashlib.sha256(encoded.encode()).hexdigest()
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
