"""Construct equal-sum complementary words and replay both physical endpoint counts."""
from dataclasses import dataclass
from itertools import count
from math import gcd

from fermat_receiver import A, D, EDGE, address_height, flipped, join, p, receiver, tagged_data
from local_word_inverse import advance, steer


def complement(word):
    return word.translate(str.maketrans("01", "10"))


def state(value, modulus):
    return value % (2 * modulus), p(value) // 2 % modulus


@dataclass(frozen=True)
class Bank:
    modulus: int
    odd: int
    dyadic: int
    initial: str
    off: str
    on: str
    width: int


def build_bank(modulus):
    dyadic = modulus & -modulus
    odd = modulus // dyadic
    length = next(i for i in count(4, 2) if pow(2, i, odd) == 1 % odd)
    initial = "0" * dyadic.bit_length()
    switches = []
    for u, v in ((2, 4), (2, 4), (1, (1 << (length - 1)) + 1)):
        first, second = (format(z, f"0{length}b") for z in (u, v))
        switches.append((second + first + initial, first + second + initial))
    off = "".join(s[0] for s in switches) * dyadic
    on = "".join(s[1] for s in switches) * dyadic
    loader_width = len(steer(modulus, (0, 0), (0, 0)))
    width = len(initial) + (odd - 1) * len(off) + loader_width
    return Bank(modulus, odd, dyadic, initial, off, on, width)


def controller_word(chart, bank, target, left, right):
    h, q, c, t, _, r, k = chart
    modulus, odd, dyadic = bank.modulus, bank.odd, bank.dyadic
    total = ((left + right + 1) << bank.width) - 1
    assert (left + right + 1) % (2 * dyadic) == 0
    divisor = gcd(h, odd)
    assert (target - k - r - t * total) % divisor == 0

    left_state, right_state = state(left, modulus), state(right, modulus)
    invariant = right_state[1] - left_state[1] + left // 2
    constant = k + r + c * total + 2 * h * invariant + h
    linear = t - c + q * total - h

    # First lift the permutation polynomial one binary digit at a time.
    dyadic_root = 0
    precision = 1
    while precision < dyadic:
        if (constant + linear * dyadic_root - q * dyadic_root**2 - target) % (2 * precision):
            dyadic_root += precision
        precision *= 2
    if dyadic == 1:
        dyadic_root = 1
    assert dyadic_root % 2 == 1

    # Join x=S mod odd to the desired odd binary root. One extra binary bit
    # chooses a specific lift for the single-address controller.
    binary_mod = 2 * dyadic
    shift = ((total - dyadic_root) * pow(binary_mod, -1, odd)) % odd
    wanted_x = dyadic_root + binary_mod * shift
    baseline = bank.initial + bank.off * (odd - 1)
    loader = steer(modulus, advance(left_state, baseline, modulus), (wanted_x, 0))
    baseline += loader
    left_read = advance(left_state, baseline, modulus)[1]
    right_read = advance(right_state, complement(baseline), modulus)[1]
    baseline_z = left_read + right_read

    # Each toggled batch changes Z by -dyadic/4 at the odd modulus.
    odd_quotient = odd // divisor
    delta = -dyadic * pow(4, -1, odd)
    error = target - k - r - t * total - 2 * h * baseline_z
    count_on = ((error // divisor) *
                pow(2 * (h // divisor) * delta, -1, odd_quotient)) % odd_quotient
    result = bank.initial + bank.on * count_on + bank.off * (odd - 1 - count_on) + loader
    assert len(result) == bank.width
    return result


def test_prefixes(chart, bank, target):
    h, _, _, t, _, r, k = chart
    divisor = gcd(h, bank.odd)
    w = ((target - k - r + t) *
         pow(t * bank.dyadic * 2 * pow(2, bank.width, divisor), -1, divisor)) % divisor
    # Long binary transitions, and either prefix allowed to vanish.
    s = 2 * bank.dyadic * (w + divisor * 181) - 1
    return ((0, s), (s, 0), (s // 3, s - s // 3))


def main():
    low = join(join(tagged_data(A, 12), tagged_data(flipped(D, *EDGE), 35)), receiver(4))
    high = join(join(tagged_data(A, 12), tagged_data(D, 35)), receiver(4))
    total_cases = 0
    for modulus in (1, 2, 4, 8, 9, 15, 20, 24, 25, 36, 45, 60):
        bank = build_bank(modulus)
        cases = 0
        for target in range(1, modulus + 1 if modulus % 2 == 0 else 2 * modulus, 2):
            for left, right in test_prefixes(low, bank, target):
                word = controller_word(low, bank, target, left, right)
                x = (left << bank.width) | int(word, 2)
                y = (right << bank.width) | int(complement(word), 2)
                assert x + y == ((left + right + 1) << bank.width) - 1
                assert x % 2 and y % 2 == 0
                # Independent exact reads of both words through the integer chart.
                lower, upper = address_height(low, x, y), address_height(high, x, y)
                assert upper == lower + 2 and (lower - target) % modulus == 0
                cases += 1
        total_cases += cases
        print("D", modulus, "common bit width", bank.width, "prefix/phase cases", cases)
    print("literal lower/upper pairs tested", total_cases)


if __name__ == "__main__":
    main()
