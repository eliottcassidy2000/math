"""Exact decimal-window controls; no assertion of the full 86 conjecture.

Default: exponent census 0..100000, 128-digit modular rejection; independent
full census 0..1000; suffix cycles through width7; digit DP through width12;
one modular-only survivor branch through width256. No computations on import.
"""
from fractions import Fraction
from hashlib import sha256
import json
import sys


def integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError("exact integer outside domain")
    return value


def period(width):
    return 4 * 5 ** (integer(width, 1) - 1)


def stable_suffix(phase, width):
    """Canonical eventual residue, even when phase itself is below width."""
    p = period(width)
    integer(phase)
    if phase >= p:
        raise ValueError("phase outside canonical period")
    n = phase + max(0, (width - phase + p - 1) // p) * p
    return pow(2, n, 10 ** width)


def nonzero_suffix(value, width):
    integer(width, 1)
    integer(value)
    if value >= 10 ** width:
        raise ValueError("residue outside decimal window")
    return "0" not in str(value).zfill(width)


def suffix_children(residue, width):
    """All legal nonzero-digit children, with the next parity register."""
    integer(width, 1)
    integer(residue)
    if (residue >= 10 ** width or residue % (1 << width)
            or residue % 5 == 0 or not nonzero_suffix(residue, width)):
        raise ValueError("not a zero-free eventual suffix")
    parity = (residue >> width) % 2
    return tuple(residue + d * 10 ** width for d in range(1, 10)
                 if d % 2 == parity)


def suffix_cycle(width):
    p = period(width)
    r = pow(2, width, 10 ** width)
    residues = [0] * p
    for offset in range(p):
        residues[(width + offset) % p] = r
        r = 2 * r % (10 ** width)
    return tuple(residues)


def digit_dp_count(width):
    """Independent count of zero-free decimal words divisible by 2**width."""
    integer(width, 1)
    modulus = 1 << width
    counts = [0] * modulus
    counts[0] = 1
    for _ in range(width):
        nxt = [0] * modulus
        for residue, count in enumerate(counts):
            if count:
                for digit in range(1, 10):
                    nxt[(10 * residue + digit) % modulus] += count
        counts = nxt
    return counts[0]


def common_decimal_places(n, gap):
    """v_10(2**(n+gap)-2**n); padded residues, not a string-length claim."""
    integer(n)
    integer(gap, 1)
    if gap % 4:
        return 0
    valuation = 0
    while gap % 5 == 0:
        valuation += 1
        gap //= 5
    return min(n, 1 + valuation)


def suffix_correlation(width, gap):
    integer(gap)
    bits = tuple(nonzero_suffix(r, width) for r in suffix_cycle(width))
    p = len(bits)
    return Fraction(sum(bits[e] and bits[(e + gap) % p]
                        for e in range(p)), p)


def sparse_count_bound(limit):
    """Exact rational upper bound; width chosen so no phase repeats below N."""
    integer(limit, 1)
    width = 1
    while period(width) <= limit:
        width += 1
    bound = 4 * width + Fraction(23, 5) ** width
    return width, min(limit, bound.numerator // bound.denominator)


def modular_survivor(depth):
    """Greedy nested zero-free suffixes, never expanding their full powers."""
    integer(depth, 1)
    phase = 1
    rows = []
    for width in range(1, depth + 1):
        if width > 1:
            old_period = period(width - 1)
            candidates = [phase + j * old_period for j in range(5)]
            phase = next(e for e in candidates
                         if nonzero_suffix(stable_suffix(e, width), width))
        p = period(width)
        # This representative has at least width actual decimal digits.
        exponent = phase + max(0, (4 * width - phase + p - 1) // p) * p
        residue = pow(2, exponent, 10 ** width)
        rows.append((width, phase, exponent, residue))
    return tuple(rows)


def scan(limit, width=128):
    """Exact complete decimal test, modular rejection before full expansion."""
    integer(limit)
    integer(width, 1)
    modulus = 10 ** width
    residue = 1
    small = 1
    zero_free = []
    suffix_rejects = full_rejects = full_expansions = 0
    max_expanded_exponent = 0
    max_expanded_digits = 0
    witness_hash = sha256()
    for n in range(limit + 1):
        # Before crossing modulus, unpadded residue is the entire number.
        text = str(residue) if small is not None else str(residue).zfill(width)
        if "0" in text:
            position = text[::-1].index("0") + 1
            suffix_rejects += 1
            witness_hash.update(f"{n}:{position}\n".encode("ascii"))
        else:
            full = str(1 << n)
            full_expansions += 1
            max_expanded_exponent = max(max_expanded_exponent, n)
            max_expanded_digits = max(max_expanded_digits, len(full))
            if "0" in full:
                full_rejects += 1
                witness_hash.update(f"{n}:{full[::-1].index('0') + 1}\n".encode("ascii"))
            else:
                zero_free.append(n)
        residue = (2 * residue) % modulus
        if small is not None:
            small *= 2
            if small >= modulus:
                small = None
    return dict(limit=limit, window=width, zero_free=zero_free,
                suffix_rejects=suffix_rejects, full_rejects=full_rejects,
                full_expansions=full_expansions,
                max_expanded_exponent=max_expanded_exponent,
                max_expanded_digits=max_expanded_digits,
                rejection_witness_sha256=witness_hash.hexdigest())


def main():
    checks = 0

    def check(condition, label):
        nonlocal checks
        checks += 1
        if not condition:
            raise RuntimeError(label)

    def rejects(call):
        try:
            call()
        except ValueError:
            return True
        return False

    # The largest possible full expansion in the declared census has 30103
    # digits. This finite setting is a resource guard, never unlimited.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(40000)

    counts = []
    for width in range(1, 13):
        count = digit_dp_count(width)
        counts.append(count)
        check(4 ** width <= count <= Fraction(23, 5) ** width,
              "all-height rational branching bounds")
        if width <= 7:
            cycle = suffix_cycle(width)
            check(len(set(cycle)) == period(width), "full suffix order")
            survivors = [r for r in cycle if nonzero_suffix(r, width)]
            check(len(survivors) == count, "independent DP/cycle count")
            check(all(r % (1 << width) == 0 and r % 5 != 0 for r in cycle),
                  "CRT cycle image")
            if width <= 6:
                lifted = set()
                for r in survivors:
                    kids = suffix_children(r, width)
                    parity = (r >> width) % 2
                    even = sum((child >> (width + 1)) % 2 == 0 for child in kids)
                    odd = len(kids) - even
                    check((even, odd) == (2, 2) if parity == 0
                          else (even, odd) in ((2, 3), (3, 2)), "parity branches")
                    check(10 * even + 13 * odd <= Fraction(23, 5) * (10 + 3 * parity),
                          "exact rational potential")
                    lifted.update(kids)
                check(lifted == {r for r in suffix_cycle(width + 1)
                                 if nonzero_suffix(r, width + 1)}, "complete lifts")

    for n in range(41):
        for gap in range(1, 241):
            difference = (1 << n) * ((1 << gap) - 1)
            actual = 0
            while difference % 10 == 0:
                difference //= 10
                actual += 1
            check(common_decimal_places(n, gap) == actual, "gap valuation")
    check(suffix_correlation(2, 20) == Fraction(9, 10), "perfect large-lag recurrence")
    check(suffix_correlation(2, 20) != Fraction(9, 10) ** 2, "independence hostile")
    check(stable_suffix(1, 2) == 52, "phase below stable threshold")
    check(not nonzero_suffix(4, 2), "padded leading zero is a real suffix digit")

    branch = modular_survivor(256)
    selected = []
    previous = None
    for width, phase, exponent, residue in branch:
        check(exponent >= 4 * width, "actual digit height")
        check(nonzero_suffix(residue, width), "nested modular survivor")
        check(residue == stable_suffix(phase, width), "phase replay")
        if previous is not None:
            check(phase % period(width - 1) == previous[1], "phase parent")
            check(residue % (10 ** (width - 1)) == previous[3], "digit parent")
        if width & (width - 1) == 0:
            selected.append(dict(width=width, exponent_bits=exponent.bit_length(),
                                 phase=phase,
                                 suffix_sha256=sha256(str(residue).encode()).hexdigest()))
        previous = (width, phase, exponent, residue)

    known = str(pow(2, 103233492954, 10 ** 250)).zfill(250)
    check(known[0] == "0" and "0" not in known[1:], "published 249-digit suffix witness")

    expected = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 18, 19,
                24, 25, 27, 28, 31, 32, 33, 34, 35, 36, 37, 39, 49, 51,
                67, 72, 76, 77, 81, 86]
    exact = scan(100000)
    check(exact["zero_free"] == expected, "finite 100000 census")
    check(sum(exact[key] for key in ("suffix_rejects", "full_rejects"))
          + len(expected) == 100001, "census partition")
    direct = [n for n in range(1001) if "0" not in str(1 << n)]
    check(scan(1000)["zero_free"] == direct, "independent full small census")
    for n in (100, 1000, 100000, 10 ** 12):
        width, bound = sparse_count_bound(n)
        check(period(width) > n and (width == 1 or period(width - 1) <= n),
              "count-bound clock selection")
        if n <= 100000:
            check(sum(e <= n for e in expected if e > 0) <= bound, "finite bound control")

    hostiles = [lambda: period(True), lambda: period(0), lambda: period(2.0),
                lambda: stable_suffix(4, 1), lambda: stable_suffix(-1, 2),
                lambda: suffix_children(50, 2), lambda: suffix_children(52.0, 2),
                lambda: suffix_children(54, 2), lambda: common_decimal_places(1, 0),
                lambda: scan(-1), lambda: nonzero_suffix(100, 2)]
    for hostile in hostiles:
        check(rejects(hostile), "typed or corrupted-state rejection")

    print("FINITE-EXACT decimal-window controls; full 86 conjecture remains OPEN")
    print("suffix_counts_width1_to12=" + json.dumps(counts))
    print("period=4*5^(m-1); 4^m<=Q_m<=((5+sqrt(17))/2)^m<(23/5)^m")
    print("positive_exponent_count_bound=" + json.dumps(
        {str(n): sparse_count_bound(n) for n in (100, 1000, 100000, 10 ** 12)}, sort_keys=True))
    print("width2_lag20_correlation=9/10; independent_product=81/100")
    print("modular_survivor_selected=" + json.dumps(selected, sort_keys=True))
    print("published_exponent103233492954: first_zero_from_right=250; only250digitscomputed")
    print("census=" + json.dumps(exact, sort_keys=True))
    print("independent_full_expansions=1001; max_exponent=1000; max_digits=302")
    print("typed_hostiles=" + str(len(hostiles)))
    print("checks=" + str(checks))


if __name__ == "__main__":
    main()
