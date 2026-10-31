"""Independent exact-integer tests; every check survives python -O."""

import hashlib
import itertools
import json
import random
import sys
from pathlib import Path

import beatty_bridge as b
import port_numeration as cg

checks = 0


def check(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def ceiling_certificate(n, result):
    """Independently certify result-1 < beta*n < result, by squaring."""
    if n == 0:
        return result == 0
    upper = 2 * result - 3 * n
    lower = upper - 2
    return lower >= 0 and lower * lower < 5 * n * n < upper * upper


def parent_certificate(n, result):
    """Independently certify beta*result <= n < beta*(result+1)."""
    if n == 0:
        return result == 0
    low = 2 * n - 3 * result
    high = 2 * n - 3 * (result + 1)
    lower_ok = result == 0 or (low > 0 and low * low > 5 * result * result)
    upper_ok = high < 0 or high * high < 5 * (result + 1) ** 2
    return result >= 0 and lower_ok and upper_ok


def state_of_word(word):
    # Last digit other than 1 determines the terminal legality state.
    return word.rstrip("1").endswith("2")


def verify_number(n, full_encoding=True):
    word = cg.encode(n)
    shifted = b.shift(n)
    p, d = b.divide(n)
    g = b.gap(n)
    check(ceiling_certificate(n, shifted), "exact ceiling certificate")
    check(parent_certificate(n, p), "exact reciprocal-floor certificate")
    check(p == 3 * n - shifted, "parent is the inverse companion coordinate")
    check(shifted == cg.value(word + "0"), "canonical weighted shift")
    check(g == (2 if state_of_word(word) else 3), "state equals Beatty gap")
    check(0 <= d < b.gap(p), "division digit allowed")
    check(n == b.shift(p) + d, "division reconstruction")
    check(cg.encode(p) == (word[:-1].lstrip("0") or "0"), "parent drops digit")
    check(d == int(word[-1]), "remainder is last CG digit")
    check(b.shift(shifted) == 3 * shifted - n, "iterated shift recurrence")
    # rho < tau iff 2*P-3*n+1 < (n+1)*sqrt(5).
    difference = (2 * shifted - 3 * n + 1) ** 2 - 5 * (n + 1) ** 2
    check(difference != 0, "strict irrational state threshold")
    check((difference < 0) == (g == 3), "rho threshold agrees with state")
    if full_encoding:
        check(b.encode(n) == word, "division encoding vs independent weight greedy")
        check(b.decode(word) == n, "ceiling-only decode")


def run():
    root = Path(__file__).parent
    old_manifest = json.loads((root / "manifest.json").read_text())
    for name, metadata in old_manifest["files"].items():
        check(hashlib.sha256((root / name).read_bytes()).hexdigest()
              == metadata["sha256"], "original verified file changed: " + name)

    for n in range(100001):
        verify_number(n)

    words = 0
    for length in range(9):
        for chars in itertools.product("012", repeat=length):
            word = "".join(chars)
            if cg.valid(word):
                n = cg.value(word)
                check(cg.valid(word[::-1]), "reversal preserves legal language")
                check(b.decode(word) == n, "legal fixed-word Beatty value")
                check(cg.value(word + "0") == b.shift(n), "all legal words shift")
                words += 1

    rng = random.Random(125_138_2618)
    for _ in range(512):
        verify_number(rng.getrandbits(512))
    for _ in range(32):
        verify_number(rng.getrandbits(4096))

    ws = cg.weights(1002)
    for i in range(1000):
        check(b.shift(ws[i]) == ws[i + 1], "1,3,8 chain is Beatty shift orbit")
        for n in (ws[i] - 1, ws[i], ws[i] + 1):
            verify_number(n, full_encoding=False)

    # Independent interval counting and exact balanced-word bound.
    interval_cases = 0
    for start in range(101):
        for length in range(101):
            actual = sum(b.gap(n) == 3 for n in range(start, start + length))
            count = b.free_count(start, length)
            check(count == actual, "telescoping Free count")
            # ceil(tau*k) = shift(k)-2*k, and tau*k is irrational for k>0.
            allowed = {0} if length == 0 else {
                b.shift(length) - 2 * length - 1,
                b.shift(length) - 2 * length,
            }
            check(count in allowed, "mechanical-word balance")
            interval_cases += 1

    # Noncanonical expressions must be normalized BEFORE applying the theorem.
    check(cg.value("22") == cg.value("100") == 8, "same value before shift")
    check(cg.value("220") == 22 and cg.value("1000") == 21,
          "raw shift and normalization do not commute")
    check(b.shift(8) == 21, "canonical restriction is essential")
    check(b.shift(b.shift(1)) == 8, "iterated shift is 8")
    check(25 < 45 < 49, "squared certificate that 6 < beta squared < 7")
    orbit_cases = 25000
    child_cases = 0
    for n in range(1, orbit_cases + 1):
        first, second = b.shift(n), b.shift(b.shift(n))
        norm = b.orbit_norm(n)
        check(norm > 0, "positive D(C) parameter")
        check(b.orbit_norm(first) == norm, "shift preserves quadratic norm")
        check(n * first + norm == (first - n) ** 2, "D(C) adjacent pair")
        check(n * second + norm == first ** 2, "D(C) outer pair")
        check(first * second + norm == (second - first) ** 2, "D(C) second pair")
        core, depth = b.orbit_core(n)
        word = cg.encode(n)
        root_word = word.rstrip("0")
        check(core == cg.value(root_word), "primitive core strips trailing zeros")
        check(depth == len(word) - len(root_word), "shift depth is zero suffix length")
        check(b.iterate_shift(core, depth) == n, "orbit decomposition reconstructs")
        check(b.divide(core)[1] in (1, 2), "primitive core has no shift predecessor")
        check(b.orbit_core(first) == (core, depth + 1), "exact orbit classification")
        check(b.orbit_norm(core) == norm, "norm attached to primitive core")
        check((norm == 1) == (root_word == "1"), "exact norm-one classification")
        check((norm == 4) == (root_word == "2"), "exact norm-four classification")
        discriminant = 2 * first - 3 * n
        check(discriminant * discriminant - 5 * n * n == 4 * norm,
              "Pell-type norm identity")
        check(discriminant % 2 == n % 2, "quadratic integer parity condition")
        for digit in range(b.gap(n)):
            child = first + digit
            increment = digit * discriminant + digit * digit
            check(b.orbit_norm(child) == norm + increment, "legal-child norm formula")
            check(digit == 0 or increment >= 4, "nonzero later digit increases norm by at least four")
            child_cases += 1
    check(b.orbit_core(0) == (0, 0) and b.orbit_norm(0) == 0, "zero orbit convention")
    check(b.orbit_norm(1) == 1 and b.orbit_norm(2) == 4, "first nonzero digit norms")
    check(b.orbit_norm(5) == b.orbit_norm(7) == 11, "norm collision")
    check(b.orbit_core(5) == (5, 0) and b.orbit_core(7) == (7, 0),
          "same norm does not classify shift orbit")
    for invalid in ("22", "212", "3", "x"):
        try:
            b.decode(invalid)
        except ValueError:
            check(True, "invalid word rejected")
        else:
            raise RuntimeError("invalid word accepted")
    for fn in (b.shift, b.parent, b.encode):
        for invalid, exception in ((-1, ValueError), (True, TypeError), (1.5, TypeError)):
            try:
                fn(invalid)
            except exception:
                check(True, "invalid number rejected")
            else:
                raise RuntimeError("invalid number accepted")

    report = {
        "status": "PASS", "optimized_python": not __debug__,
        "explicit_checks": checks,
        "consecutive_integers": 100001,
        "legal_words_through_length_8": words,
        "random_512_bit_integers": 512,
        "random_4096_bit_integers": 32,
        "large_weight_boundary_cases": 3000,
        "balanced_word_interval_cases": interval_cases,
        "D_C_orbit_classification_cases": orbit_cases,
        "legal_child_norm_cases": child_cases,
        "original_manifest_files_unchanged": len(old_manifest["files"]),
        "arithmetic": "integer isqrt, independent squared-inequality certificates; no floats",
        "boundary": "canonical numeral dynamics only, no graph Hamiltonicity inference",
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    run()
