"""Exact controls for the Collatz fibre braid; no global convergence claim.

Run from repository root:
  python3 04-computation/experiments/arithmetic_braids_20260917_collatz.py
Only the standard library is used. All checks remain active under python -O.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


def require(predicate, message):
    if not predicate:
        raise RuntimeError(message)


def valuation(n, p):
    require(n != 0, "valuation of zero")
    n, v = abs(n), 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def split(n):
    h = valuation(n, 2)
    return n // (1 << h), h


def step(n):
    return split(3 * n + 1)


def canonical_cycle(values):
    return min(tuple(values[i:] + values[:i]) for i in range(len(values)))


def bounded_words(length, total, prefix=()):
    if length == 1:
        if total >= 1:
            yield prefix + (total,)
        return
    for k in range(1, total - length + 2):
        yield from bounded_words(length - 1, total - k, prefix + (k,))


def main():
    out = {"status": "FINITE-EXACT controls; all-n proofs in companion report"}
    out["universe"] = {
        "row_prefix_length": 35,
        "forward_odd_inputs": [1, 199999],
        "inverse_odd_cores": [1, 1999],
        "inverse_exponents_per_core": 12,
        "triadic_levels": [1, 9],
        "row_valuation_block_power": 14,
        "finite_word_lengths": [1, 5],
        "finite_word_entries": [1, 4],
        "cycle_word_length_at_most": 8,
        "cycle_total_halving_exponent_at_most": 18,
        "cycle_domain": "all signed odd integer fixed points of bounded exponent words",
    }
    out["rows"] = {
        str(r): [list(split((3 * (6 * j + r) + 1) // 2)) for j in range(35)]
        for r in (1, 3, 5)
    }
    for n in range(1, 200000, 2):
        u, h = split((3 * n + 1) // 2)
        require(u % 3 != 0, "multiple of three reached")
        require((1 << (h + 1)) * u == 3 * n + 1, "address reconstruction")
        require(step(4 * n + 1)[0] == u, "braid changes target")
        require(step(4 * n + 1)[1] == h + 3, "braid height")
    inverse_checks = 0
    for u in range(1, 2000, 2):
        if u % 3 == 0:
            continue
        h0 = 1 if u % 3 == 1 else 0
        last = None
        for j in range(12):
            h = h0 + 2 * j
            n = ((1 << (h + 1)) * u - 1) // 3
            require(split((3 * n + 1) // 2) == (u, h), "inverse formula")
            require(last is None or n == 4 * last + 1, "inverse recurrence")
            last = n
            inverse_checks += 1
    out["inverse_checks"] = inverse_checks
    periods = []
    for s in range(1, 10):
        modulus = 2 * 3**s
        n, seen = 1, set()
        for _ in range(3**s):
            require(n not in seen, "short triadic period")
            seen.add(n)
            n = (4 * n + 1) % modulus
        require(n == 1 and seen == set(range(1, modulus, 2)), "not full cycle")
        periods.append({"modulus": modulus, "period": len(seen)})
    out["triadic_periods"] = periods
    for n in range(1, 200, 2):
        for t in range(1, 101):
            delta = ((4**t - 1) * (3 * n + 1)) // 3
            require(valuation(delta, 3) == valuation(t, 3), "3-adic isometry")
    H = 14
    histograms = {}
    for r in (1, 3, 5):
        hist = Counter(min(split((3 * (6 * j + r) + 1) // 2)[1], H)
                       for j in range(1 << H))
        expected = {h: 1 << (H - h - 1) for h in range(H)} | {H: 1}
        require(dict(hist) == expected, "row valuation law")
        histograms[str(r)] = dict(sorted(hist.items()))
    out["row_height_histograms_tail_capped_at_14"] = histograms
    word_checks = 0
    for length in range(1, 6):
        for word in product(range(1, 5), repeat=length):
            K, B = 0, 0
            for k in word:
                B = 3 * B + (1 << K)
                K += k
            modulus = 1 << (K + 1)
            residue = (((1 << K) - B) * pow(3**length, -1, modulus)) % modulus
            require(residue % 2 == 1, "word source parity")
            for r in (1, 3, 5):
                n = next(residue + t * modulus for t in range(3)
                         if (residue + t * modulus) % 6 == r)
                start = n
                actual = []
                for _ in word:
                    n, k = step(n)
                    actual.append(k)
                require(tuple(actual) == word, "word realization")
                require((1 << K) * n == 3**length * start + B, "affine cocycle")
                word_checks += 1
    out["finite_word_crt_checks"] = word_checks
    out["growing_prefixes"] = []
    for length in (1, 2, 3, 10, 30, 100):
        n0 = (1 << (length + 1)) - 1
        n = n0
        for _ in range(length):
            n, k = step(n)
            require(k == 1, "growth prefix")
        require(n == 2 * 3**length - 1, "growth endpoint")
        out["growing_prefixes"].append({"length": length, "start": n0,
                                       "end": n, "ratio": str(Fraction(n, n0))})
    periodic_hostiles = 0
    for modulus in range(1, 61):
        for length in range(1, 11):
            for q in range(1, 4):
                n0 = (1 << (length + 1)) * modulus * q - 1
                n = n0
                for _ in range(length):
                    n, k = step(n)
                    require(k == 1 and n % modulus == (modulus - 1) % modulus,
                            "periodic-weight hostile")
                require(n > n0, "periodic-weight growth")
                periodic_hostiles += 1
    out["periodic_weight_hostile_checks_M1to60_L1to10_q1to3"] = periodic_hostiles
    cycle_words, cycles = 0, set()
    for length in range(1, 9):
        for K in range(length, 19):
            for word in bounded_words(length, K):
                cycle_words += 1
                B, partial = 0, 0
                for k in word:
                    B = 3 * B + (1 << partial)
                    partial += k
                denominator = (1 << K) - 3**length
                if B % denominator:
                    continue
                n0 = B // denominator
                n, values = n0, []
                for k in word:
                    values.append(n)
                    n, actual_k = step(n)
                    require(actual_k == k, "cycle valuation mismatch")
                require(n == n0, "cycle closure")
                first_return = next(i for i in range(1, len(values) + 1)
                                    if values[i % len(values)] == n0)
                cycles.add(canonical_cycle(values[:first_return]))
    out["cycle_words_checked"] = cycle_words
    out["bounded_signed_cycles"] = sorted(cycles)
    require(cycles == {(1,), (-1,), (-7, -5), (-91, -17, -25, -37, -55, -41, -61)},
            "unexpected bounded cycle census")
    out["controls"] = {
        "positive": "all inverse, residue, CRT and cocycle checks pass",
        "hostile_multiple_of_three_core": "no odd precursor of core3",
        "hostile_long_growth": "100 consecutive exponent-one odd steps",
        "hostile_signed_uniqueness": "three exhibited negative cycles",
        "scope": "bounded cycle census is not a complete classification",
    }
    out["script_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    target = Path(__file__).resolve().parents[2] / "05-knowledge/results/arithmetic_braids_20260917_collatz.json"
    payload = json.dumps(out, indent=2, sort_keys=True) + "\n"
    target.write_text(payload, encoding="utf-8", newline="\n")
    print(json.dumps({"output": str(target), "sha256": sha256(payload.encode()).hexdigest(),
                      "inverse_checks": inverse_checks, "crt_checks": word_checks,
                      "cycle_words_checked": cycle_words, "cycles": sorted(cycles)}, indent=2))


if __name__ == "__main__":
    main()
