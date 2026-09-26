"""Finite controls for uniform growing-horizon Collatz multiplicative rank.

The theorem is proved in the matching note; this script audits the exact
carry-bank inequalities and the Terras/binomial tail normalization.
No floating-point comparison certifies any finite arithmetic assertion.
"""
from hashlib import sha256
from math import comb, factorial, log2, sqrt
from pathlib import Path


def compositions(total, length):
    if length == 0:
        if total == 0:
            yield ()
        return
    for k in range(1, total - length + 2):
        for tail in compositions(total - k, length - 1):
            yield (k,) + tail


def carry(word):
    c, b = 0, 1
    for k in word:
        c = 3 * c + b
        b *= 2 ** k
    return c


def prime_factors(value):
    factors = set()
    p = 2
    while p * p <= value:
        if value % p == 0:
            factors.add(p)
            while value % p == 0:
                value //= p
        p += 1
    if value > 1:
        factors.add(value)
    return factors


def factorial_inverse_strict(bound):
    a, fact = 0, 1
    assert bound > 1
    while fact * (a + 1) < bound:
        a += 1
        fact *= a
    assert factorial(a) < bound <= factorial(a + 1)
    return a


def sum_valuations(n, length):
    total = 0
    for _ in range(length):
        value = 3 * n + 1
        k = (value & -value).bit_length() - 1
        total += k
        n = value // 2 ** k
    return total


def main():
    print("SCRIPT_SHA256", sha256(Path(__file__).read_bytes()).hexdigest())
    print("EXACT CARRY UNIVERSE: all positive words lengthN,total<=3N,N=2..5")
    for n in range(2, 6):
        count = 0
        largest_bank = 0
        witness = None
        p_bound = 96 ** (n * (n - 1) // 2)
        a_bound = factorial_inverse_strict(p_bound)
        critical_cutoff = (3 ** n).bit_length() + 1
        radius = n * (n - 1) // 2
        p_critical = 2 ** ((n - 1) * critical_cutoff - radius) * 3 ** radius
        budget_counts = {"2N": 0, "critical": 0}
        for total in range(n, 3 * n + 1):
            for word in compositions(total, n):
                count += 1
                if total <= 2 * n:
                    budget_counts["2N"] += 1
                if total <= critical_cutoff:
                    budget_counts["critical"] += 1
                carries = {}
                for i in range(n):
                    for j in range(i + 1, n):
                        subword = word[i:j]
                        c = carry(subword)
                        r, k = j - i, sum(subword)
                        assert c % 2 == 1 and c % 3 != 0
                        assert c * 2 ** r < 2 ** k * 3 ** r
                        carries[i, j] = c
                for i in range(n):
                    product = 1
                    factors = set()
                    weighted_halvings = 0
                    weighted_distances = 0
                    for j in range(n):
                        if j == i:
                            continue
                        lo, hi = sorted((i, j))
                        c = carries[lo, hi]
                        product *= c
                        factors |= prime_factors(c)
                        weighted_halvings += sum(word[lo:hi])
                        weighted_distances += hi - lo
                    assert weighted_halvings <= (n - 1) * 3 * n
                    assert weighted_distances <= n * (n - 1) // 2
                    assert product < p_bound
                    if total <= 2 * n:
                        assert product < 24 ** radius
                    if total <= critical_cutoff:
                        assert product < p_critical
                    assert factorial(len(factors)) <= product
                    assert len(factors) <= a_bound
                    assert not factors.intersection({2, 3})
                    if len(factors) > largest_bank:
                        largest_bank = len(factors)
                        witness = word, i, sorted(factors)
        assert count == comb(3 * n, n)
        assert budget_counts["2N"] == comb(2 * n, n)
        assert budget_counts["critical"] == comb(critical_cutoff, n)
        print("N", n, "WORDS", count, "MAX_INCIDENT_CARRY_PRIMES", largest_bank,
              "FACTORIAL_BOUND", a_bound, "WITNESS", witness)
        print("CONDITIONED_BUDGETS", budget_counts, "CRITICAL_CUTOFF", critical_cutoff)

    print("EXACT TERRAS TAIL UNIVERSE: all odd residues modulo2^(L+1),N=1..5,L=2N,3N")
    for n in range(1, 6):
        for length in [2 * n, 3 * n]:
            modulus = 2 ** (length + 1)
            observed = sum(sum_valuations(source, n) > length
                           for source in range(1, modulus, 2))
            predicted = sum(comb(length, j) for j in range(n))
            assert observed == predicted
            if length == 3 * n:
                assert predicted * 4 ** n <= 27 ** n
            print("N", n, "L", length, "ODD_RESIDUES", modulus // 2,
                  "TAIL_OBSERVED", observed, "BINOMIAL", predicted)

    print("FACTORIAL EXCEPTION BOUNDS: log2(E_N) and log2(P_N), diagnostic decimals")
    for n in [4, 8, 16, 32, 64, 128]:
        p_bound = 96 ** (n * (n - 1) // 2)
        a_bound = factorial_inverse_strict(p_bound)
        exponent = 16 * a_bound + 40
        exception_bound = comb(3 * n, n) * comb(n, 2) * 2 ** exponent
        print("N", n, "a_N", a_bound, "log2_E", round(log2(exception_bound), 6),
              "log2_P", round(log2(p_bound), 6))
    alpha = log2(3)
    print("ASYMPTOTIC_DIAGNOSTICS_ONLY")
    print("sqrtloglog_horizon_constant", 1 / 8)
    print("exception_power_beta", (5 + alpha) / 8)
    print("largest_constant_from_this_bound", 1 / sqrt(8 * (5 + alpha)))
    print("tail_ratio", "27/32", "tail_rate_bits", 5 - 3 * alpha)
    print("conditioned_power_2N", (3 + alpha) / 8)
    print("conditioned_power_critical_cutoff", (3 * alpha - 1) / 8)
    print("UNBOUNDED-HEIGHT NONINJECTIVE HOSTILE: n_a=(4^a-1)/3 ->1->1")
    for a in [2, 3, 4, 8, 16, 64, 256]:
        source = (4 ** a - 1) // 3
        assert 3 * source + 1 == 4 ** a
        assert sum_valuations(source, 3) == 2 * a + 4
        print("a", a, "SOURCE_BITS", source.bit_length(),
              "FIRST3_VALUATION_SUM", 2 * a + 4,
              "DUPLICATE_PAIR", "(-3,4)")
    print("SCOPE: density one uniformly over source intervals; no single-orbit conclusion.")


if __name__ == "__main__":
    main()
