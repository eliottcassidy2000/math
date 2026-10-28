#!/usr/bin/env python3
"""Exact, standard-library checks for Session85 arithmetic binary codes.

No network, imported producer code, randomness, or removable assertions.
Run normally and with -O; JSON output must agree.
"""
from collections import Counter
from itertools import combinations
import json
from math import isqrt, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % q for q in range(2, isqrt(n) + 1))


def rank(columns):
    pivots = {}
    for x in columns:
        while x:
            pivot = x.bit_length()
            if pivot not in pivots:
                pivots[pivot] = x
                break
            x ^= pivots[pivot]
    return len(pivots)


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def subsets(values):
    out = [0] * (1 << len(values))
    for mask in range(1, len(out)):
        bit = mask & -mask
        out[mask] = out[mask ^ bit] ^ values[bit.bit_length() - 1]
    return out


def prime_bank(basis, columns, expected_primes):
    rows = []
    for column, p in zip(columns, expected_primes):
        d = prod(q for i, q in enumerate(basis) if column >> i & 1)
        quotient, remainder = divmod(p + 1, d)
        k = isqrt(quotient)
        require(remainder == 0 and k * k == quotient, "square class")
        require(prime(p), "composite proposed prime")
        require(all(not prime(d * j * j - 1) for j in range(1, k)),
                "proposed prime is not least in its fixed square class")
        rows.append(dict(column=column, d=d, k=k, prime=p))
    require(len(set(expected_primes)) == len(expected_primes), "duplicate prime")
    return rows


def arithmetic_check(primes, syndromes):
    """Independent path: sum every divisor via the subset zeta transform.

    Verified distinct primes mean subproducts are exactly the positive divisors.
    Integer isqrt never consults a factor-parity representation of sigma.
    """
    n = len(primes)
    products = [1] * (1 << n)
    for mask in range(1, 1 << n):
        bit = mask & -mask
        products[mask] = products[mask ^ bit] * primes[bit.bit_length() - 1]
    sums = products.copy()
    for i in range(n):
        bit = 1 << i
        for mask in range(1 << n):
            if mask & bit:
                sums[mask] += sums[mask ^ bit]
    for mask, value in enumerate(sums):
        square = isqrt(value) ** 2 == value
        require(square == (syndromes[mask] == 0), "arithmetic/parity mismatch")
    return dict(subsets_checked=len(sums), maximum_source=products[-1],
                maximum_sigma=sums[-1])


def check_code(name, basis, columns, primes, expected_weights,
               expected_circuits, expected_additions=None):
    rows = prime_bank(basis, columns, primes)
    n, r = len(columns), rank(columns)
    universe = 1 << n
    full = universe - 1
    syndrome = subsets(columns)
    code = [mask for mask, s in enumerate(syndrome) if s == 0]
    weights = dict(sorted(Counter(x.bit_count() for x in code).items()))
    require(weights == expected_weights, f"{name}: weight enumerator")
    require(len(code) == 1 << (n - r), "dimension")
    arithmetic = arithmetic_check(primes, syndrome)

    ranks = [rank(columns[i] for i in bits(E)) for E in range(universe)]
    independent = dict(sorted(Counter(E.bit_count() for E in range(universe)
                                      if ranks[E] == E.bit_count()).items()))
    circuits = Counter()
    for k in range(1, r + 2):
        for indices in combinations(range(n), k):
            E = sum(1 << i for i in indices)
            if syndrome[E] == 0 and ranks[E] == k - 1:
                circuits[k] += 1
    require(dict(sorted(circuits.items())) == expected_circuits, "circuit census")

    # Independent exhaustive erasure-fibre check: count supported codewords by
    # zeta transform; union their supports to find exactly the ambiguous bits.
    kernel_count = [int(s == 0) for s in syndrome]
    ambiguous = [E if syndrome[E] == 0 else 0 for E in range(universe)]
    for i in range(n):
        bit = 1 << i
        for E in range(universe):
            if E & bit:
                kernel_count[E] += kernel_count[E ^ bit]
                ambiguous[E] |= ambiguous[E ^ bit]
    for E in range(universe):
        require(kernel_count[E] == 1 << (E.bit_count() - ranks[E]),
                "erasure fibre count")
        for i in bits(E):
            bit = 1 << i
            recovered = not (ambiguous[E] & bit)
            require(recovered == (ranks[E] == ranks[E ^ bit] + 1),
                    "individual recovery/rank criterion")

    dual = [sum(((a & c).bit_count() % 2) << i
                for i, c in enumerate(columns)) for a in range(1 << len(basis))]
    dual = sorted(set(dual))
    require(len(dual) == 1 << r, "dual dimension")
    dual_weights = dict(sorted(Counter(x.bit_count() for x in dual).items()))
    repair_sizes = [min((d.bit_count() - 1 for d in dual if d >> i & 1),
                        default=None) for i in range(n)]
    ddual = min(d.bit_count() for d in dual if d)
    # Every coset has the same uniform projections below dual distance.
    for S in range(universe):
        if S.bit_count() < ddual:
            require(ranks[full ^ S] == r, "low-order marginal uniformity")

    result = dict(name=name, n=n, rank=r, dimension=n-r, rows=rows,
                  arithmetic=arithmetic, weights=weights,
                  circuits=dict(sorted(circuits.items())),
                  independent_erasure_sets=independent,
                  erasure_sets_checked=universe, dual_weights=dual_weights,
                  minimum_repair_sizes=repair_sizes)
    if name.startswith("hamming"):
        require(r >= 2, "marginal/sign theorem requires rank at least2")
        coset1 = [x for x, s in enumerate(syndrome) if s == 1]
        threshold = (1 << (r - 1)) - 1
        check = sum((c & 1) << i for i, c in enumerate(columns))
        selected = full ^ check
        require(selected.bit_count() == threshold, "sharp joint threshold")
        proj0 = Counter(x & selected for x in code)
        proj1 = Counter(x & selected for x in coset1)
        require(proj0 == proj1, "threshold bit marginals differ")
        sign0 = Counter(x.bit_count() % 2 for x in code)
        sign1 = Counter(x.bit_count() % 2 for x in coset1)
        require(sign0 == sign1, "separate sign laws differ")
        joint0 = Counter((x & selected, x.bit_count() % 2) for x in code)
        joint1 = Counter((x & selected, x.bit_count() % 2) for x in coset1)
        require(joint0 != joint1, "missing sharp joint hostile")
        require(min((full ^ d).bit_count() for d in dual if d) == threshold,
                "joint sign threshold")
        # One unique toggle always corrects a nonzero syndrome.
        require(columns == list(range(1, (1 << r))), "full projective bank")
        for x, s in enumerate(syndrome):
            if s:
                require(syndrome[x ^ (1 << (s - 1))] == 0, "toggle repair")
        result["marginal_threshold"] = threshold
        result["joint_sign_equal_through"] = threshold - 1
        result["joint_sign_first_hostile_columns"] = [columns[i] for i in bits(selected)]

    if expected_additions is not None:
        require(syndrome[full] == 0, "full bank not square sigma")
        additions = Counter()
        worst = []
        column_positions = {c: i for i, c in enumerate(columns)}
        for U, target in enumerate(syndrome):
            if target == 0:
                cost = 0
            elif target in column_positions and U >> column_positions[target] & 1:
                cost = 1
            else:
                cost = U.bit_count()
                T = U
                while T:
                    if syndrome[T] == target:
                        cost = min(cost, T.bit_count())
                    T = (T - 1) & U
            require(cost <= r, "addition rank bound")
            require((cost == r) == (U.bit_count() == r and ranks[U] == r),
                    "worst addition iff unused basis")
            additions[cost] += 1
            if cost == r:
                worst.append(U)
        require(dict(sorted(additions.items())) == expected_additions, "addition census")
        result["minimum_legal_addition_costs"] = dict(sorted(additions.items()))
        result["worst_addition_states"] = len(worst)
    return result


def main():
    # Scope boundary for the marginal/sign theorem: rank2 works, rank1 fails.
    small_syndromes = subsets([1, 2, 3])
    small_cosets = [[x for x, s in enumerate(small_syndromes) if s == t]
                    for t in range(4)]
    for coset in small_cosets:
        require(Counter(x.bit_count() % 2 for x in coset) == {0: 1, 1: 1},
                "rank2 separate sign balance")
        for i in range(3):
            require(Counter((x >> i) & 1 for x in coset) == {0: 1, 1: 1},
                    "rank2 one-bit marginal")
    require(Counter((x & 2, x.bit_count() % 2) for x in small_cosets[0])
            != Counter((x & 2, x.bit_count() % 2) for x in small_cosets[1]),
            "rank2 sharp joint-sign hostile")
    rank_one_cosets = [[x for x, s in enumerate(subsets([1])) if s == t]
                       for t in range(2)]
    require(Counter(x.bit_count() % 2 for x in rank_one_cosets[0])
            != Counter(x.bit_count() % 2 for x in rank_one_cosets[1]),
            "rank1 sign-law hostile")
    h7 = check_code(
        "hamming7", [2, 3, 5], list(range(1, 8)),
        [7, 2, 5, 19, 89, 59, 29],
        {0: 1, 3: 7, 4: 7, 7: 1}, {3: 7, 4: 7},
        {0: 16, 1: 56, 2: 28, 3: 28})
    h15 = check_code(
        "hamming15", [2, 3, 5, 7], list(range(1, 16)),
        [7, 2, 5, 19, 89, 59, 29, 251, 13, 83, 41, 139, 5669, 419, 839],
        {0: 1, 3: 35, 4: 105, 5: 168, 6: 280, 7: 435,
         8: 435, 9: 280, 10: 168, 11: 105, 12: 35, 15: 1},
        {3: 35, 4: 105, 5: 168},
        {0: 2048, 1: 15360, 2: 13140, 3: 1380, 4: 840})
    ext16 = check_code(
        "extended16_pure_sigma", [2, 3, 5, 7, 11],
        [1 | (v << 1) for v in range(16)],
        [7, 5, 89, 29, 13, 41, 5669, 839,
         197, 263, 109, 1319, 12473, 461, 769, 2309],
        {0: 1, 4: 140, 6: 448, 8: 870, 10: 448, 12: 140, 16: 1},
        {4: 140, 6: 448})
    # The original fifteen-prime bank + prime3 first gives C15 x F2.
    # Adding even cardinality, not square sigma alone, gives its extension.
    require(prime(3) and isqrt(3 + 1) ** 2 == 3 + 1, "zero-class prime")
    S = subsets(list(range(1, 16)) + [0])
    unguarded = [x for x, s in enumerate(S) if s == 0]
    guarded = [x for x in unguarded if x.bit_count() % 2 == 0]
    require(len(unguarded) == 4096 and len(guarded) == 2048, "extension guard")
    require(dict(sorted(Counter(x.bit_count() for x in guarded).items()))
            == ext16["weights"], "two extended-code realizations")
    # A literal Fano primal circuit toggles a codeword; it is not a parity check.
    require(7 * 2 * 5 == 70 and (7 + 1) * (2 + 1) * (5 + 1) == 144,
            "circuit hostile")
    # Overall parity is a different, one-bit target: seven Fano lines recover
    # it, despite none being a constant-zero dual check on the full code.
    S7 = subsets(list(range(1, 8)))
    C7 = [x for x, s in enumerate(S7) if s == 0]
    lines = [x for x in C7 if x.bit_count() == 3]
    access = Counter()
    for B in range(128):
        fibres = {}
        for x in C7:
            fibres.setdefault(x & B, Counter())[x.bit_count() % 2] += 1
        recovers = all(len(counts) == 1 for counts in fibres.values())
        contains_line = any(line & B == line for line in lines)
        require(recovers == contains_line, "Fano secret access structure")
        if recovers:
            access[B.bit_count()] += 1
        else:
            require(all(counts[0] == counts[1] for counts in fibres.values()),
                    "unauthorized set leaks secret parity")
    require(dict(sorted(access.items())) == {3: 7, 4: 28, 5: 21, 6: 7, 7: 1},
            "Fano parity access counts")
    print(json.dumps(dict(status="PASS", banks=[h7, h15, ext16],
                          marginal_sign_rank2_boundary="PASS",
                          marginal_sign_rank1_hostile="PASS",
                          original_bank_plus3_square_sigma_count=len(unguarded),
                          original_bank_plus3_even_guard_count=len(guarded),
                          fano_parity_authorized_sets=dict(sorted(access.items()))),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
