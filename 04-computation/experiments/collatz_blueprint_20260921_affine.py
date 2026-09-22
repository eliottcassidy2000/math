"""Exact replay of guarded affine inverse-Collatz words; no convergence claim.

Universe: all 2,047 chronological D/E words of lengths 0..10. For each,
test every positive source 1..64 plus three representatives of its legal
residue class. Exhaust all residues for lengths <=6. Assertions remain
enabled under -O because checks use an explicit exception.
"""
from fractions import Fraction as Q
from itertools import product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


IDENTITY = ((1, 0), (0, 1))
D = ((2, 0), (0, 1))
E = ((2, -1), (0, 3))
M3 = ((4, -1), (0, 3))


def matrix(word):
    a = IDENTITY
    for letter in word:
        a = matmul(D if letter == "D" else E, a)
    return a


def carry(word):
    b = m = 0
    for letter in word:
        b *= 2
        if letter == "E":
            b += 3 ** m
            m += 1
    return len(word), m, b


def legal_eval(word, n):
    require(n > 0 and isinstance(n, int), "positive integer domain")
    for letter in word:
        if letter == "D":
            n *= 2
        else:
            if n % 3 != 2:
                return None
            n = (2 * n - 1) // 3
        require(n > 0, "positive intermediate")
    return n


def rational_eval(word, n):
    n = Q(n)
    for letter in word:
        n = 2 * n if letter == "D" else (2 * n - 1) / 3
    return n


def spectral(a):
    tr = a[0][0] + a[1][1]
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    return Q(tr * tr, det)


def decode(word, endpoint):
    backwards = []
    for _ in word:
        if endpoint % 2 == 0:
            backwards.append("D")
            endpoint //= 2
        else:
            backwards.append("E")
            endpoint = (3 * endpoint + 1) // 2
    return "".join(reversed(backwards)), endpoint


def main():
    require(matmul(E, D) == M3, "M3 redundancy")
    invariants = [spectral(D), spectral(matmul(E, D)),
                  spectral(matmul(M3, D))]
    require(invariants == [Q(9, 2), Q(49, 12), Q(121, 24)],
            "normalized traces")
    records = {}
    word_count = source_checks = residue_checks = 0
    cycle_fixed_words = []
    for r in range(11):
        endpoint_residues = set()
        for letters in product("DE", repeat=r):
            word = "".join(letters)
            rr, m, b = carry(word)
            require(rr == r, "length")
            a = matrix(word)
            require(a == ((2 ** r, -b), (0, 3 ** m)), "matrix carry")
            require(a not in records, "free semigroup collision")
            records[a] = word
            modulus = 3 ** m
            residue = (b * pow(2 ** r, -1, modulus)) % modulus
            first = residue or modulus
            endpoint_first = legal_eval(word, first)
            require(1 <= endpoint_first <= 2 ** r, "least positive endpoint")
            endpoint_residue = (-b * pow(modulus, -1, 2 ** r)) % (2 ** r)
            require(endpoint_first % (2 ** r) == endpoint_residue,
                    "endpoint congruence")
            require(endpoint_residue not in endpoint_residues,
                    "disjoint endpoint classes")
            endpoint_residues.add(endpoint_residue)
            for j in range(3):
                require(legal_eval(word, first + j * modulus) ==
                        endpoint_first + j * 2 ** r, "AP correspondence")
            sources = sorted(set(range(1, 65)) |
                             {first + j * modulus for j in range(3)})
            for n in sources:
                expected_legal = (2 ** r * n - b) % modulus == 0
                actual = legal_eval(word, n)
                require((actual is not None) == expected_legal, "guard iff")
                require(rational_eval(word, n) == Q(2 ** r * n - b, modulus),
                        "independent rational evaluation")
                if expected_legal:
                    require(actual == (2 ** r * n - b) // modulus,
                            "terminal value")
                    require(decode(word, actual) == (word, n), "parity decoding")
                source_checks += 1
            if r <= 6:
                valid = [n for n in range(1, modulus + 1)
                         if legal_eval(word, n) is not None]
                require(valid == [first], "one complete source residue")
                residue_checks += modulus
            if r:
                require(spectral(a) > 4, "all nonempty positive words hyperbolic")
                delta = 2 ** r - modulus
                require(delta != 0, "unique factorization")
                fixed = Q(b, delta)
                criterion = m > 0 and delta > 0 and b % delta == 0
                require(criterion == (fixed > 0 and fixed.denominator == 1),
                        "cycle criterion")
                if criterion:
                    n = int(fixed)
                    require(legal_eval(word, n) == n, "cycle legality")
                    cycle_fixed_words.append({"word": word, "fixed": n})
            word_count += 1
        require(endpoint_residues == set(range(2 ** r)), "endpoint partition")

    good, hostile = matrix("DEDE"), matrix("DDEE")
    require(good == ((16, -7), (0, 9)), "integer cycle word")
    require(hostile == ((16, -5), (0, 9)), "rational-only cycle word")
    require(spectral(good) == spectral(hostile) == Q(625, 144), "same spectrum")
    require(legal_eval("DEDE", 1) == 1, "positive cycle control")
    require(rational_eval("DDEE", Q(5, 7)) == Q(5, 7), "rational hostile")
    require((3 * Q(5, 7) + 1) / 2 == Q(11, 7), "rational odd cycle first")
    require((3 * Q(11, 7) + 1) / 8 == Q(5, 7), "rational odd cycle second")

    d = lambda x: 2 * x
    di = lambda x: x / 2
    e = lambda x: (2 * x - 1) / 3
    ei = lambda x: (3 * x + 1) / 2
    for n in range(-10, 11):
        require(d(e(di(ei(Q(n))))) == Q(n) - Q(1, 3), "commutator")
        for j in range(13):
            require((2 ** j * Q(n) - Q(1, 3)) / 2 ** j ==
                    Q(n) - Q(1, 3 * 2 ** j), "small translations")
    for n in range(1, 100):
        if n % 3 == 2:
            require(legal_eval("E", n + 18) % 6 == legal_eval("E", n) % 6,
                    "mod18 sufficiency E")
        if n % 3 == 1:
            require(legal_eval("DE", n + 18) % 6 == legal_eval("DE", n) % 6,
                    "mod18 sufficiency M3")
    require(legal_eval("E", 5) == 3 and legal_eval("E", 11) == 7,
            "odd mod6 hostile E")
    require(legal_eval("DE", 1) == 1 and legal_eval("DE", 7) == 9,
            "odd mod6 hostile M3")
    require(all(legal_eval("DE", n) is None for n in range(3, 100, 6)),
            "M3 undefined on 3mod6")
    half_m3 = tuple(tuple(Q(x, 2) for x in row) for row in M3)
    require(spectral(half_m3) == spectral(M3), "projective scale invariance")
    require(half_m3[0][0] * half_m3[1][1] == 3, "arbitrary det3 normalization")

    print(json.dumps({
        "status": "FINITE-EXACT controls of elementary proofs; no convergence claim",
        "word_universe": {"alphabet": "DE", "length_min": 0, "length_max": 10,
                          "word_count": word_count, "distinct_matrices": len(records)},
        "source_checks": source_checks,
        "exhaustive_residue_checks_through_length6": residue_checks,
        "endpoint_partition": "all 2^r endpoint residues for every r=0..10",
        "normalized_trace_squared": [str(x) for x in invariants],
        "same_spectrum_hostile": {"invariant": "625/144", "DEDE_fixed": "1",
                                  "DDEE_fixed": "5/7"},
        "positive_integer_fixed_words_in_declared_universe": cycle_fixed_words,
        "mod6_hostiles": {"E": [[5, 3], [11, 7]], "M3": [[1, 1], [7, 9]]},
        "commutator": "x-1/3", "conjugate_translation": "-1/(3*2^j)",
        "all_checks_passed": True,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
