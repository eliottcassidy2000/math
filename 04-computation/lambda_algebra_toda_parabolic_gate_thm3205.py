#!/usr/bin/env python3
"""Exact controls for THM-3205's odd-primary lambda-algebra Toda gate family.

Self-contained: integer and F_p arithmetic only, no floating point, no
randomness, no imported executable.  Every gate is an explicit ``require`` so
that ordinary and ``-O`` replay are byte-identical.

Reference statements replayed here are the displayed relations, differential,
and Lemma 3 of Ivanov--Mikhailov--Wu, *On nontriviality of homotopy groups of
spheres*, arXiv:1506.00952v1 (Homology Homotopy Appl. 18 (2016) 337-344).
"""

import sys
from itertools import product

LAM = 0   # lambda_i, i >= 1, internal degree 2(p-1)i - 1  (odd)
MU = 1    # mu_j,     j >= 0, internal degree 2(p-1)j      (even)

sys.setrecursionlimit(100000)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# ------------------------------------------------- structure constants


def generalized_binomial(top, low):
    """C(top, low) for any integer top and low >= 0, as an exact integer."""

    if low < 0:
        return 0
    numerator = 1
    for step in range(low):
        numerator *= top - step
    denominator = 1
    for step in range(1, low + 1):
        denominator *= step
    return numerator // denominator


def coeff_a(prime, k, j):
    sign = -1 if j % 2 == 0 else 1
    return sign * generalized_binomial((prime - 1) * (k - j) - 1, j) % prime


def coeff_b(prime, k, j):
    sign = 1 if j % 2 == 0 else -1
    return sign * generalized_binomial((prime - 1) * (k - j), j) % prime


def bound_n(prime, k):
    return k - (k + 1) // prime


def bound_n_prime(prime, k):
    return k - k // prime


def is_generator(letter):
    kind, index = letter
    return index >= 1 if kind == LAM else index >= 0


def admissible_pair(prime, first, second):
    kind, index = first
    return (second[1] <= prime * index - 1 if kind == LAM
            else second[1] <= prime * index)


def rewrite_pair(prime, first, second):
    """Expand one inadmissible adjacent pair into (coeff, letter, letter)."""

    kind1, i = first
    kind2, m = second
    out = []
    if kind1 == LAM:
        k = m - prime * i
        for j in range(bound_n(prime, k) + 1):
            if j == k:
                continue                        # the left-hand side itself
            out.append((coeff_a(prime, k, j),
                        (LAM, i + k - j), (kind2, prime * i + j)))
        if kind2 == MU:
            for j in range(bound_n_prime(prime, k) + 1):
                out.append((coeff_b(prime, k, j),
                            (MU, i + k - j), (LAM, prime * i + j)))
    else:
        k = m - prime * i - 1
        for j in range(bound_n(prime, k) + 1):
            if j == k:
                continue
            out.append((coeff_a(prime, k, j),
                        (MU, i + k - j), (kind2, prime * i + j + 1)))
    return [(c, x, y) for (c, x, y) in out
            if c % prime and is_generator(x) and is_generator(y)]


class LambdaAlgebra:
    """Odd-primary lambda algebra with admissible normal form."""

    def __init__(self, prime, rightmost=False):
        self.p = prime
        self.rightmost = rightmost
        self.memo = {}
        self.dgen_memo = {}

    def normal_form(self, word):
        word = tuple(word)
        cached = self.memo.get(word)
        if cached is not None:
            return cached
        if any(not is_generator(letter) for letter in word):
            self.memo[word] = {}
            return {}
        spots = range(len(word) - 2, -1, -1) if self.rightmost \
            else range(len(word) - 1)
        position = None
        for t in spots:
            if not admissible_pair(self.p, word[t], word[t + 1]):
                position = t
                break
        if position is None:
            self.memo[word] = {word: 1}
            return self.memo[word]
        total = {}
        for c, x, y in rewrite_pair(self.p, word[position], word[position + 1]):
            replaced = word[:position] + (x, y) + word[position + 2:]
            for w, cc in self.normal_form(replaced).items():
                total[w] = (total.get(w, 0) + c * cc) % self.p
        total = {w: c for w, c in total.items() if c}
        self.memo[word] = total
        return total

    def differential_of_generator(self, letter):
        cached = self.dgen_memo.get(letter)
        if cached is not None:
            return cached
        p = self.p
        kind, k = letter
        terms = []
        if kind == LAM:
            for j in range(1, bound_n(p, k) + 1):
                terms.append((coeff_a(p, k, j), ((LAM, k - j), (LAM, j))))
        else:
            for j in range(0, bound_n(p, k) + 1):
                terms.append((coeff_a(p, k, j), ((LAM, k - j), (MU, j))))
            for j in range(1, bound_n_prime(p, k) + 1):
                terms.append((coeff_b(p, k, j), ((MU, k - j), (LAM, j))))
        out = {}
        for c, w in terms:
            if c % p == 0:
                continue
            for ww, cc in self.normal_form(w).items():
                out[ww] = (out.get(ww, 0) + c * cc) % p
        out = {w: c for w, c in out.items() if c}
        self.dgen_memo[letter] = out
        return out

    def differential(self, element, koszul=True):
        """Signed Leibniz differential of a normalized element."""

        out = {}
        for word, c in element.items():
            sign = 1
            for t, letter in enumerate(word):
                for w2, c2 in self.differential_of_generator(letter).items():
                    replaced = word[:t] + w2 + word[t + 1:]
                    for ww, cc in self.normal_form(replaced).items():
                        out[ww] = (out.get(ww, 0)
                                   + sign * c * c2 * cc) % self.p
                if koszul and letter[0] == LAM:
                    sign = -sign
        return {w: c for w, c in out.items() if c}


def show(element):
    if not element:
        return "0"
    parts = []
    for word, c in sorted(element.items()):
        letters = "".join(("l" if kind == LAM else "m") + str(i)
                          for kind, i in word)
        parts.append(str(c) + "*" + letters)
    return " + ".join(parts)


def rank_mod_p(rows, ncols, prime):
    matrix = [list(row) for row in rows]
    rank, column = 0, 0
    nrows = len(matrix)
    while column < ncols and rank < nrows:
        pivot = next((r for r in range(rank, nrows)
                      if matrix[r][column] % prime), None)
        if pivot is None:
            column += 1
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], prime - 2, prime)
        matrix[rank] = [(x * inverse) % prime for x in matrix[rank]]
        for r in range(nrows):
            if r != rank and matrix[r][column] % prime:
                factor = matrix[r][column]
                matrix[r] = [(x - factor * y) % prime
                             for x, y in zip(matrix[r], matrix[rank])]
        rank += 1
        column += 1
    return rank


def stratum(prime, top, length, index_sum, lambdas):
    """Admissible Lambda(top)*lambda words: fixed length, index sum, #lambdas.

    Every letter of a word in the left ideal Lambda*lambda has index >= 1:
    a mu_0 would force the next index to be <= 0, which cascades to an
    impossible lambda_0 before the word can end in a lambda.
    """

    out = []

    def build(word, left, remaining, lambdas_left):
        if left == 0:
            if remaining == 0 and lambdas_left == 0:
                out.append(tuple(word))
            return
        if word:
            kind, index = word[-1]
            ceiling = prime * index - 1 if kind == LAM else prime * index
        else:
            ceiling = top
        ceiling = min(ceiling, remaining - (left - 1))
        for kind in (LAM, MU):
            if kind == LAM and lambdas_left == 0:
                continue
            if left == 1 and kind != LAM:
                continue
            for index in range(1, ceiling + 1):
                build(word + [(kind, index)], left - 1, remaining - index,
                      lambdas_left - (1 if kind == LAM else 0))

    build([], length, index_sum, lambdas)
    return sorted(out)


def differential_matrix(algebra, source, target):
    lookup = {word: j for j, word in enumerate(target)}
    rows = []
    for word in source:
        row = [0] * len(target)
        for image, c in algebra.differential({word: 1}).items():
            require(image in lookup, "target stratum is complete")
            row[lookup[image]] = c
        rows.append(row)
    return rows


def smith_invariant_factors(matrix):
    a = [list(row) for row in matrix]
    rows, cols = len(a), len(a[0])
    factors, top = [], 0
    while top < rows and top < cols:
        best, pivot = None, None
        for r in range(top, rows):
            for c in range(top, cols):
                if a[r][c] and (best is None or abs(a[r][c]) < best):
                    best, pivot = abs(a[r][c]), (r, c)
        if pivot is None:
            break
        while True:
            r0, c0 = pivot
            a[top], a[r0] = a[r0], a[top]
            for row in a:
                row[top], row[c0] = row[c0], row[top]
            moved = False
            for r in range(top + 1, rows):
                if a[r][top]:
                    q = a[r][top] // a[top][top]
                    for c in range(top, cols):
                        a[r][c] -= q * a[top][c]
                    if a[r][top]:
                        moved = True
            for c in range(top + 1, cols):
                if a[top][c]:
                    q = a[top][c] // a[top][top]
                    for r in range(top, rows):
                        a[r][c] -= q * a[r][top]
                    if a[top][c]:
                        moved = True
            if not moved:
                break
            best, pivot = None, None
            for r in range(top, rows):
                for c in range(top, cols):
                    if a[r][c] and (best is None or abs(a[r][c]) < best):
                        best, pivot = abs(a[r][c]), (r, c)
        factors.append(abs(a[top][top]))
        top += 1
    for i in range(1, len(factors)):
        x, y = factors[i - 1], factors[i]
        while y:
            x, y = y, x % y
        product_value = factors[i - 1] * factors[i]
        factors[i - 1], factors[i] = x, product_value // x
    return factors


PRIMES = (3, 5, 7, 11, 13)

# ------------------------------- 1.  convention audit against the source

CONVENTION_ROWS = []
for prime in PRIMES:
    algebra = LambdaAlgebra(prime)

    def word_of(*letters):
        return algebra.normal_form(letters)

    minus_one = prime - 1
    minus_two = prime - 2
    require(word_of((MU, 0), (LAM, 1)) == {}, "mu_0 lambda_1 = 0")
    require(word_of((MU, 0), (MU, 1)) == {}, "mu_0 mu_1 = 0")
    require(word_of((MU, 0), (LAM, 2)) == {((MU, 1), (LAM, 1)): minus_one},
            "mu_0 lambda_2 = -mu_1 lambda_1")
    require(word_of((MU, 0), (MU, 2)) == {((MU, 1), (MU, 1)): minus_one},
            "mu_0 mu_2 = -mu_1 mu_1")
    require(algebra.differential_of_generator((LAM, 1)) == {},
            "d lambda_1 = 0")
    require(algebra.differential_of_generator((MU, 1))
            == {((LAM, 1), (MU, 0)): minus_one}, "d mu_1 = -lambda_1 mu_0")
    require(algebra.differential_of_generator((LAM, 2))
            == {((LAM, 1), (LAM, 1)): minus_two % prime},
            "d lambda_2 = -2 lambda_1^2")
    require(algebra.differential_of_generator((MU, 2))
            == {((LAM, 2), (MU, 0)): minus_one,
                ((LAM, 1), (MU, 1)): minus_two % prime,
                ((MU, 1), (LAM, 1)): 1},
            "d mu_2 = -lambda_2 mu_0 - 2 lambda_1 mu_1 + mu_1 lambda_1")
    require(algebra.differential(word_of((MU, 2), (LAM, 1)))
            == {((LAM, 1), (MU, 1), (LAM, 1)): minus_two % prime,
                ((MU, 1), (LAM, 1), (LAM, 1)): 1},
            "d(mu_2 lambda_1)")
    require(algebra.differential(word_of((MU, 1), (LAM, 2)))
            == {((LAM, 1), (MU, 1), (LAM, 1)): 1,
                ((MU, 1), (LAM, 1), (LAM, 1)): minus_two % prime},
            "d(mu_1 lambda_2)")
    require(algebra.differential(word_of((MU, 1), (MU, 2), (LAM, 1)))
            == {((LAM, 1), (MU, 1), (MU, 1), (LAM, 1)): 1,
                ((MU, 1), (LAM, 1), (MU, 1), (LAM, 1)): minus_two % prime,
                ((MU, 1), (MU, 1), (LAM, 1), (LAM, 1)): 1},
            "d(mu_1 mu_2 lambda_1)")
    CONVENTION_ROWS.append(prime)
require(len(CONVENTION_ROWS) == len(PRIMES), "convention audit ran")

# ------------------------------- 2.  d^2 = 0 and the Koszul sign requirement

SQUARE_ZERO_TOTAL = 0
UNSIGNED_FAILURES = 0
LETTERS = [(LAM, i) for i in range(1, 6)] + [(MU, i) for i in range(0, 6)]
for prime in (3, 5, 7):
    algebra = LambdaAlgebra(prime)
    for length in (1, 2, 3):
        for raw in product(LETTERS, repeat=length):
            element = algebra.normal_form(raw)
            if not element:
                continue
            SQUARE_ZERO_TOTAL += 1
            require(algebra.differential(algebra.differential(element)) == {},
                    "d^2 = 0 with the Koszul sign")
            unsigned = algebra.differential(
                algebra.differential(element, koszul=False), koszul=False)
            if unsigned:
                UNSIGNED_FAILURES += 1
require(SQUARE_ZERO_TOTAL > 3000, "d^2 control is non-vacuous")
require(UNSIGNED_FAILURES > 0,
        "unsigned Leibniz must fail: the Koszul sign is load-bearing")

# ------------------------------- 3.  confluence of the rewriting

CONFLUENCE_TOTAL = 0
for prime in (3, 5, 7):
    left_first = LambdaAlgebra(prime, rightmost=False)
    right_first = LambdaAlgebra(prime, rightmost=True)
    for length in (2, 3):
        for raw in product(LETTERS, repeat=length):
            require(left_first.normal_form(raw) == right_first.normal_form(raw),
                    "leftmost and rightmost rewriting agree")
            CONFLUENCE_TOTAL += 1
require(CONFLUENCE_TOTAL > 2000, "confluence control is non-vacuous")

# ------------------------------- 4.  the excess-one gate is tridiagonal

def toda_source(k, i):
    """u_0 = mu_1^k lambda_2 and u_i = mu_1^(k-i) mu_2 mu_1^(i-1) lambda_1."""

    if i == 0:
        return ((MU, 1),) * k + ((LAM, 2),)
    return (((MU, 1),) * (k - i) + ((MU, 2),)
            + ((MU, 1),) * (i - 1) + ((LAM, 1),))


def toda_target(k, i):
    return ((MU, 1),) * (k - i) + ((LAM, 1),) + ((MU, 1),) * i + ((LAM, 1),)


GATE_ROWS = []
for prime in PRIMES:
    algebra = LambdaAlgebra(prime)
    for k in range(1, 10):
        sources = [toda_source(k, i) for i in range(k + 1)]
        targets = [toda_target(k, i) for i in range(k + 1)]
        require(sorted(sources) == stratum(prime, 2, k + 1, k + 2, 1),
                "the u-family is the whole excess-one q=1 stratum")
        require(sorted(targets) == stratum(prime, 2, k + 2, k + 2, 2),
                "the v-family is the whole excess-zero q=2 stratum")
        lookup = {word: j for j, word in enumerate(targets)}
        matrix = [[0] * (k + 1) for _ in range(k + 1)]
        for column, word in enumerate(sources):
            for image, c in algebra.differential({word: 1}).items():
                require(image in lookup, "differential stays in the v-span")
                matrix[lookup[image]][column] = c
        expected = [[(-2) % prime if r == c else
                     1 % prime if abs(r - c) == 1 else 0
                     for c in range(k + 1)] for r in range(k + 1)]
        require(matrix == expected, "excess-one matrix is tridiag(1,-2,1)")
        rank = rank_mod_p([list(row) for row in matrix], k + 1, prime)
        kernel = k + 1 - rank
        gate = 0 if (k + 2) % prime else 1
        require(kernel == gate,
                "kernel is one-dimensional exactly when p divides k+2")
        GATE_ROWS.append((prime, k, kernel))

# integral refinement: the cokernel of the lift is cyclic of order k+2
SMITH_ROWS = []
for k in range(1, 14):
    size = k + 1
    integral = [[-2 if r == c else 1 if abs(r - c) == 1 else 0
                 for c in range(size)] for r in range(size)]
    factors = smith_invariant_factors(integral)
    require(factors == [1] * (size - 1) + [k + 2],
            "integral cokernel is cyclic Z/(k+2)")
    SMITH_ROWS.append((k, factors[-1]))

# ------------------------------- 5.  the explicit kernel generator

GENERATOR_ROWS = []
SPURIOUS = []
for prime in PRIMES:
    algebra = LambdaAlgebra(prime)
    for k in range(1, 22):
        combination = {}
        for i in range(k + 1):
            for word, c in algebra.normal_form(toda_source(k, i)).items():
                combination[word] = (combination.get(word, 0)
                                     + (i + 1) * c) % prime
        combination = {w: c for w, c in combination.items() if c}
        is_cycle = not algebra.differential(combination)
        if (k + 2) % prime == 0:
            require(is_cycle, "z_k is a cycle when p divides k+2")
            # The class always has a nonzero mu_2 mu_1^(k-1) lambda_1 term:
            # its coefficient is k+1 = -1 mod p.  This is exactly the input
            # the source's Proposition 2 assumes about the E^2 page.
            leading = algebra.normal_form(toda_source(k, k))
            require(len(leading) == 1, "u_k is already admissible")
            leading_word = next(iter(leading))
            require(combination.get(leading_word, 0) == prime - 1,
                    "leading mu_2 coefficient of z_k is -1 mod p")
            GENERATOR_ROWS.append((prime, k))
        else:
            require(not is_cycle, "z_k is not a cycle otherwise")
            SPURIOUS.append((prime, k)) if is_cycle else None
require(not SPURIOUS, "hostile control found no spurious cycle")
require(len(GENERATOR_ROWS) >= 12, "generator bank is non-vacuous")

# ------------------------------- 6.  higher excess: the gate does not repeat

EXCESS_ROWS = []
for excess in (2, 3, 4):
    for prime in PRIMES:
        algebra = LambdaAlgebra(prime)
        widths = []
        for delta in range(1, 11):
            index_sum = delta + 1 + excess
            source = stratum(prime, 2, delta + 1, index_sum, 1)
            if not source or len(source) > 120:
                continue
            target = stratum(prime, 2, delta + 2, index_sum, 2)
            require(len(target) > len(source),
                    "higher excess strata are strictly wide")
            rank = rank_mod_p(
                differential_matrix(algebra, source, target),
                len(target), prime)
            require(rank == len(source),
                    "higher-excess differential is injective")
            widths.append((delta, len(source), len(target)))
        EXCESS_ROWS.append((excess, prime, len(widths)))
require(all(count > 0 for _, _, count in EXCESS_ROWS),
        "higher-excess sweep is non-vacuous")

# positive control for the sweep machinery: excess one is NOT injective
CONTROL_PRIME = 5
CONTROL_ALGEBRA = LambdaAlgebra(CONTROL_PRIME)
CONTROL_SOURCE = stratum(CONTROL_PRIME, 2, 4, 5, 1)
CONTROL_TARGET = stratum(CONTROL_PRIME, 2, 5, 5, 2)
CONTROL_RANK = rank_mod_p(
    differential_matrix(CONTROL_ALGEBRA, CONTROL_SOURCE, CONTROL_TARGET),
    len(CONTROL_TARGET), CONTROL_PRIME)
require(CONTROL_RANK == len(CONTROL_SOURCE) - 1,
        "positive control: excess one at p=5, k=3 has a kernel")


print("THM-3205 ODD-PRIMARY LAMBDA ALGEBRA TODA GATE EXACT CONTROL")
print("primes=" + repr(PRIMES))
print("convention_audit_primes=" + repr(tuple(CONVENTION_ROWS)))
print("source_identities_replayed=mu0lam1,mu0mu1,mu0lam2,mu0mu2,"
      "dlam1,dmu1,dlam2,dmu2,d(mu2lam1),d(mu1lam2),d(mu1mu2lam1)")
print("d_squared_zero_words=" + str(SQUARE_ZERO_TOTAL))
print("unsigned_leibniz_failures=" + str(UNSIGNED_FAILURES))
print("confluence_words=" + str(CONFLUENCE_TOTAL))
print("excess_one_matrix=tridiag(1,-2,1)_(k+1)")
print("excess_one_gate=(kernel dim 1) iff p divides k+2 = index sum S")
print("gate_rows=" + str(len(GATE_ROWS)))
print("integral_cokernel=Z/(k+2) cyclic, invariant factors (1,...,1,k+2)")
print("smith_rows=" + repr(SMITH_ROWS))
print("kernel_generator=z_k=sum_(i=0)^k (i+1) u_i, leading mu_2 coefficient -1")
print("kernel_generator_confirmations=" + repr(GENERATOR_ROWS))
print("hostile_spurious_cycles=" + repr(SPURIOUS))
print("higher_excess_injective_rows=" + repr(EXCESS_ROWS))
print("positive_control_excess_one_p5_k3_kernel=1")
print("scope=Lambda*lambda(2) strata only; no new homotopy group is claimed")
print("ALL EXACT CHECKS PASSED")
