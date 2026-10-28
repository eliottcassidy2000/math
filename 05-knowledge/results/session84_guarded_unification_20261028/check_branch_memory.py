#!/usr/bin/env python3
"""Exact, standard-library replay for the Session84 branch-memory note.

No assertions: checks remain live under python -O. Does not import prior code.
Words are written in chronological order, so M_(uv) = M_v M_u.
"""

from collections import deque
from math import gcd, isqrt
import json


I = (1, 0, 0, 1)
B = {"U": (2, -1, 1, 0), "A": (2, 1, 1, 0), "D": (1, 2, 0, 1)}
ROOT = (2, 1)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def mul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return (a * e + b * g, a * f + b * h,
            c * e + d * g, c * f + d * h)


def act(matrix, vector):
    a, b, c, d = matrix
    x, y = vector
    return (a * x + b * y, c * x + d * y)


def det(matrix):
    a, b, c, d = matrix
    return a * d - b * c


def word_matrix(word):
    matrix = I
    for letter in word:
        matrix = mul(B[letter], matrix)
    return matrix


def phi(vector):
    m, n = vector
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def residue(vector, modulus):
    return tuple(x % modulus for x in vector)


def projective_prime(vector, prime):
    m, n = residue(vector, prime)
    require(m != 0 or n != 0, "projective zero vector")
    return "infinity" if n == 0 else m * pow(n, -1, prime) % prime


def cf_sum(m, n):
    result = 0
    while n:
        head, remainder = divmod(m, n)
        result += head
        m, n = n, remainder
    return result


def decode(m, n):
    reverse = []
    while (m, n) != ROOT:
        if m > 3 * n:
            reverse.append("D")
            m, n = m - 2 * n, n
        elif m > 2 * n:
            reverse.append("A")
            m, n = n, m - 2 * n
        else:
            reverse.append("U")
            m, n = n, 2 * n - m
        require(m > n >= 1, "inverse branch left Euclid domain")
    return "".join(reversed(reverse))


def reachable(modulus, include_bit):
    initial = residue(ROOT, modulus) + ((1,) if include_bit else ())
    seen, queue = {initial}, deque([initial])
    while queue:
        state = queue.popleft()
        for letter, matrix in B.items():
            nxt = residue(act(matrix, state[:2]), modulus)
            if include_bit:
                nxt += (state[2] * det(matrix),)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return seen


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def squarefree(n):
    return all(n % (d * d) for d in range(2, isqrt(n) + 1))


def fus_delta(n):
    proper = [d for d in divisors(n) if 1 < d < n]
    return len(proper) - sum(prime(d) for d in proper) - sum(squarefree(d) for d in proper)


def aliquot(n):
    return sum(divisors(n)) - n


def main():
    # Independent literal branch products, then a closed-form range check.
    for modulus in range(1, 101):
        word = "A" + "D" * (modulus - 1) + "UD"
        expected = (8 * modulus - 1, 4, 2 * modulus, 1)
        require(word_matrix(word) == expected, "literal universal hostile matrix")
    for modulus in range(1, 5001):
        matrix = mul(B["D"], mul(B["U"], mul((1, 2 * (modulus - 1), 0, 1), B["A"])))
        point = act(matrix, ROOT)
        require(matrix == (8 * modulus - 1, 4, 2 * modulus, 1), "hostile formula")
        require(point == (16 * modulus + 2, 4 * modulus + 1), "hostile endpoint")
        require(gcd(*point) == 1 and (sum(point) % 2) == 1, "Euclid validity")
        require(det(matrix) == -1, "hostile character")
        require(residue(point, modulus) == residue(ROOT, modulus), "endpoint collision")
        require(residue(phi(point), modulus) == residue(phi(ROOT), modulus), "triple collision")
        require(residue(act(matrix, (1, 0)), modulus) == residue((-1, 0), modulus), "pole sign loss")

    checked_words = 0
    stack = [("", I)]
    while stack:
        word, matrix = stack.pop()
        point, pole = act(matrix, ROOT), act(matrix, (1, 0))
        sign = (-1) ** word.count("A")
        require(det(matrix) == sign == pole[0] * point[1] - pole[1] * point[0], "frame identity")
        require(decode(*point) == word, "independent inverse decoder")
        require((-1) ** (point[0] + cf_sum(*point)) == sign, "continued-fraction character")
        checked_words += 1
        if len(word) < 8:
            for letter in B:
                stack.append((word + letter, mul(B[letter], matrix)))

    census = []
    for modulus in range(1, 33):
        endpoints = reachable(modulus, False)
        lifted = reachable(modulus, True)
        require(lifted == {x + (sign,) for x in endpoints for sign in (-1, 1)}, "exact two-sheet reachable product")
        census.append({"M": modulus, "endpoints": len(endpoints), "endpoint_bit_states": len(lifted)})

    small = []
    for word in ("UU", "A", "", "A" + "D" * 6 + "UD", "D" * 7):
        matrix = word_matrix(word)
        point, pole = act(matrix, ROOT), act(matrix, (1, 0))
        small.append({"word": word or "empty", "matrix": matrix,
                      "endpoint": point, "triple": phi(point), "chi": det(matrix),
                      "pole_vector": pole, "endpoint_mod7": residue(point, 7),
                      "endpoint_projective_mod7": projective_prime(point, 7),
                      "pole_projective_mod7": projective_prime(pole, 7)})
    require(small[0]["endpoint_projective_mod7"] == small[1]["endpoint_projective_mod7"] == 6, "UU/A hostile")
    require(small[2]["endpoint_mod7"] == small[3]["endpoint_mod7"], "root/frame hostile endpoint")
    require(small[2]["pole_projective_mod7"] == small[3]["pole_projective_mod7"], "root/frame hostile pole")
    require(small[2]["chi"] == -small[3]["chi"], "root/frame hostile sign")
    require(residue(small[4]["matrix"], 7) == residue(I, 7), "full-matrix finite loss")
    require(residue(word_matrix("ADUD"), 2) == residue(I, 2), "mod2 full-frame exception")

    # Show source loss under the SAME lawful fresh-prime continuation.
    require(gcd(4, 5) == gcd(12, 5) == 1, "FUS common guard")
    require((fus_delta(4), fus_delta(12), fus_delta(20), fus_delta(60)) == (-1, -1, -1, 0), "FUS source-loss hostile")
    require(gcd(2, 5) == gcd(3, 5) == 1, "aliquot common guard")
    require((aliquot(2), aliquot(3), aliquot(10), aliquot(15)) == (1, 1, 8, 9), "aliquot source-loss hostile")
    print(json.dumps({
        "status": "PASS",
        "source_pin": "ae7c073e864d05be74fed85ba0a2bc164e026a6d",
        "universal_hostile_formula_checked_M": [1, 5000],
        "literal_word_product_checked_M": [1, 100],
        "all_words_max_depth": 8,
        "all_words_checked": checked_words,
        "reachable_product_census": census,
        "small_hostiles": small,
        "FUS_hostile": {"sources": [4, 12], "delta": [-1, -1], "common_multiplier": 5, "product_delta": [-1, 0]},
        "aliquot_hostile": {"sources": [2, 3], "s": [1, 1], "common_multiplier": 5, "product_s": [8, 9]},
        "scope": "Finite replay supplements the all-modulus and all-context algebraic proofs; no spectral census rerun."
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
