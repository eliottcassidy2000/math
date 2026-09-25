"""Fano/Hamming lattice, Clifford operators, and guarded affine arithmetic.

Source note: 05-knowledge/results/creation_fano_20260925.md.
Pure integer/Fraction verification; all checks survive Python optimization.
The Cayley-Dickson basis convention is inherited from the fano_code note.
"""

from fractions import Fraction
from itertools import combinations, product
import json


def check(value, message):
    if not value:
        raise RuntimeError(message)


def dotbit(a, b):
    return (a & b).bit_count() % 2


def code():
    return {tuple(s ^ dotbit(v, x) for x in range(8))
            for s in (0, 1) for v in range(8)}


def root_vectors(words):
    roots = set()
    for i in range(8):
        for sign in (-1, 1):
            v = [0] * 8
            v[i] = 2 * sign
            roots.add(tuple(v))
    for w in words:
        support = [i for i, b in enumerate(w) if b]
        if len(support) == 4:
            for signs in product((-1, 1), repeat=4):
                v = [0] * 8
                for i, sign in zip(support, signs):
                    v[i] = sign
                roots.add(tuple(v))
    return roots


def cd_basis(i, j, dimension=8):
    # (a,b)(c,d)=(ac-db*, a*d+cb), as in the inherited code note.
    if dimension == 1:
        return 1, 0
    h = dimension // 2
    if i < h and j < h:
        return cd_basis(i, j, h)
    if i < h <= j:
        sign, k = cd_basis(i, j - h, h)
        return sign * (1 if i == 0 else -1), k + h
    if j < h <= i:
        sign, k = cd_basis(j, i - h, h)
        return sign, k + h
    sign, k = cd_basis(j - h, i - h, h)
    return -sign * (1 if i == h else -1), k


def left_unit(i):
    return tuple(sign * (k + 1) for sign, k in
                 (cd_basis(i, j) for j in range(8)))


def compose(p, q):
    """Signed permutation p after q; entries are signed one-based images."""
    return tuple((1 if x > 0 else -1) * p[abs(x) - 1] for x in q)


def apply(p, vector):
    out = [0] * len(p)
    for value, image in zip(vector, p):
        out[abs(image) - 1] = value if image > 0 else -value
    return tuple(out)


def negate(p):
    return tuple(-x for x in p)


def matrix_inner(p, q):
    return sum((1 if x * y > 0 else -1) for x, y in zip(p, q)
               if abs(x) == abs(y))


def clifford_generators():
    gammas = []
    for i in range(1, 8):
        left = left_unit(i)
        bottom = tuple((1 if x > 0 else -1) * (abs(x) + 8) for x in left)
        gammas.append(bottom + left)
    gammas.append(tuple(-(i + 9) for i in range(8)) + tuple(range(1, 9)))
    return gammas


def step(n, sigma):
    return n // 2 if n % 2 == 0 else (3 * n + sigma) // 2


def word_coeff(word, sigma):
    a, b, d = 1, 0, 1
    for bit in word:
        if bit:
            a, b = 3 * a, 3 * b + sigma * d
        d *= 2
    return a, b, d


def actual_word(n, length, sigma):
    word = []
    for _ in range(length):
        word.append(n % 2)
        n = step(n, sigma)
    return tuple(word), n


def swap15(vector):
    out = list(vector)
    out[1], out[5] = out[5], out[1]
    return tuple(out)


def pair_hadamard_doubled(z):
    # Output 2Q(z/sqrt(2)), with pairs x,x XOR 4.
    return tuple(c for x in range(4) for c in (z[x] + z[x + 4], z[x] - z[x + 4]))


def main():
    words = code()
    roots = root_vectors(words)
    check(len(words) == 16 and len(roots) == 240, "inherited code/root sizes")
    units = [left_unit(i) for i in range(1, 8)]
    check(apply(compose(units[0], units[1]), tuple(int(j == 4) for j in range(8))) !=
          apply(negate(left_unit(3)), tuple(int(j == 4) for j in range(8))),
          "octonion multiplication must not be silently made associative")
    id8 = tuple(range(1, 9))
    for i, p in enumerate(units):
        check(compose(p, p) == negate(id8), "left-unit square")
        check({tuple(x % 2 for x in apply(p, w)) for w in words} == words,
              "left units must preserve affine parity code")
        check({apply(p, v) for v in roots} == roots, "left-unit root isometry")
        for q in units[i + 1:]:
            check(compose(p, q) == negate(compose(q, p)), "left-unit anticommutation")
    gammas = clifford_generators()
    id16 = tuple(range(1, 17))
    for i, p in enumerate(gammas):
        check(compose(p, p) == negate(id16), "gamma square")
        for q in gammas[i + 1:]:
            check(compose(p, q) == negate(compose(q, p)), "gamma anticommutation")
    monomials = []
    for mask in range(256):
        p = id16
        for i, gamma in enumerate(gammas):
            if mask >> i & 1:
                p = compose(p, gamma)
        monomials.append(p)
    check(all(matrix_inner(p, p) == 16 for p in monomials), "matrix norms")
    for p, q in combinations(monomials, 2):
        check(matrix_inner(p, q) == 0, "Clifford basis orthogonality")
    D = tuple(range(1, 9)) + tuple(-i for i in range(9, 17))
    J = gammas[-1]
    S = compose(D, J)
    check(monomials[-1] in (D, negate(D)), "volume chirality")
    volume_sign = 1 if monomials[-1] == D else -1
    root16 = {v + (0,) * 8 for v in roots} | {(0,) * 8 + v for v in roots}
    for gamma in gammas:
        check({apply(gamma, v) for v in root16} == root16, "doubled lattice roots")

    # Root's three-bit parity maps, independently computed from the dynamics.
    qplus = [sum(bit << j for j, bit in enumerate(actual_word(n, 3, 1)[0]))
             for n in range(8)]
    qminus = [sum(bit << j for j, bit in enumerate(actual_word(n, 3, -1)[0]))
              for n in range(8)]
    check(qplus == [0, 5, 2, 3, 4, 1, 6, 7], "positive parity map")
    check(qminus == [0, 7, 6, 1, 4, 3, 2, 5], "negative parity map")
    for x, y in product(range(8), repeat=2):
        defect = 4 * (((x & 1) * ((y >> 1) & 1)) ^
                      (((x >> 1) & 1) * (y & 1)))
        check(qplus[x] ^ qplus[y] ^ qplus[x ^ y] == defect, "polar carry")
        check(qminus[x] ^ qminus[y] == qminus[x ^ y], "minus linearity")
    changed_words = {swap15(w) for w in words}
    common_words = words & changed_words
    check(len(common_words) == 8, "code intersection size")
    check(common_words == {tuple(s ^ dotbit(v, x) for x in range(8))
                           for s in (0, 1) for v in range(4)}, "high-bit coefficient")
    changed_roots = {swap15(v) for v in roots}
    common_roots = roots & changed_roots
    check(len(common_roots) == 112 and len(roots - changed_roots) == 128,
          "neighbor-lattice root trade")
    standard_d8_doubled = set()
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            v = [0] * 8
            v[i], v[j] = a, b
            standard_d8_doubled.add(tuple(v))
    check({pair_hadamard_doubled(v) for v in common_roots} == standard_d8_doubled,
          "explicit D8 common-root isometry")
    for v in roots:
        w = pair_hadamard_doubled(v)
        w_changed = pair_hadamard_doubled(swap15(v))
        check(w_changed == tuple(-x if i == 3 else x for i, x in enumerate(w)),
              "carry is one-coordinate reflection in D8 chart")
        check(sum(w) % 4 == 0, "original spinor glue")
        check(sum(w_changed) % 4 == (2 if w[0] % 2 else 0), "changed spinor glue")
    check({tuple(w[qminus[x]] for x in range(8)) for w in words} == words,
          "minus linear coordinate transport")

    # A weight-four affine hyperplane supplies an E8 root rho=h/sqrt(2).
    h = tuple((x >> 2) & 1 for x in range(8))
    check(h in words and swap15(h) not in words, "chosen root detects glue change")
    affine_cases = 0
    for sigma in (-1, 1):
        for n in range(1, 2049):
            for length in range(1, 9):
                word, end = actual_word(n, length, sigma)
                a, b, d = word_coeff(word, sigma)
                check(a * n + b == d * end, "word affine identity")
                # Verify the exact four-term Clifford expression on this state.
                state = tuple(n * x for x in h) + h
                moved = tuple(Fraction((a + d) * x + (a - d) * y + b * (z + t), 2 * d)
                              for x, y, z, t in zip(state, apply(D, state),
                                                   apply(S, state), apply(J, state)))
                check(moved == tuple(end * x for x in h) + h, "Clifford carry operator")
                check(sum(x * x for x in state) // 2 == 2 * (n * n + 1),
                      "lattice state norm")
                qdiff = Fraction(((d - a) * n - b) * ((d + a) * n + b), d * d)
                check(qdiff == n * n - end * end, "norm descent factorization")
                check((n * n + 1) ** 2 - (end * end + 1) ** 2 ==
                      qdiff * (n * n + end * end + 2), "quartic factorization")
                affine_cases += 1
    guard_cases = 0
    for sigma in (-1, 1):
        for length in range(1, 7):
            for word in product((0, 1), repeat=length):
                a, b, d = word_coeff(word, sigma)
                for residue in range(d):
                    check(((a * residue + b) % d == 0) ==
                          (actual_word(residue, length, sigma)[0] == word), "word guard")
                    guard_cases += 1
    expansion = []
    for sigma in (-1, 1):
        for bits in range(1, 9):
            for length in range(1, 9):
                n = 2 ** (bits + length) - sigma
                word, end = actual_word(n, length, sigma)
                check(word == (1,) * length, "all-odd expansion word")
                check(end == 3 ** length * 2 ** bits - sigma, "expansion endpoint")
                check(n % (2 ** bits) == end % (2 ** bits) and end > n,
                      "same-residue expansion")
                check(end ** 2 + 1 > n ** 2 + 1, "quadratic hostile")
                if bits == 3 and length == 3:
                    expansion.append({"sigma": sigma, "n": n, "end": end})
    minus_cycle = [5, 7, 10]
    check([step(n, -1) for n in minus_cycle] == [7, 10, 5], "minus shortcut cycle")
    print(json.dumps({"status": "EXACT; no convergence claim", "code_words": len(words),
                      "E8_roots": len(roots), "left_unit_root_checks": 7 * 240,
                      "Clifford_generators": 8, "Clifford_monomials": len(monomials),
                      "orthogonal_monomial_pairs": 32640, "volume_equals_D_times": volume_sign,
                      "doubled_lattice_root_checks": 8 * len(root16),
                      "q_plus": qplus, "q_minus": qminus,
                      "intersection_code_dimension": 3, "common_D8_roots": len(common_roots),
                      "roots_in_each_glue_difference": 128,
                      "affine_state_cases": affine_cases, "guard_cases": guard_cases,
                      "same_residue_expansion_cases": 128, "expansion_examples": expansion,
                      "minus_shortcut_cycle": minus_cycle}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
