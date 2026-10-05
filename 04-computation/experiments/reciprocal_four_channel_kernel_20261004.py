"""Exact four-channel trace and guarded affine-carrier controls.

Standard library only; the matching note contains the proofs and scope.
Explicit checks remain active under python -O.
"""

from fractions import Fraction as F
from itertools import product
from math import comb, isqrt


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


def exact_rational(x):
    if type(x) not in (int, F):
        raise ValueError("an exact integer or Fraction is required")
    return F(x)


def trace(matrix):
    return sum(matrix[i][i] for i in range(len(matrix)))


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def tensor(a, b):
    return [[a[i][j] * b[k][ell]
             for j in range(len(a[0])) for ell in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def recover_tensor(t):
    """Recover a marked positive-diagonal triangular carrier."""
    if (type(t) not in (list, tuple) or len(t) != 4
            or any(type(row) not in (list, tuple) or len(row) != 4 for row in t)):
        raise ValueError("a marked four-by-four tensor is required")
    if any(type(x) is not int for row in t for x in row):
        raise ValueError("exact integer tensor entries required")
    t = [list(row) for row in t]
    if t[0][0] <= 0 or t[3][3] <= 0:
        raise ValueError("positive diagonal carrier required")
    p, q = isqrt(t[0][0]), isqrt(t[3][3])
    if p*p != t[0][0] or q*q != t[3][3] or t[0][1] % p:
        raise ValueError("not a positive triangular tensor")
    b = t[0][1] // p
    m = [[p, b], [0, q]]
    if tensor(m, m) != t:
        raise ValueError("tensor entries do not reconstruct the carrier")
    return p, q, b


def j_value(lam):
    lam = exact_rational(lam)
    if lam <= 0:
        raise ValueError("positive multiplier required")
    return lam + 1 / lam


def stats(word):
    p, q, b = 1, 1, 0
    for a in word:
        if type(a) is not int or a < 1:
            raise ValueError("exact positive valuation letters required")
        p, q, b = 3 * p, q * 2**a, 3 * b + q
    return p, q, b


def matrix(word):
    p, q, b = stats(word)
    return [[p, b], [0, q]]


def encode(word):
    if not word:
        raise ValueError("the decoder uses nonempty positive words")
    p, q, b = stats(word)
    return j_value(F(p, q)), F(b, q - p)


def decode(j, anchor):
    """Recover a nonempty positive valuation word from typed J and its axis."""
    j, anchor = exact_rational(j), exact_rational(anchor)
    if j <= 2:
        raise ValueError("positive-word trace must exceed two")
    denominator, total, length = j.denominator, 0, 0
    while denominator % 2 == 0:
        denominator //= 2
        total += 1
    while denominator % 3 == 0:
        denominator //= 3
        length += 1
    p, q = 3**length, 2**total
    if denominator != 1 or length < 1 or total < length or j != j_value(F(p, q)):
        raise ValueError("not a positive-word trace")
    carry = (q - p) * anchor
    if carry.denominator != 1 or carry <= 0:
        raise ValueError("axis does not supply a positive integer carry")
    initial = int(carry)
    letters = []
    for remaining in range(length, 1, -1):
        difference = carry - 3**(remaining - 1)
        if difference.denominator != 1 or difference <= 0:
            raise ValueError("invalid carry prefix")
        difference = int(difference)
        a = (difference & -difference).bit_length() - 1
        if a < 1:
            raise ValueError("invalid valuation letter")
        letters.append(a)
        carry = F(difference, 2**a)
    last = total - sum(letters)
    if carry != 1 or last < 1:
        raise ValueError("invalid terminal carry or cost")
    word = tuple(letters + [last])
    if stats(word) != (p, q, initial):
        raise ValueError("decoded data do not reconstruct the input")
    return word


def source_residue(word):
    p, q, b = stats(word)
    return ((q - b) * pow(p, -1, 2 * q)) % (2 * q), 2 * q


def replay(n, word):
    if type(n) is not int or n <= 0 or n % 2 == 0:
        raise ValueError("positive odd integer source required")
    path = [n]
    for a in word:
        raw = 3 * n + 1
        actual = (raw & -raw).bit_length() - 1
        if actual != a:
            return None
        n = raw // 2**a
        path.append(n)
    return tuple(path)


def compositions(total):
    if total == 0:
        yield ()
    else:
        for first in range(1, total + 1):
            for tail in compositions(total - first):
                yield (first,) + tail


def symmetric_matrix(p, b, q, degree):
    """Polynomial action on e1^(k-i)e2^i for [[p,b],[0,q]]."""
    return [[comb(i, j) * p**(degree - i) * b**(i - j) * q**j
             if j <= i else 0 for i in range(degree + 1)]
            for j in range(degree + 1)]


def raw_representation_controls():
    # Columns: e00, e01+e10, e11, e01-e10. No diagonal eigenbasis assumed.
    basis = [[1, 0, 0, 0], [0, 1, 0, 1], [0, 1, 0, -1], [0, 0, 1, 0]]
    count = 0
    for a, b, c, d in product((-1, 0, 1), repeat=4):
        m = [[a, b], [c, d]]
        delta = a * d - b * c
        t = tensor(m, m)
        sym = [[a*a, 2*a*b, b*b], [a*c, a*d+b*c, b*d], [c*c, 2*c*d, d*d]]
        split = [row + [0] for row in sym] + [[0, 0, 0, delta]]
        need(multiply(t, basis) == multiply(basis, split), "3+1 intertwining")
        need(trace(t) == trace(multiply(m, m)) + 2 * delta, "determinant correction")
        need(trace(sym) == trace(multiply(m, m)) + delta, "symmetric character")
        count += 1
    return count


def determinant_iteration_controls():
    values = tuple(F(x, 2) for x in range(-2, 3))
    count = 0
    for a, b, c, d in product(values, repeat=4):
        m = [[a, b], [c, d]]
        t, delta = trace(m), a*d-b*c
        for _ in range(4):
            m = multiply(m, m)
            t, delta = t*t-2*delta, delta*delta
            need(trace(m) == t and m[0][0]*m[1][1]-m[0][1]*m[1][0] == delta,
                 "retained trace/determinant iteration")
            count += 1
    hostile = [[F(0), F(-1, 2)], [F(1), F(0)]]
    sequence = [trace(hostile)]
    for _ in range(2):
        hostile = multiply(hostile, hostile)
        sequence.append(trace(hostile))
    need(sequence == [F(0), F(-1), F(1, 2)], "middle quadratic is not invariant")
    need(sequence[2] != sequence[1]**2-1, "unpaid constant-determinant shortcut")
    return count


def reciprocal_controls():
    values = sorted({F(a, b) for a in range(1, 7) for b in range(1, 7)})
    count = 0
    for a, b in product(values, repeat=2):
        weights = (a*b, a/b, b/a, 1/(a*b))
        need(sum(weights) == j_value(a) * j_value(b), "ordered-pair kernel")
        need(sum(weights) == j_value(a*b) + j_value(a/b), "relative-mode retention")
        need((j_value(a/b) == 2) == (a == b), "cross-term equality boundary")
        count += 1
    for a in values:
        j = j_value(a)
        y = j + 2
        need(j_value(a*a) == j*j - 2, "Chebyshev square")
        need(j_value(a*a) + 2 == (y - 2)**2, "shifted coordinate")
        need(j == j_value(1/a), "ambient reciprocal loss")
    return len(values), count


def carrier_controls():
    words = [w for cost in range(1, 10) for w in compositions(cost)]
    characters, guards = 0, 0
    for w in words:
        p, q, b = stats(w)
        j, anchor = encode(w)
        need(decode(j, anchor) == w, "typed trace-axis decoder")
        need(recover_tensor(tensor(matrix(w), matrix(w))) == (p, q, b),
             "full marked tensor retains the carry")
        need(recover_tensor(tuple(tuple(row) for row in tensor(matrix(w), matrix(w))))
             == (p, q, b), "tuple tensor preserves the same exact carrier")
        need(F((p+q)**2, p*q) == j+2, "normalized tensor trace")
        residue, modulus = source_residue(w)
        for n in (residue, residue + modulus):
            path = replay(n, w)
            need(path is not None and path[-1] == F(p*n+b, q), "decoded source guard")
            guards += 1
        for degree in range(7):
            rep = symmetric_matrix(p, b, q, degree)
            expected = sum(p**(degree-i) * q**i for i in range(degree+1))
            need(trace(rep) == expected, "carry-blind symmetric character")
            characters += 1
    short = [w for cost in range(1, 6) for w in compositions(cost)]
    compositions_checked = 0
    for u, v in product(short, repeat=2):
        pu, qu, _ = stats(u)
        pv, qv, _ = stats(v)
        ju, au = encode(u)
        jv, av = encode(v)
        lu, lv = F(pu, qu), F(pv, qv)
        lam = lv * lu
        anchor = (lv*(1-lu)*au + (1-lv)*av)/(1-lam)
        need(decode(j_value(lam), anchor) == u+v, "lossless ordered composition")
        need(multiply(matrix(v), matrix(u)) == matrix(u+v), "matrix composition")
        need(multiply(tensor(matrix(v), matrix(v)), tensor(matrix(u), matrix(u)))
             == tensor(matrix(u+v), matrix(u+v)), "full tensor composition")
        need(ju*jv == j_value(lam)+j_value(lu/lv), "composed relative trace")
        compositions_checked += 1
    return len(words), characters, guards, compositions_checked


def hostile_controls():
    fibre = [w for w in compositions(5) if len(w) == 3]
    rows = []
    for w in fibre:
        p, q, b = stats(w)
        j, anchor = encode(w)
        residue, modulus = source_residue(w)
        need(j == F(1753, 864), "six-word common trace")
        rows.append((w, b, anchor, residue, modulus))
    need(len({row[1] for row in rows}) == 6 > 4, "four-state carry collision")
    need(len({row[3] for row in rows}) == 6, "six actual guards")
    need(stats((1,2,1,2)) == (81,64,85), "square order hostile")
    need(stats((1,1,2,2)) == (81,64,73), "nonhomomorphic square hostile")
    m = matrix((1,2))
    negated = [[-x for x in row] for row in m]
    need(tensor(m, m) == tensor(negated, negated), "untyped overall-sign loss")
    need(replay(27, (1,2)) is not None and replay(27, (1,2,1,2)) is None,
         "formal repeat does not pay its actual guard")
    need(replay(5, (4,2)) == (5,1,1), "actual word need not be first-hit word")
    for degree in range(7):
        rep = symmetric_matrix(F(1), F(-1,4), F(1), degree)
        need(trace(rep) == degree+1, "translation invisible to characters")
        if degree:
            need(rep != symmetric_matrix(F(1), F(0), F(1), degree),
                 "invisible does not mean identical")
    invalid = 0
    for args in ((F(2), F(0)), (F(5,2), F(1)), (F(13,6), F(0)),
                 (F(13,6), F(-2)), (True, F(-1)), (13/6, F(-1))):
        try:
            decode(*args)
        except ValueError:
            invalid += 1
        else:
            raise ArithmeticError("invalid typed trace-axis pair accepted")
    broken = tensor(m, m)
    broken[1][2] += 1
    tensor_invalid = 0
    for bad in (None, [[1]], [[True]*4 for _ in range(4)], broken):
        try:
            recover_tensor(bad)
        except ValueError:
            tensor_invalid += 1
        else:
            raise ArithmeticError("invalid tensor accepted")
    return rows, invalid, tensor_invalid


def main():
    raw = raw_representation_controls()
    iterations = determinant_iteration_controls()
    value_count, pair_count = reciprocal_controls()
    words, characters, guards, composition_count = carrier_controls()
    rows, invalid, tensor_invalid = hostile_controls()
    print("Four ordered channels are weighted pairs, not an oriented tournament.")
    print("Kernel: J(a)J(b)=J(ab)+J(a/b); cross term2 iff a=b>0.")
    print(f"Exact positive rational values {value_count}; ordered pair controls {pair_count}.")
    print(f"Raw matrix controls {raw}: full Sym^2 plus exterior-square intertwining and traces.")
    print("Matrix universe: all2x2 entries in{-1,0,1}, including singular boundary controls.")
    print(f"Trace/determinant iteration:625 half-integer matrices, four squarings each; {iterations} checks.")
    print("Determinant1/2 hostile: traces0,-1,1/2; the first x^2-1 step is not an invariant fiber.")
    print(f"All positive valuation words with total cost<=9: {words}.")
    print(f"Typed trace-axis round trips {words}; exact source replays {guards}.")
    print(f"Full marked tensor carrier controls {words} (list/tuple); ordered compositions {composition_count}.")
    print(f"Symmetric character controls, degrees0..6: {characters}.")
    print(f"Ordered compositions of all cost<=5 words: {composition_count}.")
    print("Six words with r3,A5,J1753/864:")
    for w, b, anchor, residue, modulus in rows:
        print(f"  {w}: carry {b}, anchor {anchor}, source {residue} mod{modulus}")
    print("Order hostile:1212 has carry85;1122 has carry73; both diagonal entries81,64.")
    print("Nontrivial translation-1/4 has all tested symmetric characters of the identity.")
    print("The full tensor retains the carry; positivity repairs its ambient overall-sign loss.")
    print("Source27 admits12 but not1212; source5 admits42 but hits1 before its last letter.")
    print(f"Invalid typed trace-axis inputs rejected {invalid}.")
    print(f"Invalid marked-tensor inputs rejected {tensor_invalid}.")
    print("No graph-size identity, intrinsic tournament orientation, or convergence claim.")


if __name__ == "__main__":
    main()
