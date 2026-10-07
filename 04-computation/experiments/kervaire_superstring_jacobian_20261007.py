"""Exact small controls for typed word storage and polynomial inverse receipts.

The superstring solver here is bounded exact subset DP, not the imported
polynomial-time approximation algorithm. No external theorem is re-audited.
Run with Python or Python -O; no assertions and no mutation on import.
"""
from itertools import permutations, product


CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def word_guard(word):
    if type(word) is not tuple or any(type(t) is not str or len(t) != 1 for t in word):
        raise ValueError("explicit tuple of one-character token labels required")


def occurs(word, tape):
    return next((i for i in range(len(tape) - len(word) + 1)
                 if tape[i:i + len(word)] == word), None)


def overlap(left, right):
    return next((i for i in range(min(len(left), len(right)), -1, -1)
                 if left[len(left) - i:] == right[:i]), 0)


def shortest_tape(words):
    if type(words) is not tuple:
        raise ValueError("tuple of explicit words required")
    for word in words:
        word_guard(word)
    distinct = sorted(set(words))
    words = tuple(w for w in distinct if w and not any(w != v and occurs(w, v) is not None for v in distinct))
    if len(words) > 10:
        raise ValueError("bounded exact solver: at most ten maximal distinct words")
    if not words:
        return ()
    states = {(1 << i, i): w for i, w in enumerate(words)}
    for mask in range(1, 1 << len(words)):
        for last in range(len(words)):
            if (mask, last) not in states:
                continue
            tape = states[mask, last]
            for nxt, word in enumerate(words):
                if mask >> nxt & 1:
                    continue
                candidate = tape + word[overlap(words[last], word):]
                key = (mask | 1 << nxt, nxt)
                if key not in states or (len(candidate), candidate) < (len(states[key]), states[key]):
                    states[key] = candidate
    full = (1 << len(words)) - 1
    return min((states[full, i] for i in range(len(words))), key=lambda w: (len(w), w))


INVERSE = {"A": "a", "a": "A", "B": "b", "b": "B"}


def inverse_word(word):
    word_guard(word)
    if any(t not in INVERSE for t in word):
        raise ValueError("unregistered polynomial-automorphism token")
    return tuple(INVERSE[t] for t in reversed(word))


def reduce_word(word):
    inverse_word(word)
    stack = []
    for t in word:
        if stack and INVERSE[t] == stack[-1]:
            stack.pop()
        else:
            stack.append(t)
    return tuple(stack)


def apply_word(word, point, degree=2):
    inverse_word(word)
    if (type(degree) is not int or degree < 2 or type(point) is not tuple
            or len(point) != 2 or any(type(c) is not int for c in point)):
        raise ValueError("exact integer point and degree at least two required")
    x, y = point
    for t in word:
        if t in ("A", "a"):
            x += (1 if t == "A" else -1) * y ** degree
        else:
            y += (1 if t == "B" else -1) * x ** degree
    return x, y


def clean(poly):
    return {m: c for m, c in poly.items() if c}


def add(p, q):
    out = dict(p)
    for monomial, coefficient in q.items():
        out[monomial] = out.get(monomial, 0) + coefficient
    return clean(out)


def scale(p, c):
    return clean({m: a * c for m, a in p.items()})


def truncate(p, cap):
    return {m: a for m, a in p.items() if sum(m) <= cap}


def multiply(p, q, cap):
    out = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            monomial = i + k, j + l
            if sum(monomial) <= cap:
                out[monomial] = out.get(monomial, 0) + a * b
    return clean(out)


def power(p, n, cap):
    out = {(0, 0): 1}
    for _ in range(n):
        out = multiply(out, p, cap)
    return out


X, Y = {(1, 0): 1}, {(0, 1): 1}
IDENTITY = (X, Y)


def compose_poly(p, pair, cap):
    powers_x = {i: power(pair[0], i, cap) for i, _ in p}
    powers_y = {j: power(pair[1], j, cap) for _, j in p}
    out = {}
    for (i, j), a in p.items():
        out = add(out, scale(multiply(powers_x[i], powers_y[j], cap), a))
    return out


def compose_map(left, right, cap):
    return tuple(compose_poly(p, right, cap) for p in left)


def shear_map(token, degree=2):
    sign = 1 if token.isupper() else -1
    if token in ("A", "a"):
        return add(X, {(0, degree): sign}), Y
    if token in ("B", "b"):
        return X, add(Y, {(degree, 0): sign})
    raise ValueError("unknown token")


def word_map(word, cap, degree=2):
    inverse_word(word)
    out = IDENTITY
    for token in word:
        out = compose_map(shear_map(token, degree), out, cap)
    return out


def polynomial_guard(pair):
    if type(pair) is not tuple or len(pair) != 2:
        raise ValueError("exact pair of polynomial dictionaries required")
    for p in pair:
        if type(p) is not dict:
            raise ValueError("polynomial dictionary required")
        for m, a in p.items():
            if (type(m) is not tuple or len(m) != 2
                    or any(type(i) is not int or i < 0 for i in m)
                    or type(a) is not int or a == 0):
                raise ValueError("canonical integer polynomial coefficients required")
    if tuple(truncate(p, 1) for p in pair) != IDENTITY:
        raise ValueError("this executable normal form requires zero constant and identity linear part")


def bounded_inverse_receipt(pair, bound):
    """Complete at a declared inverse degree bound, in the normalized form.

    Failure means only that the inverse degree is not at most `bound`.
    """
    polynomial_guard(pair)
    if type(bound) is not int or bound < 1:
        raise ValueError("positive exact integer inverse-degree bound required")
    degree = max(sum(m) for p in pair for m in p)
    cap = degree * bound
    inverse = IDENTITY
    for _ in range(cap):
        residual = compose_map(pair, inverse, cap)
        inverse = tuple(add(g, add(coord, scale(err, -1)))
                        for g, coord, err in zip(inverse, IDENTITY, residual))
    if compose_map(pair, inverse, cap) != IDENTITY:
        raise RuntimeError("formal inverse iteration did not reach declared precision")
    candidate = tuple(truncate(g, bound) for g in inverse)
    if candidate != inverse:
        return None
    if compose_map(pair, candidate, cap) != IDENTITY or compose_map(candidate, pair, cap) != IDENTITY:
        raise RuntimeError("polynomial inverse receipt failed exact composition")
    return candidate


def perm_mul(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def perm_inverse(p):
    return tuple(p.index(i) for i in range(len(p)))


def quadratic_symbolic_receipt():
    """Universal six-parameter quadratic control, over QQ; no sampling."""
    import sympy as sp

    x, y, a, b, c, d, e, f = sp.symbols("x y a b c d e f")
    h = sp.Matrix([a*x*x + b*x*y + c*y*y,
                   d*x*x + e*x*y + f*y*y])
    jac = h.jacobian((x, y))
    original = [2*a+e, b+2*f, 2*a*e-2*b*d, 4*a*f-4*c*d, 2*b*f-2*c*e]
    generators = [2*a+e, b+2*f, 2*c*d+e*f, c*e+2*f*f, 4*d*f-e*e]
    # Explicit two-way ideal membership, independent of a Groebner assertion.
    expressed = [generators[0], generators[1],
                 e*generators[0]-2*d*generators[1]+generators[4],
                 2*f*generators[0]-2*generators[2],
                 2*f*generators[1]-2*generators[3]]
    recovered = [original[0], original[1],
                 f*original[0]-original[3]/2,
                 f*original[1]-original[4]/2,
                 original[2]-e*original[0]+2*d*original[1]]
    for lhs, rhs in zip(original, expressed):
        check(sp.expand(lhs-rhs) == 0, "original ideal generator identity")
    for lhs, rhs in zip(generators, recovered):
        check(sp.expand(lhs-rhs) == 0, "reverse ideal generator identity")
    determinant = sp.expand((sp.eye(2)+jac).det()-1)
    expected = original[0]*x + original[1]*y + original[2]*x*x + original[3]*x*y + original[4]*y*y
    check(sp.expand(determinant-expected) == 0, "complete quadratic Keller equations")
    basis = sp.groebner(original, a, b, c, d, e, f, order="lex", domain=sp.QQ)
    monic = [generators[0]/2, generators[1], generators[2]/2,
             generators[3], generators[4]/4]
    check(all(sp.expand(g.as_expr()-p) == 0 for g, p in zip(basis.polys, monic))
          and len(basis.polys) == 5, "declared Groebner basis")
    cubic = jac*h
    quartic = -jac*cubic-h.subs({x:h[0], y:h[1]}, simultaneous=True)
    coefficients = []
    for degree, part in ((3, cubic), (4, quartic)):
        for coordinate, expression in enumerate(part):
            for monomial, coefficient in sp.Poly(sp.expand(expression), x, y).terms():
                quotient, remainder = basis.reduce(coefficient)
                check(remainder == 0, "universal inverse obstruction vanishes in Keller ideal")
                reconstruction = sum(q*g.as_expr() for q, g in zip(quotient, basis.polys))
                check(sp.expand(coefficient-reconstruction) == 0, "explicit coefficient ideal receipt")
                coefficients.append((degree, coordinate, monomial, remainder))
    check(len(coefficients) == 18, "all cubic and quartic inverse coefficients")
    # A separate coefficient-ring path uses the polynomial jet implementation.
    symbolic_pair = (add(X, {(2,0):a, (1,1):b, (0,2):c}),
                     add(Y, {(2,0):d, (1,1):e, (0,2):f}))
    formal_pair = []
    for coordinate in range(2):
        expression = (x, y)[coordinate]-h[coordinate]+cubic[coordinate]+quartic[coordinate]
        formal_pair.append(dict(sp.Poly(sp.expand(expression), x, y).terms()))
    for composition in (compose_map(symbolic_pair, tuple(formal_pair), 4),
                        compose_map(tuple(formal_pair), symbolic_pair, 4)):
        for actual, expected_coordinate in zip(composition, IDENTITY):
            difference = add(actual, scale(expected_coordinate, -1))
            check(all(sp.expand(value) == 0 for value in difference.values()),
                  "unrestricted six-parameter truncated inverse identity")
    print("QUADRATIC UNIVERSAL: six parameters; Keller ideal basis=(2a+e,b+2f,2cd+ef,ce+2f^2,4df-e^2)")
    print("QUADRATIC INVERSE: all18 cubic/quartic coefficient remainders zero; two-way generator and coefficient identities verified")


def main():
    words = tuple(tuple(w) for w in ("ABa", "BaB", "aBb", "BbA"))
    tape = shortest_tape(words)
    # Independent finite permutation route for this four-word population.
    trials = []
    for ordering in permutations(words):
        trial = ordering[0]
        for word in ordering[1:]:
            if occurs(word, trial) is None:
                trial += word[overlap(trial, word):]
        trials.append(trial)
    check(len(tape) == min(map(len, trials)), "independent shortest-tape length")
    positions = []
    for word in words:
        start = occurs(word, tape)
        check(start is not None, "stored occurrence")
        stop = start + len(word)
        positions.append((start, stop))
        extraction = inverse_word(tape[:start]) + tape[:stop]
        check(reduce_word(extraction) == reduce_word(word), "prefix quotient extracts occurrence")
        for point in product(range(-1, 2), repeat=2):
            check(apply_word(extraction, point) == apply_word(word, point), "exact occurrence action")
            check(apply_word(word + inverse_word(word), point) == point, "inverse token receipt")
        cap = 2 ** len(word)
        check(word_map(word + inverse_word(word), cap) == IDENTITY, "symbolic inverse composition")
    for u, v, w in product(words, repeat=3):
        cycle = inverse_word(u) + v + inverse_word(v) + w + inverse_word(w) + u
        check(reduce_word(cycle) == (), "exact chart cocycle")
    check(shortest_tape(((),)) == (), "empty tape")
    check(shortest_tape((tuple("AB"), tuple("B"), tuple("AB"))) == tuple("AB"), "contained and duplicate words")
    print("STORAGE: four explicit words, 12 input tokens; exact shortest tape", "".join(tape), "length", len(tape))
    print("OCCURRENCES:", positions, "; polynomial inverses and all64 chart cocycles checked")

    for m in range(2, 10):
        commutator = tuple("ABab")
        check(apply_word(commutator, (1, 0), m) == (0, 1), "global commutator witness")
        first = word_map(commutator, 2 * m - 1, m)
        expected = (add(X, {(m, m - 1): -m}), add(Y, {(m - 1, m): m}))
        check(first == expected, "first invisible-jet defect")
        check(tuple(truncate(p, 2 * m - 2) for p in first) == IDENTITY, "all earlier jets miss holonomy")
    print("HOLONOMY: shear degrees2..9; first commutator defect degree2m-1, yet (1,0) maps to(0,1)")

    triangular = shear_map("A")
    check(bounded_inverse_receipt(triangular, 1) is None, "too-small inverse bound")
    check(bounded_inverse_receipt(triangular, 2) == shear_map("a"), "triangular inverse bound2")
    combined = word_map(tuple("AB"), 4)
    for bound in range(1, 5):
        result = bounded_inverse_receipt(combined, bound)
        check((result is not None) == (bound == 4), "noncommuting inverse degree4")
    check(result == word_map(tuple("ba"), 4), "exact noncommuting inverse")
    # Nonconstant Jacobian boundary: local invertibility is insufficient.
    nonkeller = (add(X, {(2, 0): -1}), Y)
    for bound in range(1, 9):
        check(bounded_inverse_receipt(nonkeller, bound) is None, "Catalan inverse never passed sampled bounds")
    # Formal determinant-one, nonpolynomial boundary. Both jets compose to Id.
    for cap in range(3, 13):
        formal = ({(i, 0): 1 for i in range(1, cap + 1)}, {(0, 1): 1, (1, 1): -2, (2, 1): 1})
        inverse = ({(i, 0): (-1) ** (i - 1) for i in range(1, cap + 1)}, {(0, 1): 1, (1, 1): 2, (2, 1): 1})
        check(compose_map(formal, inverse, cap) == IDENTITY, "formal inverse jets")
        check(compose_map(inverse, formal, cap) == IDENTITY, "formal inverse jets opposite direction")
    print("INVERSE RECEIPTS: declared bound works at2/4 for triangular/noncommuting examples; support bound is not inferred")

    group = tuple(permutations(range(3)))
    a, c, identity = (0, 2, 1), (1, 0, 2), (0, 1, 2)
    values = set()
    for t in group:
        value = perm_mul(perm_mul(perm_mul(perm_mul(t, a), perm_inverse(t)), a), t)
        values.add(value)
        check(perm_mul(value, c) != identity, "unimodular equation has no S3 solution")
    check(values == {identity, a}, "full S3 value image")
    # Same two associated-graded factors, different lift orders.
    cyclic_lifts = [a for a in range(9) if a % 3 == 1]
    check(all((3 * a) % 9 != 0 for a in cyclic_lifts), "Z9 quotient class lacks order3 lift")
    split_lifts = [(a, 1) for a in range(3)]
    check(all((3 * a % 3, 3 * b % 3) == (0, 0) for a, b in split_lifts), "split extension has order3 lifts")
    check([4 * 3 ** (j + 1) - 2 for j in (0, 2, 3)] == [10, 106, 322], "imported stem indexing")
    print("GROUP HOSTILE: all6 S3 substitutions fail exponent-sum-one equation; overgroup realization remains a different target")
    print("FILTERED HOSTILE: Z9 and Z3xZ3 share two Z3 graded factors, but differ on order3 lifts")

    quadratic_symbolic_receipt()

    malformed = [lambda: shortest_tape([tuple("AB")]), lambda: inverse_word(("unknown",)),
                 lambda: apply_word(tuple("A"), (True, 0)),
                 lambda: bounded_inverse_receipt(triangular, True),
                 lambda: bounded_inverse_receipt(({(1, 0): 1.0}, Y), 1),
                 lambda: bounded_inverse_receipt((add(X, {(0, 0): 1}), Y), 1)]
    for call in malformed:
        try:
            call()
        except ValueError:
            check(True, "malformed input rejected")
        else:
            check(False, "malformed input accepted")
    print("CHECKS", CHECKS, "; all explicit controls passed")


if __name__ == "__main__":
    main()
