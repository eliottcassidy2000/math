"""Exact finite dual-field closure receipts; bounded search is not a JC proof."""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product
from math import factorial
import sympy as s

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def rational(v):
    need(type(v) in (int, Fraction, s.Integer, s.Rational) or
         isinstance(v, s.Rational) and not isinstance(v, s.Float), "exact rational")
    need(type(v) is not bool, "no boolean rational")
    return s.Rational(v)


def variables(v):
    need(type(v) is tuple and len(v) >= 1 and
         all(type(x) is s.Symbol for x in v) and len(set(v)) == len(v),
         "distinct exact source symbols")
    return v


def polynomial(h, v):
    variables(v)
    need(isinstance(h, s.Expr) or type(h) in (int, Fraction), "literal exact expression")
    need(type(h) is not bool and not isinstance(h, (float, s.Float)), "exact polynomial")
    try:
        h = s.sympify(h)
        need(not h.atoms(s.Float) and h.free_symbols <= set(v), "QQ source polynomial")
        return s.Poly(h, *v, domain=s.QQ).as_expr()
    except (s.PolynomialError, s.CoercionFailed, TypeError) as exc:
        raise ValueError("QQ source polynomial") from exc


def carrier(v, f):
    variables(v)
    need(type(f) is tuple and len(f) == len(v), "complete source map")
    f = tuple(polynomial(h, v) for h in f)
    jac = s.Matrix([[s.diff(h, x) for x in v] for h in f])
    det = s.expand(jac.det())
    need(det != 0 and not det.free_symbols, "nonzero constant Jacobian")
    dual = (jac.adjugate().T / det).applyfunc(s.expand)
    return f, dual


def derive(v, dual, j, h):
    need(type(j) is int and 0 <= j < len(v), "exact dual-field index")
    h = polynomial(h, v)
    return s.expand(sum(dual[j, i] * s.diff(h, x) for i, x in enumerate(v)))


def coefficient_matrix(v, rows):
    polynomials = [s.Poly(h, *v, domain=s.QQ) for h in rows]
    monomials = sorted(set(m for p in polynomials for m in p.monoms()))
    return s.Matrix([[p.coeff_monomial(m) for p in polynomials] for m in monomials])


def coordinates(v, basis, h):
    matrix = coefficient_matrix(v, basis + (h,))
    try:
        answer, params = matrix[:, :-1].gauss_jordan_solve(matrix[:, -1])
    except ValueError:
        return None
    need(params.rows == 0, "independent marked basis")
    return answer


@dataclass(frozen=True)
class Closure:
    variables: tuple
    source_map: tuple
    basis: tuple
    actions: tuple
    anchor: tuple


def audit(receipt):
    need(type(receipt) is Closure, "literal closure receipt")
    v = variables(receipt.variables)
    f, dual = carrier(v, receipt.source_map)
    need(type(receipt.basis) is tuple and receipt.basis, "nonempty marked basis")
    basis = tuple(polynomial(h, v) for h in receipt.basis)
    m, n = len(basis), len(v)
    need(coefficient_matrix(v, basis).rank() == m, "independent polynomial basis")
    for h in (s.Integer(1),) + v:
        need(coordinates(v, basis, h) is not None, "retain unit and every source coordinate")
    need(type(receipt.actions) is tuple and len(receipt.actions) == n, "all dual actions")
    matrices = []
    for j, raw in enumerate(receipt.actions):
        need(type(raw) is tuple and len(raw) == m and
             all(type(row) is tuple and len(row) == m for row in raw), "exact square action")
        matrix = s.Matrix([[rational(z) for z in row] for row in raw])
        for col, h in enumerate(basis):
            expected = derive(v, dual, j, h)
            need(s.expand(sum(matrix[i, col] * basis[i] for i in range(m))) == expected,
                 "untruncated action on actual source polynomials")
        need(matrix**m == s.zeros(m), "characteristic-zero nilpotence check")
        matrices.append(matrix)
    need(all(a*b == b*a for a in matrices for b in matrices), "commuting actions")
    need(type(receipt.anchor) is tuple and len(receipt.anchor) == n, "retained source point")
    anchor = tuple(rational(z) for z in receipt.anchor)
    return v, f, basis, tuple(matrices), anchor


def compile_closure(v, f, anchor=None, max_dimension=40, max_rounds=12):
    """Search only this finite budget; None means no closure receipt was found."""
    v = variables(v)
    f, dual = carrier(v, f)
    need(type(max_dimension) is int and max_dimension >= len(v)+1 and
         type(max_rounds) is int and max_rounds >= 1, "exact positive search budgets")
    if anchor is None:
        anchor = (0,) * len(v)
    need(type(anchor) is tuple and len(anchor) == len(v), "exact source anchor")
    anchor = tuple(rational(z) for z in anchor)
    basis = (s.Integer(1),) + v
    for _ in range(max_rounds):
        old_size = len(basis)
        for h in basis:
            for j in range(len(v)):
                dh = derive(v, dual, j, h)
                if coordinates(v, basis, dh) is None:
                    if len(basis) >= max_dimension:
                        return None
                    basis += (dh,)
        if len(basis) == old_size:
            matrices = []
            for j in range(len(v)):
                cols = [coordinates(v, basis, derive(v, dual, j, h)) for h in basis]
                matrix = s.Matrix.hstack(*cols)
                matrices.append(tuple(tuple(matrix[i, k] for k in range(len(basis)))
                                      for i in range(len(basis))))
            result = Closure(v, f, basis, tuple(matrices), anchor)
            audit(result)
            return result
    return None


def multiindices(total, n):
    if n == 1:
        yield (total,)
    else:
        for a in range(total+1):
            for tail in multiindices(total-a, n-1):
                yield (a,) + tail


def compose(left, right, v):
    return tuple(s.expand(h.subs(dict(zip(v, right)), simultaneous=True)) for h in left)


def reconstruct(receipt):
    v, f, basis, matrices, anchor = audit(receipt)
    m, n = len(basis), len(v)
    at = dict(zip(v, anchor))
    values = s.Matrix([[h.subs(at) for h in basis]])
    target_origin = tuple(h.subs(at) for h in f)
    result = []
    for x in v:
        column = coordinates(v, basis, x)
        terms = []
        for total in range(m):
            for alpha in multiindices(total, n):
                current = column
                term = s.Integer(1)
                for j, a in enumerate(alpha):
                    current = matrices[j]**a * current
                    term *= (v[j]-target_origin[j])**a / factorial(a)
                terms.append((values * current)[0] * term)
        result.append(s.expand(sum(terms)))
    result = tuple(result)
    need(compose(result, f, v) == v and compose(f, result, v) == v,
         "both exact inverse identities on retained source")
    return result


def modular_iterate(v, dual, j, h, count, prime):
    need(type(prime) is int and s.isprime(prime), "prime characteristic")
    need(type(count) is int and count >= 0, "exact iteration count")
    def reduce(expr):
        p = s.Poly(expr, *v, domain=s.QQ)
        out = {}
        for monomial, c in p.terms():
            a, b = map(int, s.fraction(c))
            need(b % prime != 0, "good coefficient prime")
            out[monomial] = a * pow(b, -1, prime) % prime
        return s.Poly.from_dict(out, v, modulus=prime).as_expr()
    fields = dual.applyfunc(reduce)
    h = reduce(h)
    for _ in range(count):
        h = reduce(sum(fields[j, i]*s.diff(h, x) for i, x in enumerate(v)))
    return h


def thm1300(v):
    x, y, z = v
    u = 1+x*y
    return (s.expand(u**3*z+y**2*u*(4+3*x*y)),
            s.expand(y+3*x*u**2*z+3*x*y**2*(4+3*x*y)),
            2*x-3*x**2*y-x**3*z)


def height_bound(v, field, count):
    """L1 bound for count iterates of an INTEGRAL derivation on one coordinate."""
    variables(v)
    need(type(field) is tuple and len(field) == len(v), "full integral field")
    need(type(count) is int and count >= 0, "exact derivative depth")
    polys = [s.Poly(polynomial(h, v), *v, domain=s.QQ) for h in field]
    need(all(c.q == 1 for p in polys for c in p.coeffs()), "integer coefficients")
    degree = max([1]+[p.total_degree() for p in polys if not p.is_zero])
    size = max(sum(abs(int(c)) for c in p.coeffs()) for p in polys)
    answer = 1
    for j in range(count):
        answer *= size*(1+j*(degree-1))
    return answer


def main():
    x, y, z = s.symbols('x y z')
    v = (x, y)
    examples = []
    for a, b in product(range(1, 4), repeat=2):
        u = x+y**a
        f = (u, y+u**b)
        expected = (x-(y-x**b)**a, y-x**b)
        for anchor in ((0, 0), (1, -1)):
            packet = compile_closure(v, f, anchor, 30, 8)
            check(packet is not None, "bounded positive closure")
            inv = reconstruct(packet)
            check(tuple(map(s.expand, expected)) == inv, "independent shear inverse")
            check(max(s.Poly(h, *v).total_degree() for h in inv) <= len(packet.basis)-1,
                  "finite-space inverse degree")
            for point in product(range(-2, 3), repeat=2):
                images = tuple(h.subs(dict(zip(v, point))) for h in f)
                check(tuple(h.subs(dict(zip(v, images)), simultaneous=True) for h in inv) == point,
                      "source point round trip")
            examples.append((a, b, len(packet.basis)))
    identity = compile_closure(v, v)
    check(reconstruct(identity) == v, "degree-one boundary")
    simple = compile_closure(v, (x+y*y, y))
    check(len(simple.basis) == 3 and reconstruct(simple) == (x-y*y, y), "exact initial-span closure")
    check(compile_closure(v, (x+y**3, y), max_dimension=3) is None,
          "search cap is not noninvertibility")

    # Source multiplication is essential: a linear relabelling need not be an algebra map.
    check(s.diff(2*x*x, x) != 2*x, "unmatched product relation hostile")
    # Full trace pairing detects a cubic sheet hidden by its first scalar moments.
    t = s.symbols('t')
    mul = s.Matrix(((0, 0, t), (1, 0, 0), (0, 1, 0)))
    gram = s.Matrix([[s.trace(mul**(i+j)) for j in range(3)] for i in range(3)])
    check(s.trace(mul) == 0 and s.trace(mul**2) == 0, "two scalar moments vanish")
    check(gram.det() == -27*t*t and gram.rank() == 3, "full bilinear cubic rank")

    # A coefficient bill repairs a modular zero only at a fixed depth and large prime.
    h = x
    for depth in range(1, 9):
        h = s.expand(x*x*s.diff(h, x))
        bound = height_bound((x,), (x*x,), depth)
        check(h == factorial(depth)*x**(depth+1) and bound == factorial(depth),
              "exact factorial derivative/height bill")
        prime = int(s.nextprime(bound))
        check(s.Poly(h, x, modulus=prime).as_expr() != 0, "large-prime bill prevents false zero")
    check(s.Poly(6*x**4, x, modulus=3).is_zero and height_bound((x,), (x*x,), 3) == 6,
          "prime-depth zero lacks coefficient bill")

    w = (x, y, z)
    f3, fields = carrier(w, thm1300(w))
    check(s.Matrix([[s.diff(h, q) for q in w] for h in f3]).det().expand() == -2,
          "inherited exact Keller carrier")
    for prime in (3, 5):
        for j in range(3):
            for q in w:
                check(modular_iterate(w, fields, j, q, prime, prime) == 0,
                      "automatic prime-index nilpotence")
    # Rational tests at matching depths are nonzero. The inherited theorem proves all depths.
    nonzero = x
    degrees = []
    for k in range(1, 6):
        nonzero = derive(w, fields, 2, nonzero)
        check(nonzero != 0, "rational derivative survives modular nilpotence")
        degrees.append(s.Poly(nonzero, *w).total_degree())

    # Validate exact types before comparisons; reject projected/wrong-source action packets.
    bad = [lambda: compile_closure(v, (x*x, y)),
           lambda: compile_closure(v, (x+1.0, y)),
           lambda: compile_closure(v, ('x', y)),
           lambda: compile_closure(v, v, (True, 0)),
           lambda: compile_closure(v, v, max_dimension=True),
           lambda: compile_closure(v, v, max_rounds=1.0),
           lambda: audit(replace(identity, basis=(s.Integer(1), x))),
           lambda: audit(replace(identity, basis=(s.Integer(1), x, x))),
           lambda: audit(replace(identity, source_map=(x+y*y, y))),
           lambda: audit(replace(identity, actions=(((True,0,0),(0,0,0),(0,0,0)),)*2)),
           lambda: modular_iterate(w, fields, 0, x, 2, 2),
           lambda: audit(replace(identity, anchor=(s.Float(0),0))),
           lambda: height_bound(v, (x/2, y), 3),
           lambda: height_bound(v, (x, y), True)]
    for test in bad:
        try:
            test()
        except ValueError:
            check(True, "malformed or false receipt rejected")
        else:
            check(False, "hostile must reject")
    print('Exact finite dual closure -> source-bound polynomial inverse')
    print('Positive universe: nine two-shear maps, two retained anchors each')
    print('Closure dimensions:', examples[::2])
    print('Full cubic trace Gram determinant: -27*t^2; first two moments both zero')
    print('THM1300: all three duals kill every coordinate at characteristic3/5 prime depth')
    print('Characteristic-zero V3^k(x), k1..5: nonzero; degrees', degrees)
    print('Modular lifting control: fixed derivative depth plus explicit integral coefficient bound')
    print('Budget exhaustion is inconclusive; no universal finite closure or planar JC claim')
    print('checks:', CHECKS)


if __name__ == '__main__':
    main()
