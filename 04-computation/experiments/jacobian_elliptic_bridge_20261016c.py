"""Marked Lyness elliptic presentation and its differential boundary.

Exact rational APIs; stdout is deterministic JSON.  No Keller pair, curve
isomorphism, or inverse is inferred merely from an equality of j-invariants.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product
import json
import sympy as sp


def rational(value):
    if type(value) not in (int, Fraction):
        raise ValueError("exact rational required")
    return Fraction(value)


def point(value):
    if type(value) is not tuple or len(value) != 2:
        raise ValueError("ordered rational point required")
    return tuple(rational(v) for v in value)


def invariant(source):
    x, y = point(source)
    if not x or not y:
        raise ValueError("outside the torus")
    return (x + 1) * (y + 1) * (x + y + 1) / (x * y)


def step(source):
    x, y = point(source)
    if not x:
        raise ValueError("literal Lyness denominator vanishes")
    return y, (1 + y) / x


@dataclass(frozen=True)
class Cycle:
    source: tuple
    points: tuple
    h: Fraction
    determinants: tuple


def cycle(source):
    x, y = source = point(source)
    if not x * y * (x + 1) * (y + 1) * (x + y + 1):
        raise ValueError("full literal five-step guard fails")
    states, determinants = [source], []
    for _ in range(5):
        a, b = states[-1]
        determinants.append((1 + b) / (a * a))
        states.append(step(states[-1]))
    return Cycle(source, tuple(states), invariant(source), tuple(determinants))


def audit_cycle(packet):
    if type(packet) is not Cycle:
        raise ValueError("Cycle required")
    point(packet.source)
    if type(packet.points) is not tuple or len(packet.points) != 6:
        raise ValueError("six marked states required")
    for p in packet.points:
        point(p)
    rational(packet.h)
    if type(packet.determinants) is not tuple or len(packet.determinants) != 5:
        raise ValueError("five determinants required")
    for d in packet.determinants:
        rational(d)
    if packet != cycle(packet.source):
        raise ValueError("cycle data do not match the supplied source")
    return packet


def quartic(h):
    h = rational(h)
    # Descending coefficient order; the constant coefficient is zero.
    return (Fraction(1), 2 - 2*h, h*h - 6*h + 1, -4*h, Fraction(0))


def quartic_invariants(coefficients):
    if type(coefficients) is not tuple or len(coefficients) != 5:
        raise ValueError("five ordered quartic coefficients required")
    a, b, c, d, e = (rational(x) for x in coefficients)
    I = 12*a*e - 3*b*d + c*c
    J = 72*a*c*e + 9*b*c*d - 27*a*d*d - 27*b*b*e - 2*c*c*c
    delta = 4*I**3 - J**2
    return I, J, delta


def j_invariant(h):
    I, J, delta = quartic_invariants(quartic(h))
    if not delta:
        raise ValueError("singular Lyness fibre")
    return 6912 * I**3 / delta


def short_model(h):
    I, J, delta = quartic_invariants(quartic(h))
    if not delta:
        raise ValueError("singular Lyness fibre")
    return -I/3, -J/27


def to_short_model(source):
    x, y = point(source)
    h = invariant(source)
    A, B = short_model(h)
    c, d = h*h - 6*h + 1, -4*h
    v = 2*(x + 1)*y + x*x + (3-h)*x + 2
    W, Y = d/x + c/3, d*v/(x*x)
    if Y*Y != W**3 + A*W + B:
        raise ValueError("source does not lie on its marked cubic")
    return h, (W, Y)


def from_short_model(h, target):
    h, (W, Y) = rational(h), point(target)
    A, B = short_model(h)
    if Y*Y != W**3 + A*W + B:
        raise ValueError("target is not on this source model")
    c, d = h*h - 6*h + 1, -4*h
    if W == c/3:
        raise ValueError("inverse chart denominator vanishes")
    x = d/(W - c/3)
    if x == -1:
        raise ValueError("inverse quadratic chart denominator vanishes")
    v = Y*x*x/d
    y = (v - x*x - (3-h)*x - 2)/(2*(x+1))
    if not x*y or invariant((x, y)) != h:
        raise ValueError("target lies outside the retained torus chart")
    return x, y


def mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)


def cube(a):
    return mul(mul(a, a), a)


def evaluate(coefficients, value):
    value = rational(value)
    out = Fraction(0)
    for coefficient in reversed(coefficients):
        out = out*value + rational(coefficient)
    return out


@dataclass(frozen=True)
class ModuliCondition:
    source_j: Fraction
    b_polynomial: tuple
    h_polynomial: tuple


def moduli_condition(source_j):
    source_j = rational(source_j)
    # Ascending coefficients, including all zero entries.
    b = list(cube((Fraction(5), Fraction(10), Fraction(1))))
    b[1] -= source_j
    h = list(cube(tuple(map(Fraction, (1, 12, 14, -12, 1)))))
    for exponent, coefficient in ((5, -1), (6, -11), (7, 1)):
        h[exponent] -= source_j*coefficient
    return ModuliCondition(source_j, tuple(b), tuple(h))


def audit_moduli_condition(packet):
    if type(packet) is not ModuliCondition:
        raise ValueError("ModuliCondition required")
    rational(packet.source_j)
    for values, length in ((packet.b_polynomial, 7), (packet.h_polynomial, 13)):
        if type(values) is not tuple or len(values) != length:
            raise ValueError("wrong polynomial shape")
        for coefficient in values:
            rational(coefficient)
    if packet != moduli_condition(packet.source_j):
        raise ValueError("moduli polynomial lost its source")
    return packet


def matches(packet, h):
    packet = audit_moduli_condition(packet)
    h = rational(h)
    # Reject singular values even if used as malformed candidate coordinates.
    return j_invariant(h) == packet.source_j


def canonical_quartic(parameter):
    """THM-4389 actual ordered quartic, including the retained parameter."""
    parameter = rational(parameter)
    A = tuple(map(Fraction, (7231154026500, 50541940696500,
                           6793915500000, 353642000625)))
    B = tuple(map(Fraction, (631918028977864704, 91584545734393856)))
    coefficients = [4*v for v in mul(A, B)]
    coefficients[0] += 707284001250**2 * parameter**2
    return tuple(reversed(coefficients))


def canonical_moduli_condition(parameter):
    I, J, delta = quartic_invariants(canonical_quartic(parameter))
    if not delta:
        raise ValueError("singular source quartic")
    return moduli_condition(6912 * I**3 / delta)


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and not poly[-1]:
        poly.pop()
    return tuple(poly)


def fadd(a, b, p):
    return trim(tuple(((a[i] if i < len(a) else 0) +
                       (b[i] if i < len(b) else 0)) % p
                      for i in range(max(len(a), len(b)))))


def fmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = (out[i+j] + x*y) % p
    return trim(out)


def frem(a, b, p):
    a, b = list(trim(a)), trim(b)
    if b == (0,):
        raise ValueError("zero polynomial divisor")
    while len(a) >= len(b) and any(a):
        offset = len(a)-len(b)
        coefficient = a[-1]*pow(b[-1], -1, p) % p
        for i in range(len(b)):
            a[i+offset] = (a[i+offset]-coefficient*b[i]) % p
        a = list(trim(a))
    return tuple(a)


def fgcd(a, b, p):
    while b != (0,):
        a, b = b, frem(a, b, p)
    return tuple(x*pow(a[-1], -1, p) % p for x in a)


def fpow(a, n, modulus, p):
    out = (1,)
    while n:
        if n & 1:
            out = frem(fmul(out, a, p), modulus, p)
        a = frem(fmul(a, a, p), modulus, p)
        n //= 2
    return out


def irreducible_mod3_certificate(coefficients):
    """Finite Frobenius criterion for the two declared monic degrees 6,12."""
    if type(coefficients) is not tuple or len(coefficients) not in (7, 13):
        raise ValueError("declared degree 6 or 12 required")
    if any(type(c) is not int for c in coefficients) or coefficients[-1] != 1:
        raise ValueError("monic exact integer polynomial required")
    f = trim(tuple(c % 3 for c in coefficients))
    n = len(f)-1
    exponents = (n//2, n//3)
    residues = tuple(fpow((0, 1), 3**i, f, 3) for i in exponents)
    gcds = tuple(fgcd(f, fadd(r, (0, 2), 3), 3) for r in residues)
    final = fpow((0, 1), 3**n, f, 3)
    if gcds != ((1,), (1,)) or final != (0, 1):
        raise ValueError("irreducibility certificate fails")
    return {"modulus": 3, "polynomial": f, "test_exponents": exponents,
            "frobenius_residues": residues, "gcds": gcds,
            "final_frobenius_residue": final}


def verify():
    checks = 0

    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise ArithmeticError(label)

    def reject(call, label):
        try:
            call()
        except ValueError:
            check(True, label)
        else:
            check(False, label)

    x, y, h, b = sp.symbols("x y h b")
    F = (x+1)*(y+1)*(x+y+1)
    H = F/(x*y)
    lx, ly = y, (1+y)/x
    sub = {x: lx, y: ly}
    jac = sp.det(sp.Matrix([lx, ly]).jacobian([x, y]))
    check(sp.cancel(H.subs(sub, simultaneous=True)-H) == 0, "invariant")
    check(sp.cancel(jac/(lx*ly)-1/(x*y)) == 0, "log form")
    check(sp.cancel(jac-(1+y)/x**2) == 0, "actual Jacobian")
    seq = [x, y]
    for _ in range(5):
        seq.append(sp.cancel((1+seq[-1])/seq[-2]))
    check(seq[-2:] == [x, y], "rational five period")
    C = sp.expand(F-h*x*y)
    quadratic_B = x*x+(3-h)*x+2
    v = 2*(x+1)*y+quadratic_B
    q = x**4+(2-2*h)*x**3+(h*h-6*h+1)*x*x-4*h*x
    check(sp.expand(v*v-q-4*(x+1)*C) == 0, "quartic completion")
    I = h**4-12*h**3+14*h*h+12*h+1
    J = -2*(h*h+1)*(h**4-18*h**3+74*h*h+18*h+1)
    D = h**5*(h*h-11*h-1)
    qc = sp.Poly(q, x).all_coeffs()
    qa,qb,qc2,qd,qe = qc
    check(sp.expand(12*qa*qe-3*qb*qd+qc2*qc2-I) == 0, "quartic I")
    check(sp.expand(72*qa*qc2*qe+9*qb*qc2*qd-27*qa*qd*qd-
                    27*qb*qb*qe-2*qc2**3-J) == 0, "quartic J")
    check(sp.expand(4*I**3-J**2-6912*D) == 0, "discriminant")
    check(sp.resultant(I, D, h) == 25, "no spurious singular matching root")
    check(sp.cancel(I/h**2-((h-1/h-11)**2+10*(h-1/h-11)+5)) == 0,
          "sextic quotient")
    check(sp.cancel(D/h**6-(h-1/h-11)) == 0, "sextic denominator")
    check(sp.cancel(I.subs(h, -1/h)*h**4-I) == 0, "parameter involution I")
    check(sp.cancel(D.subs(h, -1/h)*h**12-D) == 0, "parameter involution D")
    W, Y = -4*h/x+(h*h-6*h+1)/3, -4*h*v/x**2
    check(sp.cancel(Y*Y-W**3+I*W/3+J/27-16*h*h*(v*v-q)/x**4) == 0,
          "literal short model")
    check(sp.cancel(sp.diff(W,x)/(2*Y)+1/(2*v)) == 0, "differential scale")
    # Projective points at infinity have a nonzero source gradient.
    X,Z,T = sp.symbols("X Z T")
    cubic = (X+T)*(Z+T)*(X+Z+T)-h*X*Z*T
    grad = [sp.diff(cubic,z) for z in (X,Z,T)]
    for P in ((1,0,0),(0,1,0),(1,-1,0)):
        check(any(sp.expand(g.subs(dict(zip((X,Z,T),P)))) != 0 for g in grad),
              "smooth infinity")
    golden = (1+sp.sqrt(5))/2
    hg = sp.simplify(H.subs({x: golden,y:golden}))
    check(sp.simplify(hg-golden**5) == 0, "golden fixed fibre")
    check(sp.simplify(hg**2-11*hg-1) == 0, "golden singular")
    check(all(sp.simplify(g.subs({x:golden,y:golden,h:hg})) == 0
              for g in (C,sp.diff(C,x),sp.diff(C,y))), "actual singular point")
    lam = sp.symbols("lam")
    check(sp.expand((1-lam+lam*lam)**3-8*lam*lam*(1-lam)**2-
                    (lam*lam-lam-1)*(lam*lam-3*lam+1)*(lam*lam+lam-1)) == 0,
          "golden Legendre contrast")
    # Finite residue controls illustrate, but do not replace, the proof.
    check(sp.residue(sp.residue(1/(x*y), x, 0), y, 0) == 1, "log residue")
    forms = (x/y, y/x, 1/(x*(x+y)), (x+y)/(x*x*y*y), 1/(x-y))
    for A, B in product(forms, repeat=2):
        exact = sp.diff(B,x)-sp.diff(A,y)
        check(sp.residue(sp.residue(exact,x,0),y,0) == 0, "exact residue control")

    rational_sources = set(Fraction(n,d) for n in range(-7,8) for d in range(1,4))
    receipts, failed_guards = 0, 0
    for x0,y0 in product(sorted(rational_sources), repeat=2):
        if not x0*y0*(x0+1)*(y0+1)*(x0+y0+1):
            reject(lambda: cycle((x0,y0)), "native guard")
            failed_guards += 1
            continue
        packet = audit_cycle(cycle((x0,y0)))
        check(packet.points[-1] == packet.source, "actual period")
        check(all(invariant(p) == packet.h for p in packet.points), "all source fibres")
        prod = Fraction(1)
        for determinant in packet.determinants:
            prod *= determinant
        check(prod == 1, "five-step determinant")
        hh, ww = to_short_model((x0,y0))
        check(from_short_model(hh,ww) == (x0,y0), "marked curve inverse")
        condition = audit_moduli_condition(moduli_condition(j_invariant(hh)))
        check(matches(condition,hh) and matches(condition,-1/hh), "marked moduli pair")
        check(evaluate(condition.h_polynomial,hh) == 0, "degree12 equation")
        bb = hh-1/hh-11
        check(evaluate(condition.b_polynomial,bb) == 0, "degree6 equation")
        receipts += 1

    # The rational five-period identity is not an everywhere-defined action.
    hostile = (Fraction(-1),Fraction(2))
    hostile_trace = [hostile]
    for _ in range(4):
        hostile_trace.append(step(hostile_trace[-1]))
    check(hostile_trace[-1] == (0,-1), "base-point arrival")
    reject(lambda: step(hostile_trace[-1]), "literal fifth step forbidden")
    check(cycle((1,2)).determinants[0] == 3, "not a Keller map")

    source_parameters = tuple(sorted(set(Fraction(i,d) for d in (1,2,3) for i in range(-4,5))))
    canonical_rows = []
    for parameter in source_parameters:
        coeffs = canonical_quartic(parameter)
        Ii,Ji,delta = quartic_invariants(coeffs)
        check(delta != 0, "declared actual canonical fibres smooth")
        condition = audit_moduli_condition(canonical_moduli_condition(parameter))
        check(condition.h_polynomial[-1] == condition.b_polynomial[-1] == 1,
              "source compiler stays monic")
        check(condition.source_j == 6912*Ii**3/delta, "canonical source j retained")
        if parameter in (0,1):
            canonical_rows.append({"parameter": str(parameter),
                                   "j": str(condition.source_j)})
    check(canonical_moduli_condition(0).source_j != canonical_moduli_condition(1).source_j,
          "actual pencil is not collapsed to one curve")

    condition = moduli_condition(2048)
    certificates = []
    for coefficients in (condition.b_polynomial, condition.h_polynomial):
        check(all(c.denominator == 1 for c in coefficients), "integral model")
        certificate = irreducible_mod3_certificate(tuple(c.numerator for c in coefficients))
        check(certificate["gcds"] == ((1,),(1,)), "independent irreducibility")
        certificates.append(certificate)
    # Separate symbolic identities do not manufacture an isomorphism receipt.
    check(not matches(condition,1), "same word golden is not same parameter")
    reject(lambda: matches(condition,0), "singular value rejected")
    # Good reduction of E0 at 3 and its induced five-torsion Frobenius.
    e0_f3 = [(a,b) for a,b in product(range(3),repeat=2)
             if (b*b-(a*a*a+a*a-3*a-2)) % 3 == 0]
    check(len(e0_f3)+1 == 6, "E0 over F3 includes its infinity point")
    M, identity = ((0,2),(1,3)), ((1,0),(0,1))
    def mm(A,B):
        return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))%5
                           for j in range(2)) for i in range(2))
    def mv(A,v):
        return tuple(sum(A[i][k]*v[k] for k in range(2))%5 for i in range(2))
    check((M[0][0]+M[1][1])%5 == (-2)%5 and
          (M[0][0]*M[1][1]-M[0][1]*M[1][0])%5 == 3,
          "retained Frobenius trace and determinant")
    check(all((a*a+2*a+3)%5 for a in range(5)), "irreducible Frobenius polynomial")
    M2 = mm(M,M)
    check(all((M2[i][j]+2*M[i][j]+3*identity[i][j])%5 == 0
              for i,j in product(range(2),repeat=2)), "Frobenius characteristic equation")
    powers, scalar_powers = {0:identity}, []
    for k in range(1,25):
        powers[k] = mm(powers[k-1],M)
        if powers[k][0][1] == powers[k][1][0] == 0 and powers[k][0][0] == powers[k][1][1]:
            scalar_powers.append(k)
    check(scalar_powers == [6,12,18,24], "first projective return")
    check(powers[6] == ((3,0),(0,3)), "sixth Frobenius power")
    check(powers[12] == ((4,0),(0,4)), "twelfth Frobenius power")
    check(powers[24] == identity, "twenty-fourth Frobenius power")
    for vector in product(range(5),repeat=2):
        if vector == (0,0):
            continue
        orbit = [mv(powers[k],vector) for k in range(24)]
        check(len(set(orbit)) == 24 and mv(powers[24],vector) == vector,
              "full nonzero torsion-vector orbit")
        signed = [min(v,tuple(-a%5 for a in v)) for v in orbit]
        check(len(set(signed)) == 12 and signed[:12] == signed[12:],
              "torsion point modulo sign has period twelve")
    warm = cycle((1,2))
    hostiles = [lambda: cycle((True,2)), lambda: cycle((1.0,2)),
                lambda: cycle([1,2]), lambda: invariant((0,1)),
                lambda: audit_cycle(replace(warm,source=(True,2))),
                lambda: audit_cycle(replace(warm,h=float(warm.h))),
                lambda: audit_cycle(replace(warm,points=warm.points[::-1])),
                lambda: audit_cycle(replace(warm,source=(2,1))),
                lambda: moduli_condition(True), lambda: canonical_moduli_condition(0.0),
                lambda: audit_moduli_condition(replace(condition,source_j=2048.0)),
                lambda: audit_moduli_condition(replace(condition,b_polynomial=(True,)+condition.b_polynomial[1:])),
                lambda: from_short_model(1,(1,1)),
                lambda: irreducible_mod3_certificate((True,0,0,0,0,0,1)),
                lambda: irreducible_mod3_certificate((0,0,0,0,0,0,1))]
    for bad in hostiles:
        reject(bad,"malformed or false receipt")
    return {"checks": checks, "status": "PROVED scoped algebra; JC(2) OPEN",
            "rational_source_universe": {"coordinate_count": len(rational_sources),
                "native_cycles": receipts, "failed_full_guards": failed_guards},
            "canonical_parameter_controls": len(source_parameters),
            "canonical_j_samples": canonical_rows,
            "golden_lyness_parameter": "(11+5*sqrt(5))/2; singular",
            "golden_legendre_j": 2048,
            "moduli_2048": {"b_polynomial": [str(c) for c in condition.b_polynomial],
                "h_polynomial": [str(c) for c in condition.h_polynomial],
                "irreducibility_certificates": certificates},
            "frobenius_2048": {"prime":3,"affine_points":e0_f3,"point_count":6,
                "trace":-2,"torsion_prime":5,"matrix":M,
                "scalar_powers":scalar_powers,"nonzero_vector_orbit":24,
                "signed_vector_orbit":12,"projective_orbit":6},
            "literal_basepoint_hostile": [[str(a),str(b)] for a,b in hostile_trace],
            "malformed_interface_controls": len(hostiles),
            "not_claimed": ["Keller pair", "chart entry or higher-row lift",
                "curve isomorphism from j alone", "rational flattening of the marked log form",
                "obstruction to every birational conjugacy of the map"]}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
