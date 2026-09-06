#!/usr/bin/env python3
"""Three named carry controls and an independent quadratic-family audit.

No census and no general sign theorem. Exact arithmetic, with SymPy for
polynomial operations; no repository producer is imported.
"""
from hashlib import sha256
from math import factorial, gcd
import sympy as S

T, G, Z = S.symbols('tau g s')
CHECKS = 0


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def full_row(a, b, c, mass):
    channels = []
    for z in range(mass+1):
        numerator = a*mass-(a+c)*z
        if numerator % (a+b):
            continue
        y = numerator//(a+b)
        x = mass-y-z
        if min(x, y) >= 0:
            channels.append(((x, y, z), factorial(mass)//
                             (factorial(x)*factorial(y)*factorial(z))))
    return channels


def literal_rows(a, b, c, last_mass):
    state = {(0, 0): 1}
    answer = []
    for _ in range(last_mass):
        following = {}
        for (charge, z), coefficient in state.items():
            for step, dz in ((-a, 0), (b, 0), (c, 1)):
                key = charge+step, z+dz
                following[key] = following.get(key, 0)+coefficient
        state = following
        answer.append({z: coefficient for (charge, z), coefficient in state.items()
                       if charge == 0})
    return answer


def trace_matrix(p, r):
    d = p.degree()
    # Multiplication by tau, represented columnwise in the power basis.
    multiplication = S.zeros(d)
    for j in range(d):
        column = S.Poly(S.rem(T**(j+1), p.as_expr(), T), T)
        for i in range(d):
            multiplication[i, j] = column.nth(i)
    r_matrix = S.zeros(d)
    for (k,), coefficient in r.terms():
        r_matrix += coefficient*multiplication**k
    matrix = S.Matrix(d, d, lambda i, j:
                      S.trace(r_matrix*multiplication**(i+j)))

    # Independent Newton recursion, without matrix multiplication.
    coefficients = [S.Rational(x)/p.LC() for x in p.all_coeffs()]
    traces = [S.Integer(d)]
    for n in range(1, 3*d):
        value = sum(coefficients[k]*traces[n-k]
                    for k in range(1, min(n, d)+1) if n-k > 0)
        if n <= d:
            value += n*coefficients[n]
        traces.append(-value)
    newton = S.Matrix(d, d, lambda i, j:
                      sum(coefficient*traces[i+j+k] for (k,), coefficient in r.terms()))
    need(matrix == newton and matrix == matrix.T, 'multiplication/Newton trace equality')
    return matrix


def quadratic_family_audit():
    # Reduce powers recursively, independently of a Euclidean Q/P call.
    powers = [(S.Integer(1), S.Integer(0)), (S.Integer(0), S.Integer(1))]
    for _ in range(2, 5):
        constant, linear = powers[-1]
        powers.append((-linear*(G-2)*(G-3)/12, constant-linear*(G-2)))
    constant = linear = 0
    for j in range(5):
        coefficient = S.prod(2*G-k for k in range(8-j))/S.Integer(
            factorial(8-2*j)*factorial(j))
        constant += coefficient*powers[j][0]
        linear += coefficient*powers[j][1]
    constant, linear = S.factor(constant), S.factor(linear)
    h = 5141*G**4-31762*G**3+69537*G**2-62820*G+18900
    jpoly = 2269*G**3-11468*G**2+19233*G-10710
    trace = S.factor(2*constant-(G-2)*linear)
    norm = S.factor(constant**2-(G-2)*constant*linear+
                    (G-2)*(G-3)*linear**2/12)
    expected_trace = -G*(G-2)*(G-1)*(2*G-3)*(2*G-1)*jpoly/15120
    expected_norm = G**2*(G-3)*(G-2)**2*(G-1)**2*(2*G-3)**2*\
        (2*G-1)**2*(5*G-7)*h/914457600
    need(S.cancel(trace-expected_trace) == 0, 'width-four exact trace')
    need(S.cancel(norm-expected_norm) == 0, 'width-four exact norm')
    for name, polynomial in (('H', h), ('J', jpoly)):
        shifted = S.Poly(polynomial.subs(G, Z+5).expand(), Z).all_coeffs()
        need(all(coefficient > 0 for coefficient in shifted), 'shift positivity '+name)
        print('width_four_shift', name, shifted)
    print('width_four_recursive_remainder', constant, linear)
    print('width_four_trace', trace)
    print('width_four_norm', norm)
    return str(trace), str(norm)


def main():
    manifest = []
    for a, b, c in ((15, 1, 9), (21, 1, 12), (25, 1, 14)):
        need(gcd(a, gcd(b, c)) == 1, 'primitive named support')
        g = gcd(a+b, a+c)
        A, B = (a+b)//g, (a+c)//g
        native = literal_rows(a, b, c, 2*g)
        need(all(not row for row in native[:g-1]), 'earlier masses empty')
        data = [full_row(a, b, c, m) for m in (g, 2*g)]
        for m, channels in zip((g, 2*g), data):
            need(native[m-1] == {v[2]: coefficient for v, coefficient in channels},
                 'complete literal versus arithmetic row')
        p, q = [S.Poly(sum(coefficient*T**j for j, (_, coefficient) in enumerate(row)), T)
                for row in data]
        x, y, z = data[0][0][0]
        h, r = divmod(y, B)
        ez, er = (2*z)//A, (2*r)//B
        need(q.degree() == 2*h+er+ez and p.degree() == h, 'complete degree carry')
        step = (B-A, -B, A)
        need(data[1][0][0] == tuple(2*v-ez*d for v, d in zip((x, y, z), step)),
             'canonical anchor carry')
        need(p.nth(0) > 0 and p.count_roots(-S.oo, 0) == h and
             S.gcd(p, p.diff()).degree() == 0, 'simple negative first roots')
        corrected = S.Poly(S.rem(q.as_expr()*S.invert(T**ez, p.as_expr(), T),
                                 p.as_expr(), T), T)
        matrix = trace_matrix(p, corrected)
        minors = [S.factor((-1)**k*matrix[:k, :k].det()) for k in range(1, h+1)]
        need(all(value > 0 for value in minors), 'named corrected negative definiteness')
        need(S.gcd(p, q).degree() == 0, 'named coprimality')
        if a == 15:
            canonical = q.rem(p)
            need(canonical.as_expr() == -47087024-466126752*T, 'canonical sign hostile')
            rho = -S.Rational(2942939, 29132922)
            need(rho > -5 and S.cancel(p.eval(rho)/56) ==
                 S.Rational(23910838225, 848727144258084), 'hostile root location')
            need(corrected.as_expr() == 47087024*T+4743488, 'corrected hostile residue')
        print('support', (-a, b, c), 'mass', g, 'factorial_tuple', (A, B, h, r, z, x),
              'carries', (er, ez))
        print('P_Q', list(reversed(p.all_coeffs())), list(reversed(q.all_coeffs())))
        print('corrected_remainder', corrected.as_expr())
        print('signed_leading_principal_minors', minors)
        manifest.append(((-a, b, c), g, str(p), str(q), ez, str(corrected), str(minors)))
    manifest.append(quadratic_family_audit())
    print('semantic_sha256', sha256(repr(manifest).encode()).hexdigest())
    print('explicit_checks', CHECKS)
    print('PASS: three named carry controls and independent width-four symbolic audit')
    print('SCOPE: Hermite equivalence is general; negative definiteness is not proved generally')


if __name__ == '__main__':
    main()
