#!/usr/bin/env python3
"""Carried Laurent trace formula, sharp boundary jets, and recurrence hostile.

Universe: exact rows h=1..7; determinant controls h=1..4. No producer imports.
The all-height assertions are proved in the accompanying note, not inferred
from the finite root counts printed here. Checks remain active under -O.
"""
from hashlib import sha256
from math import factorial
import json
import sympy as S

y = S.Symbol('y')
ZERO = S.Poly(0, y, domain=S.QQ)
ONE = S.Poly(1, y, domain=S.QQ)
Y = S.Poly(y, y, domain=S.QQ)
GATES = 0
DIGEST = sha256()


def gate(ok, label, data=None):
    global GATES
    if not ok:
        raise RuntimeError(f'{label}: {data}')
    GATES += 1
    DIGEST.update((label + ':' + repr(data) + '\n').encode())


def falling(a, n):
    result = ONE
    for j in range(n):
        result *= S.Poly(a-j, y, domain=S.QQ)
    return result


def remainder(a, p):
    a = list(a)
    while len(a) >= len(p):
        v = a.pop()
        for j, pj in enumerate(p[:-1]):
            a[len(a)-len(p)+1+j] -= v*pj
    return a + [ZERO]*(len(p)-1-len(a))


def rows(h):
    p = [falling(y, h-j).mul_ground(S.Rational(
        factorial(2*h+1), factorial(3*h-3*j)*factorial(1+2*j)))
        for j in range(h+1)]
    q = [falling(2*y, 2*h-e).mul_ground(S.Rational(
        1, factorial(6*h-3*e)*factorial(2+2*e)))
        for e in range(-1, 2*h+1)]
    carry = q[0].exquo(p[0])
    gate(carry*p[0] == q[0], 'polynomial-carry', h)
    regular = remainder(q[1:], p)
    full = [regular[j]-carry*p[j+1] for j in range(h)]
    return p, q, carry, full, regular


def multiplication(p, r):
    h = len(r)
    cols = [remainder([ZERO]*j+r, p) for j in range(h)]
    return [[cols[j][i] for j in range(h)] for i in range(h)]


def logarithmic_trace(p, q, carry):
    h = len(p)-1
    a = p[::-1] + [ZERO]*h
    logs = [ZERO]
    for n in range(1, 2*h+1):
        logs.append(a[n]-sum(
            (logs[j]*a[n-j]).mul_ground(S.Rational(j, n))
            for j in range(1, n)))
    return carry*p[1]-q[1].mul_ground(h)+sum(
        (q[e+1]*logs[e]).mul_ground(e) for e in range(1, 2*h+1))


def rational_string(v):
    return str(S.Rational(v))


def main():
    traces = {}
    for h in range(1, 8):
        p, q, carry, full, regular = rows(h)
        matrix = multiplication(p, full)
        trace = -sum((matrix[i][i] for i in range(h)), ZERO)
        gate(trace == logarithmic_trace(p, q, carry), 'logarithmic-formula', h)
        gate(trace.degree() == 2*h, 'observed-trace-degree', h)
        gate(all(v >= 0 for v in trace.all_coeffs()), 'observed-y-coefficients', h)
        # This is a finite root-isolation observation, not an all-h theorem.
        root_part = trace.exquo(Y) if h >= 2 else trace
        gate(S.gcd(root_part, root_part.diff()).degree() == 0,
             'observed-simple-roots', h)
        roots = root_part.count_roots(S.Rational(-1, 4), 0)
        gate(roots == root_part.degree(), 'observed-quarter-interval-roots', h)
        gate(root_part.eval(S.Rational(-1, 4)) != 0 and root_part.eval(0) != 0,
             'observed-strict-root-boundaries', h)
        gate(trace.eval(h+S.Rational(1, 2)) > 0,
             'noninteger-physical-x-half-trace-control', h)
        noninteger = h+S.Rational(1, 2)
        numeric_p = S.Poly(sum(p[j].eval(noninteger)*y**j
                              for j in range(h+1)), y, domain=S.QQ)
        numeric_matrix = S.Matrix([[v.eval(noninteger) for v in row]
                                    for row in matrix])
        gate(numeric_p.count_roots(-S.oo, 0) == h
             and S.gcd(numeric_p, numeric_p.diff()).degree() == 0
             and numeric_p.eval(0) != 0,
             'noninteger-x-half-simple-negative-first-roots', h)
        gate(all(v > 0 for v in numeric_matrix.charpoly().all_coeffs()),
             'noninteger-x-half-positive-characteristic-coefficients', h)
        d = S.Rational(factorial(2*h-1), factorial(6*h))
        b = S.Rational(2*factorial(2*h), factorial(6*h+3))
        a = S.Rational((-1)**(h-1)*factorial(h-1)*factorial(2*h+1),
                       factorial(3*h))
        if h >= 2:
            alpha = S.Rational(3*h*(3*h-1)*(3*h-2), 6)
            loss = 4*alpha/((h-1)*(6*h+1)*(6*h+2)*(6*h+3))
            slope = h*d*(1-loss)
            gate(trace.eval(0) == 0 and trace.nth(1) == slope,
                 'exact-trace-slope', h)
            gate(0 < loss < S.Rational(1, 12*(h-1)),
                 'strict-carry-slope-loss-bound', h)
            no_carry_trace = trace-carry*p[1]
            gate(no_carry_trace.nth(1) == h*d,
                 'exact-no-carry-trace-slope', h)
        else:
            gate(trace == S.Poly((884*y*y+123*y+1)/90720, y),
                 'height-one-exception')
        if h <= 4:
            for label, row, order, leading in [
                ('full', full, h-1, b**h/a),
                ('no-carry', regular, h, d**h),
            ]:
                mat = multiplication(p, row)
                norm = S.Poly((-1)**h*S.Matrix(
                    [[v.as_expr() for v in rr] for rr in mat]
                ).det(method='domain-ge'), y, domain=S.QQ)
                gate(all(norm.nth(j) == 0 for j in range(order))
                     and norm.nth(order) == leading,
                     'exact-determinant-jet', (h, label))
                if label == 'full':
                    divisor = ONE
                    endpoint_factor = S.Rational(1)
                    for r in range(2, h+1):
                        divisor *= S.Poly((y+r-h)**(r-1), y, domain=S.QQ)
                        if r < h:
                            endpoint_factor *= S.Integer(r-h)**(r-1)
                    deflated = norm.exquo(divisor)
                    gate(deflated.eval(0) == b**h/(a*endpoint_factor),
                         'fully-deflated-determinant-endpoint', h)
                    gate(S.sign(deflated.eval(0)) == (-1)**(h*(h-1)//2),
                         'fully-deflated-determinant-four-periodic-sign', h)
                if h == 2 and label == 'full':
                    gate(norm.eval(S.Rational(1, 100000)) < 0
                         and trace.eval(S.Rational(1, 100000)) > 0,
                         'positive-trace-negative-norm-hostile')
                    print('NORM_HOSTILE', json.dumps({
                        'h': h, 'y': '1/100000', 'x': '-199999/100000',
                        'trace': rational_string(trace.eval(S.Rational(1, 100000))),
                        'norm': rational_string(norm.eval(S.Rational(1, 100000))),
                    }, sort_keys=True))
        traces[h] = trace.monic()
        print('TRACE_ROW', json.dumps({
            'h': h, 'degree': trace.degree(),
            'negative_roots_in_open_quarter_interval': int(roots),
            'zero_root': h >= 2,
            'coefficients_ascending': [rational_string(trace.nth(j))
                                       for j in range(2*h+1)],
        }, sort_keys=True))

    # Disprove even an unsigned quadratic-coefficient three-term recurrence.
    target = traces[3]-Y*Y*traces[2]
    columns = [Y*traces[2], traces[2], Y*Y*traces[1],
               Y*traces[1], traces[1], target]
    prime = 101
    def mod(v):
        return int(v.p % prime)*pow(int(v.q % prime), -1, prime) % prime
    matrix = S.Matrix([[mod(col.nth(j)) for col in columns] for j in range(6)])
    residue = int(matrix.det()) % prime
    gate(residue == 29, 'quadratic-three-term-recurrence-obstruction')
    print('RECURRENCE_HOSTILE', json.dumps({
        'prime': prime, 'matrix': [[int(v) for v in row] for row in matrix.tolist()],
        'augmented_determinant_mod_prime': residue,
    }, sort_keys=True))
    print('PASS', GATES, 'exact gates')
    print('SEMANTIC_SHA256', DIGEST.hexdigest())


if __name__ == '__main__':
    main()
