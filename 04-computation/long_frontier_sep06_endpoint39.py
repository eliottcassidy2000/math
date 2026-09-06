#!/usr/bin/env python3
"""Exact endpoint-39 two-rung certificate on the complete carried fibre.

Symbolic polynomial matrices over QQ[x]; independent literal count fibres
and rational principal-minor characteristic controls. No repository producer
is imported. The analytic degree bound and real-root input are in the note.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb, factorial, gcd
import json
import sympy as S

x = S.Symbol('x')
Z = S.Poly(0, x, domain=S.QQ)
I = S.Poly(1, x, domain=S.QQ)
H = 6
GATES = 0
TRACE = sha256()


def gate(ok, label, data=None):
    global GATES
    if not ok:
        raise RuntimeError(f'{label}: {data}')
    GATES += 1
    TRACE.update((label + ':' + repr(data) + '\n').encode())


def falling(a, n):
    out = I
    for j in range(n):
        out *= S.Poly(a-j, x, domain=S.QQ)
    return out


def remainder(a, p):
    a = list(a)
    while len(a) > H:
        leading = a.pop()
        offset = len(a)-H
        for j in range(H):
            a[offset+j] -= leading*p[j]
    return a + [Z]*(H-len(a))


def symbolic_rows():
    p = [falling(x+6, 6-j).mul_ground(S.Rational(
        factorial(13), factorial(18-3*j)*factorial(1+2*j)))
        for j in range(7)]
    q = [falling(2*x+12, 12-e).mul_ground(S.Rational(
        1, factorial(36-3*e)*factorial(2+2*e)))
        for e in range(-1, 13)]
    carry = S.Poly(x, x, domain=S.QQ).mul_ground(S.Rational(
        128*factorial(18), factorial(39)*factorial(13)))
    for j in range(6):
        carry *= S.Poly(2*x+2*j+1, x, domain=S.QQ)
    gate(p[-1] == I, 'first row is monic')
    gate(q[0] == carry*p[0], 'inverse carry denominator cancels identically')
    rr = remainder(q[1:], p)
    rr = [rr[j]-carry*p[j+1] for j in range(H)]
    for j, rj in enumerate(rr):
        gate(rj.degree() <= 12-j, 'response weighted degree', j)
    return p, q, rr


def matrix_multiply(a, b):
    return [[sum((a[i][j]*b[j][k] for j in range(H)), Z)
             for k in range(H)] for i in range(H)]


def symbolic_characteristic(p, rr):
    columns = [remainder([Z]*j+rr, p) for j in range(H)]
    matrix = [[columns[j][i] for j in range(H)] for i in range(H)]
    for i in range(H):
        for j in range(H):
            gate(matrix[i][j].degree() <= 12+j-i, 'multiplication entry degree', (i,j))
    power = [[I if i == j else Z for j in range(H)] for i in range(H)]
    traces = []
    for _ in range(H):
        power = matrix_multiply(power, matrix)
        traces.append(sum((power[i][i] for i in range(H)), Z))
    coefficients = [I]
    for k in range(1, H+1):
        ck = sum((coefficients[k-j]*traces[j-1] for j in range(1,k+1)), Z)
        coefficients.append(ck.mul_ground(S.Rational(-1,k)))
        gate(coefficients[k].degree() == 12*k, 'characteristic degree', k)
    # Independently check the consequence operator, not only trace identities.
    total = [[I if i == j else Z for j in range(H)] for i in range(H)]
    for k in range(1, H+1):
        total = matrix_multiply(total, matrix)
        for i in range(H):
            total[i][i] += coefficients[k]
    gate(all(value.is_zero for row in total for value in row), 'symbolic Cayley-Hamilton')
    return coefficients


def literal_fibre(g, mass):
    # Enumerate every nonnegative count triple, without a channel template.
    out = {}
    for nc in range(mass+1):
        for nb in range(mass-nc+1):
            na = mass-nb-nc
            if -39*na+(2*g-39)*nb+(3*g-39)*nc:
                continue
            out[(na,nb,nc)] = F(factorial(mass),
                               factorial(na)*factorial(nb)*factorial(nc))
    return out


def rational_remainder(a, p):
    a = list(a)
    while len(a) > H:
        leading = a.pop()/p[-1]
        offset = len(a)-H
        for j in range(H):
            a[offset+j] -= leading*p[j]
    return a + [F(0)]*(H-len(a))


def determinant(a):
    a = [row[:] for row in a]
    value = F(1)
    for j in range(len(a)):
        i = next((i for i in range(j,len(a)) if a[i][j]), None)
        if i is None:
            return F(0)
        if i != j:
            a[j],a[i] = a[i],a[j]
            value = -value
        pivot = a[j][j]
        value *= pivot
        for i in range(j+1,len(a)):
            ratio = a[i][j]/pivot
            for k in range(j+1,len(a)):
                a[i][k] -= ratio*a[j][k]
    return value


def principal_characteristic(matrix):
    result = [F(1)]
    for k in range(1,H+1):
        value = sum((determinant([[matrix[i][j] for j in subset] for i in subset])
                     for subset in combinations(range(H),k)), F(0))
        result.append((-1)**k*value)
    return result


def exact_fibre_controls(p, q, symbolic):
    sample_x = list(range(1,13)) + [100]
    primitive = 0
    for xx in sample_x:
        g = xx+19
        first, second = literal_fibre(g,g), literal_fibre(g,2*g)
        expected_first = [(xx+j,18-3*j,1+2*j) for j in range(7)]
        expected_second = [(2*xx+e,36-3*e,2+2*e) for e in range(-1,13)]
        gate(set(first) == set(expected_first), 'complete seven-channel first fibre', xx)
        gate(set(second) == set(expected_second), 'complete fourteen-channel doubled fibre', xx)
        first_scale = F(comb(g,13))
        second_scale = F(factorial(2*g),factorial(2*g-26))
        pp = [first[count]/first_scale for count in expected_first]
        qq = [second[count]/second_scale for count in expected_second]
        gate(pp == [F(pj.eval(xx)) for pj in p], 'literal first coefficient law', xx)
        gate(qq == [F(qj.eval(xx)) for qj in q], 'literal doubled coefficient law', xx)
        inv = [-pp[j+1]/pp[0] for j in range(H)]
        gate(rational_remainder([F(0)]+inv,pp) == [F(1)]+[F(0)]*5,
             'inverse phase preserves original first equation', xx)
        rr = rational_remainder(qq[1:],pp)
        without_carry = rr[:]
        rr = [rr[j]+qq[0]*inv[j] for j in range(H)]
        gate(rr != without_carry and qq[0] > 0, 'carry changes the actual quotient response', xx)
        cols = [rational_remainder([F(0)]*j+rr,pp) for j in range(H)]
        mat = [[cols[j][i] for j in range(H)] for i in range(H)]
        cc = principal_characteristic(mat)
        gate(cc == [F(ck.eval(xx)) for ck in symbolic],
             'independent principal-minor characteristic agrees', xx)
        gate(all(value > 0 for value in cc), 'literal characteristic positive', xx)
        primitive += gcd(g,39) == 1
    # At g=21 the actual support is (-39,3,24), with an earlier mass seven.
    returns = [mass for mass in range(1,22) if literal_fibre(21,mass)]
    gate(returns[0] == 7, 'omitting first-return gcd hypothesis fails', returns)
    print('FIBRE_CONTROLS '+json.dumps({'x':sample_x, 'primitive':primitive,
        'first_channels':7, 'doubled_channels':14, 'gcd_hostile_g21_first_mass':7},
        separators=(',',':')))


def main():
    p,q,rr = symbolic_rows()
    cs = symbolic_characteristic(p,rr)
    count = 0
    for k in range(1,H+1):
        shifted = S.Poly(cs[k].as_expr().subs(x,x+1), x, domain=S.QQ)
        content, primitive = shifted.primitive()
        coeffs = list(reversed(primitive.all_coeffs()))
        gate(content > 0, 'positive characteristic content', k)
        for j, value in enumerate(coeffs):
            gate(value > 0, 'positive shifted characteristic coefficient', (k,j))
        count += len(coeffs)
        print('CHAR_CERTIFICATE '+json.dumps({'k':k,'degree':12*k,
              'content':str(content),'coefficients_ascending':[int(c) for c in coeffs]},
              separators=(',',':')), flush=True)
    gate(count == 258, 'complete shifted coefficient count', count)
    exact_fibre_controls(p,q,cs)
    print('CHAR_SCHEMA v=x-1; c_k(x)=content*sum(coefficients_ascending[j]*(x-1)^j)')
    print('CLAIM endpoint39: g=x+19>=20, gcd(g,39)=1, all nonzero complex coefficients; first_nonzero_mass=g_or_2g')
    print('CANCELING_PHASES six distinct negative tau roots; complete carried response negative at every root')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())


if __name__ == '__main__':
    main()
