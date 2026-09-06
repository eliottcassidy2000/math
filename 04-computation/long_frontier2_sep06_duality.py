#!/usr/bin/env python3
"""Exact global regular-beta orbit and conserved-amplitude obstruction.

Universe: displayed formal parameter identities, three symbolic height
controls, one global H1 orbit, named primitive signs/double zero, complete
literal genuine rows, and two reflection/carry controls. No shape census.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd
import hashlib
import json

GATES = 0


def require(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


# Sparse polynomials in the formal parameter z and attenuation lambda.
def mon(i=0, j=0, c=1):
    return {(i, j): F(c)} if c else {}


def add(*ps):
    out = {}
    for p in ps:
        for k, c in p.items():
            out[k] = out.get(k, F(0))+c
    return {k: c for k, c in out.items() if c}


def scale(p, c):
    return {k: a*c for k, a in p.items() if a*c}


def mul(p, q):
    out = {}
    for (i, j), a in p.items():
        for (k, ell), b in q.items():
            key = (i+k, j+ell)
            out[key] = out.get(key, F(0))+a*b
    return {k: c for k, c in out.items() if c}


def power(p, n):
    out = mon()
    for _ in range(n):
        out = mul(out, p)
    return out


def shift(p):
    out = {}
    for (i, j), c in p.items():
        out = add(out, *(mon(k, j, c*comb(i, k)) for k in range(i+1)))
    return out


def linear(a, b):
    return add(mon(1, c=a), mon(c=b))


def rising(a, b, n):
    out = mon()
    for k in range(n):
        out = mul(out, linear(a, b+k))
    return out


def eval_poly(p, z, lam=1):
    return sum(c*z**i*lam**j for (i, j), c in p.items())


def specialize_lambda(p, lam):
    return add(*(mon(i, c=c*lam**j) for (i, j), c in p.items()))


def rows(H):
    phi = [scale(rising(1, 2*j+1, 2*H-2*j),
                 F(factorial(H), factorial(j)*factorial(3*H-3*j))) for j in range(H+1)]
    psi = [scale(rising(2, 2*j+1, 4*H-2*j),
                 F(factorial(2*H), factorial(j)*factorial(6*H-3*j))) for j in range(2*H+1)]
    return phi, psi


for H in [1, 2, 6]:
    phi, psi = rows(H)
    for j, p in enumerate(phi):
        require(mul(linear(1, 2*j+1), shift(p)) == mul(linear(1, 2*H+1), p),
                'formal regular first beta step')
    numerator = mul(linear(2, 4*H+1), linear(2, 4*H+2))
    for j, p in enumerate(psi):
        denominator = mul(linear(2, 2*j+1), linear(2, 2*j+2))
        require(mul(denominator, shift(p)) == mul(numerator, p), 'formal regular doubled beta step')

L = mon(j=1)
a = scale(mul(linear(1, 1), linear(1, 2)), F(1, 6))
A = mul(linear(2, 3), linear(2, 4))
R = mul(linear(2, 1), linear(2, 2))
c = scale(mul(R, A), F(1, 360))
N = add(mon(2, c=13), mon(1, c=31), mon(c=16))
phi, psi = rows(1)
require(phi == [a, mon()], 'genuine height-one first row')
require(psi == [c, scale(A, F(1, 3)), mon()], 'genuine height-one doubled row')
model = [c, scale(mul(L, A), F(1, 3)), mon()]
numerator = mul(linear(2, 5), linear(2, 6))
for j, p in enumerate(model):
    denominator = mul(linear(2, 2*j+1), linear(2, 2*j+2))
    require(mul(denominator, shift(p)) == mul(numerator, p), 'global attenuated orbit exact beta step')
require(model[0] == psi[0] and model[-1] == psi[-1], 'complete genuine endpoint coefficients')
response = add(power(a, 2), scale(mul(model[1], a), -1), c)
require(response == scale(mul(a, add(N, scale(mul(L, A), -10))), F(1, 30)), 'original same-root response identity')
discriminant = add(power(model[1], 2), scale(c, -4))
require(discriminant == add(scale(mul(power(L, 2), power(A, 2)), F(1, 9)), scale(mul(R, A), F(-1, 90))),
        'global root discriminant identity')
require(add(A, scale(R, -1)) == linear(8, 10), 'strict finite-parameter root boundary')
require(add(scale(power(A, 2), F(1, 90)), scale(c, -4)) == scale(mul(A, linear(8, 10)), F(1, 90)),
        'root discriminant positive at squared attenuation one tenth')
require({(0, j): v for (i, j), v in discriminant.items() if i == 4} ==
        add(mon(j=2, c=F(16, 9)), mon(c=F(-8, 45))), 'negative eventual discriminant below threshold')
require(add(scale(A, 13), scale(N, -4)) == linear(58, 92), 'exact sign-threshold gap')
require(add(mul(shift(N), A), scale(mul(N, shift(A)), -1)) ==
        add(mon(2, c=58), mon(1, c=242), mon(c=240)), 'threshold strictly increases along beta steps')
require(F(13, 40)**2-F(1, 10) == F(9, 1600), 'strict global geometry versus sign gap')

# Characteristic quotient: the incoming paired factor survives every
# attenuation, although its residual sign can change.
paired = mul(linear(1, 1), linear(1, 2))
residual = scale(add(scale(mul(L, A), 10), scale(N, -1)), F(1, 180))
require(scale(response, -1) == mul(paired, residual), 'paired characteristic factor with arbitrary attenuation')
require(specialize_lambda(residual, F(13, 40)) == scale(linear(29, 46), F(1, 360)),
        'sharp sign boundary has strictly positive linear residual')

lam = F(8, 25)
q = add(mon(2), mon(1, c=-69), mon(c=-112))
require(lam**2 > F(1, 10) and lam < F(13, 40), 'rational global orbit lies in sharp gap')
require(specialize_lambda(response, lam) == scale(mul(a, q), F(1, 150)), 'rational global orbit full response')
require(specialize_lambda(residual, lam) == scale(q, F(-1, 900)), 'paired-factor residual detects failure')
require(eval_poly(q, F(0)) < 0 and eval_poly(q, F(70)) == -42 and eval_poly(q, F(71)) == 30,
        'unique positive response transition lies between seventy and seventy-one')
require(eval_poly(response, F(70), lam) == F(-5964, 25), 'first primitive-address response')
require(eval_poly(response, F(71), lam) == F(876, 5), 'next primitive-address response')
require(gcd(70, 3) == gcd(71, 3) == 1, 'both designated masses are first support returns')
lam_zero = F(10981, 34320)
require(lam_zero == eval_poly(N, F(70))/(10*eval_poly(A, F(70))), 'exact orbit double-cancellation attenuation')
require(lam_zero**2 > F(1, 10) and lam_zero < F(13, 40), 'double-cancellation orbit retains global root geometry')
require(eval_poly(response, F(70), lam_zero) == 0 and eval_poly(a, F(70)) == 852,
        'exact original first and model second zero')


def literal(H, z, mass):
    g = 3*H+z
    found = {}
    for nc in range(mass+1):
        for nb in range(mass-nc+1):
            na = mass-nb-nc
            if -3*H*na+z*nb+(6*H+3*z)*nc == 0:
                found[(na, nb, nc)] = F(factorial(mass), factorial(na)*factorial(nb)*factorial(nc))
    return found


actual_controls = []
for z in [1, 5, 70, 71]:
    g = z+3
    pp, qq = rows(1)
    first, second = literal(1, z, g), literal(1, z, 2*g)
    first_counts = [(z+2*j, 3-3*j, j) for j in range(2)]
    second_counts = [(2*z+2*j, 6-3*j, j) for j in range(3)]
    require(set(first) == set(first_counts) and set(second) == set(second_counts), 'complete named regular charge fibres')
    require([first[n]/comb(g, 1) for n in first_counts] == [eval_poly(p, F(z)) for p in pp], 'literal first multinomial normalization')
    require([second[n]/comb(2*g, 2) for n in second_counts] == [eval_poly(p, F(z)) for p in qq], 'literal doubled multinomial normalization')
    require(all(nc >= 0 for _, _, nc in second_counts), 'regular support has no deleted inverse channel')
    require(eval_poly(response, F(z), F(1)) < 0, 'genuine interior amplitude retains actual sign')
    require(eval_poly(model[1], F(z), lam) == lam*eval_poly(psi[1], F(z)), 'model interior amplitude explicitly differs')
    actual_controls.append({'z': z, 'g': g, 'support': [-3, z, 6+3*z],
                            'first_coefficients': list(map(str, [eval_poly(p, F(z)) for p in pp])),
                            'genuine_second_coefficients': list(map(str, [eval_poly(p, F(z)) for p in qq]))})
require(sum(literal(1, 1, 4).values()) == 8, 'inherited positive literal g control')
require(sum(v*(-1)**n[2] for n, v in literal(1, 1, 8).items()) == -224,
        'inherited actual doubled literal sign control')
require(eval_poly(response, F(70), F(1)) != 0, 'model double zero is not an actual doubled cancellation')


def falling(value, n):
    out = F(1)
    for j in range(n):
        out *= value-j
    return out


# One regular complementary boundary, retaining the whole old raw carry
# range before specialization. It is separate from the singular zero block.
h, r = 3, 2
H, z = h-r, 2*r+1
old_p = [F(factorial(2*h+1), factorial(3*h-3*j)*factorial(1+2*j))*falling(F(h-r), h-j)
         for j in range(h+1)]
old_q = {e:falling(F(2*h-2*r), 2*h-e)/F(factorial(6*h-3*e)*factorial(2+2*e))
         for e in range(-1, 2*h+1)}
pp, qq = rows(H)
require(old_p == [F(0)]*r+[eval_poly(p, F(z)) for p in pp], 'complete first complementary reflection')
require(all(old_q[e] == 0 for e in range(-1, 2*r)), 'old lower carry specialized before complementary comparison')
require([old_q[e]*factorial(4*h+2) for e in range(2*r, 2*h+1)] == [eval_poly(p, F(z)) for p in qq],
        'complete doubled complementary reflection normalization')
require(gcd(z, 3*H) == 1 and z == 5, 'named reflected complement has primitive mass')

# At old h=1,x=-1, cancel the inverse-carry ratio in the generic quotient
# first. Naively evaluating the raw specialized row at t=0 is different.
x = F(-1)
old_quotient = (x+1)**2/720-(2*x+2)*(x+1)/144+(2*x+2)*(2*x+1)/1440-(2*x)*(2*x+1)/181440
require(old_quotient == F(-1, 90720), 'inherited singular inverse-carry quotient survives')
raw_specialized = [falling(F(0), 2-e)/F(factorial(6-3*e)*factorial(2+2*e)) for e in range(-1, 3)]
require(raw_specialized == [0, 0, 0, F(1, 720)], 'raw singular row is t squared over720')
require(old_quotient != 0, 'singular and nonzero complementary blocks not identified')

bank = {'universe': 'formal step/threshold identities; global H1 orbit; named literal and reflection controls',
        'symbolic_step_heights': [1, 2, 6], 'root_threshold': '1/sqrt(10)', 'sign_threshold': '13/40',
        'attenuation': str(lam), 'responses_70_71': ['-5964/25', '876/5'],
        'double_zero_attenuation': str(lam_zero), 'double_zero_phase': '-852',
        'paired_residual': '(-z^2+69z+112)/900', 'actual_controls': actual_controls,
        'singular_quotient': str(old_quotient), 'reflected_old_height_r': [h, r]}
print('STATUS exact conserved-amplitude and sharp H1 orbit thresholds; qualitative propagation REFUTED; actual all-height sign OPEN')
print(json.dumps(bank, sort_keys=True, separators=(',', ':')))
print('gates', GATES)
print('semantic_sha256', hashlib.sha256(json.dumps(bank, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
