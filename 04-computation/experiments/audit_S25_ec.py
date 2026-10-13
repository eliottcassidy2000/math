"""
audit_S25_ec.py -- independent re-derivation of the elliptic-curve identifications in section 5.4 of the S25 note.
Exact integer arithmetic (Fractions); no PARI.

(1) b/c-invariants, discriminant, j of: the fruit curve y^2 = x^3 + 109x^2 + 224x; LMFDB 910.a4 = [1,0,1,-234,1352];
    the square curve y^2 = 2x^3 + 4x^2 + 1 transformed CORRECTLY (multiply by 4, X = 2x, Y = 2y: Y^2 = X^3 + 4X^2 + 4);
    the model the note claims (Y^2 = X^3 + 4X^2 + 16); LMFDB 43.a1 = [0,1,1,0,0]; LMFDB 364.a1 = [0,1,0,-5,7].
(2) Q-isomorphism tests via (c4, c6) scaling by u^4, u^6 with u in Q.
(3) the substitution checks: what Y = 4y, X = 2x actually does to y^2 = 2x^3 + 4x^2 + 1.
(4) torsion: fruit curve, 3*(56,728) = (0,0) with the exact group law (order 6); rational 2-torsion count; point counts mod
    good primes for the torsion bound; square curve: (10,49) is a point, its image (20,98) is non-torsion (Nagell-Lutz), so rank >= 1;
    no rational 2-torsion; torsion bound from point counts.
"""
from fractions import Fraction as Fr
from math import gcd

def invariants(a):
    a1, a2, a3, a4, a6 = a
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
    D = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    assert c4 ** 3 - c6 ** 2 == 1728 * D
    return dict(b2=b2, b4=b4, b6=b6, b8=b8, c4=c4, c6=c6, D=D, j=Fr(c4 ** 3, D))

def factor(n):
    n = abs(n); out = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p; e += 1
        if e: out.append((p, e))
        p += 1
    if n > 1: out.append((n, 1))
    return out

def isomorphic_over_Q(I, J):
    """is there u in Q^x with J.c4 = u^4 I.c4 and J.c6 = u^6 I.c6 ?"""
    if I['c4'] == 0 or I['c6'] == 0 or J['c4'] == 0 or J['c6'] == 0:
        return None
    r4 = Fr(J['c4'], I['c4']); r6 = Fr(J['c6'], I['c6'])
    u2 = r6 / r4                      # u^2
    if u2 <= 0: return False
    # u^2 must be a rational square and u^4 = r4
    num, den = u2.numerator, u2.denominator
    sn, sd = int(round(num ** 0.5)), int(round(den ** 0.5))
    if sn * sn != num or sd * sd != den: return False
    return u2 * u2 == r4 and (u2 ** 3 == r6), Fr(sn, sd)

curves = {
    'fruit y^2 = x^3 + 109x^2 + 224x': (0, 109, 0, 224, 0),
    'LMFDB 910.a4 [1,0,1,-234,1352]': (1, 0, 1, -234, 1352),
    'square, correct model Y^2 = X^3 + 4X^2 + 4': (0, 4, 0, 0, 4),
    'note claim Y^2 = X^3 + 4X^2 + 16': (0, 4, 0, 0, 16),
    'LMFDB 43.a1 [0,1,1,0,0]': (0, 1, 1, 0, 0),
    'LMFDB 364.a1 [0,1,0,-5,7]': (0, 1, 0, -5, 7),
}
inv = {}
print('# (1) invariants')
for name, a in curves.items():
    I = invariants(a); inv[name] = I
    print(f'  {name}: c4 = {I["c4"]}, c6 = {I["c6"]}, Delta = {I["D"]} = {factor(I["D"])} (sign {"-" if I["D"] < 0 else "+"}), j = {I["j"]} (denominator {factor(I["j"].denominator)})')

print('\n# (2) Q-isomorphism by (c4, c6) scaling (u^4, u^6)')
pairs = [('fruit y^2 = x^3 + 109x^2 + 224x', 'LMFDB 910.a4 [1,0,1,-234,1352]'),
         ('square, correct model Y^2 = X^3 + 4X^2 + 4', 'LMFDB 364.a1 [0,1,0,-5,7]'),
         ('square, correct model Y^2 = X^3 + 4X^2 + 4', 'LMFDB 43.a1 [0,1,1,0,0]'),
         ('note claim Y^2 = X^3 + 4X^2 + 16', 'LMFDB 43.a1 [0,1,1,0,0]')]
for A, B in pairs:
    print(f'  {A}  ~  {B}: {isomorphic_over_Q(inv[A], inv[B])}')

print('\n# (3) substitution checks on y^2 = 2x^3 + 4x^2 + 1')
# X = 2x, Y = 4y  =>  Y^2 = 16 y^2 = 16(2x^3 + 4x^2 + 1) = 32 x^3 + 64 x^2 + 16 = 4 X^3 + 16 X^2 + 16   (not monic)
# X = 2x, Y = 2y  =>  Y^2 = 4 y^2 = 8x^3 + 16x^2 + 4 = X^3 + 4 X^2 + 4
import random
for _ in range(3):
    x = Fr(random.randint(-50, 50), random.randint(1, 9))
    rhs = 2 * x ** 3 + 4 * x ** 2 + 1
    X = 2 * x
    print(f'  x = {x}: 16*(2x^3+4x^2+1) = {16*rhs}; X^3+4X^2+16 = {X**3 + 4*X**2 + 16}; 4X^3+16X^2+16 = {4*X**3 + 16*X**2 + 16}; '
          f'4*(2x^3+4x^2+1) = {4*rhs}; X^3+4X^2+4 = {X**3 + 4*X**2 + 4}')

# ---- group law for y^2 = x^3 + a2 x^2 + a4 x + a6 ----
def add(P, Q, a2, a4, a6):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and y1 == -y2: return None
    if P == Q:
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - a2 - x1 - x2
    y3 = lam * (x1 - x3) - y1
    return (x3, y3)

def on_curve(P, a2, a4, a6):
    x, y = P
    return y * y == x ** 3 + a2 * x * x + a4 * x + a6

def count_points_mod_p(a, p):
    a1, a2, a3, a4, a6 = a
    cnt = 1
    for x in range(p):
        for y in range(p):
            if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % p == 0:
                cnt += 1
    return cnt

print('\n# (4) torsion and points')
a2, a4, a6 = 109, 224, 0
T = (Fr(56), Fr(728)); G = (Fr(-4), Fr(28))
print(f'  fruit: (56,728) on curve: {on_curve(T, a2, a4, a6)}; (-4,28) on curve: {on_curve(G, a2, a4, a6)}')
P = None
orders = []
for m in range(1, 13):
    P = add(P, T, a2, a4, a6)
    orders.append((m, P))
print('  multiples of (56,728):', [(m, (str(P[0]), str(P[1])) if P else 'O') for m, P in orders[:6]])
# rational roots of x^3 + 109x^2 + 224x: x = 0 and x^2 + 109x + 224 (disc 109^2 - 4*224)
disc = 109 * 109 - 4 * 224
print(f'  fruit 2-torsion: x(x^2+109x+224); disc of the quadratic = {disc} = {factor(disc)}; a square? {int(disc**0.5)**2 == disc}  -> exactly one rational 2-torsion point')
good = [p for p in (3, 11, 17, 19, 23, 29, 31) if inv['fruit y^2 = x^3 + 109x^2 + 224x']['D'] % p]
counts = [(p, count_points_mod_p(curves['fruit y^2 = x^3 + 109x^2 + 224x'], p)) for p in good]
g = 0
for p, c in counts: g = gcd(g, c)
print(f'  fruit #E(F_p) at good p: {counts}; gcd = {g}  (torsion order divides gcd; a point of order 6 exists => torsion = Z/6 iff gcd = 6)')
# square curve
print(f'  square: (10,49) on y^2 = 2x^3+4x^2+1: {49*49 == 2*1000 + 400 + 1}; image (20,98) on Y^2 = X^3+4X^2+4: {on_curve((Fr(20), Fr(98)), 4, 0, 4)}')
D = inv['square, correct model Y^2 = X^3 + 4X^2 + 4']['D']
print(f'  Nagell-Lutz: 98^2 = {98*98} divides |Delta| = {abs(D)}? {abs(D) % (98*98) == 0}  -> (20,98) is NOT torsion -> rank >= 1 (PROVED, elementary)')
roots = [r for r in range(-10, 11) if r ** 3 + 4 * r * r + 4 == 0]
print(f'  square rational 2-torsion: integer roots of X^3+4X^2+4 among divisors of 4: {roots} (none -> no 2-torsion)')
good = [p for p in (3, 5, 11, 17, 19, 23) if D % p]
counts = [(p, count_points_mod_p((0, 4, 0, 0, 4), p)) for p in good]
g = 0
for p, c in counts: g = gcd(g, c)
print(f'  square #E(F_p) at good p: {counts}; gcd = {g} -> torsion trivial iff gcd = 1')
print('done')
