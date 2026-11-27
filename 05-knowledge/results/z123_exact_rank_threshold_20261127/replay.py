"""Exact certificate, Python 3 standard library. No matroid enumeration required."""
from functools import cache
from itertools import combinations, product
from math import comb

KEYS = ((), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3))
FANO = {(): 1, (1,): 7, (2,): 7, (1, 2): 21}
COLOOP = {(): 1}

def packing(q):
    m = q*q + q + 1
    P, L, B = q**4, q**3 * (q**3 + q*q + q + 1), q*q * (q*q + 1) * m
    D = P * comb(m, 2) + comb(m, 3)
    return dict(zip(KEYS, (
        1, P + m, L + P*m + comb(m, 2), B + D,
        q*L + 2*P*m + 2*comb(m, 2),
        q*q*B + m*q*q*(q*q+1) + 3*D,
        (q*q+q)*B + m*q**4*(q*q+1) + 3*D,
        q*(q*q+q)*B + 2*m*q**4*(q*q+1) + 6*D,
    )))

def rank_splits(H, r, s):
    for I in product(range(r+1), repeat=len(H)):
        J = tuple(h-i for h, i in zip(H, I))
        if tuple(sorted(I)) != I or tuple(sorted(J)) != J:
            continue
        if any(j < 0 or j > s for j in J):
            continue
        yield (tuple(sorted({i for i in I if 0 < i < r})),
               tuple(sorted({j for j in J if 0 < j < s})))

@cache
def terms(H, r, s):
    return tuple(rank_splits(H, r, s))

def total(f, g, r, s):
    def count(*H):
        return sum(f[I]*g[J] for I, J in terms(H, r, s))
    return count

def flags(f, r):
    return {I: f(*I) for l in range(r) for I in combinations(range(1, r), l)}

def from_flags(f, r):
    a = f(r-1)
    b = f(r-2) + f(r-3, r-1) - f(r-3, r-2)
    if r == 5:
        return a, b
    c = (f(3) + f(2, 5) - f(2, 3) + f(1, 4) - f(1, 3)
         + f(1, 3, 5) - f(1, 3, 4) - f(1, 2, 5) + f(1, 2, 3))
    return a, b, c

def direct(f, col):
    _, w, y, z, u, v, s, t = [f[K] for K in KEYS]
    i = f[()]
    if col:
        return (w + 8*y + 14*z + 8*i,
                21*u + 14*v + 27*s - 13*w + 30*y + 21*z - 6*i,
                49*u + 14*t + 34*v + 40*s - 42*w + 37*y - 12*z - 5*i)
    return (y + 7*z + 7*i, -7*w + y + 13*z + 7*v + 14*s + i,
            -19*w + 8*y - 6*z - 8*u + 19*v + 20*s + 8*t + 8*i)

@cache
def coeff(q, col=0):
    return direct(packing(q), col)

def reciprocal(a, b, c):
    return a, b-3, c-2*a

def disc(a, b, c):
    a,b,c = reciprocal(a,b,c)
    return a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c

def translated(a,b,c):
    a,b,c = reciprocal(a,b,c)
    return (1, a-6, b-4*a+12, c-2*b+4*a-8)

def newton(f, q, degree):
    vals = [f(q+i) for i in range(degree+2)]
    coeffs = []
    while vals:
        coeffs.append(vals[0])
        vals = [y-x for x, y in zip(vals, vals[1:])]
    assert coeffs.pop() == 0
    return coeffs

def nonzero_sign(vals, sign):
    assert all(sign*v > 0 for v in vals)

# Compare the closed formulas to independent coordinate splits on a full
# integral basis of block-flag values (no assumption about q).
for K in KEYS:
    f = {J: int(J == K) for J in KEYS}
    g = total(f, FANO, 4, 3)
    assert from_flags(g, 6) == direct(f, 0)
    assert from_flags(total(flags(g, 7), COLOOP, 7, 1), 6) == direct(f, 1)
    _, w, y, z, u, v, s, t = [f[J] for J in KEYS]
    i = f[()]
    assert from_flags(g, 5) == (w+7*y+7*z+i, -7*w+15*y+z+14*u+7*v+6*s-6*i)

def family(col):
    def seq(q):
        return coeff(q, col)
    def discriminant(q):
        return disc(*seq(q))
    def edge(q):
        return translated(*seq(q))[-1]
    def gap(q):
        a,b,c = seq(q)
        return b*b-a*c
    def gap1(q):
        a,b,c = seq(q)
        return a*a-b
    def middle(q):
        a,b,c = seq(q)
        return c-b
    root_switch = next(q for q in range(2, 70) if discriminant(q) < 0)
    gap_switch = next(q for q in range(2, 300) if gap(q) < 0)
    assert all(discriminant(q) > 0 and min(translated(*seq(q))) > 0
               for q in range(2, root_switch))
    assert all(gap(q) > 0 for q in range(2, gap_switch))
    for f, shift, degree, sign in (
        (discriminant, root_switch, 33, -1), (edge, root_switch, 9, 1),
        (gap, gap_switch, 17, -1), (gap1, 2, 16, 1), (middle, 2, 9, 1),
    ):
        nonzero_sign(newton(f, shift, degree), sign)
    print(f'coloop={col}: integer root switch {root_switch}; integer LC switch {gap_switch}')
    for q in ((29, 173) if not col else (31, 181)):
        a,b,c = seq(q)
        print('  q', q, 'size', q**4+q*q+q+8+col, 'a,b,c', a,b,c, 'G', gap(q))
    from fractions import Fraction
    from math import factorial
    leading = tuple(Fraction(newton(lambda q: seq(q)[i], 2, d)[-1],factorial(d))
                    for i,d in enumerate((8,8,9)))
    print('  leading', *leading)

family(0)
family(1)

def ctrl(q):
    a,b = from_flags(total(packing(q), FANO, 4, 3), 5)
    return a-5, b-3*a+5, (a+1)**2 + 4 - 4*b

for i,d in enumerate((8,8,16)):
    nonzero_sign(newton(lambda q: ctrl(q)[i], 2, d), 1)
print('rank-five affine/Fano control: positive gamma coefficients and discriminant')
