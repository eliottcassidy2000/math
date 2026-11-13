from math import comb
from fractions import Fraction as F
for q in range(5):
    for b in range(1, 6):
        x, z = (1, 0)
        for n in range(30):
            expected = sum((comb(n + j, n - j) * q ** (n - j) * b ** j for j in range(n + 1)))
            if not x == expected:
                raise RuntimeError((q, b, n))
            x, z = ((q + b) * x + q * b * z, x + q * z)
D = 645

def add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def neg(x):
    return (-x[0], -x[1])

def mul(x, y):
    return (x[0] * y[0] + D * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def inv(x):
    d = x[0] ** 2 - D * x[1] ** 2
    return (x[0] / d, -x[1] / d)

def scale(a, x):
    return (a * x[0], a * x[1])
one = (F(1), F(0))
U = (F(127, 2), F(5, 2))
Up = (U[0], -U[1])
if not mul(U, Up) == one:
    raise RuntimeError('mul(U, Up) == one')
U3 = mul(mul(U, U), U)
if not add(U3, inv(U3)) == (F(2048002), F(0)):
    raise RuntimeError('add(U3, inv(U3)) == (F(2048002), F(0))')
t = scale(-432, U3)
y = inv(t)
j = add(add(t, (F(864), F(0))), scale(432 ** 2, inv(t)))
if not j == (-F(960 ** 3), F(0)):
    raise RuntimeError('j == (-F(960 ** 3), F(0))')
f = mul(add(one, scale(432, y)), inv(add(one, scale(-432, y))))
if not f == (F(0), F(320, 8127)):
    raise RuntimeError('f == (F(0), F(320, 8127))')
if not scale(36 * 5418, f) == (F(0), F(24 * 320)):
    raise RuntimeError('scale(36 * 5418, f) == (F(0), F(24 * 320))')
if not 36 * (F(263) - F(5418, 2)) == 24 * -3669:
    raise RuntimeError('36 * (F(263) - F(5418, 2)) == 24 * -3669')
a = lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k)
b = [sum((a(j) * comb(n + j, n - j) * (-432) ** (n - j) for j in range(n + 1))) for n in range(4)]
if not b == [1, -312, 114264, -44196288]:
    raise RuntimeError('b == [1, -312, 114264, -44196288]')
print('PASS:750 positive transfer checks; exact quadratic norm, trace, j and marked-observer identities; initial companion coefficients.')
