"""54 = 9 * 6 = (sqrt-3)^4 * 6 with sqrt-3 = 1 + 2w: scale the representation of 6 by 1 + 2w and verify exactly."""
import sympy as sp

w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2


def mul(x, y):
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


six = [(12, -149), (150, 541), (123, 427), (297, -309)]
r3 = (1, 2)                                   # 1 + 2w = sqrt(-3)
s = mul(r3, r3); print('(1 + 2w)^2 =', s, ' (1 + 2w)^4 =', mul(s, s))
f54 = [mul(r3, a) for a in six]
print('bases for 54:', f54, ' max |coord| =', max(max(abs(x), abs(y)) for x, y in f54))
val = sp.nsimplify(sp.expand(sum((a + b * w) ** 4 for a, b in f54)))
sx = sy = 0
for a in f54:
    q = mul(a, a); f = mul(q, q); sx += f[0]; sy += f[1]
print('sum (basis 1, w):', (sx, sy), ' symbolic:', val)
