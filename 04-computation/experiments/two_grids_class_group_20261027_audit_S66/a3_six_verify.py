"""Exact verification of the representation of 6 found by the audit's B = 800 search (two independent evaluations)."""
import sympy as sp

reps = {
    'rep A': [(12, -149), (150, 541), (123, 427), (297, -309)],
    'rep B (complex conjugate of A)': [(304, 427), (606, 309), (161, 149), (391, 541)],
}
w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2          # omega = exp(2 pi i / 3)


def mul(x, y):                                          # (a + b w)(c + d w), w^2 = -1 - w
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


for name, R in reps.items():
    # 1. integer arithmetic in the basis (1, w)
    sx = sy = 0
    terms = []
    for (a, b) in R:
        s = mul((a, b), (a, b))
        f = mul(s, s)
        terms.append(f)
        sx += f[0]
        sy += f[1]
    # 2. symbolic complex arithmetic with w = e^(2 pi i/3)
    val = sp.nsimplify(sp.expand(sum((a + b * w) ** 4 for a, b in R)))
    print(f'{name}: bases {R}')
    print(f'   fourth powers (x, y) = x + y w: {terms}')
    print(f'   sum in basis (1, w): ({sx}, {sy});  symbolic value: {val};  max |a|,|b| = {max(max(abs(a), abs(b)) for a, b in R)}')
    print(f'   norms N(alpha) = a^2 - ab + b^2: {[a * a - a * b + b * b for a, b in R]}')
# conjugation check: conj(a + b w) = (a - b) - b w
A = reps['rep A']
conjA = sorted(tuple(abs(c) for c in ((a - b), -b)) if False else ((a - b), -b) for a, b in A)
print('conj of rep A (up to sign):', [(a - b, -b) for a, b in A])
# residues mod 2 and mod 16 of the four bases (class pattern for n = 6 = 6 mod 16)
for a, b in A:
    s = mul((a, b), (a, b)); f = mul(s, s)
    print(f'   base ({a},{b}) = ({a % 2},{b % 2}) mod 2; fourth power mod 16 = ({f[0] % 16},{f[1] % 16})')
