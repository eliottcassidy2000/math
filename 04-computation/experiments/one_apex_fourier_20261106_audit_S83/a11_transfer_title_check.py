# Audit S83, step 11: exact S_n(r/s) for large n by the (A mod s, last bit) automaton, to test THM-4886's TITLE clause
# "12 * 2^(n-1) (q/2^n)^(1/4) for q <= 2^(n-6)" (the body and the note restrict it to 2^((n-6)/3) <= q <= 2^(n-6)).
# Automaton (Session 193 section 2): bit 1: (H, A, p) -> (H, 2A + 1, 1); bit 0: (H + A + p, 2A, 0).  Words in X_n start
# with 1, i.e. state (1, 1, 1) after the first letter.  V[a, p] = sum of e(theta H) over prefixes with A = a mod s, last bit p.
import math, cmath
from fractions import Fraction
import numpy as np

def S_transfer(theta, n):
    r, s = theta.numerator, theta.denominator
    w = [cmath.exp(2j * math.pi * r * k / s) for k in range(s)]     # e(theta k), k mod s
    V = np.zeros((s, 2), dtype=complex)
    V[1 % s, 1] = w[1 % s]                                         # H = 1 after the leading 1
    for _ in range(n - 1):
        W = np.zeros_like(V)
        for a in range(s):
            for p in (0, 1):
                x = V[a, p]
                if x == 0:
                    continue
                W[(2 * a + 1) % s, 1] += x                         # append 1
                W[(2 * a) % s, 0] += x * w[(a + p) % s]            # append 0: H += A + p
        V = W
    return V.sum()

# cross-check against brute force for small n
from oa_engine import g_block, S_exact
for th in (Fraction(1, 7), Fraction(1, 3), Fraction(3, 16)):
    for n in (10, 16, 20):
        a = abs(S_transfer(th, n)); b = abs(S_exact(g_block(n), th))
        assert abs(a - b) < 1e-6 * 2 ** n, (th, n, a, b)
print('transfer evaluation agrees with brute force (theta = 1/7, 1/3, 3/16; n = 10, 16, 20)')
print(' theta  q    n   exact |S_n|/2^(n-1)   title bound 12 (q/2^n)^(1/4)   body (A1) bound 6.5 q^(-1/2)   title clause holds?')
for th in (Fraction(1, 7), Fraction(1, 3)):
    q = (2 * th).denominator
    first_fail = None
    for n in range(10, 161):
        ex = abs(S_transfer(th, n)) / 2 ** (n - 1)
        tb = 12 * (q / 2 ** n) ** 0.25
        if ex > tb and first_fail is None:
            first_fail = n
        if n in (26, 40, 60, 80, 100, 120, 140, 160) or n == first_fail:
            print(f' {str(th):>5s} {q:2d}  {n:3d}   {ex:.3e}             {tb:.3e}                    {6.5 / math.sqrt(q):.3f}                     {ex <= tb}')
    print(f'   first n with exact > 12 (q/2^n)^(1/4) at theta = {th}: {first_fail}  (q = {q} <= 2^(n-6) for every n >= 9)')
