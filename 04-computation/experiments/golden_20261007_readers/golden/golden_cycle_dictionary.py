#!/usr/bin/env python3
"""Golden reading of Collatz T-periodic points (T(x) = x/2 or 3x+1, THM-4528 convention).

For a periodic parity word w = b_0..b_{n-1} (no '11' cyclically) the periodic point is the rational
x_w = B/(2^A - 3^r) and Theta(x_w) = P_w(phi)/(phi^n - 1) with P_w(phi) = sum_j b_j phi^(n-1-j).
kappa_n(x) := class of P_w(phi) in R_n = Z[phi]/(phi^n - 1).  kappa_n(T x) = phi * kappa_n(x).

Checks: equivariance, the fibre structure of kappa_n (bijective except {0,-1,-2} -> 0 for even n),
the n = 5 labels in F_11 (paper's H / -H), the n = 10 picture (paper's Prop 5.4 / Thm 6.1),
and the -17 cycle in R_18 = Z[phi]/(76).
"""
from fractions import Fraction
from itertools import product

def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a

def phi_pow(k):
    """phi^k = F_{k-1} + F_k phi for k >= 1; phi^0 = 1."""
    if k == 0:
        return (1, 0)
    return (fib(k - 1), fib(k))

def mul(x, y):
    a, b = x; c, d = y
    return (a * c + b * d, a * d + b * c + b * d)

def golden_words(n):
    """All cyclic words of length n over {0,1} with no two cyclically adjacent 1s."""
    out = []
    for w in product((0, 1), repeat=n):
        if n == 1:
            if w == (1,):
                continue   # word '1' then '1' adjacent cyclically
        ok = all(not (w[i] == 1 and w[(i + 1) % n] == 1) for i in range(n))
        if ok:
            out.append(w)
    return out

def periodic_point(w):
    """Rational x with T-parity word w^infty (T(x)=x/2 if even, 3x+1 if odd). Returns Fraction."""
    # compose affine maps: x -> 3x+1 (b=1) , x -> x/2 (b=0)
    a, c = Fraction(1), Fraction(0)   # map x -> a x + c
    for b in w:
        if b:
            a, c = 3 * a, 3 * c + 1
        else:
            a, c = a / 2, c / 2
    # fixed point x = a x + c
    return c / (1 - a)

def Tmap(x):
    # 2-adic parity of a rational with odd denominator
    num, den = x.numerator, x.denominator
    assert den % 2 == 1
    if num % 2 == 0:
        return x / 2
    return 3 * x + 1

def parity(x):
    return x.numerator % 2

def Pw(w):
    n = len(w)
    s = (0, 0)
    for j, b in enumerate(w):
        if b:
            p = phi_pow(n - 1 - j)
            s = (s[0] + p[0], s[1] + p[1])
    return s

class Rn:
    def __init__(self, n):
        self.n = n
        g = phi_pow(n)
        self.a, self.b = g[0] - 1, g[1]           # phi^n - 1 = a + b phi
        a, b = self.a, self.b
        self.D = a * (a + b) - b * b              # signed norm
    def red(self, x):
        m, k = x; a, b, D = self.a, self.b, self.D
        s = ((a + b) * m - b * k) // D
        t = (-b * m + a * k) // D
        return (m - s * a - t * b, k - s * b - t * (a + b))

def main():
    print('n  #words  |R_n|  #classes  max-fibre  collapsed-fibre(points)                 equivariant')
    for n in range(1, 23):
        R = Rn(n)
        W = golden_words(n)
        cls = {}
        equi = True
        for w in W:
            x = periodic_point(w)
            # verify the word really is x's parity word
            y = x
            for j in range(n):
                assert parity(y) == w[j], (w, x)
                y = Tmap(y)
            assert y == x
            c = R.red(Pw(w))
            cls.setdefault(c, []).append(x)
            # equivariance: class of rotated word == phi * class
            w2 = w[1:] + w[:1]
            c2 = R.red(Pw(w2))
            if c2 != R.red(mul((0, 1), c)):
                equi = False
        fibres = [v for v in cls.values() if len(v) > 1]
        mx = max(len(v) for v in cls.values())
        note = '; '.join(str(sorted(v)) for v in fibres)
        print(f'{n:2d} {len(W):6d} {abs(R.D):6d} {len(cls):8d} {mx:6d}   {note[:48]:48s} {equi}')
        assert len(cls) == abs(R.D)

    # n = 5 labels: R_5 = O/(phi^5-1) = O/(4-phi) -> F_11, m + k phi -> m + 4k
    print('\nn = 5: period-5 rational T-points and their F_11 labels (paper: H={1,3,4,5,9} gold, -H teal)')
    H = {1, 3, 4, 5, 9}
    R5 = Rn(5)
    seen = set()
    for w in golden_words(5):
        x = periodic_point(w)
        c = R5.red(Pw(w))
        lab = (c[0] + 4 * c[1]) % 11
        cl = 'H' if lab in H else ('-H' if lab else '0')
        print(f'   word {"".join(map(str, w))}  x = {str(x):7s}  kappa = {c}  label {lab:2d}  class {cl}')
    # near-perfect matching M_0 in Collatz terms
    labels = {}
    for w in golden_words(5):
        x = periodic_point(w); c = R5.red(Pw(w)); labels[(c[0] + 4 * c[1]) % 11] = x
    print('   M_0 = {{-mu, mu} : mu in H}:', [(str(labels[m]), str(labels[(-m) % 11])) for m in sorted(H)])
    for alpha in range(11):
        Ma = [tuple(sorted(((alpha - m) % 11, (alpha + m) % 11))) for m in H]
        if alpha in (0, 8):
            print(f'   M_{alpha} (missing {labels[alpha]}):', [(str(labels[u]), str(labels[v])) for u, v in Ma])

    # n = 10: (alpha, beta) coordinates; period-5 points have beta = 0; T^5 acts as (alpha,beta)->(alpha,-beta)
    print('\nn = 10: R_10 = O/(11) = F_11 x F_11 via (alpha,beta) = (m+4k, m+8k)')
    R10 = Rn(10)
    pts = {}
    for w in golden_words(10):
        x = periodic_point(w)
        c = R10.red(Pw(w))
        al, be = (c[0] + 4 * c[1]) % 11, (c[0] + 8 * c[1]) % 11
        pts.setdefault((al, be), []).append(x)
    beta0 = sorted((k, [str(v) for v in vs]) for k, vs in pts.items() if k[1] == 0)
    print('   beta = 0 classes:', len(beta0))
    for k, vs in beta0:
        print('     ', k, vs)
    # T^5 swap check on primitive period-10 points
    ok = True; npairs = 0
    for (al, be), vs in pts.items():
        if be == 0:
            continue
        for x in vs:
            y = x
            for _ in range(5):
                y = Tmap(y)
            # find class of y
            found = [k for k, v in pts.items() if y in v]
            if found != [(al, (-be) % 11)]:
                ok = False
        npairs += 1
    print('   T^5 sends (alpha,beta) -> (alpha,-beta) on all 110 primitive points:', ok, '; classes with beta != 0:', npairs)
    # Psi_- correspondence: [alpha,beta] -> {alpha + 1/beta, alpha - 1/beta}; labels of period-5 points: alpha = 2*label
    inv2 = pow(2, -1, 11)
    print('   example Psi_- images (pairs of period-5 points):')
    shown = 0
    for (al, be), vs in sorted(pts.items()):
        if be == 0 or be > 5:
            continue
        bi = pow(be, -1, 11)
        e1, e2 = (al + bi) % 11, (al - bi) % 11
        # convert R_10 alpha-labels to R_5 labels: alpha = 2*label
        l1, l2 = e1 * inv2 % 11, e2 * inv2 % 11
        if shown < 6:
            print(f'     {{{str(vs[0])}, T^5}} -> {{{labels[l1]}, {labels[l2]}}}')
            shown += 1

    # n = 18: the -17 cycle
    print('\nn = 18: R_18 = O/(phi^18 - 1) = O/(76) since phi^18 - 1 = L_9 phi^9 = 76 phi^9')
    g = phi_pow(18); g9 = phi_pow(9)
    print('   phi^18 - 1 =', (g[0] - 1, g[1]), ' 76*phi^9 =', (76 * g9[0], 76 * g9[1]))
    R18 = Rn(18)
    x0 = Fraction(-17)
    w = []
    y = x0
    for _ in range(18):
        w.append(parity(y)); y = Tmap(y)
    assert y == x0
    c = R18.red(Pw(tuple(w)))
    print('   -17 word', ''.join(map(str, w)), ' kappa =', c, ' |R_18| =', abs(R18.D))
    # 19-component: phi -> roots 5 and 15 mod 19
    for r in (5, 15):
        vals = []
        ww = tuple(w)
        for i in range(18):
            wr = ww[i:] + ww[:i]
            cc = Pw(wr)
            vals.append((cc[0] + r * cc[1]) % 19)
        print(f'   19-component (phi -> {r}, order {next(k for k in range(1,19) if pow(r,k,19)==1)}):', vals)
    # Theta(-17) exact
    num = Pw(tuple(w))
    print('   Theta(-17) = P(phi)/(phi^18-1) = P(phi)/(76 phi^9); P =', num)

if __name__ == '__main__':
    main()
