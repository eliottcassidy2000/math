# Direct check of the two inequalities behind the bilinear bound.
# (1) Lemma: sum_{0<|d|<D} min(K, 1/(2||beta d||))^2 <= 2 (D/q + 1)(3K^2 + 2qK) when |beta - r/q| <= 1/q^2.
# (2) Core: sum_{u in U_p, w in W_eps} |F(u,w)|^2 <= 2^L K_u K_w (1 + R), R = 2(2^L/q+1) sqrt((3+2q/K_u)(3+2q/K_w)).
import math, random, cmath
from fractions import Fraction
def norm(x): return abs(x - round(x))
random.seed(7)
bad = 0; worst = 0.0
for t in range(3000):
    q = random.randint(1, 400); r = random.randint(0, q - 1)
    while math.gcd(r, q) != 1: r = random.randint(0, q - 1)
    beta = r / q + random.uniform(-1, 1) / q**2
    K = random.randint(1, 2000); D = random.randint(2, 3000)
    lhs = sum(min(K, 1 / (2 * norm(beta * d)) if norm(beta * d) > 0 else float('inf'))**2 for d in range(1, D)) * 2  # +-d symmetric
    rhs = 2 * (D / q + 1) * (3 * K * K + 2 * q * K)
    worst = max(worst, lhs / rhs); bad += lhs > rhs
print(f'lemma (1): violations {bad} of 3000; max lhs/rhs = {worst:.3f}')
def Hw(w):
    H, A, p = 1, 0, 0
    for ch in w:
        if ch == '1': A = 2 * A + 1; p = 1
        else: H += A + p; A = 2 * A; p = 0
    return H
def Aw(w): return int(w, 2) if w else 0
def Cw(w): return sum((1 - int(w[j])) << j for j in range(len(w)))
def pw(w): return int(w[-1]) if w else 0
def ew(w): return 1 if (w and w[0] == '0') else 0
def words(m): return [format(i, '0%db' % m) for i in range(1 << m)] if m else ['']
bad2 = 0; out = []
for (a, L, l) in [(5, 4, 5), (6, 5, 6), (7, 4, 7), (6, 6, 6), (8, 3, 8)]:
    V = words(L)
    for theta in (Fraction(1, 7), Fraction(3, 37), Fraction(5, 101), (math.sqrt(5) - 1) / 2, Fraction(11, 257)):
        th = float(theta)
        q = (2 * theta).denominator if isinstance(theta, Fraction) else Fraction(2 * th).limit_denominator(1 << ((a + L + l) // 2)).denominator
        for p in (0, 1):
            for eps in (0, 1):
                Us = ['1' + format(j, '0%db' % (a - 2)) + str(p) for j in range(1 << (a - 2))]
                Ws = [str(1 - eps) + rest for rest in words(l - 1)]
                tot = 0.0
                for u in Us:
                    for w in Ws:
                        F = sum(cmath.exp(2j * math.pi * th * (Hw(v) + p * ew(v) + pw(v) * eps + Aw(u) * Cw(v) + Aw(v) * Cw(w))) for v in V)
                        tot += abs(F)**2
                Ku, Kw = len(Us), len(Ws)
                R = 2 * (2**L / q + 1) * math.sqrt((3 + 2 * q / Ku) * (3 + 2 * q / Kw))
                bnd = 2**L * Ku * Kw * (1 + R)
                bad2 += tot > bnd
                out.append((a, L, l, q, tot / (2**L * Ku * Kw), 1 + R))
for o in out[::4]: print(f'a={o[0]} L={o[1]} l={o[2]} q={o[3]:>4}: sum|F|^2 / diag = {o[4]:.3f}  <= 1+R = {o[5]:.2f}')
print('core (2): violations', bad2, 'of', len(out))
