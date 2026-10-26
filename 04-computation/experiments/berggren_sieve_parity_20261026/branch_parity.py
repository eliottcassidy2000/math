# The F8 branch parity (parity of the number of middle steps A in the Berggren word) at the integer level.
#   (P1) companion formula: with M_w the 2x2 word matrix (U=[2 -1;1 0], A=[2 1;1 0], D=[1 2;0 1]) and
#        (m,n) = M_w (2,1), (al,ga) = M_w (1,0) (the image of the pole t = infinity), (-1)^#A = al*n - ga*m,
#        0 <= al <= m, 0 <= ga <= n;
#   (P1') regular-continued-fraction automaton: #A = number of even heads in the digit recoding of m/n;
#   (P2) every-modulus obstruction: (m,n) = (16M+2, 4M+1) is congruent to (2,1) mod M, has word A D^(M-1) U D
#        and odd branch parity; so no function of the triple modulo any fixed M determines the parity;
#   (P3) sign(a-b) sees only the trailing A-run: Q(child) = -Q(parent) under A, = (a+b)(parent) under D,
#        = -(a+b)(parent) under U, Q(m,n) = m^2-2mn-n^2 = a-b = N((m-n)+n*sqrt2); on the A ray (-1)^#A = -(a-b).
import sys
from math import gcd
from fractions import Fraction
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1500

def descent(m, n):
    w = []
    while (m, n) != (2, 1):
        if m > 3 * n: w.append('D'); m -= 2 * n
        elif m > 2 * n: w.append('A'); m, n = n, m - 2 * n
        else: w.append('U'); m, n = n, 2 * n - m
    return ''.join(reversed(w))            # root-first application order

MAT = {'U': (2, -1, 1, 0), 'A': (2, 1, 1, 0), 'D': (1, 2, 0, 1)}
def apply(word, v):
    x, y = v
    for ch in word:
        a, b, c, d = MAT[ch]; x, y = a * x + b * y, c * x + d * y
    return x, y

def cf(m, n):
    d = []
    while n:
        d.append(m // n); m, n = n, m % n
    return d

def automaton_A(m, n):
    d = cf(m, n)
    if d[-1] == 1 and len(d) > 1: d = d[:-2] + [d[-2] + 1]      # canonical last digit >= 2
    count = 0
    while True:
        if d == [2]: return count
        h = d[0]
        if h % 2 == 0:
            if len(d) == 1: return count                       # t = h even integer: only D steps remain
            count += 1; d = d[1:]
        else:
            if len(d) == 1: raise ValueError('odd integer slope cannot occur for opposite-parity coprime m, n')
            if len(d) == 2: return count                       # [h; c1] -> ... -> [1; 1] = 2
            d = [d[2] + 1] + d[3:]

ok1 = ok1b = ok2 = ok3 = tot = 0; box = True
for m in range(2, MMAX + 1):
    for n in range(1, m):
        if (m - n) % 2 == 0 or gcd(m, n) != 1: continue
        tot += 1
        w = descent(m, n); nA = w.count('A')
        assert apply(w, (2, 1)) == (m, n)
        al, ga = apply(w, (1, 0))
        ok1 += (al * n - ga * m) == (-1) ** nA
        box &= 0 <= al <= m and 0 <= ga <= n
        ok1b += automaton_A(m, n) == nA
        # (P3) trailing A-run
        a, b = m * m - n * n, 2 * m * n
        r = len(w) - len(w.rstrip('A')); rest = w[:len(w) - r]
        base = (+1 if rest and rest[-1] == 'D' else -1)        # D child: a-b = (a+b)(parent) > 0; U child / root: negative
        ok3 += (1 if a > b else -1) == (-1) ** r * base
print(f"(P1) companion determinant = (-1)^#A on {ok1}/{tot} primitive triples with m <= {MMAX}; companion in box: {box}")
print(f"(P1') continued-fraction automaton = #A on {ok1b}/{tot}")
print(f"(P3) sign(a - b) = (-1)^(trailing A-run) * (+1 after D, -1 after U or at the root) on {ok3}/{tot}")
# (P2)
bad = 0
for M in range(1, 2001):
    m, n = 16 * M + 2, 4 * M + 1
    w = descent(m, n)
    good = (m % M, n % M) == (2 % M, 1 % M) and w == 'A' + 'D' * (M - 1) + 'UD' and w.count('A') == 1
    a, b, c = m * m - n * n, 2 * m * n, m * m + n * n
    good &= (a % M, b % M, c % M) == (3 % M, 4 % M, 5 % M)
    bad += not good
print(f"(P2) M = 1..2000: (16M+2, 4M+1) == (2,1) mod M, word A D^(M-1) U D, one middle step, triple == (3,4,5) mod M; failures {bad}")
print("    M = 7 example:", (16 * 7 + 2, 4 * 7 + 1), [(16 * 7 + 2) ** 2 - (4 * 7 + 1) ** 2, 2 * (16 * 7 + 2) * (4 * 7 + 1), (16 * 7 + 2) ** 2 + (4 * 7 + 1) ** 2])
# A ray: Pell sign
x = (2, 1); out = []
for k in range(8):
    m, n = x; out.append((m, n, m * m - n * n - 2 * m * n, (-1) ** k)); x = apply('A', x)
print("(P3) A ray (m, n, a - b, (-1)^#A):", out)
