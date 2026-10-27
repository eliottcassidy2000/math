"""A000568(n), strong counts, intransitivity (C3-structure) classes c(n), strong c_s(n), and labelled
intransitivity counts l(n) for n <= NMAX from the prime cycle-index TYPE census (files t3.txt .. t{K}.txt).
TYPE count=C k=K aut=G AUT: ct:mult ... SC=s ANTI: ct:mult ...   (ct = cycle type like 1^4,3^2, ; mult = elements)
Species (session56 / THM-4803): A = SEQ(S), S = X + sum_Q Z_{Aut Q}(A-1);  c = MSET(c_s), c_s = X + sum_{Q/dual} Z_{Aut+-Q}(c).
Labelled: l_s(x) = x + (1/2) sum_k P_lab(k) l(x)^k / k!  (P_lab(k) = labelled primes of order k), 1 + l = exp(l_s)."""
import sys, re, os
from fractions import Fraction
from math import factorial
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 11
DATA = os.environ.get("DATA", ".")
types = []   # (k, count, G, [(ctdict, mult)], sc, [(ctdict, mult)])
def parse_ms(s):
    res = []
    for tok in s.split():
        ct, mult = tok.rsplit(':', 1); d = {}
        for part in ct.strip(',').split(','):
            if part: l, e = part.split('^'); d[int(l)] = int(e)
        res.append((d, int(mult)))
    return res
for k in range(3, NMAX):
    fn = os.path.join(DATA, f't{k}.txt')
    if not os.path.exists(fn): continue
    for line in open(fn):
        m = re.match(r'TYPE count=(\d+) k=(\d+) aut=(\d+) AUT:(.*?) SC=(\d) ANTI:(.*)$', line.strip())
        if not m: continue
        cnt, kk, G, A, sc, anti = int(m.group(1)), int(m.group(2)), int(m.group(3)), parse_ms(m.group(4)), int(m.group(5)), parse_ms(m.group(6))
        assert kk == k and sum(mu for _, mu in A) == G
        types.append((k, cnt, G, A, sc, anti))
print("types by order:", {k: sum(1 for t in types if t[0] == k) for k in range(3, NMAX)})
print("primes by order:", {k: sum(t[1] for t in types if t[0] == k) for k in range(3, NMAX)})
def polya(ms, G, series, n):
    """[x^n] (1/G) sum_{(ct,mult)} mult * prod_l (sum_j series[j] x^{lj})^{c_l}"""
    tot = Fraction(0)
    for ct, mult in ms:
        poly = [Fraction(0)] * (n + 1); poly[0] = Fraction(1)
        for l, c in ct.items():
            base = [Fraction(0)] * (n + 1)
            for j in range(1, n + 1):
                if l * j <= n: base[l * j] = Fraction(series[j])
            for _ in range(c):
                new = [Fraction(0)] * (n + 1)
                for a in range(n + 1):
                    if poly[a] == 0: continue
                    pa = poly[a]
                    for b in range(n + 1 - a):
                        if base[b]: new[a + b] += pa * base[b]
                poly = new
        tot += mult * poly[n]
    return tot / G
a = [Fraction(0)] * (NMAX + 1); s = [Fraction(0)] * (NMAX + 1); a[1] = s[1] = Fraction(1)
def seqpow(s, n):
    f = [Fraction(0)] * (n + 1); f[0] = Fraction(1)
    for m in range(1, n + 1): f[m] = sum(s[j] * f[m - j] for j in range(1, m + 1))
    return f[n]
for n in range(2, NMAX + 1):
    s[n] = sum(cnt * polya(A, G, a, n) for (k, cnt, G, A, sc, anti) in types if k <= n)
    a[n] = seqpow(s, n)
print("strong S(n), n=1..%d:" % NMAX, [int(s[n]) for n in range(1, NMAX + 1)])
print("A000568(n), n=1..%d:" % NMAX, [int(a[n]) for n in range(1, NMAX + 1)])
print("known A000568: 1,1,2,4,12,56,456,6880,191536,9733056,903753248")
c = [Fraction(0)] * (NMAX + 1); cs = [Fraction(0)] * (NMAX + 1); c[1] = cs[1] = Fraction(1)
def euler_mset(cs, n):
    # Euler transform: prod_j (1 - x^j)^(-cs[j]); c(n) = (1/n) sum_{k=1}^{n} b(k) c(n-k), b(k) = sum_{d | k} d cs[d]
    b = [Fraction(0)] * (n + 1)
    for k in range(1, n + 1): b[k] = sum(Fraction(d) * cs[d] for d in range(1, k + 1) if k % d == 0)
    f = [Fraction(0)] * (n + 1); f[0] = Fraction(1)
    for m in range(1, n + 1): f[m] = sum(b[k] * f[m - k] for k in range(1, m + 1)) / m
    return f[n]
for n in range(2, NMAX + 1):
    tot = Fraction(0)
    for (k, cnt, G, A, sc, anti) in types:
        if k > n: continue
        if sc: tot += cnt * polya(A + anti, 2 * G, c, n)        # Aut+- has order 2G, one class per prime
        else: tot += Fraction(cnt, 2) * polya(A, G, c, n)      # converse pairs: Aut+- = Aut, two primes per class
    cs[n] = tot; c[n] = euler_mset(cs, n)
print("strong intransitivity classes c_s(n), n=1..%d:" % NMAX, [int(cs[n]) for n in range(1, NMAX + 1)])
print("intransitivity classes c(n), n=1..%d:" % NMAX, [int(c[n]) for n in range(1, NMAX + 1)])
print("known c (THM-4803, n<=10): 1,1,2,3,9,32,247,3328,93612,4774417")
# labelled counts: EGF.  l_s = x + (1/2) sum_k Plab(k) l^k / k!,  1 + l = exp(l_s)
Plab = {k: sum(Fraction(cnt * factorial(k), G) for (kk, cnt, G, A, sc, anti) in types if kk == k) for k in range(3, NMAX)}
print("labelled primes P_lab(k):", {k: int(v) for k, v in Plab.items()})
# work with EGF coefficient lists e[n] = l(n)/n!
L = [Fraction(0)] * (NMAX + 1); Ls = [Fraction(0)] * (NMAX + 1)
def mul(p, q, n):
    r = [Fraction(0)] * (n + 1)
    for i in range(n + 1):
        if p[i] == 0: continue
        for j in range(n + 1 - i): r[i + j] += p[i] * q[j]
    return r
def powser(p, k, n):
    r = [Fraction(0)] * (n + 1); r[0] = Fraction(1)
    for _ in range(k): r = mul(r, p, n)
    return r
def expser(p, n):  # exp of a series with p[0]=0
    r = [Fraction(0)] * (n + 1); r[0] = Fraction(1); term = r[:]
    for m in range(1, n + 1):
        term = mul(term, p, n); term = [t / m for t in term]
        r = [r[i] + term[i] for i in range(n + 1)]
    return r
for n in range(1, NMAX + 1):
    # iterate: Ls depends on L (orders < n only), then L[n] from exp
    Ls_n = Fraction(1) if n == 1 else Fraction(0)
    for k in range(3, n + 1):
        if k in Plab: Ls_n += Fraction(1, 2) * Plab[k] / factorial(k) * powser(L, k, n)[n]
    Ls[n] = Ls_n
    L[n] = expser(Ls, n)[n]
print("labelled intransitivity structures l(n), n=1..%d:" % NMAX, [int(L[n] * factorial(n)) for n in range(1, NMAX + 1)])
print("known l (THM-4803, n<=10): 1,1,2,11,208,8563,680426,101747483,28904318404,15824793243139")
