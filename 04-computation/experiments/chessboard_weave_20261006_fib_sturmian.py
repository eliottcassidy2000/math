#!/usr/bin/env python3
"""Sturmian reading of the monotile boundary substitution (chessboard-weave, 2026-10-06).

Source: PingYou Ltd 2026, "A chiral aperiodic polygon with Fibonacci monodromy",
Thm 5.2: on F(a,c) the boundary substitution is sigma = iota*rho with
    sigma(a) = cac,  sigma(c) = cacac,  sigma = Inn(cac) o q^3,  q(a) = c, q(c) = ac.
sigma is positive, so it is a monoid morphism of {a,c}*.

Conventions.  a = one horizontal rook step (the path crosses a vertical grid line),
c = one vertical rook step (crosses a horizontal grid line).
cut(p,q) = cutting sequence (sequence of rook steps) of the straight segment joining
the CENTRES of two squares with offset (p,q) = p columns right, q rows up.
Lothaire generators (letter b renamed c): E: a<->c,  phi: a->ac, c->a,  phit: a->ca, c->a.

Sections
  A. sigma = E phi phi phit E; the 7 Sturmian morphisms with matrix Q^3; sigma is the
     unique one with palindromic images; sigma = q' q' q with q' = mirror of q.
  B. sigma^n(a), sigma^n(c), n <= 8: Parikh vectors vs columns of Q^{3n}, palindromes,
     balance (brute force), Christoffel class + palindromic conjugate, central tests,
     equality with cut(.,.) computed geometrically.
  C. the fixed point x = sigma^oo(c) equals the rounding word s_{1/phi,1/2}
     (cutting sequence of the golden ray from the centre of a square); its odd
     Fibonacci prefixes are the centre-to-centre flights of the Fibonacci leapers.
Run:  python3 chessboard_weave_20261006_fib_sturmian.py
"""
from math import gcd, isqrt
from itertools import product
import time
import numpy as np


def fib(k):
    if k < 0:
        return (-1) ** (k + 1) * fib(-k)
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def app(m, w):
    return ''.join(m[x] for x in w)


def comp(*fs):
    """comp(f1, f2, ..., fk) = f1 o f2 o ... o fk."""
    res = {'a': 'a', 'c': 'c'}
    for f in reversed(fs):
        res = {x: app(f, res[x]) for x in 'ac'}
    return res


def parikh(w):
    return (w.count('a'), w.count('c'))


def matrix(f):
    # columns = Parikh vectors of the images; rows = (#a, #c)
    pa, pc = parikh(f['a']), parikh(f['c'])
    return ((pa[0], pc[0]), (pa[1], pc[1]))


def is_pal(w):
    return w == w[::-1]


def is_balanced_fast(w):
    """Brute force: every pair of equal-length factors differ by <= 1 in #c."""
    x = (np.frombuffer(w.encode(), dtype=np.uint8) == ord('c')).astype(np.int32)
    P = np.concatenate(([0], np.cumsum(x, dtype=np.int32)))
    L = len(w)
    buf = np.empty(L + 1, dtype=np.int32)
    for l in range(1, L):
        d = np.subtract(P[l:], P[:L + 1 - l], out=buf[:L + 1 - l])
        if int(d.max()) - int(d.min()) > 1:
            return False, l
    return True, None


def christoffel_lower(p, q):
    """lower Christoffel word with p letters a and q letters c (gcd(p,q)=1)."""
    n = p + q
    return ''.join('c' if ((k + 1) * q) // n - (k * q) // n else 'a' for k in range(n))


def rounding_word(p, q):
    """centred digitisation: letter k is c iff floor((k+1)q/n+1/2) - floor(kq/n+1/2) = 1."""
    n = p + q
    r = lambda k: (2 * k * q + n) // (2 * n)
    return ''.join('c' if r(k + 1) - r(k) else 'a' for k in range(n))


def cut_centres(p, q):
    """Geometric cutting sequence of the segment (1/2,1/2) -> (1/2+p, 1/2+q).
    Crossing of x=i at relative time (2i-1)/(2p), of y=j at (2j-1)/(2q); integer keys
    after multiplying by 2pq.  Returns (word, number of simultaneous crossings)."""
    if p == 0:
        return 'c' * q, 0
    if q == 0:
        return 'a' * p, 0
    ev = [((2 * i - 1) * q, 0, 'a') for i in range(1, p + 1)] + \
         [((2 * j - 1) * p, 1, 'c') for j in range(1, q + 1)]
    ev.sort()
    keys = [e[0] for e in ev]
    ties = len(keys) - len(set(keys))
    return ''.join(e[2] for e in ev), ties


def r_golden(k):
    """floor(k/phi + 1/2) exactly: = floor((k*sqrt5 - k + 1)/2), floor(k sqrt5) = isqrt(5k^2)."""
    return (isqrt(5 * k * k) - k + 1) // 2


def x_prefix(N):
    r = [r_golden(k) for k in range(N + 1)]
    return ''.join('c' if r[k + 1] - r[k] else 'a' for k in range(N))


def periods(w):
    n = len(w)
    pi = [0] * n
    for i in range(1, n):
        k = pi[i - 1]
        while k and w[i] != w[k]:
            k = pi[k - 1]
        if w[i] == w[k]:
            k += 1
        pi[i] = k
    per = {n}
    b = pi[-1] if n else 0
    while b:
        per.add(n - b)
        b = pi[b - 1]
    return per


def is_central(w):
    """central = palindrome with coprime periods p,q and |w| = p+q-2 (incl. unary words)."""
    n = len(w)
    if n == 0 or len(set(w)) == 1:
        return True
    if not is_pal(w):
        return False
    per = periods(w)
    for p in per:
        q = n + 2 - p
        if q >= 1 and gcd(p, q) == 1 and (q in per or q > n):
            return True
    return False


def manacher_odd(s):
    n = len(s)
    d1 = [0] * n
    l, r = 0, -1
    for i in range(n):
        k = 1 if i > r else min(d1[l + r - i], r - i + 1)
        while i - k >= 0 and i + k < n and s[i - k] == s[i + k]:
            k += 1
        d1[i] = k
        if i + k - 1 > r:
            l, r = i - k + 1, i + k - 1
    return d1


def palindromic_rotations(w):
    """offsets k such that rotation w[k:]+w[:k] is a palindrome (|w| odd)."""
    L = len(w)
    assert L % 2 == 1
    d1 = manacher_odd(w + w)
    h = (L - 1) // 2
    return [k for k in range(L) if d1[k + h] >= h + 1]


t0 = time.time()
sigma = {'a': 'cac', 'c': 'cacac'}
E = {'a': 'c', 'c': 'a'}
phi = {'a': 'ac', 'c': 'a'}
phit = {'a': 'ca', 'c': 'a'}
q = {'a': 'c', 'c': 'ac'}          # the paper's Fibonacci automorphism
qm = {'a': 'c', 'c': 'ca'}         # its mirror: qm(x) = reverse(q(x))

print("== A. sigma inside the Sturmian monoid <E, phi, phit> ==")
dec = comp(E, phi, phi, phit, E)
print("E o phi o phi o phit o E =", dec, " equals sigma:", dec == sigma)
print("q = E o phit o E:", comp(E, phit, E) == q, "; q' := E o phi o E =", comp(E, phi, E),
      "= mirror of q:", all(comp(E, phi, E)[x] == q[x][::-1] for x in 'ac'))
print("sigma = q' o q' o q:", comp(qm, qm, q) == sigma)
q3 = comp(q, q, q)
print("q^3 =", q3, "; sigma(x).cac == cac.q^3(x) for x in {a,c}:",
      all(sigma[x] + 'cac' == 'cac' + q3[x] for x in 'ac'))
print("matrix(sigma) =", matrix(sigma), "(rows #a,#c; columns images of a,c) = Q^3 = [[F2,F3],[F3,F4]]")

# all morphisms in <E,phi,phit> of word length <= 7 with matrix Q^3
gens = {'E': E, 'phi': phi, 'phit': phit}
found = {}
for L in range(1, 8):
    for word in product(gens, repeat=L):
        f = comp(*[gens[g] for g in word])
        if matrix(f) == ((1, 2), (2, 3)):
            key = (f['a'], f['c'])
            if key not in found or len(word) < len(found[key]):
                found[key] = word
print("distinct morphisms in <E,phi,phit> (words of length <= 7) with matrix Q^3:", len(found),
      " (Seebold: |f(a)|+|f(c)|-1 = 7)")
# order them as conjugates f_0..f_6 starting from the standard one (different last letters)
f0 = [k for k in found if k[0][-1] != k[1][-1]]
assert len(f0) == 1
chain = [f0[0]]
while chain[-1][0][0] == chain[-1][1][0]:
    fa, fc = chain[-1]
    z = fa[0]
    chain.append((fa[1:] + z, fc[1:] + z))
print("conjugate chain f_0..f_%d (f_{i+1}(x) = f_i(x) with first letter moved to the end):" % (len(chain) - 1))
for i, (fa, fc) in enumerate(chain):
    tags = []
    if (fa, fc) == (sigma['a'], sigma['c']):
        tags.append('= sigma')
    if (fa, fc) == (q3['a'], q3['c']):
        tags.append('= q^3')
    if is_pal(fa) and is_pal(fc):
        tags.append('both images palindromes')
    print("   f_%d: a->%-4s c->%-6s in St: %s  %s" % (i, fa, fc, (fa, fc) in found, ' '.join(tags)))
print("set of chain == set found:", set(chain) == set(found))
print("the 8 words of length 3 in {q, q'} and the conjugate they give:")
for word in product(['q', "q'"], repeat=3):
    f = comp(*[q if g == 'q' else qm for g in word])
    print("   %-12s -> f_%d" % (' o '.join(word), chain.index((f['a'], f['c']))))
# sanity: each of the 7 maps a long Sturmian prefix to a balanced word
xs = x_prefix(3000)
print("each f_i maps the length-3000 prefix of a Sturmian word to a balanced word:",
      all(is_balanced_fast(app({'a': fa, 'c': fc}, xs))[0] for fa, fc in chain))

print()
print("== sanity of the central-word test: #central words of length n == totient(n+2) ==")
def totient(m):
    return sum(1 for k in range(1, m + 1) if gcd(k, m) == 1)
ok = True
for n in range(0, 15):
    cnt = sum(1 for t in product('ac', repeat=n) if is_central(''.join(t)))
    ok &= (cnt == totient(n + 2))
print("n = 0..14 all match:", ok)

print()
print("== B. the boundary words sigma^n(a), sigma^n(c), n = 0..8 ==")
W = {'a': 'a', 'c': 'c'}
rows = []
for n in range(0, 9):
    if n > 0:
        W = {x: app(sigma, W[x]) for x in 'ac'}
    for x in 'ac':
        w = W[x]
        L = len(w)
        pa, pc = parikh(w)
        if x == 'a':
            exp = (fib(3 * n - 1), fib(3 * n))
        else:
            exp = (fib(3 * n), fib(3 * n + 1))
        pal = is_pal(w)
        centre = w[(L - 1) // 2] if L % 2 else None
        tb = time.time()
        bal, bad = is_balanced_fast(w) if L > 1 else (True, None)
        tbal = time.time() - tb
        C = christoffel_lower(pa, pc)
        P = C[1:-1]
        off = (C + C).find(w)
        pred_off = None
        if L % 2 == 1 and pc > 0 and pa > 0:
            pred_off = ((L - 1) // 2) * pow(pc, -1, L) % L
        prots = palindromic_rotations(C) if L % 2 == 1 else None
        cw, ties = cut_centres(pa, pc)
        rw = rounding_word(pa, pc) if L % 2 == 1 else None
        xp = x_prefix(L)
        print("sigma^%d(%s): len %6d = F_%d; Parikh (#a,#c) = (%d,%d) = Q^{3n} column %s: %s"
              % (n, x, L, 3 * n + (1 if x == 'a' else 2), pa, pc, '1' if x == 'a' else '2', (pa, pc) == exp))
        print("     palindrome %s, centre letter %s; balanced (brute force, %.1fs) %s; gcd(#a,#c) = %d"
              % (pal, centre, tbal, bal, gcd(pa, pc)))
        print("     Christoffel a.P.c of slope #c/#a = %d/%d: P palindrome %s, P central %s;"
              " w = rotation of it by %d (predicted %s); palindromic rotations of the class: %s"
              % (pc, pa, is_pal(P), is_central(P), off, pred_off, prots))
        print("     w central: %s; w[1:-1] central: %s; w == rounding word: %s; w == cut(%d,%d): %s"
              " (corner ties %d); w == prefix of x: %s"
              % (is_central(w), is_central(w[1:-1]) if L > 2 else None, w == rw, pa, pc, w == cw, ties, w == xp))
        if L <= 34:
            print("     w =", w)
        rows.append((n, x, L, (pa, pc) == exp, pal, bal, off >= 0, w == cw, ties, w == xp))

print("summary: all Parikh == Q^{3n} columns: %s; all palindromes: %s; all balanced: %s; all in Christoffel class: %s;"
      " all == centre-to-centre cut with no corner: %s; all prefixes of x: %s"
      % (all(r[3] for r in rows), all(r[4] for r in rows), all(r[5] for r in rows), all(r[6] for r in rows),
         all(r[7] and r[8] == 0 for r in rows), all(r[9] for r in rows if r[2] > 1 or r[1] == 'c')))

print()
print("== C. the fixed point x = lim sigma^n(c) ==")
X = W['c']   # sigma^8(c), length F_26
print("sigma^8(c) (len %d) == prefix of s_{1/phi,1/2}, x(k) = round((k+1)/phi) - round(k/phi): %s"
      % (len(X), X == x_prefix(len(X))))
print("sigma(x) == x on the first %d letters: %s" % (len(X), app(sigma, X)[:len(X)] == X))
fa = X.count('a') / len(X)
print("letter frequencies in sigma^8(c): a %.9f (1/phi^2 = %.9f), c %.9f (1/phi = %.9f)"
      % (fa, ((5 ** .5 - 1) / 2) ** 2, 1 - fa, (5 ** .5 - 1) / 2))
# the Fibonacci (characteristic) word of the same slope, for contrast: c_alpha(k) = floor((k+2)alpha) - floor((k+1)alpha)
def fl_alpha(k):   # floor(k/phi) = floor((k sqrt5 - k)/2)
    return (isqrt(5 * k * k) - k) // 2
cfib = ''.join('c' if fl_alpha(k + 2) - fl_alpha(k + 1) else 'a' for k in range(60))
print("x       =", X[:60])
print("c_alpha =", cfib, "(Fibonacci/characteristic word, same slope: different word)")
def manacher_even(s):
    """d2[i] = number of even palindromes centred between s[i-1] and s[i]."""
    n = len(s)
    d2 = [0] * n
    l, r = 0, -1
    for i in range(n):
        k = 0 if i > r else min(d2[l + r - i + 1], r - i + 1)
        while i - k - 1 >= 0 and i + k < n and s[i - k - 1] == s[i + k]:
            k += 1
        d2[i] = k
        if i + k - 1 > r:
            l, r = i - k, i + k - 1
    return d2
N = len(X)
d1, d2 = manacher_odd(X), manacher_even(X)
pref_pal = [l for l in range(1, N + 1)
            if (l % 2 == 1 and d1[(l - 1) // 2] >= (l + 1) // 2) or (l % 2 == 0 and d2[l // 2] >= l // 2)]
print("lengths l <= %d with x[0:l] a palindrome (Manacher, odd and even):" % N, pref_pal)
oddfib = [fib(m) for m in range(2, 27) if m % 3 != 0 and fib(m) <= N]
print("odd Fibonacci numbers F_m (m>=2, 3 !| m) <= %d:" % N, sorted(set(oddfib)), " equal:", sorted(set(oddfib)) == pref_pal)
print("brute-force cross-check of the Manacher list for l <= 3000:",
      [l for l in range(1, 3001) if is_pal(X[:l])] == [l for l in pref_pal if l <= 3000])
print("prefixes x[0:F_m] vs centre-to-centre flights cut(F_{m-2},F_{m-1}), m = 2..26:")
for m in range(2, 27):
    L = fib(m)
    p_, q_ = fib(m - 2), fib(m - 1)
    cw, ties = cut_centres(p_, q_)
    pre = X[:L]
    if ties == 0:
        st = "equal: %s" % (pre == cw)
    else:
        # the flight meets exactly one corner (its midpoint); the two resolutions differ by ac <-> ca there
        mid = (L - 1) // 2   # crossings 'mid' and 'mid+1' (0-based) are simultaneous
        alt = cw[:mid] + cw[mid + 1] + cw[mid] + cw[mid + 2:]
        which = 'cw' if pre == cw else ('alt' if pre == alt else 'neither')
        st = "flight meets a corner (%d tie); x-prefix = resolution '%s' at the corner: %s" % (
            ties, pre[mid:mid + 2], which != 'neither')
    print("   m=%2d len F_m=%7d leaper (%6d,%6d) %s: %s" % (m, L, p_, q_, 'colour-preserving' if (p_ + q_) % 2 == 0 else 'colour-changing  ', st))
print("bi-infinite fixed point rev(x).x: unique seed c.c (sigma(y) ends with c, sigma(z) starts with c); 'cc' occurs in x:",
      'cc' in X)
print("elapsed %.1fs" % (time.time() - t0))
