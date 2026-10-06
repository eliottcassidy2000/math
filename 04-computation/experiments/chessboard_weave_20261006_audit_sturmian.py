#!/usr/bin/env python3
"""Independent audit (written from scratch) of Theorem 4.1 / Corollary 4.2 of
05-knowledge/results/chessboard_weave_20261006.md.
 * re-derive sigma = iota o rho from rho(a) = c^-1 a^-1 c^-1, rho(c) = c^-1 a^-1 c^-1 a^-1 c^-1 (free group);
 * sigma = Inn(cac) q^3, sigma = E phi phi phit E, sigma = q' q' q;
 * cut(p,q) (rook-step cutting sequence of the centre-to-centre segment, a = horizontal, c = vertical), computed
   geometrically with exact integer comparisons, vs sigma^n(a), sigma^n(c) for n = 0..8;
 * the fixed point x vs the centred rounding word round((k+1)/phi) - round(k/phi) (exact via isqrt);
 * step (iii) rounding agreement; palindromes; balance (n <= 6); corner/tie claims of Cor 4.2;
 * the Sturmian morphisms with matrix Q^3 (enumerated from the generators E, phi, phit).
"""
from math import isqrt, gcd
from fractions import Fraction as Fr
from collections import deque

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


# ---------- free group ----------
def red(w):
    out = []
    for x in w:
        if out and out[-1][0] == x[0] and out[-1][1] == -x[1]:
            out.pop()
        else:
            out.append(x)
    return out


def inv(w):
    return [(l, -e) for (l, e) in reversed(w)]


def W(s):  # 'cAc' style: lowercase = generator, uppercase = inverse
    return [(ch.lower(), 1 if ch.islower() else -1) for ch in s]


def S(w):
    return ''.join(l if e == 1 else l.upper() for (l, e) in w)


def apply(m, w):
    out = []
    for (l, e) in w:
        img = m[l]
        out += img if e == 1 else inv(img)
    return red(out)


rho = {'a': W('CAC'), 'c': W('CACAC')}
iota_aut = {'a': W('A'), 'c': W('C')}  # letter inversion automorphism
sigma_from_rho = {x: apply(iota_aut, rho[x]) for x in 'ac'}
sigma_from_rho_anti = {x: inv(rho[x]) for x in 'ac'}  # iota as w -> w^-1
print("iota o rho (iota = letter inversion):", {x: S(sigma_from_rho[x]) for x in 'ac'})
print("iota o rho (iota = w -> w^-1):       ", {x: S(sigma_from_rho_anti[x]) for x in 'ac'})
sigma = {'a': W('cac'), 'c': W('cacac')}
check(all(sigma_from_rho[x] == sigma[x] for x in 'ac'), "sigma = iota rho")
q = {'a': W('c'), 'c': W('ac')}
q3 = {x: apply(q, apply(q, apply(q, W(x)))) for x in 'ac'}
g = W('cac')
inn = {x: red(g + q3[x] + inv(g)) for x in 'ac'}
print("q^3:", {x: S(q3[x]) for x in 'ac'}, " Inn(cac) q^3 (w -> g w g^-1):", {x: S(inn[x]) for x in 'ac'})
check(all(inn[x] == sigma[x] for x in 'ac'), "sigma = Inn(cac) q^3")
E = {'a': W('c'), 'c': W('a')}
phi = {'a': W('ac'), 'c': W('a')}
phit = {'a': W('ca'), 'c': W('a')}
comp = {x: apply(E, apply(phi, apply(phi, apply(phit, apply(E, W(x)))))) for x in 'ac'}
check(all(comp[x] == sigma[x] for x in 'ac'), "sigma = E phi phi phit E")
qp = {'a': W('c'), 'c': W('ca')}
comp2 = {x: apply(qp, apply(qp, apply(q, W(x)))) for x in 'ac'}
check(all(comp2[x] == sigma[x] for x in 'ac'), "sigma = q' q' q")
print("E phi phi phit E:", {x: S(comp[x]) for x in 'ac'}, " q'q'q:", {x: S(comp2[x]) for x in 'ac'})

# ---------- string substitution ----------
F = [0, 1]
for _ in range(40):
    F.append(F[-1] + F[-2])


def Fm(n):  # allow F_{-1} = 1
    return 1 if n == -1 else F[n]


def sig(w):
    return ''.join('cac' if ch == 'a' else 'cacac' for ch in w)


def cut(p, q, tie=None):
    """rook-step cutting sequence of segment (1/2,1/2)->(p+1/2,q+1/2); a = crossing a vertical line (horizontal step),
    c = crossing a horizontal line (vertical step). Returns (word, number of corner crossings)."""
    out = []
    i, j = 1, 1
    corners = 0
    while i <= p or j <= q:
        if i > p:
            out.append('c'); j += 1; continue
        if j > q:
            out.append('a'); i += 1; continue
        si, sj = (2 * i - 1) * q, (2 * j - 1) * p  # compare (2i-1)/(2p) vs (2j-1)/(2q)
        if si < sj:
            out.append('a'); i += 1
        elif sj < si:
            out.append('c'); j += 1
        else:
            corners += 1
            out.append(tie if tie else 'ac')
            i += 1; j += 1
    return ''.join(out), corners


A, C = 'a', 'c'
x_long = None
for n in range(0, 9):
    ca, ka = cut(Fm(3 * n - 1), F[3 * n])
    cc, kc = cut(F[3 * n], F[3 * n + 1])
    check(A == ca and kc == 0 and ka == 0, f"sigma^{n}(a) = cut(F_{3*n-1},F_{3*n})")
    check(C == cc, f"sigma^{n}(c) = cut(F_{3*n},F_{3*n+1})")
    check(A == A[::-1] and C == C[::-1], f"palindromes n={n}")
    check((A.count('a'), A.count('c')) == (Fm(3 * n - 1), F[3 * n]) and (C.count('a'), C.count('c')) == (F[3 * n], F[3 * n + 1]), f"letter counts n={n}")
    if n >= 1:
        check(C.startswith(A), f"sigma^{n}(a) prefix of sigma^{n}(c)")
    print(f"n={n}: |sigma^n(a)|={len(A)}, |sigma^n(c)|={len(C)}, cut identities OK={A == ca and C == cc}, corners {ka},{kc}")
    x_long = C
    A, C = sig(A), sig(C)


# ---------- fixed point vs centred rounding word ----------
def R(k):  # round(k/phi) = floor(k/phi + 1/2), exact
    if k == 0:
        return 0
    X = isqrt(5 * k * k)  # floor(k sqrt5), k sqrt5 irrational
    return (X - k + 1) // 2


Lx = len(x_long)
rw = ''.join('c' if R(k + 1) - R(k) == 1 else 'a' for k in range(Lx))
check(rw == x_long, "prefix of x = centred rounding word")
print(f"x prefix (length {Lx}) equals the centred rounding word round((k+1)/phi)-round(k/phi):", rw == x_long)
# Fibonacci word (corner-start, intercept 0, c_alpha with alpha=1/phi) differs:
def Fl(k):  # floor(k/phi), exact
    return (isqrt(5 * k * k) - k) // 2


fw = ''.join('c' if Fl(k + 2) - Fl(k + 1) == 1 else 'a' for k in range(40))
print("first 40 letters of x:", x_long[:40])
print("characteristic word c_{1/phi} (first 40):", fw)

# ---------- step (iii): rounding agreement ----------
for m in range(3, 27):
    Fm_ = F[m]
    dis = [k for k in range(0, Fm_ + 1) if R(k) != ((2 * k * F[m - 1] + Fm_) // (2 * Fm_))]
    if Fm_ % 2 == 1:
        check(dis == [], f"rounding agreement odd F_{m}")
    else:
        # rational value is a half-integer at k = F_m/2: report
        if m <= 12:
            print(f"   even F_{m}={Fm_}: disagreements (round-half-up convention) at k = {dis}")
print("step (iii): round(k/phi) = round(k F_{m-1}/F_m) for 0<=k<=F_m whenever F_m odd, m<=26: OK" if not any('rounding agreement' in f for f in fails) else "step (iii) FAILED")

# ---------- step (iv): cut(p,q) = rational rounding word, p+q odd ----------
for p in range(0, 40):
    for qq in range(0, 40):
        if p + qq == 0 or gcd(p, qq) != 1 or (p + qq) % 2 == 0:
            continue
        s = p + qq
        w = ''.join('c' if (2 * (k + 1) * qq + s) // (2 * s) - (2 * k * qq + s) // (2 * s) == 1 else 'a' for k in range(s))
        check(w == cut(p, qq)[0], f"cut = rounding word ({p},{qq})")
print("step (iv): cut(p,q) = rounding word with slope q/(p+q) for coprime p+q odd, p,q<40: checked")

# ---------- corners iff p,q odd ----------
for p in range(1, 60):
    for qq in range(1, 60):
        if gcd(p, qq) != 1:
            continue
        k = cut(p, qq)[1]
        check((k > 0) == (p % 2 == 1 and qq % 2 == 1) and k <= 1, f"corner rule ({p},{qq})")
print("corner iff p,q both odd (single corner, at the midpoint): checked coprime p,q<60")

# ---------- tie resolution at even prefixes F_{3k} ----------
ties = []
for k in range(1, 9):
    Lk = F[3 * k]
    pref = x_long[:Lk]
    w_ac, _ = cut(F[3 * k - 2], F[3 * k - 1], tie='ac')
    w_ca, _ = cut(F[3 * k - 2], F[3 * k - 1], tie='ca')
    res = 'ac' if pref == w_ac else ('ca' if pref == w_ca else 'neither')
    ties.append(res)
print("tie resolution of x at prefixes of length F_{3k}, k=1..8:", ties)
check(ties == ['ca', 'ac'] * 4, "ties alternate ca, ac, ...")

# ---------- balance (brute force) for n <= 6 ----------
try:
    import numpy as np
    A6 = 'a'
    C6 = 'c'
    for n in range(6):
        A6, C6 = sig(A6), sig(C6)
    for name, wd in (("sigma^6(a)", A6), ("sigma^6(c)", C6)):
        arr = np.frombuffer(wd.encode(), dtype=np.uint8) == ord('c')
        cs = np.concatenate([[0], np.cumsum(arr)])
        ok = True
        for l in range(1, len(wd) + 1):
            win = cs[l:] - cs[:-l]
            if win.max() - win.min() > 1:
                ok = False
                break
        check(ok, f"balanced {name}")
        print(f"{name} (length {len(wd)}) balanced:", ok)
except ImportError:
    pass

# ---------- Sturmian morphisms with matrix Q^3 ----------
gens = {'E': {'a': 'c', 'c': 'a'}, 'phi': {'a': 'ac', 'c': 'a'}, 'phit': {'a': 'ca', 'c': 'a'}}


def app_s(m, w):
    return ''.join(m[ch] for ch in w)


start = ('a', 'c')
seen = {start}
dq = deque([start])
while dq:
    fa, fc = dq.popleft()
    for gname, gm in gens.items():
        new = (app_s(gm, fa), app_s(gm, fc))
        if len(new[0]) + len(new[1]) <= 8 and new not in seen:
            seen.add(new)
            dq.append(new)
target = [f for f in seen if (f[0].count('a'), f[0].count('c'), f[1].count('a'), f[1].count('c')) == (1, 2, 2, 3)]
print("Sturmian morphisms with matrix Q^3 (a->1a2c, c->2a3c):", len(target), sorted(target))
check(len(target) == 7, "seven Sturmian morphisms")
pal = [f for f in target if f[0] == f[0][::-1] and f[1] == f[1][::-1]]
print("with both images palindromes:", pal)
check(pal == [('cac', 'cacac')], "sigma unique palindromic")
# conjugation chain: f -> (u^-1 f(a) u, u^-1 f(c) u) when both images start with u
nxt = {}
for f in target:
    if f[0][0] == f[1][0]:
        h = (f[0][1:] + f[0][0], f[1][1:] + f[1][0])
        nxt[f] = h
firsts = [f for f in target if f not in nxt.values()]
chain = []
if len(firsts) == 1:
    f = firsts[0]
    chain.append(f)
    while f in nxt:
        f = nxt[f]
        chain.append(f)
print("conjugation chain:", chain)
q3s = ('cac', 'accac')
qp3s = (app_s({'a': 'c', 'c': 'ca'}, app_s({'a': 'c', 'c': 'ca'}, 'c')), app_s({'a': 'c', 'c': 'ca'}, app_s({'a': 'c', 'c': 'ca'}, 'ca')))
print("q'^3 =", qp3s, " q^3 =", q3s, " sigma position in chain:", chain.index(('cac', 'cacac')) if ('cac', 'cacac') in chain else None)
check(len(chain) == 7 and chain[3] == ('cac', 'cacac') and chain[0] == qp3s and chain[-1] == q3s, "sigma middle of chain q'^3 .. q^3")

# ---------- parity classes, delta along the Fibonacci walk ----------
cls = [(1, 0)]
for _ in range(3):
    u, v = cls[-1]
    cls.append((v % 2, (u + v) % 2))
print("Q on nonzero classes mod 2:", cls)
check(cls == [(1, 0), (0, 1), (1, 1), (1, 0)], "Q cycles parity classes")


def nd(x):
    f = x - (x.numerator // x.denominator)
    return min(f, 1 - f)


def delta2(a, b):
    if a == 0 or b == 0:
        return Fr(0)
    cands = set()
    for v in (a, b):
        for k in range(2 * v + 1):
            cands.add(Fr(k, 2 * v))
    for d in (a + b, abs(a - b)):
        if d:
            for k in range(d + 1):
                cands.add(Fr(k, d))
    return max(min(nd(t * a), nd(t * b)) for t in cands)


for n in range(0, 31):
    a, b = F[n], F[n + 1]
    claim = Fr(1, 2) if n % 3 == 1 else Fr(1, 2) - Fr(1, 2 * F[n + 2])
    if n <= 13:
        check(delta2(a, b) == claim, f"delta(F_{n},F_{n+1}) exact")
    else:
        s = a + b
        check(Fr(s // 2, s) == claim, f"delta(F_{n},F_{n+1}) via Thm 2.1")
print("delta(F_n,F_{n+1}) claim checked n<=30 (exact candidate evaluation for n<=13):",
      [str(delta2(F[n], F[n + 1])) for n in range(0, 6)])

# ---------- step (ii) general claim: do unbounded palindromic prefixes force intercept 0 or 1/2? ----------
# lower mechanical word s_{alpha,rho}(n) = floor((n+1)alpha+rho) - floor(n alpha+rho), alpha = 1/phi (letter c = 1).
# rho = j*alpha (j >= 1) gives the characteristic word c_alpha (j = 1) and its shifts T^{j-1} c_alpha.
def floor_k_over_phi_plus_j(k):  # floor(k/phi) exact
    return (isqrt(5 * k * k) - k) // 2


def shifted_char(jshift, length):
    # T^jshift c_alpha, c_alpha(n) = floor((n+2)/phi) - floor((n+1)/phi)
    return ''.join('c' if floor_k_over_phi_plus_j(n + 2 + jshift) - floor_k_over_phi_plus_j(n + 1 + jshift) == 1 else 'a'
                   for n in range(length))


for jshift in (0, 1, 2, 3):
    wd = shifted_char(jshift, 20000)
    pals = [L for L in range(2, 20001) if wd[:L] == wd[:L][::-1]]
    print(f"T^{jshift}(c_alpha) [intercept {jshift+1}*alpha mod 1]: palindromic prefix lengths (>=2, <=20000):", pals[-6:], "count", len(pals))
    if jshift >= 1:
        check(len(pals) >= 5 and pals[-1] > 1000, f"shifted characteristic word has long palindromic prefixes j={jshift}")

print()
print("TOTAL FAILURES:", len(fails))
