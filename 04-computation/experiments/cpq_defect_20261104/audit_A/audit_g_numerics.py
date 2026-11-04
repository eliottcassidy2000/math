#!/usr/bin/env python3
"""audit_A (g): the note's section 5 (folds, mirrors, the ball, -7/4) recomputed independently.
mpmath: MRB constant C, B = sum (-1)^n (r_n - 1) with r_n = exp(-W_0(-ln n/n)) (r_1 = 1, r_2 = 4), B = 3 - length(E),
C = sqrt2 - 1 - int_E f', B - C, m'(1.01), ball volumes, x*, mirror dimensions, the alternating ball sum against
e^pi erfc(sqrt pi), monotonicity of V(n)/V(n-1); sympy: the (f^3 - x)/(f - x) factorisations at c = -7/4, -29/16
and the multipliers; gp: the cubic fields; nauty: the twelve 5-tournaments and the self-converse counts n <= 7.
"""
import subprocess, itertools
from mpmath import mp, mpf, lambertw, log, exp, sqrt, pi, gamma, erfc, nsum, inf, findroot, quad, diff
import sympy as sp
import pynauty

mp.dps = 40

def r(n):
    n = mpf(n)
    if n == 1: return mpf(1)
    if n == 2: return mpf(4)
    return exp(-lambertw(-log(n) / n, 0))

def f(x): return x ** (1 / x)

print("== mirror bases ==")
for n in (2, 3, 4, 5, 10):
    print("  r_%d = %s   f(r_n) - f(n) = %s" % (n, mp.nstr(r(n), 20), mp.nstr(f(r(n)) - f(n), 5)))
# Cayley expansion check at n = 3
z = log(3) / 3
cay = sum(mpf(m + 1) ** (m - 1) / mp.factorial(m) * z ** m for m in range(0, 61))
print("  Cayley sum (60 terms) at n = 3: %s vs r_3 = %s, gap %s" % (mp.nstr(cay, 10), mp.nstr(r(3), 10), mp.nstr(r(3) - cay, 4)))
# m'(x) at x = 1.01
x0 = mpf("1.01")
mval = exp(-lambertw(-log(x0) / x0, -1))   # mirror of x in (1, e) lies above e: branch -1
fp = lambda u: f(u) * (1 - log(u)) / u ** 2
print("  m(1.01) = %s, m'(1.01) = f'(x)/f'(m) = %s, asymptotic -m^2/(ln m - 1) = %s" % (mp.nstr(mval, 8), mp.nstr(fp(x0) / fp(mval), 8), mp.nstr(-mval ** 2 / (log(mval) - 1), 8)))

print("\n== the two constants ==")
C = nsum(lambda n: (-1) ** n * (n ** (1 / n) - 1), [1, inf], method="alternating")
print("  MRB C = %s" % mp.nstr(C, 32))
print("  claimed 0.18785964246206712024851793405...")
B = nsum(lambda n: (-1) ** n * (r(n) - 1), [1, inf], method="alternating")
print("  B     = %s" % mp.nstr(B, 32))
print("  claimed 2.09679112439152612572503067508...")
print("  B - C = %s (claimed 1.90893148192945900547651274103...)" % mp.nstr(B - C, 32))
BC = 4 - sqrt(2) + nsum(lambda n: (-1) ** n * (r(n) - n ** (1 / n)), [3, inf], method="alternating")
print("  4 - sqrt2 + sum_{n>=3} (-1)^n (r_n - n^(1/n)) = %s" % mp.nstr(BC, 32))
# B = 3 - length(E) with E = union [r_{2k+2}, r_{2k+1}], k >= 1, and C = sqrt2 - 1 - int_E f': both are the pairing
# of consecutive terms of a convergent alternating series (an exact rearrangement: grouping consecutive terms of a
# convergent series does not change its sum); numerically: paired partial sums to K and the tail bound (next term).
K = 4000
lenE_K = sum(r(2 * k + 1) - r(2 * k + 2) for k in range(1, K + 1))
intE_K = sum(f(mpf(2 * k + 1)) - f(mpf(2 * k + 2)) for k in range(1, K + 1))   # = int over [r_{2k+2}, r_{2k+1}] of f'
altB_K = sum((-1) ** n * (r(n) - 1) for n in range(1, 2 * K + 3))
altC_K = sum((-1) ** n * (mpf(n) ** (1 / mpf(n)) - 1) for n in range(1, 2 * K + 3))
print("  paired partial sums K=%d: 3 - length_K(E) = %s = alternating partial sum of B to 2K+2 = %s (tail < r_(2K+3) - 1 = %s)" % (
    K, mp.nstr(3 - lenE_K, 12), mp.nstr(altB_K, 12), mp.nstr(r(2 * K + 3) - 1, 5)))
print("  paired partial sums K=%d: sqrt2 - 1 - int_K = %s = alternating partial sum of C to 2K+2 = %s (tail < %s); |B - (3 - length_K)| = %s, |C - (...)| = %s" % (
    K, mp.nstr(sqrt(2) - 1 - intE_K, 12), mp.nstr(altC_K, 12), mp.nstr(mpf(2 * K + 3) ** (1 / mpf(2 * K + 3)) - 1, 5),
    mp.nstr(abs(B - (3 - lenE_K)), 5), mp.nstr(abs(C - (sqrt(2) - 1 - intE_K)), 5)))
print("  check one interval: int_[r_4, r_3] f' = %s, f(3) - f(4) = %s" % (mp.nstr(quad(fp, [r(4), r(3)]), 15), mp.nstr(f(mpf(3)) - f(mpf(4)), 15)))
for n in (10 ** 3, 10 ** 4, 10 ** 5):
    zz = log(n) / mpf(n)
    print("  (r_n - n^(1/n)) / (ln n/n)^2 at n = %d: %s" % (n, mp.nstr((r(n) - mpf(n) ** (1 / mpf(n))) / zz ** 2, 8)))

print("\n== unit ball ==")
V = lambda x: pi ** (x / 2) / gamma(x / 2 + 1)
print("  V(n), n = 0..8:", [mp.nstr(V(mpf(n)), 6) for n in range(9)])
xs = findroot(lambda x: diff(V, x), 5.25)
print("  x* = %s, V(x*) = %s" % (mp.nstr(xs, 8), mp.nstr(V(xs), 8)))
for n in range(6):
    m = findroot(lambda x: V(x) - V(mpf(n)), 13 - 1.3 * n)
    print("  mirror of %d: %s (V = %s vs %s)" % (n, mp.nstr(m, 8), mp.nstr(V(m), 8), mp.nstr(V(mpf(n)), 8)))
S = nsum(lambda n: (-1) ** n * V(n), [0, inf])
print("  sum (-1)^n V(n) = %s" % mp.nstr(S, 30))
print("  e^pi erfc(sqrt pi) = %s (claimed 0.28205917617568264518411963...)" % mp.nstr(exp(pi) * erfc(sqrt(pi)), 30))
rat = [sqrt(pi) * gamma((mpf(n) + 1) / 2) / gamma(mpf(n) / 2 + 1) for n in range(1, 15)]
print("  V(n)/V(n-1), n=1..14:", [mp.nstr(x, 6) for x in rat], " strictly decreasing:", all(rat[i] > rat[i + 1] for i in range(len(rat) - 1)))
print("  n=5: 16/15 = %s, n=6: 5pi/16 = %s" % (mp.nstr(rat[4], 10), mp.nstr(rat[5], 10)))

print("\n== -7/4 and -29/16 (sympy) ==")
x = sp.symbols('x')
for cval in (sp.Rational(-7, 4), sp.Rational(-29, 16)):
    fx = x ** 2 + cval
    f3 = fx.subs(x, fx.subs(x, fx))
    quo = sp.cancel((f3 - x) / (fx - x))
    print("  c = %s: (f^3 - x)/(f - x) = %s" % (cval, sp.factor(quo)))
    for fac, mult in sp.factor_list(quo)[1]:
        if sp.degree(fac, x) == 3:
            roots = sp.Poly(fac, x).all_roots()
            multiplier = sp.simplify(sp.prod([2 * rt for rt in roots]))
            print("     cubic factor %s (mult %d): roots %s, multiplier prod 2x_i = %s" % (fac, mult, [sp.N(rt, 8) for rt in roots], multiplier))
print("  c(c+1)^2 at c=-7/4:", sp.Rational(-7, 4) * sp.Rational(-3, 4) ** 2, " f^3(0) =", (sp.Rational(21, 16)) ** 2 - sp.Rational(7, 4))
cyc = [sp.Rational(-7, 4)]
for _ in range(3): cyc.append(cyc[-1] ** 2 - sp.Rational(29, 16))
print("  x^2 - 29/16 orbit of -7/4:", cyc)

print("\n== cubic fields (gp) ==")
gpin = r"""
default(parisize, "64M");
P1 = 8*x^3 + 4*x^2 - 18*x - 1;
print("disc(P1) = ", poldisc(P1), "  field disc = ", nfdisc(P1), "  isisom with x^3+x^2-2x-1 (Q(zeta7)+): ", nfisisom(P1, x^3+x^2-2*x-1) != 0);
print("Galois group of P1: ", polgalois(P1));
P2 = 64*x^3 + 16*x^2 - 164*x + 23;
print("disc(P2) = ", poldisc(P2), " = ", factor(poldisc(P2)), "  field disc = ", nfdisc(P2), " = ", factor(nfdisc(P2)), "  galois: ", polgalois(P2));
print("P2 multiplier 8*prod roots = ", 8*(-23/64));
"""
out = subprocess.run(["gp", "-q"], input=gpin, capture_output=True, text=True)
print(out.stdout.strip()); print(out.stderr.strip())

print("\n== five-vertex tournaments ==")
def digraph6_to_adj(line):
    data = [ord(ch) - 63 for ch in line[1:]]; n = data[0]; bits = []
    for d in data[1:]:
        for b in range(5, -1, -1): bits.append((d >> b) & 1)
    adj = {v: [] for v in range(n)}; idx = 0
    for u in range(n):
        for v in range(n):
            if bits[idx]: adj[u].append(v)
            idx += 1
    return n, adj
def strong(adj, n):
    for s in range(n):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for v in adj[u]:
                if v not in seen: seen.add(v); st.append(v)
        if len(seen) < n: return False
    return True
for n in range(1, 8):
    lines = subprocess.run(["gentourng", "-z", "-q", str(n)], capture_output=True, text=True).stdout.split()
    sc = 0; auts = []; st = 0
    for line in lines:
        nn, adj = digraph6_to_adj(line)
        conv = {v: [] for v in range(nn)}
        for u in adj:
            for v in adj[u]: conv[v].append(u)
        c1 = pynauty.certificate(pynauty.Graph(nn, directed=True, adjacency_dict=adj))
        c2 = pynauty.certificate(pynauty.Graph(nn, directed=True, adjacency_dict=conv))
        if c1 == c2: sc += 1
        if n == 5:
            _, g1, g2, _, _ = pynauty.autgrp(pynauty.Graph(nn, directed=True, adjacency_dict=adj))
            auts.append(int(round(g1 * 10 ** g2)))
            if strong(adj, nn): st += 1
    extra = ""
    if n == 5:
        from collections import Counter
        extra = "  |Aut| distribution %s, strong %d" % (dict(Counter(auts)), st)
    print("  n=%d: %d tournaments, %d self-converse%s" % (n, len(lines), sc, extra))
