"""Self-converse decomposition and labelled prime counts from the TYPE census (orders 3..10).
(1) sc(n) = SC_nonstrong(n) + SC_strong_nonprime(n) + SCP(n): validated against the census for n<=10, applied at n=11.
(2) labelled strong s_lab(n) from 2^C(n,2); labelled strong non-prime at n=11 from P_lab(k), k<=10; hence P_lab(11) exactly."""
import re, os, sys
from fractions import Fraction
from math import factorial, comb
NMAX = 11
exec(open('species11.py').read().split('print("types by order:"')[0])   # loads types, parse_ms
def polya_ms(ms, series_fn, n):
    """[x^n] sum_{(ct,mult)} mult * prod_l (series_fn(l) evaluated at x^l)^{c_l}; series_fn(l) returns list coeffs s[0..n]"""
    tot = Fraction(0)
    for ct, mult in ms:
        poly = [Fraction(0)] * (n + 1); poly[0] = Fraction(1)
        for l, cnum in ct.items():
            ser = series_fn(l)
            base = [Fraction(0)] * (n + 1)
            for j in range(1, n + 1):
                if l * j <= n: base[l * j] = Fraction(ser[j])
            for _ in range(cnum):
                new = [Fraction(0)] * (n + 1)
                for a_ in range(n + 1):
                    if poly[a_] == 0: continue
                    for b in range(n + 1 - a_):
                        if base[b]: new[a_ + b] += poly[a_] * base[b]
                poly = new
        tot += mult * poly[n]
    return tot
A = [0, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 0]          # all classes, n<=10 (verified by species)
S = [0, 1, 0, 1, 1, 6, 35, 353, 6008, 178133, 9355949, 0]            # strong classes
census_sc = [0, 1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784]
SCP = {k: sum(t[1] for t in types if t[0] == k and t[4] == 1) for k in range(3, 11)}
print("SC primes by order:", SCP)
sc = [0] * 13; scs = [0] * 13; sc[1] = scs[1] = 1
def seqp(p, t):  # ordered p-tuples of strong classes with total size t
    f = [Fraction(0)] * (t + 1); f[0] = Fraction(1)
    for _ in range(p):
        g = [Fraction(0)] * (t + 1)
        for i in range(t + 1):
            if f[i] == 0: continue
            for j in range(1, t + 1 - i): g[i + j] += f[i] * S[j]
        f = g
    return f[t]
def sc_nonstrong(m):
    tot = Fraction(0)
    for p in range(1, m // 2 + 1):
        if m % 2 == 0: tot += seqp(p, m // 2)
        for j in range(1, m - 2 * p + 1):
            if (m - j) % 2 == 0 and (m - j) // 2 >= p: tot += seqp(p, (m - j) // 2) * scs[j]
    return tot
def sc_strong_nonprime(m):
    tot = Fraction(0)
    for (k, cnt, G, Aut, isc, anti) in types:
        if k >= m or not isc: continue
        orb = polya_ms(Aut, lambda l: A, m) / G
        fixed_anti = polya_ms(anti, lambda l: (sc if l % 2 == 1 else A), m)
        orbpm = (polya_ms(Aut, lambda l: A, m) + fixed_anti) / (2 * G)
        tot += cnt * (2 * orbpm - orb)
    return tot
for m in range(2, 12):
    ns = sc_nonstrong(m); snp = sc_strong_nonprime(m)
    if m <= 10:
        sc[m] = ns + snp + SCP.get(m, 0); scs[m] = sc[m] - ns
        print(f"n={m}: SC nonstrong={ns} SC strong nonprime={snp} SC prime={SCP.get(m,0)} total={sc[m]} census={census_sc[m]} {'OK' if sc[m]==census_sc[m] else 'MISMATCH'}")
    else:
        print(f"n=11: SC nonstrong={ns} SC strong nonprime={snp}  => SCP(11) = SC(11) - {ns + snp}")
# labelled
tl = [2 ** comb(m, 2) for m in range(0, 12)]
T = [Fraction(tl[m], factorial(m)) for m in range(12)]  # EGF, T[0]=1
# S_lab = 1 - 1/T : compute 1/T series
inv = [Fraction(0)] * 12; inv[0] = Fraction(1)
for m in range(1, 12): inv[m] = -sum(T[j] * inv[m - j] for j in range(1, m + 1))
slab = [int(-inv[m] * factorial(m)) if m >= 1 else 0 for m in range(12)]
print("labelled strong s_lab(n), n=1..11:", slab[1:])
Plab = {k: sum(Fraction(cnt * factorial(k), G) for (kk, cnt, G, Aut, isc, anti) in types if kk == k) for k in range(3, 11)}
def powser(p, k, n):
    r = [Fraction(0)] * (n + 1); r[0] = Fraction(1)
    for _ in range(k):
        q = [Fraction(0)] * (n + 1)
        for i in range(n + 1):
            if r[i] == 0: continue
            for j in range(n + 1 - i): q[i + j] += r[i] * p[j]
        r = q
    return r
Tp = [Fraction(0)] + T[1:]
for m in range(3, 12):
    nonprime = sum(Plab[k] / factorial(k) * powser(Tp, k, m)[m] for k in range(3, min(m, 11)) if k in Plab) * factorial(m)
    plab_m = slab[m] - nonprime
    print(f"n={m}: labelled strong={slab[m]} labelled strong non-prime={int(nonprime)} => P_lab({m})={int(plab_m)}" + (f" census={int(Plab[m])} {'OK' if plab_m==Plab[m] else 'MISMATCH'}" if m in Plab else ""))
    if m == 11: print("ell(11) = 997585691239358 + P_lab(11)/2 =", 997585691239358 + plab_m / 2)
