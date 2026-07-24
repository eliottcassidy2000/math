#!/usr/bin/env python3
"""
snippet_eq27_fingerprint_opus_S4.py   opus-2026-07-23-S4

Given the object is CERTAIN:  (2457/6592) log(8847357/2974400) - log(1285/896) > 1/25,
hunt for the CONSTRUCTION that produced arguments A=1285/896, B=8847357/2974400.
Diagnostics:
 (a) are A,B near a NAMED constant or its convergents? (irrationality-measure reading)
 (b) are 2974400,8847357 (and 896,1285) hypergeometric = ratios of binomials/factorials?
 (c) do log_A,log_B,LHS match THM-2000 figurate-mass closed forms / differences?
 (d) fingerprint primes 2949119, 3001, 381347, 727, 257 -- values of simple polynomials?
"""
import mpmath as mp
from sympy import isprime, factorint, binomial, factorial
mp.mp.dps = 60

A = mp.mpf(1285) / 896
B = mp.mpf(8847357) / 2974400
tA = mp.mpf(389) / 2181
tB = mp.mpf(5872957) / 11821757
logA, logB = mp.log(A), mp.log(B)
LHS = mp.mpf(2457) / 6592 * logB - logA

print("A =", mp.nstr(A, 18), "  B =", mp.nstr(B, 18))
print("logA =", mp.nstr(logA, 18), " logB =", mp.nstr(logB, 18), " logB/logA =", mp.nstr(logB / logA, 15))
print("tB =", mp.nstr(tB, 18), " (near 1/2? ->", mp.nstr(mp.mpf(1)/2 - tB, 6), ")")
print("LHS =", mp.nstr(LHS, 25), "  margin over 1/25 =", mp.nstr(LHS - mp.mpf(1)/25, 12))

# (a) A,B vs named constants
print("\n(a) A,B vs named constants (value, and A/const):")
named = {'sqrt2': mp.sqrt(2), '2^(1/3)': mp.mpf(2)**(mp.mpf(1)/3), '3^(1/3)': mp.mpf(3)**(mp.mpf(1)/3),
         'e^(1/3)': mp.e**(mp.mpf(1)/3), 'pi/sqrt(5)': mp.pi/mp.sqrt(5), 'golden': (1+mp.sqrt(5))/2,
         '3': mp.mpf(3), 'e': mp.e, 'sqrt(e)': mp.sqrt(mp.e), 'pi-1/6': mp.pi-mp.mpf(1)/6}
for nm, cv in named.items():
    print(f"    A/{nm:9s} = {mp.nstr(A/cv,10):12s}   B/{nm:9s} = {mp.nstr(B/cv,10)}")

# convergents of B (is B a convergent of a nice number? show CF)
print("\n    CF(A):", mp.taylor and [], "->", end=" ")
def cf(x, n=15):
    out = []
    for _ in range(n):
        a = int(mp.floor(x)); out.append(a); x -= a
        if x == 0: break
        x = 1/x
    return out
print("CF(A)=", cf(A), "  CF(B)=", cf(B))
print("    CF(logB/logA)=", cf(logB/logA), "  CF(2457/6592)=", cf(mp.mpf(2457)/6592))
print("    CF(tB)=", cf(tB))

# (b) hypergeometric test: ratios of binomials / factorials near A,B and the integers
print("\n(b) hypergeometric / binomial fingerprint:")
for target in [896, 1285, 2974400, 8847357, 5872957]:
    hits = []
    for n in range(1, 40):
        for k in range(0, n + 1):
            b = int(binomial(n, k))
            if b == target: hits.append(f"C({n},{k})")
            if b and target % b == 0 and 1 < target // b <= 200000 and target // b in [int(binomial(nn,kk)) for nn in range(1,30) for kk in range(nn+1)]:
                pass
    print(f"    {target}: binomial hits = {hits if hits else 'none direct'}")
# central binomials & lcm(1..n) scale
print("    C(2n,n) scale near 2974400:", [(n, int(binomial(2*n, n))) for n in range(8, 13)])
from sympy import ilcm
lc = 1
for n in range(1, 25):
    lc = ilcm(lc, n)
    if 100000 < lc < 10**8:
        print(f"    lcm(1..{n}) = {lc}")

# (c) THM-2000 figurate masses (closed forms) vs logA, logB, LHS
print("\n(c) figurate-mass closed forms vs the logs:")
masses = {
    'M(3,2)=2': mp.mpf(2), 'M(4,2)=pi^2/6': mp.pi**2/6, 'M(5,2)=3log3-pi/sqrt3': 3*mp.log(3)-mp.pi/mp.sqrt(3),
    'M(6,2)=2log2': 2*mp.log(2), 'M(4,3)=18-24log2': 18-24*mp.log(2), 'M(5,3)=pi^2/3-2': mp.pi**2/3-2,
    'M(4,4)=21-2pi^2': 21-2*mp.pi**2, 'sum1/F3=4pi^2/3-12': 4*mp.pi**2/3-12,
}
for nm, mv in masses.items():
    print(f"    {nm:26s}={mp.nstr(mv,12):14s}  -logA:{mp.nstr(mv-logA,6):10s} -logB:{mp.nstr(mv-logB,6):10s}")
# pairwise differences vs logA, logB
print("    pairwise mass differences hitting logA=0.36057 or logB=1.09008:")
mk = list(masses.items())
for i in range(len(mk)):
    for j in range(len(mk)):
        if i == j: continue
        d = mk[i][1] - mk[j][1]
        if abs(d - logA) < 1e-4 or abs(d - logB) < 1e-4 or abs(d - LHS) < 1e-4:
            print(f"      {mk[i][0]} - {mk[j][0]} = {mp.nstr(d,10)}")

# (d) fingerprint primes as polynomial values
print("\n(d) fingerprint primes -- simple closed forms?")
for p in [257, 727, 3001, 381347, 2949119, 103, 389]:
    print(f"    {p}: prime={isprime(p)}", end="  ")
    # near a power?
    for base in [2, 3, 5, 6, 10]:
        e = mp.log(p) / mp.log(base)
        if abs(e - round(e)) < 0.02: print(f"~{base}^{round(e)}", end=" ")
    # near k^2, k^3, k(k+1)/2, 2^k+1 ...
    r2 = mp.sqrt(p); r3 = p**(mp.mpf(1)/3)
    if abs(r2 - round(r2)) < 0.02: print(f"={round(r2)}^2?", end=" ")
    if abs(r3 - round(r3)) < 0.02: print(f"={round(r3)}^3?", end=" ")
    # 2^k+1 / 2^k-1
    lg = mp.log(p - 1) / mp.log(2)
    if abs(lg - round(lg)) < 0.01: print(f"=2^{round(lg)}+1", end=" ")
    lg = mp.log(p + 1) / mp.log(2)
    if abs(lg - round(lg)) < 0.01: print(f"=2^{round(lg)}-1", end=" ")
    print()

# (e) relations among the ratios: is B ~ A^3, B ~ 3, etc.
print("\n(e) structural relations:")
print("    logB - 3*logA =", mp.nstr(logB - 3*logA, 12), "  (B vs A^3)")
print("    log3 - logB   =", mp.nstr(mp.log(3) - logB, 12), "  (B vs 3)")
print("    A^3 =", mp.nstr(A**3, 12), "  B =", mp.nstr(B, 12))
print("    is (2457/6592) ~ logA/logB + margin?  logA/logB=", mp.nstr(logA/logB, 12),
      "  2457/6592=", mp.nstr(mp.mpf(2457)/6592, 12))
