# opus-2026-07-15-S320 -- HYP-6960: THE mu6 SPINE OF LRC(14).
# [a] chi13-grading of leave-one-out families {1..13}\{k}: exact M, argmax t*,
#     denominator factorization, Eisenstein-norm test.
# [b] the tight locus of {1..13}: all t with min_v ||vt|| = 1/14 (exact scan):
#     conjecture = {p/14 : gcd(p,14) = 1} = mu6 mod 14 = <3>.
# [c] mediant-family denominators q = 3N+2 through the Eisenstein lens.
from fractions import Fraction
from math import gcd

def exact_M(S, qmax=400):
    # M(S) = max_t min_v ||vt|| over t in (0,1): scan all Farey points with
    # q <= qmax (walls have bounded denominators for small speed sets; qmax
    # generous) -- exact rational arithmetic.
    best = (Fraction(0), None)
    for q in range(2, qmax + 1):
        for p in range(1, q):
            if gcd(p, q) != 1: continue
            t = Fraction(p, q)
            m = min(min((v*t) % 1, 1 - (v*t) % 1) for v in S)
            if m > best[0]: best = (m, t)
    return best

def eisenstein_norm(q):
    # q representable as a^2 - a b + b^2 <=> every prime = 2 mod 3 appears to
    # an even power
    n, p, ok = q, 2, True
    while p*p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0: n //= p; e += 1
            if p % 3 == 2 and e % 2 == 1: ok = False
        p += 1
    if n > 1 and n % 3 == 2: ok = False
    return ok

QR13 = {1, 3, 4, 9, 10, 12}
print("[a] leave-one-out families {1..13}\\{k}: exact M, argmax, chi13(k):")
rows = []
for k in range(1, 14):
    S = [v for v in range(1, 14) if v != k]
    M, t = exact_M(S)
    chi = 'QR' if k in QR13 else ('0' if k == 13 else 'QNR')
    q = t.denominator
    rows.append((k, chi, M, t, q, eisenstein_norm(q)))
    print(f"   k={k:2d} ({chi:3s}): M = {M} = {float(M):.6f}  t* = {t}  "
          f"q = {q} (Eisenstein norm: {eisenstein_norm(q)})")
# grading check
byM = {}
for k, chi, M, t, q, e in rows:
    byM.setdefault(M, []).append((k, chi))
print("   M-value classes:")
for M in sorted(byM): print(f"      M = {M}: {byM[M]}")

print("\n[b] the tight locus of the AP {1..13}: all t with min = 1/14 (q <= 200):")
S = list(range(1, 14))
tight = []
for q in range(2, 201):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        t = Fraction(p, q)
        m = min(min((v*t) % 1, 1 - (v*t) % 1) for v in S)
        if m == Fraction(1, 14): tight.append(t)
print(f"   tight times: {tight}")
print(f"   = p/14 with p in (Z/14)* = <3> = mu6: "
      f"{sorted(t.numerator for t in tight) == [1, 3, 5, 9, 11, 13]}")
pw = [pow(3, i, 14) for i in range(6)]
print(f"   powers of 3 mod 14: {pw} (cyclic order 6: {sorted(pw) == [1,3,5,9,11,13]})")

print("\n[c] mediant denominators q = 3N+2 and Phi_6/Eisenstein:")
for N in (7, 13, 19, 25, 31):
    q = 3*N + 2
    print(f"   N={N}: q = {q} = {'Eisenstein norm' if eisenstein_norm(q) else 'NOT Eisenstein'}"
          f"; N mod 6 = {N % 6}")
print(f"   deep well: 183 = Phi_6(14) = {14**2 - 14 + 1} "
      f"(Eisenstein norm: {eisenstein_norm(183)}); 14^k mod 183: "
      f"{[pow(14, i, 183) for i in range(1, 7)]}")
