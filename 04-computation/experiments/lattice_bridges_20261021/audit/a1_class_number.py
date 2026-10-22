#!/usr/bin/env python3
"""AUDIT of THM-4751 (1), (2), (5): Dirichlet's formula A_q - B_q = (2 - (2/q)) h(-q) for primes q = 3 mod 4, q > 3,
with h(-q) computed by COUNTING REDUCED BINARY QUADRATIC FORMS (no L-function, no PARI), B_q by Euler's criterion.
Also: the q with B_q <= 1; the Heegner primes; which q attain the largest possible B_q = (m-3)/2 for q > 11;
the number of primes in the range; the flip-distance minimum over multipliers; class number one list < 20000.
Output: a1_class_number.out (and a1_hlist.txt, a list 'q h' for the PARI cross-check a1_class_number.gp).
"""
import os, math
HERE = os.path.dirname(os.path.abspath(__file__))
LIM = 20000

def primes_upto(N):
    s = bytearray([1]) * (N + 1); s[0] = s[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]

def h_reduced_forms(D):
    """class number of the (negative, fundamental) discriminant D = -q by counting reduced forms (a,b,c):
    b^2 - 4ac = D, |b| <= a <= c, and b >= 0 whenever |b| == a or a == c."""
    q = -D; h = 0; a = 1
    while 3 * a * a <= q:
        for b in range(-a + 1, a + 1):
            if (b - D) % 2: continue            # b = D mod 2
            num = b * b - D
            if num % (4 * a): continue
            c = num // (4 * a)
            if c < a: continue
            if b < 0 and (a == c): continue      # (a,-b,a) ~ (a,b,a)
            h += 1
        a += 1
    return h

def legendre(a, q):
    t = pow(a % q, (q - 1) // 2, q)
    return 0 if t == 0 else (1 if t == 1 else -1)

out = open(os.path.join(HERE, 'a1_class_number.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

P = [q for q in primes_upto(LIM) if q % 4 == 3]
say(f"primes q = 3 mod 4 below {LIM}: {len(P)} (including q = 3); with 7 <= q < {LIM}: {len([q for q in P if q >= 7])}")
bad = []; small = []; maxB = []; h1 = []; hlist = []
for q in P:
    m = (q - 1) // 2
    h = h_reduced_forms(-q); hlist.append((q, h))
    if h == 1: h1.append(q)
    if q == 3: continue
    chi2 = legendre(2, q)
    A = sum(1 for a in range(1, m + 1) if legendre(a, q) == 1); B = m - A
    if A - B != (2 - chi2) * h: bad.append(q)
    # multiplier-minimised flip distance from R_q = {1..m} to P_q = QR: min over a of #{k <= m : -k in aQR}
    QR = {x * x % q for x in range(1, q)}
    d_res = sum(1 for k in range(1, m + 1) if (q - k) in QR)                 # a in QR
    d_non = sum(1 for k in range(1, m + 1) if (q - k) not in QR)             # a in NQR (aQR = NQR)
    if min(d_res, d_non) != B or d_res != B or d_non != A: bad.append(('dist', q))
    if B <= 1: small.append((q, B, h))
    if q > 11 and B == (m - 3) // 2 and (m - 3) % 2 == 0: maxB.append((q, h, chi2))
    if q > 11 and B > (m - 3) / 2: bad.append(('Bbound', q))
say(f"Dirichlet formula A_q - B_q = (2 - (2/q)) h(-q) with h from reduced forms, flip distance = B_q = d(a in QR) < d(a in NQR) = A_q: "
    f"mismatches = {bad}")
say(f"B_q <= 1 at (q, B, h): {small}")
say(f"class number one among these primes (q < {LIM}): {h1}")
say(f"Heegner primes B_q: {[(q, (q - 1)//2 - sum(1 for a in range(1, (q - 1)//2 + 1) if legendre(a, q) == 1)) for q in (7, 11, 19, 43, 67, 163)]}")
say(f"q > 11 attaining the maximal possible B_q = (m-3)/2 [(2-(2/q))h = 3]: (q, h, (2/q)) = {maxB}")
say("  -> 'least correlated with R_q' is shared by the h = 1 primes q = 3 mod 8 and the h = 3 primes q = 7 mod 8")
# arc-overlap fraction
q = 23; m = 11; QR = {x * x % q for x in range(1, q)}
agree = sum(1 for i in range(q) for k in range(1, m + 1) if k in QR)          # arcs i -> i+k common to R and P
say(f"arc overlap check q=23: common arcs {agree} of {q*m} = {agree/(q*m):.6f}; formula 1/2 + h/(q-1) with h=3: {0.5 + 3/(q-1):.6f}")
open(os.path.join(HERE, 'a1_hlist.txt'), 'w').write(''.join(f"{q} {h}\n" for q, h in hlist))
out.close()
