#!/usr/bin/env python3
"""Lane nt (mac-mini-2026-10-07-oaimath3): structure of the 101 negative discriminants with
class group of exponent <= 2 (output of exp2_sieve.c), the conductor lemma, genus counts,
and the Hilbert-class-field / Hadwiger-Nelson field-plane table.  Our own code; exact arithmetic."""
from math import gcd, isqrt, sqrt, log, pi, e
from itertools import combinations
import sys

vals = [int(x) for x in open(sys.argv[1] if len(sys.argv) > 1 else 'exp2_1e13.out') if not x.startswith('#')]

def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def squarefree(n): return all(v == 1 for v in factor(n).values())

def fundamental(N):  # D = -N
    if N % 4 == 3: return squarefree(N)
    if N % 4 == 0:
        n = N // 4  # D = -4n, m = -n must be 2,3 mod 4  <=>  n = 2,1 mod 4
        return n % 4 in (1, 2) and squarefree(n)
    return False

def fund_part(N):
    """write D=-N = D0 f^2 with D0 fundamental; return (|D0|, f)"""
    best = None
    for f in range(1, isqrt(N) + 1):
        if N % (f * f) == 0 and (N // (f * f)) % 4 in (0, 3) and fundamental(N // (f * f)):
            best = (N // (f * f), f)
    return best

def classno(N):
    h = 0; a = 1
    while 3 * a * a <= N:
        for b in range(-a + 1, a + 1):
            if (b * b + N) % (4 * a): continue
            c = (b * b + N) // (4 * a)
            if c < a or (b < 0 and a == c): continue
            if gcd(gcd(a, abs(b)), c) != 1: continue
            h += 1
        a += 1
    return h

def mu(N):  # number of assigned genus characters (Cox Prop 3.11 / Thm 3.15)
    r = len([p for p in factor(N) if p != 2])
    if N % 4 == 3: return r
    n = N // 4
    if n % 4 == 3: return r
    if n % 4 in (1, 2): return r + 1
    if n % 8 == 4: return r + 1
    return r + 2

fund = [N for N in vals if fundamental(N)]
print(f"total {len(vals)}; fundamental {len(fund)}; non-fundamental {len(vals)-len(fund)}")
ok = all(classno(N) == 2 ** (mu(N) - 1) for N in vals)
print("h(D) = 2^(mu-1) (one class per genus) for all:", ok)
hs = {}
for N in fund: hs.setdefault(classno(N), []).append(N)
print("fundamental by class number:", {h: len(v) for h, v in sorted(hs.items())})
print("  h=1:", hs[1]); print("  largest fundamental:", max(fund), "h =", classno(max(fund)))
# conductor structure
cond = {}
for N in vals:
    D0, f = fund_part(N)
    cond.setdefault(f, []).append((N, D0))
print("conductors f occurring:", sorted(cond), " all divide 840:", all(840 % f == 0 for f in cond))
for f in sorted(cond):
    if f > 1: print(f"  f={f}: (N, |D0|) =", cond[f])
# brute-force: for every fundamental |D0| in list and every f <= 2000, exponent<=2 only if in list
def has_nonamb(N):
    a = 2
    while 3 * a * a <= N:
        for b in range(1 if N & 1 else 2, a, 2):
            if (N + b * b) % (4 * a) == 0 and (N + b * b) // (4 * a) > a: return True
        a += 1
    return False
extra = []
for D0 in fund:
    for f in range(2, 2001):
        N = D0 * f * f
        if N in vals: continue
        # quick: exponent<=2 for D0 f^2 implies f | 840 by the lemma; test all f<=2000 anyway
        if f > 60 and 840 % f: continue  # covered by the lemma (prime >= 11 or p^2 | f or 16 | f)
        if not has_nonamb(N): extra.append((D0, f))
print("extra (D0,f) with exponent<=2 not in list (f<=60 or f|840):", extra)

# idoneal numbers and the odd 'x^2+xy+ny^2' analogue
ido = sorted(N // 4 for N in vals if N % 4 == 0)
odd = sorted((N + 1) // 4 for N in vals if N % 4 == 3)
print("idoneal (65):", len(ido)); print("odd-discriminant analogue n=(N+1)/4 (36):", odd)

# Tatuzawa / EKN Lemma 3 lower bound sanity check on all fundamental D with 73131<=|D|<=3e5? (cheap range)
c = 0.655 / (pi * e)
print(f"EKN Lemma 3 constant 0.655/(pi e) = {c:.5f}")

# ---------------- field planes ----------------
def sqclass_local(d, l):
    """return (is_square_in_Q_l, Q_l(sqrt d)/Q_l unramified) for nonzero integer d"""
    v = 0; u = d
    while u % l == 0: u //= l; v += 1
    if l == 2:
        if v % 2: return (False, False)
        return (u % 8 == 1, u % 4 == 1)
    if v % 2: return (False, False)
    return (pow(u % l, (l - 1) // 2, l) == 1, True)

kappa = {2: '4', 3: '3', 5: '4', 7: '4', 11: '5', 13: '6*', 17: '5-6', 19: '5*'}
def primes_upto(n): return [p for p in range(2, n + 1) if all(p % q for q in range(2, isqrt(p) + 1))]

sq_ido = [m for m in ido if m % 4 == 1 and squarefree(m)]
print("\nsquarefree idoneal m = 1 mod 4 (19 expected):", sq_ido, len(sq_ido))
print("plane F^2, F=Q(sqrt p: p|m); L=F(i)=H(Q(sqrt-m)); Lemma R data per small prime l")
for m in sq_ido:
    ps = sorted(factor(m)) if m > 1 else []
    gens = [-1] + ps
    V = []
    for k in range(len(gens) + 1):
        for S in combinations(gens, k):
            d = 1
            for g in S: d *= g
            V.append(d)
    h = classno(4 * m)
    unram2 = any(p % 4 == 3 for p in ps)
    notes = []
    best = None
    for l in primes_upto(60):
        Vsq = [d for d in V if sqclass_local(d, l)[0]]
        Vun = [d for d in V if sqclass_local(d, l)[1]]
        cD = all(d > 0 for d in Vsq)
        cI = all(d > 0 for d in Vun)
        if not cD: continue
        if cI:  # c in inertia: ramified in L/L+
            bound = 2 if l == 2 else 3
            notes.append(f"{l}:ram->{bound}")
        else:
            bound = kappa.get(l, '?')
            notes.append(f"{l}:inert,q={l}->k({l})={bound}")
        if best is None: best = notes[-1]
    tri = 3 in ps
    lower = []
    if tri: lower.append('triangle>=3')
    if 3 in ps and 11 in ps: lower.append('Moser>=4')
    if all(p in ps for p in (3, 5, 7)): lower.append('L_9 spindle>=4')
    if all(p in ps for p in (3, 5, 11)): lower.append('Heule Z[w1,w3,w4]>=5')
    print(f"m={m:5d} h(-4m)={h:2d} F=Q({','.join('sqrt'+str(p) for p in ps) or '1'}) "
          f"L/L+ unram at 2: {unram2}; c-fixed primes l<60: {notes[:4]}; lower: {lower}")
