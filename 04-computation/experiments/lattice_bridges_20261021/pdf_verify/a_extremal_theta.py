#!/usr/bin/env python3
"""(a) Extremal theta series of even unimodular lattices in dims 24, 48, 72, 96 (and 80).
Exact integer arithmetic. Basis E4^(k/4-3j) Delta^j of M_k(SL2(Z)), k = n/2.
For n = 24m the extremal form is the unique f in M_{12m} with f = 1 + O(q^{m+1}).
For n = 80 (k = 40, dim M_40 = 4) the 'extremal' form is 1 + O(q^{4}) (min 8 = 2*floor(80/24)+2)."""
from fractions import Fraction
N = 40  # number of q-coefficients

def sigma(n, k):
    return sum(d**k for d in range(1, n+1) if n % d == 0)

def mul(a, b):
    c = [0]*N
    for i, x in enumerate(a):
        if x == 0: continue
        for j in range(N - i):
            c[i+j] += x*b[j]
    return c

E4 = [1] + [240*sigma(n, 3) for n in range(1, N)]
# Delta = q prod (1-q^n)^24
P = [1] + [0]*(N-1)
for n in range(1, N):
    for _ in range(24):
        P = [P[i] - (P[i-n] if i >= n else 0) for i in range(N)]
Delta = [0] + P[:N-1]
assert Delta[:6] == [0, 1, -24, 252, -1472, 4830]

def power(f, e):
    r = [1] + [0]*(N-1)
    for _ in range(e):
        r = mul(r, f)
    return r

def extremal(k, zeros):
    """unique f in M_k with f = 1 + 0 q + ... + 0 q^zeros + ...; basis E4^a Delta^j with 4a+12j = k
    (k multiple of 4 here)."""
    J = k // 12
    basis = []
    for j in range(J+1):
        a4 = (k - 12*j)
        assert a4 % 4 == 0
        basis.append(mul(power(E4, a4//4), power(Delta, j)))
    assert zeros == J
    f = [0]*N
    coeffs = []
    for j in range(J+1):
        target = 1 if j == 0 else 0
        a = target - f[j]
        assert basis[j][j] == 1
        coeffs.append(a)
        f = [f[i] + a*basis[j][i] for i in range(N)]
    return f, coeffs

def bern(n):
    B = [Fraction(0)]*(n+1)
    B[0] = Fraction(1)
    from math import comb
    for m in range(1, n+1):
        B[m] = -sum(comb(m+1, k)*B[k] for k in range(m))/(m+1)
    return B

B = bern(60)
out = []
for n in (24, 48, 72, 96):
    k = n//2
    f, coeffs = extremal(k, k//12)
    out.append(f"dim {n}: weight {k}, dim M_k = {k//12+1}; coefficients a_j of E4^(.)Delta^j: {coeffs}")
    out.append(f"  extremal theta: " + ", ".join(f"q^{i}:{f[i]}" for i in range(0, 9)))
    # Eisenstein comparison
    Ek = [Fraction(1)] + [Fraction(-2*k, 1)/B[k]*sigma(i, k-1) for i in range(1, N)]
    out.append(f"  Siegel average E_{k}: " + ", ".join(f"q^{i}:{float(Ek[i]):.6g}" for i in range(1, 8)))
    neg = [i for i in range(N) if f[i] < 0]
    out.append(f"  negative coefficients among q^0..q^{N-1}: {neg}")
print("\n".join(out))
# explicit answer for dim 96
f96, _ = extremal(48, 4)
print("\nDIM 96 extremal: q^5 =", f"{f96[5]:,}", " q^6 =", f"{f96[6]:,}", " q^7 =", f"{f96[7]:,}")
print("paper 3 q^6 = 2,972,108,280,960,000 ->", f96[6] == 2972108280960000)
print("paper 4 q^6 = 184,139,089,776,640 ->", f96[6] == 184139089776640)
print("both papers q^5 = 565,866,362,880 ->", f96[5] == 565866362880)
# what could paper 4's number be? test a few possibilities
cands = {}
import itertools
for i in range(1, 12):
    if f96[i] == 184139089776640: cands['coef'] = i
print("paper4 q^6 value equals some coefficient of the extremal form?", cands)
# Siegel average numbers of short vectors in dimension 96 (mass-weighted average over the genus II_96)
E48 = [Fraction(1)] + [Fraction(-96, 1)/B[48]*sigma(i, 47) for i in range(1, 8)]
print("B_48 =", B[48])
for i in range(1, 7):
    print(f"  E_48 coefficient q^{i} (avg # vectors of norm {2*i}): {float(E48[i]):.6e}")
