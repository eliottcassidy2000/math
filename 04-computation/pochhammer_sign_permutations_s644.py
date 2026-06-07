#!/usr/bin/env python3
"""
S644 / HYP-2322 — Rising & falling factorials carry the Sn / An (sign) structure.

Continuing the Galois arc (S642 commutator depth, S643 Aut solvable / Aut(T) <= An)
through the rising/falling factorial lens.

  rising  x^(n) = x(x+1)...(x+n-1)  [ascPochhammer]:  coeffs = UNSIGNED Stirling 1st
          c(n,k) = #{perms with k cycles}; at x=1 -> sum c(n,k) = n! = |Sn| (ALL perms).
  falling (x)_n = x(x-1)...(x-n+1)   [descPochhammer]: coeffs = SIGNED s(n,k)=(-1)^{n-k}c(n,k);
          sign(perm with k cycles) = (-1)^{n-k}; at x=1 -> sum s(n,k) = #even - #odd
          = (1)_n = 0 for n>=2  ==>  |An| = n!/2.
  antipode sigma:x->-x exchanges them: (x)_n = (-1)^n (-x)^(n).

Verifies the dictionary against ACTUAL permutation counts, and ties to S643.
No external libs.
"""
from itertools import permutations
from math import factorial

# ---------- Stirling numbers of the first kind via the polynomials ----------
def poly_mul(p, q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): r[i+j]+=a*b
    return r
def falling_coeffs(n):       # (x)_n = prod_{i=0}^{n-1} (x - i) ; coeffs low->high
    p=[1]
    for i in range(n): p=poly_mul(p, [-i, 1])   # (x - i)
    return p
def rising_coeffs(n):        # x^(n) = prod_{i=0}^{n-1} (x + i)
    p=[1]
    for i in range(n): p=poly_mul(p, [i, 1])    # (x + i)
    return p

# ---------- actual permutation cycle counts ----------
def cycles(p):
    n=len(p); seen=[False]*n; c=0
    for i in range(n):
        if not seen[i]:
            c+=1; j=i
            while not seen[j]: seen[j]=True; j=p[j]
    return c
def perm_sign(p):
    return (-1)**(len(p)-cycles(p))

print("="*68)
print("(A) falling-factorial coeffs = signed cycle counts; rising = unsigned")
print("="*68)
print("  n | (x)_n coeffs (signed s)        | x^(n) coeffs (unsigned c) | check vs perms")
for n in range(0,7):
    fc=falling_coeffs(n); rc=rising_coeffs(n)
    # actual cycle-count histogram
    hist=[0]*(n+1)
    for p in permutations(range(n)):
        hist[cycles(p)]+=1
    # unsigned rising coeff k should equal hist[k]; signed falling coeff k = (-1)^{n-k} hist[k]
    ok_r = all(rc[k]==hist[k] for k in range(n+1))
    ok_f = all(fc[k]==((-1)**(n-k))*hist[k] for k in range(n+1))
    print(f"  {n} | {str(fc):30s} | {str(rc):24s} | rising=c:{ok_r} falling=s:{ok_f}")

print("\n" + "="*68)
print("(B) evaluate at x=1:  rising -> n! = |Sn|;  falling -> 0 = #even-#odd")
print("="*68)
print("  n | x^(n)|_1 (=sum c) | n!=|Sn| | (x)_n|_1 (=#even-#odd) | #even | #odd | |An|")
for n in range(0,8):
    r1=sum(rising_coeffs(n)); f1=sum(falling_coeffs(n))
    perms=list(permutations(range(n)))
    ev=sum(1 for p in perms if perm_sign(p)==1)
    od=sum(1 for p in perms if perm_sign(p)==-1)
    print(f"  {n} | {r1:16d} | {factorial(n):7d} | {f1:22d} | {ev:5d} | {od:4d} | {ev}")
print("  -> rising@1 = n! = |Sn| (all perms, the generic Galois group);")
print("     falling@1 = 0 for n>=2  <=>  #even=#odd  <=>  |An|=n!/2 (the sign/index-2).")
print("     (n=0,1: falling@1 = 1, no balance yet -- below the An threshold, like S2,S3")
print("      being solvable; the '0' switches on exactly when the sign map is onto.)")

print("\n" + "="*68)
print("(C) the antipode sigma:x->-x exchanges rising<->falling")
print("="*68)
for n in range(0,6):
    fc=falling_coeffs(n); rc=rising_coeffs(n)
    # (x)_n should equal (-1)^n * (-x)^(n): coeff_k of (-1)^n*rising(-x) = (-1)^n*(-1)^k*rc[k]
    ok=all(fc[k]==((-1)**n)*((-1)**k)*rc[k] for k in range(n+1))
    print(f"  n={n}: (x)_n = (-1)^n (-x)^(n) ? {ok}")
print("  -> the arc's antipodal involution sigma (S638/S643) IS the rising<->falling")
print("     duality. Sn(all=rising) and An(sign=falling) are exchanged by sigma:x->-x;")
print("     on tournaments sigma is the exiled anti-automorphism (S643) -- same sigma.")

print("\n" + "="*68)
print("(D) tie to S643: |An| = n!/2 is the index of the sign map; Aut(T) <= An always")
print("="*68)
print("  falling@1=0 (n>=2)  <=>  sign:Sn->{+-1} is ONTO  <=>  [Sn:An]=2.")
print("  S643: every tournament Aut(T) lands in An (odd-order => even perms). So the")
print("  perspective group sits in the kernel of the sign map = the falling-factorial's")
print("  'even half'. The discriminant (square <=> Galois<=An) is the sign map's invariant;")
print("  the falling factorial vanishing at 1 is its global count. Rising = the whole Sn")
print("  (Abel-Ruffini side), falling = the An/sign half (Feit-Thompson side).")
