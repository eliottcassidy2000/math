#!/usr/bin/env python3
"""
klein-2026-07-20-S341 -- INDEPENDENT VERIFICATION of the outside GAUSSIAN MOMENT CONJECTURE
counterexample, and a search for improvements.

CLAIM SUPPLIED BY THE OWNER (outside source):
    P = (1+Z2)(W1(1-Z1) + W2),  Q = Z2,  in 4 real Gaussians via complex Z_j, W_j
    E[P^m] = 0 for all m >= 1,  but  E[Q P^m] = m! != 0 for all m >= 1
    => GMC(4) is false.  (cubic, 6 terms)

GMC(n): if E[P^m] = 0 for all m >= 1 then E[Q P^m] = 0 for m >> 0, for every Q.
(the Mathieu-Zhao subspace property of ker E for the Gaussian measure on R^n)

EXACT WICK CALCULUS.  For independent standard complex Gaussians Z_j (circularly symmetric,
E[Z]=E[Z^2]=0, E[Z conj(Z)] = 1):
        E[ prod_j Z_j^{a_j} conj(Z_j)^{b_j} ] = prod_j delta(a_j,b_j) * a_j!
Every expectation below is an EXACT INTEGER computed from that rule -- no sampling.

The statement "4 real Gaussians" forces 2 independent complex Gaussians, so Z1,Z2,W1,W2
cannot all be independent.  Several conjugate-pairings are consistent with the notation;
all are tested and the one that reproduces the claim is identified.
"""
from fractions import Fraction
from itertools import product as iproduct

# ---------------------------------------------------------------- polynomials
# monomial: tuple of 2K ints (a_1,b_1,...,a_K,b_K) = exponents of Z_k and conj(Z_k)
def pmul(A, B):
    C = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = tuple(x + y for x, y in zip(ma, mb))
            C[m] = C.get(m, 0) + ca * cb
    return {m: c for m, c in C.items() if c != 0}

def padd(A, B):
    C = dict(A)
    for m, c in B.items():
        C[m] = C.get(m, 0) + c
        if C[m] == 0: del C[m]
    return C

def pscale(A, s):
    return {m: c * s for m, c in A.items() if c * s != 0}

def ppow(A, k):
    R = {(0,) * len(next(iter(A))): 1} if A else {}
    for _ in range(k): R = pmul(R, A)
    return R

def fact(n):
    r = 1
    for i in range(2, n + 1): r *= i
    return r

def E(A):
    """exact Gaussian expectation via the Wick rule"""
    tot = 0
    for m, c in A.items():
        ok, val = True, 1
        for k in range(0, len(m), 2):
            a, b = m[k], m[k + 1]
            if a != b: ok = False; break
            val *= fact(a)
        if ok: tot += c * val
    return tot

def var(K, k, conj=False):
    """the monomial Z_k  (or conj(Z_k))  in K complex Gaussians"""
    e = [0] * (2 * K)
    e[2 * k + (1 if conj else 0)] = 1
    return {tuple(e): 1}

def const(K, c):
    return {(0,) * (2 * K): c} if c else {}

# ---------------------------------------------------------------- hypotheses
# 2 independent complex Gaussians (= 4 REAL Gaussians).  Name them A and B.
K = 2
A_, Ab = var(K, 0), var(K, 0, True)
B_, Bb = var(K, 1), var(K, 1, True)
one = const(K, 1)

HYPS = {
 "H2: Z1=A, Z2=conj(A), W1=B, W2=conj(B)":      (A_, Ab, B_, Bb),
 "H3: Z1=A, Z2=B,       W1=conj(A), W2=conj(B)":(A_, B_, Ab, Bb),
 "H4: Z1=conj(A), Z2=A, W1=B, W2=conj(B)":      (Ab, A_, B_, Bb),
 "H5: Z1=A, Z2=conj(B), W1=B, W2=conj(A)":      (A_, Bb, B_, Ab),
 "H6: Z1=conj(A),Z2=conj(B),W1=A, W2=B":        (Ab, Bb, A_, B_),
}

print("=" * 80)
print("INDEPENDENT VERIFICATION -- which conjugate pairing realises the claim?")
print("=" * 80)
print(f"{'hypothesis':<48} {'E[P^m], m=1..5':<22} {'E[Q P^m], m=1..5'}")
winner = None
for name, (Z1, Z2, W1, W2) in HYPS.items():
    # P = (1 + Z2) * ( W1*(1 - Z1) + W2 )
    inner = padd(pmul(W1, padd(one, pscale(Z1, -1))), W2)
    P = pmul(padd(one, Z2), inner)
    Q = Z2
    ep = [E(ppow(P, m)) for m in range(1, 6)]
    eq = [E(pmul(Q, ppow(P, m))) for m in range(1, 6)]
    tag = ""
    if all(v == 0 for v in ep) and all(eq[m - 1] == fact(m) for m in range(1, 6)):
        tag = "   <== MATCHES THE CLAIM EXACTLY"; winner = (name, P, Q, Z1, Z2, W1, W2)
    elif all(v == 0 for v in ep) and all(v != 0 for v in eq):
        tag = "   <== disproves GMC (different constant)"
        if winner is None: winner = (name, P, Q, Z1, Z2, W1, W2)
    print(f"{name:<48} {str(ep):<22} {eq}{tag}")

if winner is None:
    print("\nNO PAIRING REPRODUCED THE CLAIM -- reporting that rather than forcing it.")
else:
    name, P, Q, Z1, Z2, W1, W2 = winner
    print("\n" + "=" * 80)
    print(f"CONFIRMED under {name}")
    print("=" * 80)
    print(f"  P has {len(P)} monomials; total degree = {max(sum(m) for m in P)}")
    print(f"  P = {sorted(P.items())}")
    ep = [E(ppow(P, m)) for m in range(1, 11)]
    eq = [E(pmul(Q, ppow(P, m))) for m in range(1, 11)]
    print(f"  E[P^m],   m=1..10 : {ep}")
    print(f"  E[Q P^m], m=1..10 : {eq}")
    print(f"  E[Q P^m] == m!    : {[eq[m-1] == fact(m) for m in range(1,11)]}")
    print("\n  => E[P^m] = 0 for all m tested, E[Q P^m] = m! != 0 for all m tested.")
    print("     GMC IS FALSE IN 4 REAL GAUSSIAN VARIABLES.  VERIFIED EXACTLY, NOT SAMPLED.")

# ================================================================ THE MECHANISM
print("\n" + "=" * 80)
print("THE MECHANISM -- a CLOSED FORM for all m, not a check up to m=10")
print("=" * 80)
print("""
 Write A, B for the two independent standard complex Gaussians (= 4 real).  The verified
 reading is W_j = conj(Z_j), so with Z1=A, Z2=B:

     P = (1 + Bb)( A(1 - Ab) + B ),   Q = Bb          [ Xb := conj(X) ]

 A and B are INDEPENDENT, so expanding U = A(1-Ab) + B by the binomial theorem,
 U^m = sum_k C(m,k) (A(1-Ab))^k B^{m-k}, the expectation FACTORISES:

     E[P^m] = sum_k C(m,k) * E_A[ A^k (1-Ab)^k ] * E_B[ (1+Bb)^m B^{m-k} ]

 Both factors are one-line Wick evaluations, because E[A^p Ab^q] = delta_pq p!:

     E_A[ A^k (1-Ab)^k ]        = (-1)^k k!        (only the Ab^k term survives)
     E_B[ (1+Bb)^m B^{m-k} ]    = C(m, m-k) (m-k)! = C(m,k) (m-k)!

 Hence, using the collapse  C(m,k)^2 k! (m-k)! = m! C(m,k):

     E[P^m]   = sum_k C(m,k)^2 (-1)^k k! (m-k)!  =  m! * sum_k (-1)^k C(m,k)  =  m! (1-1)^m

 and with the Q = Bb insertion shifting the B-index by one, C(m,m-k-1) = C(m,k+1):

     E[Q P^m] = sum_k C(m,k)(-1)^k k! * C(m,k+1)(m-k)!  =  m! * sum_k (-1)^k C(m,k+1)
              = m! * ( 1 - (1-1)^m )

 THE WHOLE COUNTEREXAMPLE IS THE BINOMIAL THEOREM:

     E[P^m]   = m! (1-1)^m       = 0     for all m >= 1
     E[Q P^m] = m! (1 - (1-1)^m) = m!    for all m >= 1

 Same computation, same m!; the Q-insertion just deletes the k = -1 term that made the
 alternating sum vanish.  PROVED FOR ALL m, not verified to m = 10.
""")

# ---- machine-check the closed form against brute Wick, and the two simplifications
def binom(n, k):
    if k < 0 or k > n: return 0
    r = 1
    for i in range(k): r = r * (n - i) // (i + 1)
    return r

print(" machine-check of the closed form vs brute Wick expansion:")
Z1, Z2, W1, W2 = Ab, Bb, A_, B_          # the confirmed H6 pairing
inner = padd(pmul(W1, padd(one, pscale(Z1, -1))), W2)
P = pmul(padd(one, Z2), inner); Q = Z2
for m in range(1, 9):
    lhs, rhs = E(ppow(P, m)), fact(m) * (0 ** m if m > 0 else 1)
    lhs2, rhs2 = E(pmul(Q, ppow(P, m))), fact(m)
    assert lhs == 0 == rhs and lhs2 == rhs2, (m, lhs, lhs2)
print("   m = 1..8: E[P^m] = m!(1-1)^m = 0 and E[QP^m] = m! -- both match exactly.")

print("\n" + "=" * 80)
print("IMPROVEMENT 1 -- THE CONSTANT IN (1 - Z1) IS INERT: 6 terms collapse to 4")
print("=" * 80)
print(" E_A[A^k (c - Ab)^k] = (-1)^k k! for EVERY c, since only the Ab^k term survives.")
print(" So c = 0 is legal and A(0 - Ab) = -|A|^2.  The counterexample shortens to")
print("     P' = (1 + Bb)(B - |A|^2),   Q = Bb      -- 4 monomials, still cubic.")
P2 = pmul(padd(one, Bb), padd(B_, pscale(pmul(A_, Ab), -1)))
ep2 = [E(ppow(P2, m)) for m in range(1, 11)]
eq2 = [E(pmul(Bb, ppow(P2, m))) for m in range(1, 11)]
print(f"   monomials: {len(P2)}   degree: {max(sum(mm) for mm in P2)}")
print(f"   E[P'^m],   m=1..10 : {ep2}")
print(f"   E[Q P'^m], m=1..10 : {eq2}   == m! : {all(eq2[m-1]==fact(m) for m in range(1,11))}")

print("\n" + "=" * 80)
print("IMPROVEMENT 2 -- A ONE-PARAMETER FAMILY (the original is mu = 1)")
print("=" * 80)
print(" P_mu = (1 + mu*Bb)(A(1 - lam*Ab) + B) gives E[P^m] = m!(mu - lam)^m by the same")
print(" computation, so ANY lam = mu kills every moment; and then E[Q P^m] = m! mu^{m-1}.")
for mu in (1, 2, 3, -1):
    Pm = pmul(padd(one, pscale(Bb, mu)), padd(pmul(A_, padd(one, pscale(Ab, -mu))), B_))
    ep = [E(ppow(Pm, m)) for m in range(1, 7)]
    eq = [E(pmul(Bb, ppow(Pm, m))) for m in range(1, 7)]
    pred = [fact(m) * mu ** (m - 1) for m in range(1, 7)]
    print(f"   mu=lam={mu:>2}:  E[P^m]={ep}   E[QP^m]={eq}   == m! mu^(m-1): {eq==pred}")
print("\n so the counterexample is not isolated -- it is a LINE of counterexamples through")
print(" the diagonal lam = mu, degenerating only at mu = 0.")
