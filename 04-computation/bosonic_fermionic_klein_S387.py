#!/usr/bin/env python3
"""
klein-2026-07-20-S387 -- THE BOSONIC/FERMIONIC DICHOTOMY: the determinant side (Vandermonde,
alternating, TRANSITIVITY, cancels -- THM-1805 + boxeph THM-1800) has a symmetric partner, the
PERMANENT/HAFNIAN side, and THAT is the Gaussian moment engine.  The dichotomy explains WHY
GMC(2) is hard while transitivity is a clean classical fact.

Owner: "think more about discriminants and determinants and other similar concepts."

  FERMIONIC / alternating:  det, Pfaffian, Vandermonde, discriminant.  SIGNED sum -> intransitivity
     CANCELS by the 3-cycle involution (THM-1805); only n! transitive tournaments survive.
  BOSONIC / symmetric:      permanent, hafnian.  UNSIGNED sum -> NOTHING cancels; every tournament
     contributes with +.  E[Z^a Zb^b] = delta_ab a! = per(J_a) is the BOSONIC (Gaussian) moment.
  => the Laplace moment engine lives on the PERMANENT side.  No sign, no cancellation -- which is
     exactly why it has no finite-cutoff detection (the EMP floor grows with degree, THM-1790),
     unlike the determinant/transitivity side that collapses in one alternating step.
"""
import itertools
from math import factorial
from fractions import Fraction as Fr

def perm(M):
    n=len(M); tot=0
    for s in itertools.permutations(range(n)):
        p=1
        for i in range(n): p*=M[i][s[i]]
        tot+=p
    return tot
def det(M):
    n=len(M); tot=0
    for s in itertools.permutations(range(n)):
        sg=1; seen=[False]*n
        for i in range(n):
            if seen[i]: continue
            j=i;L=0
            while not seen[j]: seen[j]=True;j=s[j];L+=1
            if L%2==0: sg=-sg
        p=sg
        for i in range(n): p*=M[i][s[i]]
        tot+=p
    return tot

print("="*82)
print("(1) THE COMPLEX GAUSSIAN MOMENT IS THE PERMANENT (bosonic Wick)")
print("="*82)
print(" E[|Z|^{2a}] = a! = per(J_a) (all-ones a x a permanent = #matchings of a Z's to a Zb's):")
for a in range(1,7):
    Ja=[[1]*a for _ in range(a)]
    print(f"   a={a}: E[|Z|^{2*a}] = {factorial(a)}   per(J_{a}) = {perm(Ja)}   equal: {factorial(a)==perm(Ja)}")
print("   (contrast det(J_a)=0 for a>=2: the FERMIONIC moment vanishes -- Pauli exclusion.)")
print(f"   det(J_2)={det([[1,1],[1,1]])}, det(J_3)={det([[1,1,1]]*1+[[1,1,1]]*2)}")

print("\n"+"="*82)
print("(2) det vs per OF THE VANDERMONDE MATRIX = signed vs unsigned tournament sum")
print("="*82)
print(" V_ij = x_i^{j} (0-indexed powers).  det V = Vandermonde = sum_T sgn(T) x^score (transitivity);")
print(" per V = sum_T x^score (ALL tournaments, unsigned).  Verify the tournament readings:")
from collections import defaultdict
def tour_sums(n):
    P=[(i,j) for i in range(n) for j in range(i+1,n)]
    signed=defaultdict(int); unsigned=defaultdict(int)
    for bits in range(1<<len(P)):
        wins=[0]*n; up=0
        for k,(i,j) in enumerate(P):
            if bits>>k&1: wins[i]+=1
            else: wins[j]+=1; up+=1
        signed[tuple(wins)]+=(-1)**up; unsigned[tuple(wins)]+=1
    return signed,unsigned
for n in (3,4):
    signed,unsigned=tour_sums(n)
    surv_signed=sum(1 for c in signed.values() if c!=0)
    total_tour=2**(n*(n-1)//2)
    n_fact=factorial(n)
    print(f" n={n}: det side -- nonzero-coeff scores (transitive survivors) = {surv_signed} = n! = {n_fact}: {surv_signed==n_fact}")
    print(f"        per side -- sum of all unsigned coeffs = {sum(unsigned.values())} = #tournaments = {total_tour}: {sum(unsigned.values())==total_tour}")
    print(f"        per V(1,1,..,1) = {perm([[1]*n for _ in range(n)])}  (=n! only because all x_i equal; the GRADED per is the score enumerator)")

print("\n"+"="*82)
print("(3) THE LOAD-BEARING DICHOTOMY")
print("="*82)
print("""
  FERMIONIC (det / Pfaffian / Vandermonde / discriminant):  a SIGN twists the sum, and the
  3-cycle-reversing involution CANCELS all intransitivity in ONE step -> only the n! transitive
  tournaments survive (THM-1805).  Transitivity is therefore a clean, finite, classical fact.

  BOSONIC (permanent / hafnian / the GAUSSIAN moment E[Z^a Zb^b]=delta a!):  NO sign.  Nothing
  cancels; every tournament / every Wick matching contributes with +.  The moment engine
  E[P^m]=L_s(CT_u[Lambda_s^m]) is this permanent/hafnian functional.  With no alternating
  collapse available, detecting its nullcone cannot be done in one step or at a finite moment
  cutoff -- which is EXACTLY the EMP floor: detection depth >= d+1, growing with radial degree
  (THM-1790).  The permanent has no cheap vanishing (it is #P-hard; the determinant is easy).

  SO: GMC(2) is hard for the SAME structural reason the permanent is hard and the determinant is
  easy.  boxeph's transitivity/discriminant results live on the tractable FERMIONIC side; GMC(2)
  lives on the intractable BOSONIC side.  The bridge that couples them (angular DvdK = fermionic
  GIT nullcone; radial EMP = bosonic Laplace) must cross the det/per divide -- which is why it
  cannot be finite-degree-uniform (THM-1770/1790).
""")
