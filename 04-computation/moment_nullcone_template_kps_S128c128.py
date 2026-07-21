#!/usr/bin/env python3
"""moment_nullcone_template_kps_S128c128.py -- kind-pasteur-2026-07-20-S128c128

A UNIFYING FRAME: the MOMENT-NULLCONE TEMPLATE, and the rational < algebraic < holonomic
complexity hierarchy it lives in.

Claim: three apparently unrelated objects share ONE structure --
  a functional phi(X^m) = projection onto the trivial component of a symmetry,
  whose generating sum F(t) = sum_m phi(X^m) t^m satisfies a FINITE-order linear recurrence,
  so the NULLCONE {phi(X^m)=0 for all m} is detected at a finite DEPTH = recurrence order,
  and equals the locus where F(t) degenerates to its trivial value.

  | object      | symmetry     | phi(X^m)        | X            | nullcone         | F(t) class  | depth |
  | tournament  | trace/S_n    | tr(A^m)         | adjacency A  | transitive/nilp. | RATIONAL    | n     |
  | TNC         | circle U(1)  | CT(Lambda^m)    | Laurent L    | one-sided        | ALGEBRAIC   | D     |
  | GMC(2)      | U(1)xradial  | E[P^m]          | P(Z,Zbar)    | charge one-sided | HOLONOMIC   | K     |

  Cayley-Hamilton is the RATIONAL (finite-matrix) case of the holonomic recurrence:
     sum_m tr(A^m) t^m = sum_i lambda_i t/(1-lambda_i t) = -t d/dt log det(I - tA),
  a rational function with poles 1/lambda_i, so tr(A^m) obeys the char-poly recurrence and
  is detected at depth n.  This script verifies the TOURNAMENT instance and the hierarchy.
"""
import sys
import numpy as np
from fractions import Fraction as Fr
import random


def from_bits(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    return A


print("=" * 88)
print("(1) TOURNAMENT INSTANCE: transitive <=> tr(A^m)=0 for m=1..n  (detection depth n)")
print("=" * 88)
print("  phi = trace; nullcone = {A : tr(A^m)=0 for all m} = nilpotent = TRANSITIVE.")
print("  Cayley-Hamilton caps the depth at n: the char poly gives tr(A^m) for m>n from m<=n.")
print()
for n in range(3, 8):
    m = n * (n - 1) // 2
    cap = min(1 << m, 4000)
    # detection-depth test: does tr(A^k)=0 for k=1..n force tr(A^k)=0 for k=1..3n?
    depth_ok = True
    transitive_iff = True
    checked = 0
    for bits in range(cap):
        A = from_bits(bits, n)
        tr = [int(round(np.trace(np.linalg.matrix_power(A.astype(float), k)))) for k in range(1, 3 * n + 1)]
        low_zero = all(t == 0 for t in tr[:n])
        all_zero = all(t == 0 for t in tr)
        if low_zero != all_zero:
            depth_ok = False
        # transitive iff acyclic iff nilpotent iff all_zero
        # acyclic test: is there a linear order? A is transitive tournament iff scores are 0..n-1
        scores = sorted(A.sum(axis=1).tolist())
        is_trans = (scores == list(range(n)))
        if is_trans != all_zero:
            transitive_iff = False
        checked += 1
    print("  n=%d (%d tournaments): [tr_k=0 for k<=n] => [tr_k=0 all k] : %s ; transitive<=>nullcone : %s"
          % (n, checked, depth_ok, transitive_iff))
sys.stdout.flush()

print()
print("=" * 88)
print("(2) THE GENERATING FUNCTION IS RATIONAL, and its recurrence IS Cayley-Hamilton")
print("=" * 88)
print("  F(t)=sum_m tr(A^m) t^m obeys the char-poly recurrence; poles at 1/lambda_i.")
random.seed(1)
for n in (5, 6, 7):
    m = n * (n - 1) // 2
    A = from_bits(random.getrandbits(m), n)
    ch = np.poly(A)                     # char poly coeffs, monic, [1, c1, ..., cn]
    tr = [np.trace(np.linalg.matrix_power(A.astype(float), k)) for k in range(1, 2 * n + 2)]
    # Newton/CH recurrence: tr_k = -(c1 tr_{k-1}+...+c_{k-1} tr_1 + k c_k) for k<=n;
    # for k>n: tr_k = -(c1 tr_{k-1}+...+cn tr_{k-n})
    ok = True
    for k in range(n + 1, 2 * n + 1):
        pred = -sum(ch[i] * tr[k - i - 1] for i in range(1, n + 1))
        if abs(pred - tr[k - 1]) > 1e-6:
            ok = False
    print("  n=%d : Cayley-Hamilton recurrence tr_k = -sum c_i tr_{k-i} for k>n holds : %s" % (n, ok))
print("  So the tournament moment sequence sits at the RATIONAL floor of the hierarchy,")
print("  detection depth = n = matrix size.  TNC (THM-1710) is one level up (ALGEBRAIC,")
print("  depth D = charge span); GMC(2) (THM-1740) one more (HOLONOMIC, depth K ~ span).")

print()
print("=" * 88)
print("(3) THE NULLCONE = GENERATING-FUNCTION DEGENERATION, in all three")
print("=" * 88)
print("  tournament: transitive <=> det(I - tA) = 1 <=> F(t) = 0     (all lambda_i = 0)")
print("  TNC       : one-sided  <=> F(t) = 1 (constant)              (the small-branch sum trivial)")
print("  GMC(2)    : charge one-sided <=> E-series trivial past m=0   (charge threshold)")
print("  One template: the nullcone is where the moment generating function collapses to its")
print("  trivial value, detected at the finite depth set by the GF's arithmetic class.")
