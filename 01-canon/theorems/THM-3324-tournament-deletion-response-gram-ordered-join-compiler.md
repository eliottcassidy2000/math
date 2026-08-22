---
id: THM-3324
title: "Tournament deletion-response Gram ordered-join compiler"
status: >
  PROVED by exact algebra + VERIFIED-EXACT in the declared finite universes.
  For s_T=(P_T,zN_T)^T, marked deletion response r_Tv=s_(T-v), and
  Gamma_T=sum_v r_Tv(z)r_Tv(w)^T, ordered join obeys an exact two-coordinate
  split-product law and a closed 2-by-2 Gram congruence law.  Gamma consists
  of three kernels D=sum pp, E=sum p(z)w nu(w), and
  F=sum z nu(z)w nu(w).  F is load-bearing: order-six masks 73 and 83 agree
  on s,D,E,E^T but have a nonzero rank-one F difference, exposed in D after
  joining a singleton; no such reduced-interface collision occurs through
  order five.  In the P+/-zN basis, every Gram channel diagonalizes and has an
  explicit product-sum formula for an arbitrary repeated join; identical
  factors give k times one Gram channel times a channel product to power k-1.
  The full (s,Gamma) interface is commutative and still forgets
  SCC order: K1 ordered-join C3 and the reverse have equal interfaces and
  different score sequences.  No isomorphism, Hamiltonian, path-cover,
  substitution, or chronology consequence follows.
audit: >
  The exact companion checks all 2553 ordered labelled factor pairs of total
  order at most six.  It compares direct joins with both response laws, every
  marked deletion, every Gram block, and the reverse join.  Newton-trace
  characteristic polynomials are independently checked by permutation
  determinants for joins and deletions.  It exhausts the reduced-interface
  minimality boundary through order five and the reverse-order boundary
  through total order three.  Normal and optimized transcripts byte-match;
  the source has zero assertion nodes and zero floating literals.
source: root/creative-synthesis-next/2026-08-03
depends_on:
  - THM-3322-tournament-switching-second-moment-deletion-gram-and-order-join-law
related:
  - THM-3315-tournament-cut-switching-centered-coronal-walk-compiler
  - THM-1862-order-join-reduction-principle
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
script: 04-computation/tournament_deletion_response_gram_order_join_compiler_20260803.py
output: 05-knowledge/results/tournament_deletion_response_gram_order_join_compiler_20260803.out
script_sha256: eb6a315b732aff147a178b5840724df949288e9dee20e84d2bb9c49f9a88816f
output_sha256: 08b6e2dc752830943c699517aa600610cbc084f72f7d5eb53b8601e8218556fe
semantic_sha256: e4e973e83afebe69d28f7da009a4cb5e75c87c53d2d7668b0b24b71666f42fe8
hash_basis: LF-normalized bytes
---

# THM-3324 -- tournament deletion-response Gram ordered-join compiler

**PROVED by exact algebra + VERIFIED-EXACT in the stated finite universes.**

Let `T` be a tournament and retain THM-3322's centered characteristic
polynomial and all-ones adjugate numerator

```text
P_T(z)=det(I-zC_T),
N_T(z)=1^T adj(I-zC_T) 1.                                 (1)
```

Use `P_empty=1,N_empty=0`.

## 1. Split response algebra

Define

```text
s_T(z)=(P_T(z),zN_T(z))^T,
r_Tv(z)=s_(T-v)(z)=(p_v(z),z nu_v(z))^T.                  (2)
```

For pairs, put

```text
(a,b) star (c,d)=(ac+bd,ad+bc),
L_(a,b)=((a,b),(b,a)).                                    (3)
```

The change of coordinates `(a,b)->(a+b,a-b)` diagonalizes `star`, so it is
associative and commutative.  THM-3322's ordered-join law is

```text
s_(X triangleright Y)=s_X star s_Y.                       (4)
```

Deletion respects the factor in which its marked vertex lies:

```text
r_(X triangleright Y),v = r_Xv star s_Y   for v in X,
r_(X triangleright Y),v = s_X star r_Yv   for v in Y.     (5)
```

In particular the new deletion polynomial consumes both old coordinates:

```text
p^J_v=p^X_v P_Y+z^2 nu^X_v N_Y       for v in X.          (6)
```

Thus the scalar deletion Gram cannot be expected to close by itself.

## 2. Closed Gram law

Define the marked-response Gram kernel

```text
Gamma_T(z,w)=sum_v r_Tv(z) r_Tv(w)^T
            =((D_T(z,w),E_T(z,w)),
              (E_T(w,z),F_T(z,w))),                       (7)
```

where

```text
D_T=sum_v p_v(z)p_v(w),
E_T=sum_v p_v(z) w nu_v(w),
F_T=sum_v z nu_v(z) w nu_v(w).                            (8)
```

Put

```text
H_T(z)=((P_T(z),zN_T(z)),(zN_T(z),P_T(z))).               (9)
```

Summing the outer products in `(5)` proves

```text
Gamma_(X triangleright Y)(z,w)
 =H_Y(z) Gamma_X(z,w) H_Y(w)^T
 +H_X(z) Gamma_Y(z,w) H_X(w)^T.                          (10)
```

Equations `(4)` and `(10)` give an exact associative compiler for repeated
ordered joins.  The requested scalar `D_(X triangleright Y)` is the
upper-left entry.  Its first-factor contribution is

```text
P_Y(z)P_Y(w)D_X
+P_Y(z)wN_Y(w)E_X(z,w)
+zN_Y(z)P_Y(w)E_X(w,z)
+zN_Y(z)wN_Y(w)F_X,                                      (11)
```

with the analogous `X`-weighted contribution from `Gamma_Y`.

## 3. Diagonal product-sum closed form

Let

```text
S=((1,1),(1,-1)),
u_(T,+)(z)=P_T(z)+zN_T(z),
u_(T,-)(z)=P_T(z)-zN_T(z),
Gamma_hat_T=S Gamma_T S^T.                                (12)
```

Since `S H_T S^(-1)=diag(u_(T,+),u_(T,-))`, equations `(4)` and `(10)`
become, for channel signs `a,b` in `{+,-}`,

```text
u_(X triangleright Y,a)=u_(X,a)u_(Y,a),

Gamma_hat_(X triangleright Y,ab)(z,w)
 =u_(Y,a)(z)u_(Y,b)(w) Gamma_hat_(X,ab)(z,w)
 +u_(X,a)(z)u_(X,b)(w) Gamma_hat_(Y,ab)(z,w).             (13)
```

Induction gives a nonrecursive formula for
`J=T_1 triangleright ... triangleright T_k`:

```text
u_(J,a)(z)=product_i u_(T_i,a)(z),

Gamma_hat_(J,ab)(z,w)
 =sum_i Gamma_hat_(T_i,ab)(z,w)
        product_(j!=i) u_(T_j,a)(z)u_(T_j,b)(w).          (14)
```

For `k` identical factors this collapses to

```text
Gamma_hat_(T^triangleright k,ab)(z,w)
 =k Gamma_hat_(T,ab)(z,w)
    (u_(T,a)(z)u_(T,b)(w))^(k-1).                        (15)
```

Thus every second-response channel of a repeated join has a closed
product-power expression in the repetition count.  This is an exact
algebraic compiler, not a bit-complexity bound or an injective tournament
encoding.  The original Gram is recovered from
`Gamma_T=(1/4)S Gamma_hat_T S^T`.

## 4. The final Gram coordinate is necessary

In THM-3322's labelled-mask convention, the order-six tournaments with masks
`73` and `83` have score sequences

```text
(5,3,2,2,2,1),             (4,4,3,2,1,1).                (16)
```

They have the same response

```text
P =1+6z+30z^2+80z^3+168z^4+192z^5+128z^6,
zN=6z+30z^2+112z^3+216z^4+256z^5+128z^6,                 (17)
```

and agree on `D,E,E^T`.  Their last Gram blocks differ by

```text
F_73-F_83
 =128 z^3 w^3(1+2z+4z^2)(1+2w+4w^2).                    (18)
```

Joining a singleton on either side exposes this omitted coordinate in the
upper-left block:

```text
D_(73 triangleright K1)-D_(83 triangleright K1)
 =128 z^4 w^4(1+2z+4z^2)(1+2w+4w^2).                    (19)
```

Therefore `(P,N,D)` is insufficient to compose `D`; even adjoining the mixed
kernel `E` is insufficient.  Exhaustion of all labelled tournaments through
order five finds no equal `(s,D,E,E^T)` profile with unequal `F`.  Hence order
six is the sharp first failure of this reduced interface, without classifying
all order-six collisions.

## 5. The full interface still loses SCC order

Commutativity of `star` and exchange of the two summands in `(10)` give

```text
(s_(X triangleright Y),Gamma_(X triangleright Y))
 =(s_(Y triangleright X),Gamma_(Y triangleright X)).       (20)
```

But

```text
K1 triangleright C3: score sequence (3,1,1,1),
C3 triangleright K1: score sequence (2,2,2,0).             (21)
```

The two tournaments are nonisomorphic and have different ordered SCC
condensations, while their full `(s,Gamma)` interfaces agree.  All reverse
joins of total order at most three are isomorphic; a nontrivial strong
tournament needs at least three vertices, so the `1+3` hostile is sharp.

## 6. Audit and scope

The companion checks all `2553` ordered labelled factor pairs of total order
at most six.  For every pair it compares direct calculation with `(4)`, all
marked responses in `(5)`, all four blocks in `(10)`, and the reverse-order
interface.  Newton traces and independent permutation determinants agree for
every joined tournament and marked deletion.  The semantic digests are

```text
join     d78fb46977680cb7b8dba9d0d1c0f77614c866c3768a084321f011e23cc3e8b8
hostile  24654c5a8dc4d41ec641cbd074c8afad1c51a03a6cfdf35ad7e22b9045b8a9e3
semantic e4e973e83afebe69d28f7da009a4cb5e75c87c53d2d7668b0b24b71666f42fe8
```

The interface retains the second tensor moment of marked deletion responses.
It forgets vertex ownership, ordered SCC condensation, higher deletion
moments, individual switching cuts, Hamiltonian and path-cover data,
substitution runs, and chronology.  No reconstruction or classification
claim follows.

Run

```text
python3 04-computation/tournament_deletion_response_gram_order_join_compiler_20260803.py
python3 -O 04-computation/tournament_deletion_response_gram_order_join_compiler_20260803.py
```

and compare LF-normalized bytes with the declared output.

**QED.**
