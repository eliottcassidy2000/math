---
id: THM-3322
title: "Tournament switching second moment, deletion Gram, and order-join law"
status: >
  PROVED by exact algebra + VERIFIED-EXACT in the declared finite universes.
  For centered tournament adjacency C, P(z)=det(I-zC),
  B(z)=adj(I-zC), and N_d(z)=d^T B(z)d, the uniform cut-cube second moment is
  determined exactly by P plus the vertex-deletion Gram kernel
  D_T(z,w)=sum_i B_ii(z)B_ii(w).  Both remaining trace kernels have closed
  divided-difference formulas in P.  Equal P first fails to determine D_T at
  order seven; an explicit pair has a rank-one Gram difference.  For selected
  switched targets, order join composes by
  P_join=P1P2+z^2N1N2 and N_join=N1P2+P1N2, equivalently P+/-zN multiply.
  This law is commutative and therefore forgets SCC order; K1 join C3 and C3
  join K1 are a sharp hostile.  No Hamiltonian, path-cover, substitution, or
  chronology consequence follows.
audit: >
  The companion exhausts all 1099 labelled tournaments through order five and
  all 16933 cut signs modulo global negation, comparing direct switched
  determinants with the algebraic formulas and relabel controls.  It exhausts
  every labelled tournament through order six for the minimal P-only boundary,
  verifies the explicit order-seven collision, and checks 1369 marked order-
  join pairs of total order at most five.  Characteristic polynomials use both
  Newton traces and independent permutation determinants.  Normal and
  optimized outputs byte-match; the source has zero assertion nodes and zero
  floating literals.
source: root/creative-synthesis-next/2026-08-03
depends_on:
  - THM-3315-tournament-cut-switching-centered-coronal-walk-compiler
related:
  - THM-1862-order-join-reduction-principle
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
  - THM-2501-switching-fourth-moment-signed-c4-and-gram-energy
script: 04-computation/tournament_switching_coronal_second_moment_join_compiler_20260803.py
output: 05-knowledge/results/tournament_switching_coronal_second_moment_join_compiler_20260803.out
script_sha256: c33ad0ee6b3333401e24ad94e59ea7d851b283af1849b1cbdc4a14c645fc5f4a
output_sha256: 6747c637bf1bc6b1ff17162b83b53699175def74b28467968b4f3f2cba578c11
semantic_sha256: 40427da1f4e5a2d18a971f25c3a81ccaac48314b81d16b2e1cc4f409142ba345
hash_basis: LF-normalized bytes
---

# THM-3322 -- tournament switching second moment, deletion Gram, and order-join law

**PROVED by exact algebra + VERIFIED-EXACT in the stated finite universes.**

Let `T` be a labelled tournament of order `n`, with row-to-column adjacency
`A`, centered matrix `C=2A-J`, and a cut sign `d~-d`.  Put

```text
P(z)=det(I-zC),       B(z)=adj(I-zC),       N_d(z)=d^T B(z)d. (1)
```

Then `C^T=-2I-C`.

## 1. Exact cut-cube second moment

Let `M(z)=tr B(z)`.  Since `d_i^2=1`, the Walsh expansion is

```text
N_d(z)=M(z)+sum_(i<j)(B_ij(z)+B_ji(z))*d_i*d_j.            (2)
```

The degree-two characters are orthogonal on the cut cube, so

```text
Cov(N_d(z),N_d(w))
 =sum_(i<j)(B_ij(z)+B_ji(z))(B_ij(w)+B_ji(w)).             (3)
```

Define

```text
K_-(z,w)=tr(B(z)B(w)),
K_+(z,w)=tr(B(z)B(w)^T),
D_T(z,w)=sum_i B_ii(z)B_ii(w).                             (4)
```

Entrywise expansion of `(4)` turns `(3)` into

```text
E[N_d(z)N_d(w)]
 =M(z)M(w)+K_-(z,w)+K_+(z,w)-2D_T(z,w).                   (5)
```

The diagonal cofactor is the centered characteristic polynomial of the
vertex-deleted tournament:

```text
B_ii(z)=det(I-zC[V\{i}]).                                  (6)
```

Thus `D_T` is literally the Gram kernel of the vertex-deletion polynomial
deck.  It is invariant under switching and relabelling.

## 2. Everything except the deletion Gram is P-determined

Jacobi's identity gives

```text
M(z)=nP(z)-zP'(z).                                         (7)
```

For `R_z=(I-zC)^(-1)`, the ordinary resolvent identity and
`C^T=-2I-C` give

```text
R_z R_w=(zR_z-wR_w)/(z-w),
R_z R_w^T=(zR_z+wR_w^T)/(z+w+2zw).                         (8)
```

Taking traces after multiplying by `P(z)P(w)` proves the polynomial
identities

```text
(z-w)K_-(z,w)=zM(z)P(w)-wP(z)M(w),                         (9)
(z+w+2zw)K_+(z,w)=zM(z)P(w)+wP(z)M(w).                    (10)
```

Consequently `(5)` is determined exactly by `P` and `D_T`; conversely, `P`
and the second moment recover `D_T`.  For the switched adjacency determinant
`Q_d=P-zN_d` from THM-3315,

```text
Cov(Q_d(z),Q_d(w))=zw*Cov(N_d(z),N_d(w)).                  (11)
```

## 3. Sharp P-only failure at order seven

For the standard bit-mask convention, the order-seven masks `4164` and
`5122` have the same centered characteristic polynomial

```text
P=1+7z+42z^2+140z^3+360z^4+576z^5+576z^6+256z^7,         (12)
```

but different score sequences

```text
(6,5,3,2,2,2,1),          (5,5,3,3,2,2,1).                (13)
```

Their deletion-Gram difference factors as

```text
D_4164-D_5122
 =128*z^4*w^4*(1+2z+4z^2)*(1+2w+4w^2).                   (14)
```

Therefore their second moments differ by `-2` times `(14)`.  Exhaustion of
all labelled tournaments through order six finds no equal-`P`, unequal-`D_T`
pair.  Hence order seven is the sharp first failure of `P` alone, though this
does not classify all order-seven collisions.

## 4. Order-join compiler

Let `S_i=T_i^(d_i)` be two selected switched targets, and let `(P_i,N_i)` be
their centered denominator and all-ones numerator.  Form the order join
`S=S_1 triangle S_2`, with every cross arc directed from `S_1` to `S_2`.
Write

```text
Q_i=P_i-zN_i=det(I-2zA(S_i)).                              (15)
```

Block upper triangularity gives `Q_S=Q_1Q_2`.  The inverse block matrix gives
the total-walk coronal identity

```text
G_S(2z)=G_1(2z)+G_2(2z)+2zG_1(2z)G_2(2z).                 (16)
```

Using `G_i=N_i/Q_i` and `P_S=Q_S+zN_S` yields

```text
N_S=N_1P_2+P_1N_2,
P_S=P_1P_2+z^2N_1N_2.                                     (17)
```

Equivalently,

```text
P_S +/- zN_S=(P_1 +/- zN_1)(P_2 +/- zN_2).                (18)
```

This is an associative two-coordinate compiler for repeated order joins of
the already selected switched targets.  It does not assert that joining first
and then applying an arbitrary combined cut obeys `(17)`.

## 5. Exact loss of SCC order

The law `(17)` is commutative although order join is not.  Indeed,

```text
K1 triangle C3: score sequence (3,1,1,1),
C3 triangle K1: score sequence (2,2,2,0),                  (19)
```

but both have

```text
P=1+4z+12z^2+16z^3+16z^4,
N=4+12z+24z^2+16z^3.                                      (20)
```

Thus `(P,N)` composes exactly while forgetting the ordered condensation of
strong components, including source-versus-sink placement.  Neither the
second moment nor the join pair restores individual cut labels, higher
cut-cube moments, Hamiltonian/path-cover data, substitution runs, or
chronology.

QED.
