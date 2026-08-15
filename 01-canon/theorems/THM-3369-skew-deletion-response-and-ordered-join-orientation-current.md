---
id: THM-3369
title: "Skew deletion response and ordered-join orientation current"
status: >
  PROVED by exact algebra + VERIFIED-EXACT in the declared finite universes.
  Contracting the marked two-coordinate deletion responses of THM-3324
  against the intrinsic skew adjacency K=A-A^T gives a relabel-invariant
  bivariate 2-by-2 current Omega.  Its first response is not a new sidecar:
  q_T=sum_v r_Tv=n*s_T-z*s_T'.  Under ordered join, Omega is the sum of two
  familiar H-conjugates and one antisymmetric cross-factor outer product.
  The resulting (s,Omega) interface is exact and noncommutative, with a
  closed ordered-pair formula for arbitrary repeated joins.  It distinguishes
  K1 ordered-join C3 from the reverse, the sharp SCC-order hostile invisible
  to the full symmetric (s,Gamma) interface.  This recovers one orientation
  bit, not the full SCC order, vertex ownership, chronology, or isomorphism.
source: root/repository-archaeology-2026-08-14
depends_on:
  - THM-3324-tournament-deletion-response-gram-ordered-join-compiler
related:
  - THM-3315-tournament-cut-switching-centered-coronal-walk-compiler
  - THM-3322-tournament-switching-second-moment-deletion-gram-and-order-join-law
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
script: 04-computation/tournament_skew_deletion_response_order_join_thm3369.py
output: 05-knowledge/results/tournament_skew_deletion_response_order_join_thm3369.out
script_sha256: a8472c3189542e172a48802b2663cee02c434e2e1089248e7a19a3c59813e6c0
output_sha256: bf009893390dd865c751ad8713c3c22f16eaede24b2092146b87fd8bf9ee83b6
semantic_sha256: cd70e45bc435acb281084bdee703dcf36e7b53b6877003b8260a5ff40671f12e
hash_basis: LF-normalized bytes
---

# THM-3369 -- skew deletion response and ordered-join orientation current

**PROVED by exact algebra + VERIFIED-EXACT in the stated finite universes.**

THM-3324 closes the second tensor moment of marked deletion responses under
ordered join.  Its Gram law is necessarily commutative and therefore forgets
which strong component came first.  The missing operation column is not a
higher symmetric moment.  Retaining the intrinsic skew relation before
contracting the same response deck produces an exact orientation current.

## 1. The skew deletion response

Let `T` be a tournament on `n` vertices, with row-to-column adjacency `A`.
Keep THM-3324's notation

```text
C_T=2A-J,
P_T(z)=det(I-zC_T),
N_T(z)=1^T adj(I-zC_T)1,
s_T(z)=(P_T(z),zN_T(z))^T,
r_Tv(z)=s_(T-v)(z).                                      (1)
```

Put

```text
K_T=A-A^T=C_T+I.                                         (2)
```

Thus `K_T` has zero diagonal and its off-diagonal entries are the intrinsic
arc signs.  Define the bivariate matrix polynomial

```text
Omega_T(z,w)
 =sum_(u,v in V(T)) (K_T)_(u,v) r_Tu(z) r_Tv(w)^T.       (3)
```

Equivalently, if `R_T(z)` is the `2`-by-`n` matrix with columns `r_Tv(z)`,
then

```text
Omega_T(z,w)=R_T(z) K_T R_T(w)^T.                        (4)
```

Formula `(4)` makes three immediate properties exact.

1. Simultaneous relabelling conjugates `K_T` and permutes the columns of
   `R_T`, so `Omega_T` is a labelled-isomorphism invariant.
2. Skewness gives

   ```text
   Omega_T(w,z)^T=-Omega_T(z,w).                         (5)
   ```

   In particular `Omega_T(z,z)` is a skew `2`-by-`2` matrix and hence one
   scalar exterior channel.
3. Reversing every arc transposes every centered deletion matrix, preserving
   `s_T` and the deletion responses while negating `K_T`.  Therefore

   ```text
   Omega_(T^op)=-Omega_T.                                (6)
   ```

This is an antisymmetric relation contraction, not a tournament manufactured
from ties and not a new time parameter.

## 2. The first deletion response is an Euler derivative

Set

```text
q_T(z)=sum_v r_Tv(z).                                    (7)
```

No new sidecar is needed to compose `q_T`.  Coordinatewise,

```text
q_T(z)=n s_T(z)-z s_T'(z).                               (8)
```

For the first coordinate, introduce a scalar `t`:

```text
det(tI-zC_T)=t^n P_T(z/t).                               (9)
```

Differentiating at `t=1` sums the principal cofactors and gives

```text
sum_v P_(T-v)(z)=nP_T(z)-zP_T'(z).                       (10)
```

For the second coordinate, use the rank-one deformation

```text
F(t,u)=det(tI-zC_T-uJ)
      =t^n P_T(z/t)-u t^(n-1)N_T(z/t).                  (11)
```

The mixed coefficient which chooses one `J` entry and one remaining scalar
diagonal entry has the cofactor sign

```text
partial_t partial_u F(1,0)=-sum_v N_(T-v)(z)
  =-((n-1)N_T(z)-zN_T'(z)).                            (11a)
```

The two minus signs therefore cancel, giving

```text
sum_v N_(T-v)(z)=(n-1)N_T(z)-zN_T'(z).                  (12)
```

Multiplying `(12)` by `z` is the second coordinate of `(8)`, since

```text
n(zN_T)-z(zN_T)'=z((n-1)N_T-zN_T').                     (13)
```

Thus the first deletion moment is the Euler derivative of the split response.

## 3. Exact two-factor ordered-join law

For `s=(a,b)^T`, write

```text
H_s(z)=((a(z),b(z)),(b(z),a(z))).                        (14)
```

For a tournament `T`, abbreviate this as `H_T=H_(s_T)`.  THM-3324 proves

```text
s_(X triangleright Y)=s_X star s_Y,
r_(X triangleright Y),x=H_Y r_Xx,       x in X,
r_(X triangleright Y),y=H_X r_Yy,       y in Y,          (15)
```

where `(a,b) star (c,d)=(ac+bd,ad+bc)`.  The skew adjacency of the ordered
join has the block form

```text
K_(X triangleright Y)=((K_X,J),(-J,K_Y)).                (16)
```

Put

```text
a_(X|Y)(z)=H_Y(z)q_X(z),
a_(Y|X)(z)=H_X(z)q_Y(z).                                (17)
```

Substituting `(15)--(17)` into `(4)` gives the exact law

```text
Omega_(X triangleright Y)(z,w)
 =H_Y(z)Omega_X(z,w)H_Y(w)^T
  +H_X(z)Omega_Y(z,w)H_X(w)^T
  +a_(X|Y)(z)a_(Y|X)(w)^T
  -a_(Y|X)(z)a_(X|Y)(w)^T.                              (18)
```

The first line is the symmetric THM-3324 grammar.  The second line is an
exterior cross-factor current and changes sign when the two factors are
reversed.  Since `(8)` recovers each `q` from `s`, the pair

```text
(s_T,Omega_T)                                            (19)
```

is a closed exact interface for ordered join.  Its noncommutativity is not an
extra label attached by hand; it is the sign of the literal cross arcs in
`(16)`.

## 4. Arbitrary repeated joins

Let

```text
J=T_1 triangleright ... triangleright T_k.                (20)
```

Because all matrices `H_T` have the split form `(14)`, they commute.  For
each factor set

```text
G_i(z)=product_(ell != i) H_(T_ell)(z),
q_tilde_i(z)=G_i(z)q_(T_i)(z).                           (21)
```

Every vertex in factor `i` has its deletion response multiplied by `G_i`.
Every cross block from factor `i` to factor `j` has sign `+1` exactly when
`i<j`.  Therefore `(18)` iterates without recursion to

```text
Omega_J(z,w)
 =sum_i G_i(z)Omega_(T_i)(z,w)G_i(w)^T
  +sum_(i<j) [
      q_tilde_i(z)q_tilde_j(w)^T
      -q_tilde_j(z)q_tilde_i(w)^T ].                    (22)
```

Thus the full order dependence is a signed sum over factor pairs.  Repeated
identical factors have equal transformed first responses, so their pairwise
wedge vanishes; `(22)` deliberately does not invent an order distinction
where the factors themselves provide none.

## 5. The sharp source/sink hostile is detected

Let `C3` be the directed triangle.  THM-3324 proves that

```text
K1 triangleright C3,
C3 triangleright K1                                      (23)
```

have the same full symmetric `(s,Gamma)` interface, although their score
sequences are respectively

```text
(3,1,1,1),                 (2,2,2,0).                    (24)
```

For both factors the internal current vanishes.  The cross term in `(18)` is
nonzero, and reversing the factors reverses its sign:

```text
Omega_(C3 triangleright K1)=-Omega_(K1 triangleright C3). (25)
```

For example, with coordinates ordered as `(P,zN)`, the upper-right entry is

```text
[Omega_(K1 triangleright C3)]_(0,1)
 =-24 w^3 (z+1)(4z^2+2z+1).                             (26)
```

This is nonzero, so `(s,Omega)` restores the source-versus-sink placement
lost by `(s,Gamma)` on the canonical hostile.  Total order four is sharp for
this phenomenon: below it, two nonempty factors have no nontrivial strong
component, and the reverse joins are isomorphic transitive tournaments.

## 6. Preserved predicate, loss, and boundary

The connection contract is exact:

```text
source: marked deletion responses plus intrinsic arc signs
map:    R K R^T
target: ordered-join composition and one SCC-order hostile
kept:   the antisymmetric response channel
lost:   the kernel of B -> R_T(z) B R_T(w)^T, vertex ownership,
        individual arcs, and full component order
sidecar:s_T, which supplies H_T and q_T=n*s_T-z*s_T'
test:   K1 triangleright C3 versus its reverse.            (27)
```

The current can vanish even when order exists.  In particular, the pairwise
term in `(22)` vanishes for identical transformed first responses.  Therefore
`Omega` is neither an injective tournament invariant nor a reconstruction of
the SCC condensation.  Adjacency powers still count relation walks rather
than chronology.  No Hamiltonian, path-cover, substitution-run, switching,
isomorphism, LRC, or time-evolution consequence follows.

## 7. Exact audit

The exact integer companion, importing only the already audited matrix
primitives of THM-3324, checks:

- all `1,099` labelled tournaments of orders one through five for `(8)`,
  `(5)`, and converse negation `(6)`;
- `1,589` relabellings through order four;
- all `2,553` ordered labelled factor pairs of total order at most six for
  the direct response and current law `(18)`;
- all `339` ordered labelled factor triples of total order at most six for
  the nonrecursive formula `(22)`;
- the five reverse joins of total order at most three; and
- the literal `K1/C3` hostile, including an independent permutation-
  determinant calculation and the factorization `(26)`.

The source contains zero assertion nodes and zero floating literals.  Normal
and optimized modes byte-match the stored fourteen-line transcript.

Reproduce with

```text
python3 04-computation/tournament_skew_deletion_response_order_join_thm3369.py
python3 -O 04-computation/tournament_skew_deletion_response_order_join_thm3369.py
```

**QED.**
