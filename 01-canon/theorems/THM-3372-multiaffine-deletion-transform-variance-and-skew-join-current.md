---
id: THM-3372
title: "Multiaffine deletion transform, variance certificate, and skew join current"
status: >
  PROVED + VERIFIED-EXACT.  The whole Hamiltonian deletion deck has a
  multiaffine OCF transform whose monomials are uncovered vertex sets.  Its
  diagonal D has nonnegative coefficients, all higher deletion sums as Taylor
  derivatives, and exact ordered-join product.  The first two derivatives at
  one give a nonnegative variance which vanishes exactly on transitive
  tournaments.  Contracting the marked first deletions against A-A^T gives a
  scalar skew current xi with a closed ordered-join law and no extra sidecar;
  it distinguishes K1 ordered-join C3 from the reverse although D does not.
  Neither D nor xi is asserted injective or chronological.
source: root/repository-archaeology-2026-08-14
depends_on:
  - THM-002-ocf
  - THM-026-deletion-ratio-formula
related:
  - THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect
  - THM-506-permanental-companion
  - THM-3166-falling-factorial-order-join-path-colour-transform
  - THM-3369-skew-deletion-response-and-ordered-join-orientation-current
script: 04-computation/tournament_multiaffine_deletion_transform_skew_current_thm3372.py
output: 05-knowledge/results/tournament_multiaffine_deletion_transform_skew_current_thm3372.out
script_sha256: 8cbd74312c199af4d86fb3ae13823fa5f9e99776ce4ebc8c425845ceb32e439f
output_sha256: b800a2b6bdcb25d6c75ebfe9e921cb08533695fad30850260cdcf2c6687a1d96
semantic_sha256: c66025a0fa9e9e36bc8f9e2576f95998ff6c3e6035eaa871f97c8997dbea7688
hash_basis: LF-normalized bytes
---

# THM-3372 -- the deletion deck has a multiplicative exterior response

**PROVED + VERIFIED-EXACT.**

THM-026 records only the first deletion sum.  The full deck has a natural
transform in which every derivative is another deletion level, ordered join
is multiplication, and the first two levels carry a sharp variance.  A skew
contraction then restores the source/sink bit that the commutative transform
must forget.

## 1. Multiaffine deletion transform

Let `T` be a tournament on vertex set `V`, `|V|=n`.  Write `H(empty)=1`,
consistent with the OCF and ordered-join multiplication.  Define

```text
calD_T((y_v))
  = sum_(X subset V) H(T-X) product_(v in X)(y_v-1).     (1)
```

If `S` is an independent set in the odd-cycle conflict graph `Omega(T)`, let
`U(S)` be the union of the vertices of its cycles.  Then

```text
calD_T((y_v))
  = sum_(S in Ind Omega(T)) 2^|S| product_(v notin U(S)) y_v.  (2)
```

*Proof.*  Apply the OCF (THM-002) to every `H(T-X)` in `(1)` and interchange
the two finite sums.  A fixed packing `S` occurs exactly when
`X subset V\U(S)`, and its inner contribution is

```text
sum_(X subset V\U(S)) product_(v in X)(y_v-1)
  = product_(v notin U(S))(1+(y_v-1)).                 (3)
```

This is `(2)`.  In particular, `calD_T` retains the uncovered vertex set of a
packing.  It does not retain how the covered union decomposes into cycles. ∎

On the diagonal put

```text
D_T(y)=calD_T(y,...,y)
      =sum_S 2^|S| y^(n-|U(S)|).                       (4)
```

Thus `D_T` is monic with nonnegative integer coefficients,
`D_T(1)=H(T)`, and

```text
ord_(y=0) D_T
  =n-max_S |U(S)|.                                     (5)
```

The constant term is the weighted number of full-cover odd-cycle packings.
In the THM-505/506 master packing polynomial `Phi`, formula `(4)` is the
odd-cycle specialization

```text
D_T(y)=y^n Phi_T(x_ell=2y^(-ell) for odd ell,
                 x_ell=0 for even ell).                (6)
```

The apparent Laurent powers cancel termwise.  The new coordinate here is the
deletion/Taylor interpretation, not a claim that `Phi` was absent.

## 2. Every deletion level is a derivative

For `v in V`, differentiating `(2)` with respect to `y_v` retains precisely
the odd-cycle packings avoiding `v`.  Hence

```text
partial_(y_v) calD_T=calD_(T-v).                       (7)
```

Repeated diagonal differentiation sums over ordered distinct vertices, so for
every `0<=k<=n`,

```text
D_T^(k)(y)=k! sum_(X subset V, |X|=k) D_(T-X)(y).      (8)
```

At `y=1`, this recovers the complete symmetric deletion deck:

```text
D_T^(k)(1)=k! sum_(|X|=k) H(T-X).                      (9)
```

In particular, THM-026 is exactly the logarithmic derivative
`D_T'(1)/D_T(1)` together with the OCF expectation in `(4)`.

## 3. Ordered join is multiplication

Every directed cycle of `X triangleright Y` lies inside one factor, since all
cross arcs point from `X` to `Y`.  Its odd-cycle packings therefore split
uniquely between the two factors.  Formula `(4)` gives

```text
D_(X triangleright Y)(y)=D_X(y)D_Y(y).                (10)
```

So `D` is a commutative operation interface.  As expected, `(10)` alone
forgets factor order.

## 4. A sharp variance certificate for transitivity

Set

```text
B_0=H(T),
B_1=sum_v H(T-v),
B_2=sum_(u<v) H(T-u-v).                                (11)
```

Give each odd-cycle packing probability `2^|S|/H(T)`, and put
`W(S)=n-|U(S)|`.  Equations `(8)--(9)` say

```text
B_1=B_0 E[W],                 2B_2=B_0 E[W(W-1)].      (12)
```

Consequently

```text
Delta_T:=B_0(2B_2+B_1)-B_1^2
        =B_0^2 Var(W)>=0.                              (13)
```

Moreover,

```text
Delta_T=0  iff  T is transitive.                       (14)
```

Indeed, a transitive tournament has only the empty packing, so `W=n` is
constant.  Every nontransitive tournament contains a directed triangle.  Its
empty packing has uncovered count `n`, while the singleton triangle packing
has uncovered count `n-3`; both have positive weight.  Thus the variance is
strictly positive.  The smallest hostile is `C3`, where `Delta=18`, while the
transitive `T3` has `Delta=0`.

## 5. Exteriorizing the marked deletion response

Let `A_T` be the `0/1` adjacency matrix and

```text
K_T=A_T-A_T^T.
```

Write `d_T(z)=(D_(T-v)(z))_(v in V)` as a row.  Define the scalar bivariate
current

```text
xi_T(z,w)=d_T(z) K_T d_T(w)^T
         =sum_(u,v) K_(u,v)D_(T-u)(z)D_(T-v)(w).       (15)
```

It is invariant under relabelling and satisfies

```text
xi_T(w,z)=-xi_T(z,w),          xi_(T^op)=-xi_T.        (16)
```

The first marked response is already closed by `(8)`:

```text
sum_v D_(T-v)(z)=D_T'(z).                              (17)
```

Let `[fg](z)` mean `f(z)g(z)`.  Splitting the skew adjacency of an ordered
join into its two diagonal blocks and cross blocks `+J,-J` gives

```text
xi_(X triangleright Y)(z,w)
 =D_Y(z)D_Y(w)xi_X(z,w)+D_X(z)D_X(w)xi_Y(z,w)
 +[D_Y D_X'](z)[D_X D_Y'](w)
 -[D_X D_Y'](z)[D_Y D_X'](w).                         (18)
```

This is an exact noncommutative operation interface `(D,xi)` with no
independent first-deletion sidecar.  The proof uses only
`D_((X-x) triangleright Y)=D_(X-x)D_Y`, its `Y` analogue, and summation of the
all-ones cross blocks via `(17)`.

## 6. Sharp order hostile and information loss

For the singleton `K1` and directed triangle `C3`,

```text
D_(K1 triangleright C3)(y)
 =D_(C3 triangleright K1)(y)=y^4+2y,                  (19)

xi_(K1 triangleright C3)(z,w)=6(w^3-z^3),
xi_(C3 triangleright K1)(z,w)=6(z^3-w^3).             (20)
```

Thus the commutative deck transform sees the same three order-four types as
the OCF covered-vertex grading, while `xi(0,1)=+6` versus `-6` restores the
source/sink bit.  This is complementary to THM-3369: that theorem
exteriorizes a two-coordinate spectral response, whereas `(15)` exteriorizes
the entire scalar Hamiltonian deletion transform.

There is also a structural kernel.  If `T` is self-converse, relabelling
invariance and `(16)` give `xi_T=-xi_T`, hence `xi_T=0`.  In the pair-order bit
convention of the companion, the nonisomorphic order-five masks `8` and `10`
are both self-converse and have

```text
D_T(y)=y^5+6y^2+2,                  xi_T(z,w)=0.        (20a)
```

Their multisets of marked deletion polynomials differ: mask `8` has two
copies of `y^4` and three of `y^4+4y`, while mask `10` has one copy of `y^4`,
two of `y^4+2y`, and two of `y^4+4y`.  Thus even `(D,xi)` loses information
which the multiaffine or marked deck retains.  Exhaustion shows this is the
first cross-isomorphism collision through order five.

The typed contract is

```text
source:  multiaffine Hamiltonian deletion deck + intrinsic arc signs
map:     calD -> diagonal D; then d K d^T
target:  exact ordered-join composition and a sharp order witness
kept:    weighted uncovered-set data; one antisymmetric deletion channel
lost:    cycle decomposition inside U(S), kernel of b -> d(z)b d(w)^T,
         vertex ownership, individual arcs, full SCC order, chronology
sidecar: none beyond D, since sum_v D_(T-v)=D'
test:    K1 triangleright C3 versus C3 triangleright K1.              (21)
```

No injectivity, reconstruction, Hamiltonian optimization, switching,
isomorphism, LRC, dynamical, or chronological consequence is claimed.

## 7. Exact audit

The standard-library-only integer companion checks:

- all `1,099` labelled tournaments of orders one through five against the
  literal multiaffine deletion expansion and an independently enumerated OCF
  packing expansion;
- all `6,504` derivative/deletion identities `(8)`;
- all `1,099` variance cells, with exactly the `153=sum_(n=1)^5 n!`
  transitive labelled tournaments at equality;
- converse skewness and nonzero-current counts `0,0,0,16,320` by order;
- injectivity of `(D,xi)` on isomorphism classes through order four and the
  explicit first loss `(20a)` at order five;
- all `2,553` ordered labelled factor pairs of total order at most six for
  `(10)` and `(18)`; and
- the literal `K1/C3` formulas `(19)--(20)`.

Normal and optimized modes byte-match the stored thirteen-line transcript.
Reproduce with

```text
python3 04-computation/tournament_multiaffine_deletion_transform_skew_current_thm3372.py
python3 -O 04-computation/tournament_multiaffine_deletion_transform_skew_current_thm3372.py
```

All hashes are pinned in the frontmatter.

**QED.**
