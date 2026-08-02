---
id: THM-3105
title: "Weighted Jordan repair Newton attachment spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The finite
  determinant Newton support gives the exact robust criterion for a unit
  repair of a nilpotent matrix under competing weighted errors.  For a
  generic one-parameter error and Jordan type lambda_1>=...>=lambda_r, the
  lower attachment faces are s^k z^(d-Lambda_k), where Lambda_k is the sum
  of the k largest blocks.  Along s=t^a,z=t^b the generic determinant order
  is min_k[ka+(d-Lambda_k)b], so pure repair has the sharp threshold
  a>lambda_1 b.  Special directions and identical unshifted Smith barcodes
  can have different shifted spectra; endpoint minors or characteristic
  coefficient jets are load-bearing.  Conjugate local pairs turn every
  nonzero leading attachment into a positive real norm, but this theorem
  does not prove that a physical GMC direction hits a nonzero face.
source: root-weighted-jordan-attachment-2026-08-02
audit: >
  An independent hostile audit rederived the unrestricted cycle-cover lower
  bound, endpoint-minor coefficient, simultaneous generic Zariski-open
  realization, arc valuation and sharp longest-block threshold, both
  special-direction hostiles, conjugate-pair norm sign, and the precise
  THM-3101 application.  Normal, optimized, and stored exact companions
  agree byte-for-byte after LF normalization; both declared hashes and the
  documentation checker pass.  The audit retains the finite nilpotent
  shifted-determinant scope and does not infer genericity of a physical GMC
  deformation.
depends_on:
  - THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
related:
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3101-zero-primary-next-moment-unit-suspension-and-conditional-codimension-five-spectrum
script: 04-computation/gmc_weighted_jordan_repair_newton_spectrum_thm3105.py
output: 05-knowledge/results/gmc_weighted_jordan_repair_newton_spectrum_thm3105.out
script_sha256: 6fb60635ccc6e4a12006b0e81a147b8636c9a20f998ba55c3bc824a63a40f4ca
output_sha256: 513f5df4e6e590d1294e242732ab8414f6494f617894cc4af771accffdf6f9ae
hash_basis: LF-normalized bytes
---

# THM-3105 -- weighted Jordan repair Newton attachment spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2985 identifies the first ambient normal determinant and THM-2960 turns a
chosen arc into Smith-bar lifetimes.  A unit repair introduces a second,
distinguished direction.  The relevant object is no longer the entrywise
size of the old error or the unshifted Smith multiset.  It is the finite
Newton support of the shifted determinant.

This distinction exactly repairs the failed nilpotent version of THM-3101.
A length-two Jordan chain can consume two repair factors with one error
edge, so an `o(z)` error need not be small relative to the determinant's
`z^2` repair.

## 1. The exact finite repair polygon

Let `K` be a characteristic-zero field and let

```text
N,P,E_1,...,E_s in Mat_d(K),             P in GL_d(K). (1)
```

Assume that

```text
Q=P^(-1)N                                                   (2)
```

is nilpotent.  This holds, in particular, when `N` and `P` are commuting
multiplication operators, `N` is nilpotent, and `P` is a unit.  Replace each
error by `P^(-1)E_i` and define

```text
D(u_1,...,u_s,z)
 =det(Q+sum_i u_i E_i+z I)
 =sum_((alpha,q) in S) c_(alpha,q) u^alpha z^q.        (3)
```

The finite set

```text
S subset N^s x {0,...,d}                               (4)
```

is the **repair attachment spectrum**.  Nilpotence gives

```text
c_(0,d)=1,                  c_(0,q)=0 for q<d.          (5)
```

For arbitrary positive scales `u_i(C),z(C)->0`, the exact identity is

```text
det(N+sum_i u_i(C) P E_i+z(C)P)
 =det(P) D(u(C),z(C)).                                  (6)
```

Consequently the pure repair has the robust asymptotic

```text
det(...)=det(P) z(C)^d(1+o(1))                          (7)
```

whenever every nonpure support point obeys

```text
u(C)^alpha z(C)^q=o(z(C)^d)
 for (alpha,q)!=(0,d), c_(alpha,q)!=0.                 (8)
```

If several support points have the smallest weight, their displayed finite
sum is the leading face polynomial.  It may vanish on a special direction;
one then passes to the next face.  Thus `(8)` is a coefficientwise robust
sufficient condition, while the actual face polynomial is the exact
necessary datum on a specified arc.

For a monomial arc

```text
u_i=c_i t^a_i+o(t^a_i),       z=c_z t^b+o(t^b),
c_i c_z!=0,                                              (9)
```

put

```text
w(alpha,q)=sum_i alpha_i a_i+qb.                       (10)
```

If the initial face in `(3)` does not vanish at `(c_1,...,c_s,c_z)`, then

```text
ord_t det(...)=min_((alpha,q) in S) w(alpha,q).         (11)
```

This is the exact Newton/arc interface.  Over `K[[t]]`, the right side is
also the sum of the Smith exponents of the repaired matrix.  The full
support `(3)`, not an analogy with a polygon, proves `(11)`.

## 2. One generic error and the Jordan attachment law

Now use one error parameter `s`, put `P=I`, and let the nilpotent Jordan type
of `Q` over an algebraic closure be

```text
lambda_1>=lambda_2>=...>=lambda_r>=1,
sum_i lambda_i=d.                                       (12)
```

Write

```text
Lambda_0=0,
Lambda_k=sum_(i=1)^k lambda_i             (1<=k<=r),
Lambda_k=d                                (k>=r).       (13)
```

Every monomial `s^k z^q` in

```text
det(Q+sE+zI)                                             (14)
```

obeys the sharp combinatorial bound

```text
q>=d-Lambda_k.                                          (15)
```

Indeed, view a determinant term as a cycle cover.  The nonzero entries of
`Q` are the forward edges in the disjoint Jordan chains.  After deleting the
`q` vertices carrying `z`-loops, every remaining path of `Q`-edges must be
closed by an `E`-edge.  With `k` error edges, at most `k` chains can be closed,
and their total number of vertices is at most `Lambda_k`.  This proves
`(15)`.

The bound is attained generically for every `0<=k<=r`.  Let `f_i,l_i` be the
first and last vertices of block `i`.  For a set `T` of `k` blocks, the cycle
covers using every internal Jordan edge of those blocks contribute, up to a
fixed nonzero sign,

```text
det[E_(l_i,f_j)]_(i,j in T) s^k z^(d-sum_(i in T)lambda_i). (16)
```

If `lambda_k>lambda_(k+1)`, the lower coefficient at
`s^k z^(d-Lambda_k)` contains the unique top-`k` endpoint determinant.  At a
tie it is the corresponding finite sum over maximizing subsets.  These are
nonzero polynomials in the entries of `E`: setting all non-endpoint entries
to zero and choosing signed block closures gives the exact positive model

```text
det(Q+sE+zI)=product_(i=1)^r (z^lambda_i+s).            (17)
```

Therefore one Zariski-open set of error directions realizes every lower
face

```text
s^k z^(d-Lambda_k),                 0<=k<=r.            (18)
```

This is a generic-direction statement.  It does not say that every physical
or coordinate error direction hits every endpoint minor.

## 3. Sharp arc threshold and Smith attachment lifetime

Put

```text
s=t^a,                         z=t^b,        a,b>0.      (19)
```

For a generic error direction avoiding the finitely many initial-face
cancellations, `(15)--(18)` give

```text
ord_t det(Q+t^aE+t^bI)
 =min_(0<=k<=r) [ka+(d-Lambda_k)b].                    (20)
```

The `k=0` pure repair has order `db`.  Since

```text
Lambda_k<=k lambda_1,                                  (21)
```

it is the unique lowest face exactly when

```text
a>lambda_1 b.                                          (22)
```

At equality the longest-chain attachment ties it; below equality it arrives
first.  In scale notation, the robust sufficient condition is

```text
s=o(z^lambda_1).                                       (23)
```

It is sharp already for one block:

```text
det(zI+J_lambda+s E_(lambda,1))
 =z^lambda+(-1)^(lambda-1)s.                           (24)
```

Changing the sign of the closure entry reverses the attachment sign.  This
is the smallest real-residue hostile to a positivity claim.

Equation `(20)` is a Smith-bar identity only **after** the repair direction
has been inserted.  The Smith bars of the unshifted error family alone do
not determine it.  For every `A>=1`, the two arcs

```text
Q_1(t)=[0    1],             Q_2(t)=[t^A  1]
       [t^A  0]                    [t^A  0]             (25)
```

have the same special fibre `J_2` and the same Smith exponents `(0,A)`: each
has a unit entry and determinant `-t^A`.  Nevertheless,

```text
det(Q_1+zI)=z^2-t^A,
det(Q_2+zI)=z^2+t^A z-t^A.                             (26)
```

Thus even the pair “special Jordan type + unshifted Smith multiset” loses a
shifted characteristic coefficient.  One must retain `(3)`, the endpoint
jets `(16)`, or equivalently the principal exterior/characteristic
coefficients.  At the opposite special direction `E=I`,

```text
det(Q+sI+zI)=(s+z)^d,                                  (27)
```

so the long-chain faces disappear and `s=o(z)` suffices.  Genericity in
`(18)` is therefore load-bearing.

## 4. Conjugate norm positivity

The sign hostile `(24)` is real.  In a real finite algebra with no real
residue field, every local factor occurs with its complex conjugate.  If the
first nonzero attachment coefficient on one member is `c`, the conjugate
factor contributes `conjugate(c)`.  Their product is

```text
|c|^2>0.                                                (28)
```

Nilpotent thickness does not change this norm identity.  Hence on a
conjugate-paired physical resultant, finding any nonzero leading attachment
face settles its sign positively.  What remains difficult is proving that a
single face is nonzero uniformly over every physical composition/gap
direction.  The present theorem separates that nonvanishing debt from the
sign debt.

## 5. Application and exact boundary for THM-3101

On THM-3101's scheme-zero primary factor, `Q=0`.  Its Jordan type is

```text
(1,1,...,1),                                            (29)
```

so `(23)` reduces exactly to the proved entrywise condition `s=o(z)`.  This
recovers the repaired proof without loss.

If the old vanishing moment is nilpotent but nonzero, then
`lambda_1>=2`.  The slowest prior defining-row error and the physical repair
have exponential bases

```text
alpha_m=((m-1)/(m+1))^(m-1),
rho_m=(m/(m+1))^m,
alpha_m>rho_m^2.                                        (30)
```

Thus the available scalar estimate `s=o(z)` cannot cross even the generic
length-two attachment threshold.  It may still succeed in a special
direction such as `(27)`, or a different mixed face may be nonzero and hence
positive by `(28)`.  The next physical theorem must compute the endpoint or
characteristic jets of the actual lower moment deformations on every
zero-primary component.  Neither nilpotence, Smith bars alone, nor one norm
order can substitute for that calculation.

This theorem does not prove the physical direction generic, the endpoint
minors nonzero, arbitrary-gap uniformity on the nilpotent stratum, SFC(4),
or any new Gaussian Moment Conjecture case.  It supplies the exact finite
attachment object and the cheapest decisive tests for that frontier.

## 6. Exact evidence

The dependency-free companion verifies:

1. every Jordan partition through dimension eight and all `217` generic
   sparse attachment faces;
2. every combinatorially possible dense-error support point through
   dimension seven, with no monomial below `(15)`;
3. `4,224` monomial arcs and the exact valuation formula `(20)`;
4. the longest-block dominance inequality through dimension twenty;
5. eight same-special-fibre, same-Smith-bar hostiles `(25)--(26)`; and
6. both signs of the one-block closure wall through dimension twelve.

Run

```text
python 04-computation/gmc_weighted_jordan_repair_newton_spectrum_thm3105.py
python -O 04-computation/gmc_weighted_jordan_repair_newton_spectrum_thm3105.py
```

Both modes must equal the stored transcript after LF normalization.

**QED.**
