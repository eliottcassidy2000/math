# Primitive lattice directions carry a projective Hermite precision law

**Status: PROVED in the stated scope; exact audit PASS; independently proof-audited by the root
agent; see the [independent root audit](creative_20260906_root_audit.md),
including 60 literal SymPy Smith checks.** The proof extends the integer-node precision mechanism to
homogeneous polynomial observations on primitive lattice directions, including
points at infinity. This result has no theorem ID and makes no external
novelty or priority claim. The complete partition remains nonmetric in general.

## Inheritance and research choice

The closest proved mechanism is
[THM-4435 / four-node metric blindness and universal Hermite precision](../../01-canon/theorems/THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md):
construct inverse columns and prove that their denominators are attained.
[THM-4439 / all-node two-jet metric precision by terminal clusters](../../01-canon/theorems/THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md)
is the metric input. The canonical hostile is THM-4435's pair with identical
metric and largest exponent but unequal intermediate Smith factors. The
corrected near miss is MISTAKE-547: a change of observer can expose cancellation
that the metric forgot. The least-used relevant sidecar is the **integral
coefficient lattice under a change of projective chart**.

The initial concept board was: complete local Hasse observations; primitive
lattice directions; homogeneous binary forms; projective change of chart;
terminal residue clusters; and exact inverse denominators. Searches of current
canon and result prose for projective/Hermite/confluent/Smith combinations did
not find this observer. Existing tensor-grid transport and normalized reciprocal
jet updates were recovered and excluded as already developed routes. The cheap
hostile was to rescale one lattice representative without retaining primitivity;
it changes the Smith module, as Section 6 records. Retaining primitive vectors
and their unimodular tangents repairs that failed quotient.

## 1. The integral observer and exact denominator theorem

Let `n>=1`. Let `v_i=(a_i,b_i)` be primitive integer vectors, with distinct
rational directions: `det(v_i,v_j)!=0` whenever `i!=j`. Choose any integer
vector `w_i` satisfying `det(v_i,w_i)=1`. Such a choice exists by Bezout.
On the free rank `2n` lattice of homogeneous binary forms of degree `2n-1`,
observe

```text
E(F) = ([T^0] F(v_i+T w_i), [T^1] F(v_i+T w_i))_i.
```

These are the complete first two Hasse jets in the declared local coordinate.
Define, with empty products equal to one and empty sums zero,

```text
Delta_ij = det(v_i,v_j),
d_i      = product_(j!=i) Delta_ij,
c_i      = 2 sum_(j!=i) det(w_i,v_j)
                    product_(k!=i,j) Delta_ik.
```

Then `E` is nonsingular and its largest integer Smith factor is exactly

```text
s_max(E) = lcm_i [ |d_i|^3 / gcd(|d_i|,c_i) ].        (1)
```

For each prime `p`, the sharp uniform extra observation precision is

```text
L_p = max_i { 2v_p(d_i), 3v_p(d_i)-v_p(c_i) },        (2)
```

where `v_p(0)=infinity`. For every `N>=1`, all observations modulo
`p^(N+L_p)` determine all homogeneous coefficients modulo `p^N`; if `L_p>0`,
one fewer digit fails uniformly. For `n=1`, `s_max=1` and every `L_p=0`.

### Proof and denominator attainment

For a variable vector `Z=(X,Y)`, put

```text
t_i(Z)=det(v_i,Z),     s_i(Z)=det(Z,w_i),
Q_i(Z)=product_(j!=i) det(Z,v_j)^2.
```

The pair `(s_i,t_i)` is an integral unimodular linear coordinate system:
`s_i(v_i+T w_i)=1` and `t_i(v_i+T w_i)=T`. Also

```text
Q_i(v_i)=d_i^2,
[T] Q_i(v_i+T w_i)=d_i*c_i.
```

The two cardinal forms of degree `2n-1` are

```text
Phi_(i,1) = t_i Q_i/d_i^2,
Phi_(i,0) = (d_i s_i-c_i t_i) Q_i/d_i^3.             (3)
```

At every other direction `v_j`, the factor `det(Z,v_j)^2` kills both
observations. At `v_i`, direct substitution in `(3)` gives the prescribed
value/derivative pair `(0,1)` or `(1,0)`. Thus `(3)` is the complete inverse
matrix, proving nonsingularity.

Each linear form `det(Z,v_j)` is primitive because `v_j` is primitive.
Gauss's elementary content lemma therefore makes `Q_i` primitive. The form
`t_i` is also primitive. Hence the exact least coefficient denominator of
`Phi_(i,1)` is `d_i^2`. Since `(s_i,t_i)` is an integral unimodular pair,
the coefficient gcd of `d_i s_i-c_i t_i` is exactly `gcd(d_i,c_i)`.
Multiplication by the primitive `Q_i` preserves that content. Therefore the
exact denominator of `Phi_(i,0)` is `|d_i|^3/gcd(|d_i|,c_i)`. This quantity
is a multiple of `d_i^2`, so the lcm of these value-column denominators is
the least denominator of the entire inverse.

For any nonsingular integer matrix, its least inverse denominator is its
largest Smith factor: apply the integral unimodular Smith changes to the
inverse diagonal matrix. This proves `(1)` and `(2)`, including sharpness
for actual integral inputs in Smith coordinates. It is not merely a bound
on the determinant or on arbitrary rational input noise.

## 2. Projective covariance preserves the whole Smith module

The observer is independent of the auxiliary integer tangents up to integral
unimodular target changes. Indeed any other tangent is `w_i+k_i v_i`, and
homogeneity of degree `D=2n-1` gives

```text
V_i -> V_i,             D_i -> D_i+D*k_i V_i.        (4)
```

Correspondingly `c_i` changes to `c_i+2(n-1)k_i d_i`, leaving the gcd in `(1)`
unchanged. Changing the sign of a primitive representative gives unit row
changes after the tangent is adjusted. Thus rational directions, represented
primitively, are the intrinsic data.

Every `g in GL_2(Z)` acts integrally and unimodularly on homogeneous forms by
substitution: its inverse also has integer entries. If `det(g)=1`, take
`v_i'=g v_i`, `w_i'=g w_i`. The transformed observations equal observations
of the substituted form. If `det(g)=-1`, take `w_i'=-g w_i`, adding a sign to
the derivative row. Thus **the full integer Smith module**, not only `(1)`,
is invariant under this projective action.

The same proof works over `Z_p` with `g in GL_2(Z_p)`, any local primitive
representatives, and unit rescalings of those representatives. The latter
and local-coordinate changes act by triangular matrices with unit diagonal
on each complete jet block. This also proves full local Smith covariance
for arbitrary positive multiplicities on homogeneous degree `sum m_i-1`:
under `w -> w+k v`, substitute

```text
F(v+T(w+k v)) = (1+kT)^(sum m_i-1)
               F(v+(T/(1+kT))w).
```

Truncation through each complete order `m_i-1` gives a unitriangular block.
No higher-jet metric-only precision theorem is inferred from covariance.

## 3. A projective determinant and residue-class splitting

The complete two-jet determinant satisfies

```text
|det E| = product_(i<j) |Delta_ij|^4.                (5)
```

One direct proof uses the ordinary confluent Vandermonde determinant on the
open chart `b_i!=0`. Write `x_i=a_i/b_i`, `D=2n-1`, and
`P(X)=F(X,1)`. The value row is `b_i^D P(x_i)` and the derivative row is

```text
D*(w_i)_2*b_i^(D-1) P(x_i) - b_i^(D-2) P'(x_i),
```

because `det(v_i,w_i)=1`. Their two-by-two row change has determinant
`-b_i^(2D-2)`. Substituting
`x_i-x_j=Delta_ij/(b_i*b_j)` into the confluent Vandermonde cancels all
powers of `b_i`, giving `(5)` up to its irrelevant sign. To include directions
at infinity, first apply the integer determinant-one shear `(X,Y) ->
(X,kX+Y)`, choosing an integer `k` outside the finitely many forbidden values
for which `k*a_i+b_i=0`. Section 2 preserves the observer determinant up to a
unit, and the shear preserves every bracket. All new second coordinates are
nonzero, so the already proved affine calculation applies. No limiting
argument or chart with unit denominators is needed for this determinant
identity over the rationals.

Fix a prime `p` and partition the directions by their reductions in
`P^1(F_p)`. Directions in different classes have unit bracket. The full
`p`-Smith module is the direct sum of the complete two-jet modules of those
classes, each on its own homogeneous degree `2|C|-1` source lattice.

Here is an integral proof of that splitting. For a class `C`, multiply its
source form by

```text
R_C(Z)=product_(j not in C) det(Z,v_j)^2.
```

This gives a homomorphism `B` from the direct sum of the class source
lattices to the full source lattice. Its observations vanish outside `C`.
Within `C`, multiplication by `R_C` acts on the two jets by a triangular
matrix whose two diagonal entries are `R_C(v_i)`, both units. Thus `E B`
is a block sum of the class observers, up to target unit matrices. Formula
`(5)` says that the valuation of `det E` is the sum of the valuations of
the class determinants, since every cross-class bracket is a unit.
Taking determinants of `E B` therefore gives `v_p(det B)=0`. Hence `B` is
unimodular over `Z_p`, proving the complete module splitting.

This proof works even when all `p+1` projective residue classes occur. One
cannot then put every direction into a single integral affine chart; the
splitting, rather than an illegal chart choice, is the missing operation.

## 4. The terminal-tree precision law for lattice directions

For each prime `p`, define the projective pair depths and mass

```text
h_ij = v_p(Delta_ij),        S_i = sum_(j!=i) h_ij.
```

These are nonnegative and define the usual rooted projective residue tree.
A terminal cluster `C` of depth `f` is a nonsingleton ball intersection
whose distinct points all have pair depth exactly `f`. Its children are
leaves and `S_i` is constant on `C`, written `S_C`. At positive depth,
`2<=|C|<=p`. A depth-zero terminal root may instead have `p+1` leaves.
Then

```text
L_p = max_(terminal C) [2S_C+max(0,f_C-[|C|=p])].    (6)
```

Use zero for a singleton observer. Formula `(6)` says that the sharp
coefficient-recovery precision depends only on the bracket valuation tree
for **all primitive lattice directions**, including infinity.

To prove it, use Section 3 to split the projective residue classes. In each
class choose `g in GL_2(Z_p)` taking its common reduction to the finite chart,
and normalize its representatives to `(x_i,1)` with `x_i in Z_p`. All
normalizing denominators are units; Section 2 preserves the full local Smith
module, and bracket valuations become `v_p(x_i-x_j)`. Apply
[THM-4439 / terminal-cluster precision](../../01-canon/theorems/THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md)
to each affine class and take the maximum.

The inherited theorem is stated for integer nodes; its extension here to
`Z_p` nodes needs no extra completeness assumption. Approximate finitely many
nodes by integers closely enough to preserve every pair depth and all matrix
entries modulo `p^(T+1)`, where `T` is the known determinant valuation.
A full-rank integral matrix with determinant valuation `T` has its Smith
partition determined modulo `p^(T+1)`, since every determinantal divisor has
valuation at most `T`. The approximation transfers both sides exactly.

Cross-class terms contribute zero to every `S_i`. If every class is a
singleton, `(5)` gives `L_p=0`, and the depth-zero terminal root in `(6)`
also contributes zero, including its possible `p+1` children. This proves
all boundary cases.

For an explicit complete-partition application, fix a prime `p` and `e>=1`
and take the `p+2` directions

```text
(1,0), (1,p^e), and (a,1) for 0<=a<p.
```

Every one of the `p+1` projective residue classes occurs. Only the first two
directions collide, at depth `e`. Residue splitting and the inherited two-node
law give the **full** `p`-Smith exponent list

```text
p odd:  0 repeated 2p+2 times, e, 3e;
p=2:    0 repeated 6 times, e+1, 3e-1.              (7)
```

In particular the sharp loss is `3e` at odd primes and `3e-1` at two. This
is a concrete family whose integral precision cannot be justified by placing
all directions in one affine chart with unit denominators.

## 5. The actual connection and its limitations

The source is complete affine Hermite observation on integer polynomial
coefficients. Homogenization `P(X) -> Y^(2n-1)P(X/Y)` and local coordinates
map it to the homogeneous observer on primitive rational directions. The
map retains the integral coefficient lattice precisely when the chart
changes and representative normalizations are units at the prime under
study. It preserves the entire Smith module under those changes. The
projective quotient forgets representative content; primitive normalization
and a unimodular tangent restore exactly that coordinate.

The resulting transfer is concrete: the exact worst precision for polynomial
jets on primitive lattice directions is computed from pairwise determinant
valuations alone, using `(6)`. Formula `(1)` also gives its exact global
integer value without constructing the full observer matrix. A common
projective chart may fail; the residue-block multiplication map in Section 3
supplies the missing integral gluing argument.

This enlarges the domain of the existing precision theorem. It does not
resolve LRC(14), arbitrary-jet metric precision, or full metric-only Smith
classification. The latter is still false because all affine hostile
examples embed as primitive vectors `(x_i,1)`.

## 6. Hostiles, exact checks, and next question

**Nonprimitive quotient failure.** Cubic forms observed at `(1,0),(0,1)` with
unimodular tangents have Smith factors `(1,1,1,1)`. Replacing only `(1,0)` by
`(2,0)` and retaining transverse vector `(0,1)` keeps the rational directions
but scales its value and derivative rows by `8` and `4`. The dyadic exponents
become `(0,0,2,3)`. No integer unimodular tangent exists at `(2,0)`. The first
failed implication is that projective equality preserves an arbitrary
integral representative's observer. The repair is primitivity; for general
vectors one must retain content and its derivative-order weights.

**Full-partition boundary.** THM-4435's affine vectors from
`8*(0,1,2,5)` and `8*(0,1,3,4)` still give respective dyadic lists

```text
(0,0,4,7,12,16,19,26),       (0,0,4,7,12,17,18,26).
```

They have the same bracket metric under the inherited relabeling and the
same sharp loss `26`. Thus `(6)` does not imply that every Smith factor is
metric-only, even after projective transport.

The companion audit checks exact cardinal columns against an independently
constructed Gauss-Jordan inverse; exact full Smith lists by rational DVR
elimination; every determinant against `(5)`; the largest exponent against
both `(1)` and `(6)`; and all minors in two small controls. It tests infinity,
singletons, every projective residue class at four primes, deep infinity
clusters, random signed primitive directions, tangent shifts, determinant
one and determinant minus one projective maps, and both hostiles above.
The explicit finite universe is printed in the output. These finite checks
support the analytic proof; no sampling claim is promoted to a universal
quantifier.

Reproduce from the repository root:

```text
python3 -B 04-computation/creative_20260906_smith_bridge.py
python3 -B -O 04-computation/creative_20260906_smith_bridge.py
```

Source: [creative_20260906_smith_bridge.py](../../04-computation/creative_20260906_smith_bridge.py).
Matching output: [creative_20260906_smith_bridge.out](creative_20260906_smith_bridge.out).
The [independent root audit](creative_20260906_root_audit.md) also compares
60 fresh symbolic-jet matrices with full integer Smith forms using SymPy;
all exact denominators and bracket determinants agree.
Normal and optimized outputs agree byte for byte. The completed audit has 70
direction configurations, 280 inverse/Smith/tree prime comparisons, 360 full
projective/tangent Smith symmetry rows, 992 nonempty all-minor determinants,
and 69,828 exact gates. The script uses no floating point or Python assertion
statements.

```text
source_sha256:
1f13a63bf47c25d57427b89e4072ce2e31b093f0c3ab6a6be2590db26dbb1509
output_sha256:
d99bddb5475833d62b2bf14900222554d231cb519c0a66341d3a76930b8c37fc
semantic_sha256:
1238cff122b60e35b0dfb12a91029ca55215a04c5874dd80820f96e2ad4db9ee
```

**OPEN next question.** For arbitrary multiplicities on primitive directions,
projective covariance and residue splitting still hold. Can the exceptional
intermediate minor ideals already isolated in the affine three-node banks
be written intrinsically in determinant brackets, so that their unit sidecar
transports across projective charts without rebuilding all minors? The
present two-jet precision theorem supplies the lawful coordinate operation,
not that higher-jet invariant classification.
