---
id: THM-2770
title: "Tree-incidence A/D Weyl clutch and the four-vertex fan dichotomy"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every rooted tree, reduced integral incidence is unimodular and
  identifies the graceful obstruction with an A_(n-1) vertex discriminant
  clutched to a D_(n-1) edge discriminant.  At n=4, D3=A3: the star supports
  the B3 fan with 48 chambers, while the path adds one diagonal which makes
  twelve unimodular Farey splits and 60 chambers.  Nevertheless the unique
  balanced coefficient is 120 for the star and zero for the path, although
  P4 is graceful.  Fan unimodularity is therefore not a signed coefficient
  selector.  No Graceful Tree, Keller, modular-action, or LRC closure follows.
source: a4-resolvent-next-gate/tree-incidence-weyl-clutch-2026-07-28
depends_on:
  - THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge
  - THM-2765-rooted-nullstellensatz-linear-range-distinct-edge-labeling
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
related:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
script: 04-computation/tree_incidence_a_d_weyl_clutch_thm2770.py
output: 05-knowledge/results/tree_incidence_a_d_weyl_clutch_thm2770.out
script_sha256: e61c3b20b3eb8fbef3872e9840a1b6cb7f04bf3b0469f5cd46b26e4e272c654b
output_sha256: 5d8c651cb7e0b27ed60f5da9fe2fb7546199fa84e136711c09d39b9f48bc0b01
hash_basis: LF-normalized bytes
---

# THM-2770 -- the graceful polynomial is an `A/D` Weyl clutch

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The factor three in THM-2765 is not merely a count of convenient polynomial
factors.  After quotienting translation, a tree gives an integral change of
coordinates between vertex potentials and edge gradients.  Under that
change, the graceful obstruction is the product of two reflection
discriminants:

```text
vertex injectivity       = type A_(n-1),
unsigned edge separation = type D_(n-1).                (1)
```

This makes the special role of binary signs and ternary permutations exact
at four vertices, where `D3=A3`.  It also exposes a sharp failure boundary:
even a unimodular common fan does not prevent signed coefficient
cancellation.

## 1. Reduced tree incidence is an integral isomorphism

Let `T` be a tree on `V={0,...,n-1}`.  Choose a root, orient every edge away
from it, and order nonroot vertices parent before child.  For the edge whose
child is `v`, write

```text
y_v=x_v-x_parent(v).                                    (2)
```

Translation acts freely on vertex potentials, so use either the quotient

```text
L_V=Z^V / Z(1,...,1)                                    (3)
```

or its unique root-zero section.  The reduced incidence map is

```text
partial_T:L_V -> Z^E,             x |-> (y_e)_e.        (4)
```

In parent-before-child order its matrix is unit lower triangular: the row
for child `v` has `+1` in column `v` and, unless the parent is the root,
`-1` in an earlier column.  Therefore

```text
det(partial_T)=1.                                        (5)
```

The inverse is the path-sum formula

```text
x_v=sum_(e in path(root,v)) y_e.                         (6)
```

Thus `(4)` is a unimodular lattice isomorphism, not only a rational change
of variables.  Changing the root, edge orientations, or edge order changes
the edge coordinates by a signed permutation and changes the displayed
discriminants below by at most a global sign.

## 2. The exact `A_(n-1)/D_(n-1)` factorization

Fix vertex and edge orders and put

```text
Delta_A(x)=product_(u<v)(x_v-x_u),
Delta_D(y)=product_(e<f)(y_f-y_e)(y_f+y_e)
          =product_(e<f)(y_f^2-y_e^2).                  (7)
```

The THM-2765 graceful obstruction is exactly

```text
Phi_T(x)=Delta_A(x) Delta_D(partial_T x).                (8)
```

The first factor is the positive-root discriminant of `A_(n-1)`.  The
second is the positive-root discriminant of `D_(n-1)` in the unimodular edge
coordinates.  Consequently

```text
deg Phi_T
 = |A_(n-1)^+|+|D_(n-1)^+|
 = binom(n,2)+(n-1)(n-2)
 = (n-1)(3n-4)/2.                                       (9)
```

For an integral point, `Delta_A!=0` says all vertex labels are distinct.
It also forces every edge gradient to be nonzero, because parent and child
are vertices.  Conditional on that, `Delta_D!=0` says no two gradients are
equal or opposite.  Hence

```text
Phi_T(x)!=0
iff x is injective and {|y_e|:e in E} is pairwise distinct.  (10)
```

The `1+2` exponent increment in THM-2765 is therefore one earlier
`A`-root plus the two sign mirrors of one earlier `D`-root.  It is not a
hidden order-three action.

## 3. Why rank three is the exact binary/ternary coincidence

For three edge coordinates, define four linear forms by

```text
2z_0= y_1+y_2+y_3,       2z_1= y_1-y_2-y_3,
2z_2=-y_1+y_2-y_3,       2z_3=-y_1-y_2+y_3.             (11)
```

Their sum is zero.  Their six pair differences are exactly, up to sign,

```text
y_i-y_j,                 y_i+y_j        (1<=i<j<=3).    (12)
```

Thus the half-Hadamard map `(11)` identifies the `D3` root arrangement with
the `A3` root arrangement.  Its four rows are the even sign vectors of
THM-2766.  The ambient signed-edge group is

```text
W(B3)=C2^3 semidirect S3,                                (13)
```

while the product-oriented subgroup is

```text
W(D3)=V4 semidirect S3 isomorphic_to S4.                 (14)
```

This is the exact coexistence of binary sign choices and a ternary edge
permutation in the tree polynomial.  It is a finite semidirect product,
not the modular free product `C2*C3`; THM-2768 treats that distinct carrier.

## 4. The two four-vertex tree clutches

There are two unlabelled four-vertex trees.  Root the star at its centre and
write `y_i=x_i-x_0` on its three edges.  Equations `(7)--(8)` give

```text
Phi_star
 =y_1y_2y_3 product_(i<j)(y_j-y_i)^2(y_j+y_i).          (15)
```

Ignoring multiplicity, its hyperplanes are

```text
y_i=0,                         y_i-y_j=0,
y_i+y_j=0                      (1<=i<j<=3),             (16)
```

which is exactly the `B3` reflection arrangement.  The three difference
hyperplanes occur twice in `(15)`.

For the path rooted at an endpoint, let its successive gradients be
`y_1,y_2,y_3`.  The six vertex differences are

```text
y_1,y_2,y_3,y_1+y_2,y_2+y_3,y_1+y_2+y_3.               (17)
```

Therefore

```text
Phi_path
 =y_1y_2y_3 (y_1+y_2)^2 (y_2+y_3)^2 (y_1+y_3)
  (y_2-y_1)(y_3-y_1)(y_3-y_2)(y_1+y_2+y_3).             (18)
```

Its distinct support is the same `B3` arrangement plus the single diagonal

```text
H:y_1+y_2+y_3=0.                                        (19)
```

The difference is not cosmetic: it is the complete four-vertex fan
dichotomy.

## 5. Twelve exact Farey flanks

A `B3` chamber fixes the three coordinate signs and a strict order of their
absolute values.  Hence there are

```text
2^3 3! = 48                                              (20)
```

chambers.  The diagonal `(19)` meets a chamber precisely when the unique
largest absolute coordinate has sign opposite to the other two.  There are
six absolute orders and two global signs, so exactly twelve chambers are
cut.  Every cut adds one chamber, and therefore

```text
number_of_chambers(path support)=48+12=60.               (21)
```

The local split is the rank-two Farey move of THM-2056 in one higher-rank
orthant.  In a cut chamber write the ordered magnitudes as

```text
0<u=a<v=a+b<w=a+b+c,             a,b,c>0.               (22)
```

After a global sign, `(19)` is `u+v-w=0`, hence simply

```text
a=c.                                                     (23)
```

The positive `(a,b,c)` orthant splits into the cones with ray columns

```text
(e_a,e_b,e_a+e_c),             (e_c,e_b,e_a+e_c).       (24)
```

Their determinants are `+1` and `-1`.  Thus all twelve new flanks are
unimodular primitive-ray subdivisions.  This is a genuine Farey-fan
connection, but it preserves only the hyperplane support.  It forgets
factor multiplicities and coefficient signs.

## 6. The balanced coefficient separates star and path

Both four-vertex polynomials have degree twelve.  Therefore the only
full-degree monomial with every vertex exponent at most three is

```text
x_0^3x_1^3x_2^3x_3^3.                                   (25)
```

Let `P_T(y)` denote `(15)` or `(18)`.  Since `y_i=x_i-x_0`, direct
coefficient extraction gives

```text
[x_0^3x_1^3x_2^3x_3^3] Phi_T
 =-sum_(a_1+a_2+a_3=12, a_i>=3)
    [y_1^a1 y_2^a2 y_3^a3]P_T product_i binom(a_i,3).   (26)
```

The sign is common because the `x_0` excess is three.  The complete nonzero
contribution tables collapse to

```text
star:  3*(-40)+6*(+40)=120,
path:  20+40+80-40-80-20=0.                            (27)
```

Under another labelling or factor order the star coefficient may change
sign, but its absolute value remains `120`; the path coefficient remains
zero.  Thus the star reaches the average exponent floor by the ordinary
coefficient-grid theorem, while the path does not.

This is not a graceful obstruction.  The path labelling

```text
(x_0,x_1,x_2,x_3)=(0,3,1,2)                             (28)
```

has edge differences `(3,2,1)` and is graceful.  Consequently

```text
unimodular fan + an actual grid survivor
does not imply a balanced nonzero coefficient.           (29)
```

The common-refinement geometry is unsigned; the interpolation coefficient
is a signed sum.  A cancellation-aware orientation, involution, or reference
phase is therefore genuinely load-bearing in any attempt to improve
THM-2765 by balanced monomials.  Raw characteristic-two parity also misses
the successful star witness because `120=0 mod 2`.

## 7. Exact verification and scope

Run

```bash
python 04-computation/tree_incidence_a_d_weyl_clutch_thm2770.py
python -O 04-computation/tree_incidence_a_d_weyl_clutch_thm2770.py
```

The exact companion uses explicit exceptions and no truth-bearing Python
assertions.  It checks the incidence/path-sum inverse on all `5,913`
recursive parent arrays through eight vertices, the positive-root degree
identity through twelve vertices, the half-Hadamard `D3=A3` map, the
`B3`/diagonal supports, all `48` signed-order chambers, the twelve cut
chambers and their two unimodular determinants, both independent central
coefficient extractions, all six recursive four-vertex presentations, and
the graceful `P4` hostile `(28)`.  Normal and optimized runs byte-match the
stored transcript.

```text
PROVED HERE (candidate):  integral tree-incidence isomorphism;
                          exact A_(n-1)/D_(n-1) clutch;
                          graceful nonvanishing equivalence;
                          D3=A3 half-Hadamard realization;
                          star B3 and path B3+diagonal supports;
                          chamber counts 48 and 60;
                          twelve unimodular Farey splits;
                          central coefficients 120 and zero;
                          sharp P4 coefficient-method no-go.

NOT PROVED:               a balanced coefficient for larger trees;
                          an improved uniform range beyond THM-2765;
                          the Graceful Tree Conjecture;
                          a PSL2(Z) action or tournament equivalence;
                          a Keller/JC(2)/DC(2) exclusion;
                          a physical LRC(14) carrier or closure.         (30)
```

QED (candidate).
