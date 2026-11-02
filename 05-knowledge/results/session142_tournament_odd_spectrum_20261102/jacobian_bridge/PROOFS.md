# Session 142: Jacobian-to-spectrum bridge audit

**Status:** PROVED elementary bridge/obstruction lemmas; FINITE-EXACT controls;
no proof of tournament spectrum completeness. Base repository commit:
`91cc1210ed2bc94ca52686b8091c0eec2c57a441`.

## 1. A genuine integral-Jacobian obstruction from Redei parity

For `n>=1`, index the unordered pairs `e=(i,j), i<j`, and write
`a_ij=x_e`, `a_ji=1-x_e`, `a_ii=0`. Define

```
H_n(x) = sum_(pi in S_n) product_(r=1)^(n-1) a_(pi_r,pi_(r+1)).
```

At a Boolean point this is the Hamilton-path count of the encoded labeled
tournament, with bit one meaning `i -> j`. Each path uses an unordered pair
at most once, so `H_n` is multi-affine, with integer coefficients and degree
at most `n-1`.

**Proposition 1 (coefficientwise parity and determinant obstruction).**
There is a unique multi-affine `Q_n in Z[x]` such that

```
H_n = 1 + 2 Q_n.
```

Consequently every entry of `grad H_n` lies in `2 Z[x]`. Every square integer
polynomial response map containing `H_n` as one coordinate has Jacobian
determinant in `2 Z[x]`; it cannot have determinant `+1` or `-1`. With `k`
such count coordinates its determinant is divisible by `2^k`.

**Proof.** Write `H_n=sum_S c_S x^S`. Boolean Mobius inversion gives

```
c_S = sum_(A subset S) (-1)^(|S|-|A|) H_n(1_A).
```

Here coordinates outside `S` are zero. Redei's oddness theorem gives
`H_n(1_A)=1 mod 2`. For nonempty `S` the displayed sum is zero mod two.
The all-zero orientation is transitive, so `c_empty=1`. Dividing the
nonconstant coefficients by two proves the first claim. Differentiation
preserves the factor two; multilinearity of the determinant proves the rest.
QED.

The characteristic-zero distinction is consequential: multiplying the response
by `1/2` is harmless over Q or C but is not a unimodular change of the original
integer target lattice. The correct target for a spectrum lifting attempt is
`q=(h-1)/2`, with forbidden values `3,10`, not an unnormalized `h` coordinate.
This is an exact algebraic bridge, not an analogy with the Jacobian conjecture.

## 2. Continuous relaxation erases every spectral gap

**Proposition 2 (weighted image is an interval).** The image of `[0,1]^E`
under `H_n` is exactly `[min_T h(T), max_T h(T)]`.

**Proof.** Independently orient edge `e` according to Bernoulli parameter
`x_e`. A path's indicator has expectation the corresponding product in the
definition of `H_n`. Thus `H_n(x)=E[h(T)]`, a convex combination of the vertex
values, and lies between their minimum and maximum. Conversely the continuous
image of the connected cube is connected and contains the minimum and maximum
vertex values, hence contains their entire interval. QED.

There are especially small **regular** weighted preimages of both forbidden
counts. In lexicographic pair order, the following pairs differ in one bit:

| n | endpoint masks | varied edge | H(t), 0<=t<=1 | target at t=1/2 |
|---|---|---|---|---|
| 5 | 4, 260 | (2,4) | 5+4t | 7 |
| 6 | 20, 28 | (0,4) | 19+4t | 21 |

The endpoint counts are verified by independent subset-DP and direct
permutation enumeration. Multi-affineness proves the affine formulas, not
merely interpolation at tested sample points. The derivative is constantly
four, so even the extension from this parameter line to R is a globally
bijective proper affine map. The obstruction is exactly that the unique
preimage of the target is `t=1/2`, not in the admissible source `{0,1}`.
After Proposition 1's normalization the responses are `2+2t` and `9+2t`.

An even stronger guard failure survives **unimodularity**. Parameterize the
act of flipping one edge by `t`, with the original orientation at zero and
its reversal at one:

| n | endpoint masks | flipped edge | H(t) | normalized Q(t) |
|---|---|---|---|---|
| 5 | 4, 5 | (0,1) | 5-2t | 2-t |
| 6 | 20, 16 | (0,3) | 19-2t | 9-t |

In either row the normalized response is an affine automorphism of `Z`
with derivative `-1`, also a proper globally invertible map over R or C.
The respective forbidden targets `q=3,10` have the unique **integral**
preimage `t=-1`, outside the legal set `{0,1}` (and even outside `[0,1]`).
Thus an actual tournament arc template can have an integral, unimodular,
proper inverse yet fail the required target solely at the source guard.
Signed polynomial evaluations at `t=-1` are algebraic controls, not claimed
weighted tournaments or new realizations of 7 and 21.

The checks also enumerate every labeled tournament through order six
(`1+2+8+64+1024+32768=33867` inputs), confirm oddness and omission of 7 and 21
in that finite universe, and recover the coefficientwise parity by exact
Mobius transforms. This census is **not** an all-order exclusion proof.
The all-order omissions are separately recorded by
`THM-1370-h-spectrum-omits-7-21-all-n.md`; they do not establish completeness.

## 3. What really upgrades an affine response to an integral legal lift

Let `M` be an integer `r x d` matrix, and suppose a proposed construction has
response `b+M u`, where **every** `u in N^d` is known to be a legal object.
Define `L=M Z^d`, `S=M N^d`, and `C=M R_(>=0)^d`.

**Proposition 3 (lattice and effectivity gates).**

1. Rational rank determines only the rational linear span. For an integer
   target `y`, an unrestricted integral lift exists iff `y-b in L`.
2. A legal nonnegative integral lift exists iff `y-b in S`, not merely in
   `C intersect L`.
3. If `r=d` and `det M=+1 or -1`, then `S=C intersect Z^d`. Hence in this
   special square case lattice lifting is automatic but the cone guard is
   still essential.

**Proof.** The first two equivalences are their exact image definitions
(Smith normal form decides membership in `L`). In the square unimodular case
`M^-1` is integral by the adjugate formula; if `y in C intersect Z^d`, its
unique inverse is both integral and nonnegative. The reverse inclusion is
immediate. QED.

Two hostile controls separate the gates.

* `M=[2,3]` has full rational row rank and primitive image lattice `Z`.
  Its positive real image is `R_(>=0)`, but 1 is absent from `M N^2`:
  `2u+3v=1` has no nonnegative integer solution, although `(-1,1)` is an
  integer solution. Primitive maximal minors do not imply monoid saturation.
* `M=[[1,1],[1,2]]` has determinant one. The positive integer target `(1,3)`
  has the unique inverse `(-1,2)`, outside the source cone. Unimodularity alone
  does not make every positive target legal.

These are parameter-map statements, not claims that either displayed matrix
has already been realized by a tournament gadget. For tournament templates,
the legal parameter domain can be stricter than `N^d`; all additional source,
orientation, endpoint, disjointness, and gluing constraints must also be kept.

## 4. Positive weights can supply a polynomial integral inverse, not legality

**Proposition 4 (integral positive-graded inverse).** Let positive integer
weights `w_1,...,w_d` be fixed. Suppose `F_i in Z[x_1,...,x_d]` is weighted
homogeneous of weight `w_i`, and `det JF=+1 or -1`. Then `F` has an inverse
in `Z[y_1,...,y_d]`.

**Proof.** Group variables by increasing weight. In a weight-`w` output any
occurrence of a weight-`w` variable must be a linear monomial, while every
nonlinear monomial uses only smaller-weight variables. Thus

```
F_[w] = A_w x_[w] + P_w(x_[<w]).
```

The constant Jacobian at zero is block diagonal with integer blocks `A_w`.
Its determinant is their product, `+1` or `-1`, so each integer determinant
is `+1` or `-1`. Each block inverse is integral. Recursively substitute
`x_[w]=A_w^-1(y_[w]-P_w(x_[<w]))`; this yields an integer polynomial inverse.
QED.

This is a genuine sufficient replacement for local differential invertibility
under the stated grading. Even here the legal source guard remains separate:

```
F(x,y)=(x+y^2,y), weights (2,1), det JF=1,
F^-1(0,1)=(-1,1).
```

So a positive-weight, polynomially invertible, integral and proper response
can still fail to lift a nonnegative or Boolean target into a nonnegative or
Boolean source. A future tournament argument must prove that its inverse
lands in the actual construction domain, not just in affine space.

**Source correction needed.** The main theorem of
`THM-1370-elliptic-jacobian-conjecture-all-dimensions.md` survives, but its last
scope paragraph says that `F=X+H` with ordinary homogeneous `H` of degree >=2
is never positively weighted. The displayed shear has `H=(y^2,0)` and refutes
that sentence. The repaired scope is: positive weighted homogeneity is an
additional constraint, not one automatically obeyed by arbitrary
identity-plus-homogeneous maps; some triangular ones do obey it.
Also its phrase that Jacobian failure "requires an indefinite grading"
must be restricted to counterexamples that admit the asserted diagonal
weighted-homogeneous grading. It does not prove that every counterexample
admits such a grading. These prose corrections do not refute its main theorem.

## 5. Fresh exact audit of the three-variable Jacobian map

Let `u=1+xy` and

```
F=(u^3 z+y^2 u(4+3xy),
   y+3x u^2 z+3xy^2(4+3xy),
   2x-3x^2y-x^3z).
```

These are genuine integer-coefficient polynomials on all affine three-space.
The standard-library checker independently expands their coefficients,
differentiates, expands all six determinant terms and obtains exactly `-2`.
It also verifies

```
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
```

Here is a compact independent hand certificate. Only for computing on the
dense open `x!=0`, set

```
w=1+xy, c=2-3xy-x^2z,
A=w^2+w-w^3 c, B=4w+2-3w^2 c.
```

Then `F=(A/x^2,B/x,xc)` and `det d(x,w,c)/d(x,y,z)=-x^3`.
Factoring the powers of `x` from the other Jacobian gives

```
det dF/d(x,w,c) = x^-3 det [[-2A, A_w, A_c],
                           [-B,  B_w, B_c],
                           [ c,    0,   1]]
                       = 2 x^-3.
```

The last equality expands to
`-2A B_w+B A_w+c(A_w B_c-A_c B_w)=2`.
The chain rule gives `det JF=-2` on `x!=0`; as a polynomial identity it
holds everywhere. The rational chart used in the calculation does not impose
any restriction on the original polynomial map or on the collision at `x=0`.

These identities rule out a global inverse. They are verification, not a new
discovery. Current primary sources independently document the same example:
[AFP entry, 2026-07-20](https://isa-afp.org/entries/Jacobian_Counterexample.html)
and [van den Essen, 2026-09-15](https://arxiv.org/abs/2609.17795).
The AFP formalization is **CITED, not locally built**. The older attribution
warning in THM-1300 is historical; current primary sources, not the warning
alone, should govern present provenance claims. No attribution claim is needed
for any proposition above.
Noninjectivity by itself does **not** imply nonsurjectivity, and this packet
makes no such inference. The tournament problem asks for surjectivity onto
the permitted odd integers, not for uniqueness or a canonical inverse. The
one-edge source-guard witnesses above address that realizability issue
directly; the Keller collision only blocks an unsupported inverse shortcut.

## 6. Recovered mechanisms and the exact missing hypothesis

All source paths and blob/byte hashes are in `source_manifest.json`.

| Source | Actual proved/imported mechanism | What it does not supply here |
|---|---|---|
| [THM-2230](sources/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md), planar response fiber | Given a unit mate over a characteristic-zero **field**, `ker D_P=k[P]`; quotient identifies each nonempty fiber exactly and maps onto **im D_P** | No surjectivity onto all polynomials; no integer-domain theorem; no tournament construction |
| [THM-3225](sources/THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape.md), affine clutch/no escape | A specific two-accessory-field family retains critical-resultant multiplicity after explicit boundary accounting | No all-Keller properness theorem, no integer lift, no general spectrum coverage |
| [THM-3477](sources/THM-3477-jacobian-effectivity-moving-numerator-jordan-prufer-response-normal-form.md), moving-numerator effectivity | Exact additive-module debt `B/E_H`; zero-root finite blocks versus nonzero-root divisible rays | Not a quotient ring and not a proof of legal global chart entry; formal/local solvability need not pay effectivity debt |
| [THM-3652](sources/THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof.md), Wright criterion | In dimension two over a characteristic-zero field, **elementary/diagonal factorability** of the Jacobian yields a tame inverse | Determinant alone does not imply that factorability; field inverse need not preserve integer or legal parameter sources |
| [THM-4626](sources/THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md), dual-field equivalence | Under its Keller dual-field setting and imported generic-fiber theorem, line fibers force local nilpotence and a polynomial slice splitting | It is not a theorem about arbitrary response fields or graph counts; using an equivalent open planar-JC assertion would be circular |
| [THM-4727](sources/THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md), mixed-weight quotient | Torus quotient Jacobian has a divisor factor; reconstruction retains orbit scale and polynomial-entry divisibility conditions | Quotient surjectivity alone does not reconstruct a global legal source point; no count-spectrum theorem |
| [THM-1370](sources/THM-1370-elliptic-jacobian-conjecture-all-dimensions.md), elliptic restriction | Contracting positive weights add global control; Proposition 4 gives a constructive integral matching-weight variant | No arbitrary tournament family is shown to satisfy that grading or to preserve its legal domain |

The useful synthesis is therefore a **proof-obligation chain**, not a claim
that the conjectures are equivalent:

```
exact normalized response identity
 -> enough global image/inverse information
 -> integer lattice lift
 -> effective/legal source lift
 -> a tournament with the prescribed count.
```

A constructive sufficient finish would be a finite legal gadget atlas with
exact response laws and a well-founded inverse descent covering every
`q>=0` except `3,10`, terminating in verified base tournaments. An affine
atlas must prove its congruence/lattice conditions and guard inequalities;
polynomial templates must additionally prove polynomial inverse/entry where
used. Neither local rank, an exact response quotient, a continuous image
interval, nor a bounded table discharges those global coverage obligations.

The known all-order omission of 7 and 21 is not the missing step. The remaining
step is to construct every other odd value. The source literature states that
as Conjecture 6.2 in [Chan--Pak, Computational complexity of counting
coincidences](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).
No proposition in this packet closes that conjecture.

## 7. Reproduction and audit boundaries

```
python3 verify_bridge.py > verification.json
python3 -O verify_bridge.py > verification_optimized.json
cmp verification.json verification_optimized.json
```

No third-party packages or network are needed. Failures use explicit
exceptions and remain active under `-O`. The exact proof of Proposition 1
uses Redei's theorem; the finite coefficient check is only a consistency
control. Propositions 2--4 are elementary proofs. The repository mechanisms
are read and scope-audited here, not globally re-proved or fully replayed.
No Git edits, commits, or pushes were made by this lane.
