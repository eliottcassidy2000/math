---
id: THM-4009
title: "Euclidean covering-transference short-relation compression"
status: >
  CITED + PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A hypothetical primitive LRC(14)
  counterexample has a nonzero integer Graver relation a with ||a||_2<14,
  sum a_i^2<=195, ||a||_1<=50, and |a_i|<=13. Its support-two branch has
  exactly 47 reduced ratios, 3,666 unoriented ratio/support packets, and
  7,332 oriented labelled assignments. The projected
  half-lattice centre also forces some odd-sum relation, not proved to be the
  short one. This sharply supersedes THM-3743's numerical cap but remains a
  necessary reduction, not a proof of LRC(14).
source: root / 2026-08-24
audit: >
  PASS. The primary standard-library companion recomputes the
  strict Euclidean threshold, all integer conversions, the 47 pair ratios,
  55,459 primitive absolute histograms, 5,030,161 signed-multiset types,
  labelled Euclidean-ball counts, the rank-eleven Hadamard terminal, the Pell
  intersection, parity scope, and AP hostile. Normal and optimized runs agree
  byte-for-byte with the frozen output. A separately implemented referee uses
  a multiplicity generating function, polynomial convolution, Moebius sieve,
  balanced-squares optimization, and direct hostile proofs. It independently
  reproduces every implication and finite count. It caught the distinction
  between 3,666 unoriented packets and 7,332 oriented labelled assignments;
  the theorem and primary companion now carry the repair. The cited input is
  only Banaszczyk's
  Euclidean covering-radius/dual-minimum transference inequality; the
  projected-zonotope equivalence and quotient-dual algebra are inherited and
  reproved in scope below.
depends_on:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
related:
  - THM-1008-lrc13-descent-floor
  - THM-1043-the-spread-ladder-and-the-refined-covering-map
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
  - THM-2169-bounded-relation-on-every-lrc-deletion
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-3825-prime-colour-valuation-two-cube-decoder
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-3995-scale-two-parity-hole-support-and-integer-variance-tariff
  - THM-4002-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family
  - THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression
  - THM-4015-first-kind-character-sensitive-foster-transference
script: 04-computation/lrc14_euclidean_covering_transference_audit_thm4009.py
output: 05-knowledge/results/lrc14_euclidean_covering_transference_audit_thm4009.out
script_sha256: 4aedcfaebde66295d07b600f063e880a6fb4a81f11f5d6a597f9f35835394df5
output_sha256: 0924afeded809ac23b65cb4f6a635e39e4b24804041970e83ab8fcd826348952
independent_script: 04-computation/lrc14_euclidean_covering_transference_independent_referee_thm4009.py
independent_output: 05-knowledge/results/lrc14_euclidean_covering_transference_independent_referee_thm4009.out
independent_script_sha256: 3830131aaed502bde29e933d49c9c740c53718900ab55ae8bf736cac7c743511
independent_output_sha256: 0967cff265a8de3c8bafd3745eac946c6ce0fdfca42bcef9de890db907c7da77
hash_basis: raw working-tree bytes
---

# THM-4009 -- the lonely zonotope contains a large Euclidean ball

**CITED + PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
The only new external input is
Banaszczyk's Euclidean transference theorem from
[W. Banaszczyk, *New bounds in some transference theorems in the geometry of
numbers*, Math. Ann. 296 (1993), 625--635](https://doi.org/10.1007/BF01445125):
Theorem (2.2), journal page 632, states in the dualized notation used here
for a `d`-dimensional Euclidean lattice `L`,

```text
2 mu(L) lambda_1(L*) <= d.                              (1)
```

Here `mu` is the Euclidean covering radius and `lambda_1` the Euclidean
length of a shortest nonzero vector. The factor `2` is load-bearing. The
lonely-runner zonotope equivalence of Beck--Hosten--Schymura and the exact
projected dual from THM-3743 are the inherited inputs.

Let `n=(n_1,...,n_13)` be a primitive vector of distinct positive integer
speeds. Put

```text
V=n^perp,                         dim V=12,
pi:R^13 -> V                      orthogonal projection,
Lambda=pi(Z^13),
Z(n)=pi([1/14,13/14]^13).                              (2)
```

If `n` were a counterexample to LRC(14), there would be a nonzero
`a in Z^13` with

```text
a dot n=0,
||a||_2<14,
sum_i a_i^2<=195,
||a||_1<=50,
max_i |a_i|<=13.                                      (3)
```

One may choose `a` to be a Graver relation. Independently, a counterexample
forces some integer relation `b dot n=0` with odd coefficient sum; this
possibly different parity relation is not bounded by `(3)`. LRC(14) remains
open.

## 1. Exact Euclidean inball

Write

```text
[1/14,13/14]^13=(1/2)1+(3/7)[-1,1]^13,
c=(1/2)pi(1).                                           (4)
```

For every `x in V` with `||x||_2<=1`, each coordinate satisfies
`|x_i|<=1`, and `pi(x)=x`. Therefore

```text
B_V(0,1) subset pi([-1,1]^13),
B_V(c,3/7) subset Z(n).                                 (5)
```

This elementary inclusion is the missing scale coordinate. It does not use
a volume estimate, a Gaussian envelope, or the false lattice-tail
factorization excluded elsewhere in the repository.

There is also an exact speed-dependent version. For
`K=pi([-1,1]^13)` and a unit vector `u in V`, the support function is

```text
h_K(u)=sum_i |u_i|.                                    (5a)
```

Hence the centred Euclidean inradius of `K` is
`rho(n)=min_(u in V,||u||_2=1)||u||_1`. Normalize instead by
`||u||_1=1` and maximize the convex function `||u||_2` on each signed
cross-polytope slice. An extreme point has exactly two nonzero coordinates:
the equations `sum |u_i|=1` and `n dot u=0` leave at most two positive
absolute-coordinate variables, and balance requires opposite signs. The
pair `(i,j)` gives

```text
||u||_1/||u||_2=(n_i+n_j)/sqrt(n_i^2+n_j^2).           (5b)
```

The right side decreases with the ratio of the larger speed to the smaller.
Writing `m=min_i n_i` and `M=max_i n_i`, one obtains the exact formula

```text
rho(n)=(m+M)/sqrt(m^2+M^2)>1,
B_V(c,(3/7)rho(n)) subset Z(n).                         (5c)
```

The universal radius `3/7` is the unbounded-spread limit, not the exact
radius of any fixed positive speed row.

A hypothetical counterexample makes the **closed** zonotope `Z(n)` disjoint
from `Lambda`. Hence the closed ball in `(5)` is disjoint from `Lambda`.
Because a Euclidean lattice is discrete, the nearest lattice point to `c`
exists, and

```text
dist(c,Lambda)>3/7,
mu(Lambda)>3/7.                                         (6)
```

The strict inequality is justified by closed-ball disjointness; replacing it
by a merely interior-free assertion would lose the final integer unit.

## 2. Banaszczyk transference gives the short row

THM-3743 proves the exact dual identity

```text
Lambda*=Z^13 intersect n^perp.                          (7)
```

Apply `(1)` to the twelve-dimensional lattice `Lambda`. From `(6)`,

```text
lambda_1(Lambda*)
 <=6/mu(Lambda)
 <14/rho(n)
 =14 sqrt(m^2+M^2)/(m+M)
 <14.                                                   (8)
```

A shortest nonzero vector `a in Lambda*` is therefore an ordinary integer
speed relation satisfying the first two lines of `(3)`. Since its squared
Euclidean norm is an integer strictly below `196`, it is at most `195`, and
every coordinate has absolute value at most `13`.

Cauchy--Schwarz gives

```text
||a||_1^2<=13 sum_i a_i^2<=2535<51^2,                  (9)
```

so `||a||_1<=50`. This integer cap is arithmetically sharp under the norm
constraint: the unique sorted absolute vector of sum `50` with minimum
squared norm is

```text
(3,3,4,4,4,4,4,4,4,4,4,4,4),
```

whose squared norm is `194`. This says nothing about its realizability as a
speed relation or counterexample certificate.

Formula `(8)` strengthens every bounded-spread branch. If `M/m<=R`, the
integer square-norm and `l1` caps for several useful exact thresholds are

```text
R                 2    3    5    13    21    50
sum a_i^2 <=     108  122  141  169   178   188
||a||_1 <=        37   39   42   46    47    49.        (9a)
```

These are conditional refinements, not a global spread bound. In particular,
the `R=13` line gives the strict square threshold `<170`, hence the integer
cap `169`. It is a calibration rather than a live counterexample branch:
THM-1043's `n=14` spread rung already supplies the explicit lonely time
`1/(14m)` whenever `M/m<=13`. Conversely, THM-1008 forces the largest-to-
second-largest ratio below `13` in every counterexample. Thus any live row is
globally spread but cannot consist of one dominant top speed.

Choose `a` shortest in Euclidean norm. If `a=u+v` were a nontrivial conformal
decomposition by integer kernel vectors, then every coordinate of `u` is
bounded in absolute value by the corresponding coordinate of `a`, with a
strict inequality somewhere. Thus `||u||_2<||a||_2`, a contradiction. Hence
`a` is conformally indecomposable: it is a Graver element of the one-row
speed matrix.

## 3. Exact finite atlases

If `a` has support two on speeds `n_i,n_j`, primitiveness forces, up to
global sign,

```text
(a_i,a_j)=(n_j/g,-n_i/g),        g=gcd(n_i,n_j).        (10)
```

Consequently its reduced positive coefficients `p<q` satisfy

```text
p^2+q^2<=195.                                           (11)
```

There are exactly `47` such coprime ratios and `3,666` choices after selecting
an **unoriented** coordinate support. Assigning which labelled coordinate
carries `p` and which carries `q` gives `7,332` oriented assignments modulo
global sign. Every ratio has `p+q<=19`. The
old square-triangular/Pell selector packet `(1,5),(5,29),(29,169)` leaves only
`(1,5)` under `(11)`.

If one discards the Euclidean norm and keeps only the derived `l1<=50` cap,
the corresponding support-two envelope has `386` coprime ratios with
`p+q<=50`.  Exactly `339` are coarse false positives for `(11)`, leaving the
native `47`.  The `386` count must not be quoted as the THM-4009 atlas.

For higher support, the finite exact companion gives a practical orbit-level
upper atlas before coordinate labelling:

```text
support:                 2    3    4    5     6     7     8
absolute histograms:    47  209  566  1177  2057  3180  4490

support:                 9    10    11     12     13
absolute histograms:  5911  7374  8805  10188  11455.   (12)
```

These are primitive nondecreasing positive coefficient multisets with square
sum at most `195`; the impossible support-two multiset `(1,1)` is removed.
There are `55,459` total absolute histograms and `5,030,161` nonconstant sign-
multiset types modulo equal-entry permutations and global sign. They are an
ambient atlas, not a claim that every type annihilates a distinct positive
speed row.

For comparison, the exact labelled Euclidean ball contains

```text
711,202,814,025,242
```

nonzero integer vectors. This is still too large for blind enumeration, but
it is more than `2,782,564,628,359` times smaller than THM-3743's ambient
`l1<=356` universe. The histogram atlas `(12)`, not the labelled ball, is the
right first traversal.

## 4. Join with the rank-eleven relation code

Let `W` be THM-2052's rational span of support-at-most-three relations of
height `Q=91^6`. If `dim W=12`, the inherited finite terminal already fires.
If `dim W=11`, then

```text
a notin W  => dim span(W,a)=12,
a in W     => W contains a nonzero row with sum a_i^2<=195. (13)
```

In the first branch, eleven independent rows of `W` have Euclidean norm at
most `sqrt(3)Q`; the last row has squared norm at most `195`. Cofactors and
Hadamard now give the explicit terminal

```text
max_i n_i <= floor(sqrt(195*3^11)*Q^11)

 =11639011567946276516330452125265832396450210671398535269626998205847761612881723251682730932524351641613494572161619231675129410858204. (14)
```

It has `134` digits. This cap is not computationally attractive by itself;
its value is the exact rank fork. Every row from `(3)` of support at most
three has height at most `13<Q`, so it is already one of THM-2052's
generators and lies in `W`. Only support at least four can supply the outside-
`W` rank increment.

The unresolved branch is therefore much narrower than before:

```text
rank-eleven star space contains a primitive Graver row
with coefficient square sum at most 195.                (15)
```

This does not by itself force a crossing row in the two-component incidence
model: a component may absorb a vector satisfying `(15)`. The component and
owner sidecars remain necessary.

## 5. The half-lattice parity sidecar

Equation `(4)` gives

```text
2c=pi(1) in Lambda.                                    (16)
```

A counterexample has `c notin Lambda`. Since `Lambda**=Lambda`, some
`b in Lambda*` pairs nonintegrally with `c`. But for every such integer
relation,

```text
<b,c>=(1/2) sum_i b_i.                                 (17)
```

Thus some speed relation has odd coefficient sum. This is a genuine parity
sidecar coming from the centre's nontrivial order-two character. It must not
be silently attached to the shortest row: `(1)` supplies an unphased minimum,
and no argument here bounds the shortest **odd-character** relation. A
character-sensitive covering transference theorem would be a strictly new
ingredient.

## 6. Hostile, loss ledger, and stopping boundary

The arithmetic progression `(1,2,...,13)` has the relation
`(1,-2,1,0,...)`, with squared norm `6`, and is nevertheless lonely at the
valid boundary time `t=1/14`. Hence even the much sharper resonance is only
necessary.

```text
source:      lattice-point-free projected LRC(14) zonotope
target:      Euclidean shortest vector in the integer relation lattice
map:         inball depth -> covering radius -> dual minimum
preserved:   counterexample => one Graver row with square norm <=195
destroyed:   support, sign partition, endpoint owner, phase, arrival, parity
sidecars:    unbounded odd character; THM-2052/2169; two-component owner data
hostile:     AP13 has square norm 6 and is safely lonely
next test:   traverse the 55,459 absolute histograms inside rank-eleven stars,
             retaining component and owner words.                         (18)
```

The previous Banaszczyk-related lattice-tail no-go in the repository is not a
counterexample to this argument: it refuted a claimed uniform absolute
Fourier/theta summation. Here `(1)` is used once as an exact covering-radius
transference inequality, with no tail factorization.

## 7. Reproduction

```bash
python3 -B 04-computation/lrc14_euclidean_covering_transference_audit_thm4009.py
python3 -B -O 04-computation/lrc14_euclidean_covering_transference_audit_thm4009.py
python3 -B 04-computation/lrc14_euclidean_covering_transference_independent_referee_thm4009.py
python3 -B -O 04-computation/lrc14_euclidean_covering_transference_independent_referee_thm4009.py
```

Each primary/referee pair matches its frozen output byte-for-byte. The two
scripts use different counting mechanisms. They verify the finite arithmetic
and hostile controls; Sections 1--6 carry the mathematical implications. The
theorem does not recover simultaneous loneliness or prove LRC(14). **QED.**
