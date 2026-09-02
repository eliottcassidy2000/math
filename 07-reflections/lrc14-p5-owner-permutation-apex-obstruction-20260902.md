# LRC(14) unit-tail owner permutations and the `p=5` apex obstruction

**Status:** PROVED elementary owner/component theorem + VERIFIED-EXACT one-row
control.  This is not yet canonical, does not close the minority-anchor branch,
and does not prove LRC(14).

## Inheritance and portfolio

The anchor is the degree-twelve minority row from
[THM-4330](../01-canon/theorems/THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve.md),

```text
S_h={2h} union W,       |W|=12,       W positive odd,       420|h.
```

The closest proved mechanisms are THM-4330 `MC2a--MC2d` (two-sheet nearest-
integer owners), its unit equality-wall safe-band lemma `MC9a--MC9f`, and
[THM-771](../01-canon/theorems/THM-771-sheet-endpoint-defect-and-reduced-winding-pierce.md)
(strict endpoint slack and owner overlap).  The prime-seven token polynomial
in [THM-773](../01-canon/theorems/THM-773-prime-seven-sheet-monodromy-and-tournament-fibre.md)
is the direct model for the `p=5` checksum below.  The canonical hostile is
THM-4330's `P=9009` row: marginal capacity is exactly one, but endpoint slack
and repeated owners leave three free sheets.  The corrected near miss is to
read capacity equality as coverage.  The least-used sidecar is the affine
intercept `n_w=round(w theta)`, not only the tail speed or its capacity.

The concept board was:

```text
unit singleton mask | owner permutation | endpoint disappearance
tooth-centre determinant | core apex margin | finite-field moments
```

The niche was the componentwise determinant left after the capacity quotient.
The wildcard was the `F_5` polynomial checksum.  Both returned positive signal
and meet in the theorem below.

## The unit-owner permutation theorem

Let `2<=d<=7`.  Let `C` be a nonempty finite set of positive integers and let
`T` consist of exactly `d` distinct positive integers, each coprime to `d`.
Put

```text
S=dC union T,
G_C={theta in R/Z:min_(c in C)||c theta||>=1/14}.
```

For a chosen real representative of `theta`, define the strict bad-sheet mask

```text
B_w(theta)={j in Z/dZ:||w(theta+j)/d||<1/14}.          (1)
```

### A. Exact singleton and owner formula

For every `w in T`,

```text
B_w(theta)=empty,                                      if ||w theta||>=d/14,
B_w(theta)={kappa_w(theta)},                           if ||w theta||< d/14,

kappa_w(theta)=-n_w(theta) w^(-1) mod d,
n_w(theta)=round(w theta).                             (2)
```

The nearest integer in the second line is unique.  In particular, equality
in the first line gives an empty **strict** mask.  When `d=7` and
`w theta` is a half-integer, nearest-integer ownership is not defined and is
not needed: the mask is empty.

Changing the representative from `theta` to `theta+1` sends

```text
n_w -> n_w+w,                 kappa_w -> kappa_w-1.    (3)
```

Thus individual labels are gauge-dependent, while coverage, bijectivity, and
the determinants below are gauge-invariant.

### B. Strict counterexample iff owner permutation everywhere

The following are equivalent:

1. `M(S)<1/14`;
2. for every `theta in G_C`, every tail is active,
   `||w theta||<d/14`, and

   ```text
   w -> kappa_w(theta)
   ```

   is a bijection from `T` to `Z/dZ`.

This is an exact characterization, not merely a sufficient test.  It is the
information retained when the equality `d` in the scalar capacity bound is
resolved into locations.

### C. Component address and determinant box

Assume the equivalent strict-counterexample conditions.  Let `I` be a
positive-length connected component of `G_C`, lifted to a closed interval in
`R`.  For each `w in T` there is one integer `n_w`, constant on `I`, such that

```text
I subset (n_w/w-d/(14w), n_w/w+d/(14w)).              (4)
```

The owner labels `-n_w w^(-1)` form one fixed permutation on the component.
For distinct `u,v in T`, set

```text
Delta_(u,v)=n_u v-n_v u,       D_(u,v)=|Delta_(u,v)|. (5)
```

Then

```text
d does not divide Delta_(u,v),
gcd(u,v) divides Delta_(u,v),                              (6)

|I| < min(
  d/(7u),
  d/(7v),
  d/(14u)+d/(14v)-D_(u,v)/(uv)
).                                                        (7)
```

The third expression is positive.  It is not an estimate inserted by hand:
the minimum on the right of `(7)` is the exact length of the intersection of
the two open teeth in `(4)`.

### D. Apex-intercept obstruction

Put

```text
R=max C,       mu=M(C)>1/14,       rho=(mu-1/14)/R,    (8)
```

and choose any maximizer `theta_0` of the core envelope.  A strict
counterexample necessarily satisfies, for every `w in T`,

```text
||w theta_0||+w rho < d/14,                            (9)
```

and its owners at `theta_0` form a permutation.  For every pair `u!=v`, with
the constant integers from the core-safe radius-`rho` interval,

```text
D_(u,v)/(uv)+2rho < d/(14u)+d/(14v).                  (10)
```

Consequently either reversed nonstrict inequality is a safety certificate.
This retains the intercept `||w theta_0||`; dropping it from `(9)` recovers
THM-4330's coarser height certificate.

When `|C|=13-d`, cited LRC through thirteen total runners gives

```text
mu>=1/(14-d),       rho>=d/[14(14-d)R].               (11)
```

Thus every strict counterexample on this equality branch obeys

```text
w<(14-d)R                                             (12)
```

for every tail and, pairwise,

```text
2uv < (14-d)R (u+v-14D_(u,v)/d).                     (13)
```

Equation `(12)` is the known scale-free safe-band boundary.  Equation `(13)`
is its owner-address strengthening: `D` is nonzero modulo `d`, so it cannot be
discarded as an arbitrary nonnegative error.

### E. The `p=5` moment checksum

For `d=5`, at any phase where all five owners are defined, they form an exact
sheet permutation if and only if

```text
product_(w in T)(X-kappa_w)=X^5-X              in F_5[X],       (14)
```

or equivalently

```text
(sum kappa, sum kappa^2, sum kappa^3, sum kappa^4)
  =(0,0,0,-1)                                  in F_5^4.       (15)
```

This is the five-sheet analogue of THM-773's seven-token checksum.  The new
content for the minority-anchor lane is not the finite-field identity alone;
it is the coupling of `(14)` to the physical affine integers `n_w`, the
strict endpoint disappearance in `(2)`, and the metric bounds `(7)--(13)`.

## Proof

Suppose first that sheet `j` lies in `B_w(theta)`.  There is an integer `m`
with

```text
|w theta+w j-dm|<d/14<=1/2.                           (16)
```

Hence `n=dm-wj` is a nearest integer to `w theta`; the strict inequality makes
it unique.  It satisfies `n+wj=0 mod d`, whose unique solution is
`j=-n w^(-1)`.  Conversely, if `||w theta||<d/14`, choose the unique nearest
integer `n` and that residue `j`.  Dividing

```text
|w theta-n|<d/14
```

by `d` proves that `j` is bad.  This proves `(2)`.  At equality, any bad lift
would produce an integer strictly closer than `d/14`, which is impossible.
This also handles the `d=7` half-integer without assigning it an owner.
Equation `(3)` follows from `round(w(theta+1))=round(w theta)+w` whenever the
owner is defined.

For `theta in G_C`, all `d` lifts `(theta+j)/d` preserve the core because

```text
||dc(theta+j)/d||=||c theta||.                         (17)
```

By `(2)`, the `d` tails kill all `d` lifts exactly when all masks are
singletons and their labels are distinct.  Conversely, every physical phase
maps to one quotient phase under multiplication by `d`; if that quotient is
outside `G_C`, the core already kills it.  This proves the iff in B, including
both directions and all boundary points.

Under strict failure, `(2)` puts a connected closed component `I` inside the
open active set of every tail.  The components of that active set are exactly
the teeth in `(4)`, so one integer address `n_w` is constant on all of `I`.
For two owners,

```text
Delta_(u,v)=uv(kappa_v-kappa_u) mod d.                 (18)
```

The owner difference is nonzero and `uv` is a unit, proving the first part of
`(6)`; the second is immediate from `(5)`.  The tooth centres are separated by
`D/(uv)` and their radii are `d/(14u),d/(14v)`.  Direct interval intersection
therefore gives exact length

```text
min(d/(7u),d/(7v),d/(14u)+d/(14v)-D/(uv)).            (19)
```

Because the closed interval `I` is strictly inside that open intersection,
its length is strictly smaller.  This proves `(7)` and also shows that the
third term is positive.

The core envelope is `R`-Lipschitz.  Hence the closed circular interval of
radius `rho` around `theta_0` lies in `G_C`.  Under strict failure it must be
strictly contained in every tail tooth.  Comparing its centre and radius with
that tooth gives `(9)`; comparing its length with the two-tooth intersection
gives `(10)`.  Substitution of `(11)` yields `(12)--(13)`.

Finally, a permutation of `F_5` plainly gives `(14)--(15)`.  Conversely,
Newton identities in degrees one through four are invertible in `F_5` and
turn `(15)` into

```text
e_1=e_2=e_3=0,       e_4=-1.
```

Thus the owner polynomial is `X^5-X-e_5`.  Evaluating it at any one of its
`F_5` roots and using `a^5=a` forces `e_5=0`, proving `(14)` and the converse.
QED.

## A simultaneous marginal hostile closed by the apex intercept

Take

```text
h=420,
C=(13,168,349,375,711,737,1073,1099),
T=(1,21,23,327,689),
S=5C union T.
```

Then `5*168=840=2h`; the other twelve speeds are distinct positive odds, and
the row is primitive.  Exactly five odd relatives are nonmultiples of five,
so this is the `r_5=5` equality branch.

The core maximum is exact:

```text
M(C)=90/181,                 theta_0=7/181.            (20)
```

Indeed every member of `C` is congruent to `13` or `168` modulo `181`, and
all eight distances at `theta_0` are `90/181`.  For the upper bound, suppose
both `||13 theta||` and `||168 theta||` exceeded

```text
90/181=1/2-1/362.
```

Write

```text
13theta=A+1/2+e,       168theta=B+1/2+f,
|e|,|f|<1/362.
```

The identity `168(13theta)-13(168theta)=0` would make a half-integer cancel
`168e-13f`, but

```text
|168e-13f|<181/362=1/2,
```

which is impossible.  Thus the pair `{13,168}` already supplies the upper
bound in `(20)`.

At `theta_0`, all five tails are active and their exact data are

```text
w       1       21       23       327       689
n_w     0        1        1        13        27
error  7/181   34/181   20/181   64/181    64/181
kappa   0        4        3         1         2.       (21)
```

Thus the central owner deck is a perfect permutation; an instantaneous
collision check does not close it.  Nevertheless

```text
rho=1079/2784866,
5R/[14(mu-1/14)]=994595/1079>689,                     (22)
```

so the old height-only safe-band certificate also does not fire, while the
left sides of `(9)` already exceed `5/14` for tails `327` and `689`.

Move by only

```text
theta_0-theta_1=3/276206<rho,
theta_1=59/1526.                                      (23)
```

Tail `327` reaches equality and has no strict owner.  In fact tail `689` is
also inactive, and the five masks become

```text
1:{0},       21:{4},       23:{3},       327:empty,       689:empty. (24)
```

The free sheet `j=2` gives

```text
t=(theta_1+2)/5=3111/7630,
min_(s in S)||st||=551/7630=1/14+3/3815.              (25)
```

The same row is complete for every THM-366 denominator `2,...,14`, represents
every required sign class in each useful doubled-denominator unit bank, and
has no adaptive-capacity closure at any represented divisor.  Its exact
capacity minimum is `1`, uniquely at `d=5`.  Both THM-4330 half-turn clocks
are killed with clearance only `11/336`.  These controls make the affine
intercept in `(9)` the decisive added coordinate among those marginal tests.

The row does have `1,416` witnesses on the full integer half-turn grid.  It is
therefore not claimed as a hostile to that larger grid, and `(25)` is a new
certificate for a safe row, not a counterexample.

## Connection contract and scope

```text
source:       equality-capacity row dC union T with d unit tails
target:       labelled singleton tokens, then addressed open-tail teeth
map:          theta=d t; n_w=round(w theta); kappa_w=-n_w w^(-1) mod d
preserved:    strict full-sheet cover and existence of a physical safe lift
destroyed:    a token at one phase forgets distance to its next endpoint
sidecar:      tooth integer n_w, core component/apex, R, mu, and rho
test:         permutation polynomial plus inequalities (9)--(10)
```

The symbolic theorem holds for `2<=d<=7`, but the scale-free corollary uses
the cited lower-dimensional LRC input and the cardinality `|C|=13-d`.  The
permutation conclusion uses the equality `|T|=d`; with more tails, full cover
is polynomial divisibility rather than a bijection and collisions need not
free a sheet.  The `p=5` branch is the new application here.  The `p=7`
tokenization overlaps THM-771/773, while the `d=2,3` height shadows overlap
older affine cones.  No arbitrary `r_5>5` branch, entry theorem, or proof of
LRC(14) follows.

## Exact audit

The companion checks `710,910` exact requirements.  It verifies the singleton
formula at `72,044` rational `(d,w,theta)` cases, including both equality
boundaries and every `d=7` half-integer tie; it verifies `47,716` gauge shifts;
and it independently rebuilds `40,690` pairwise tooth intersections and
determinant residues.  It then audits every displayed quantity of the anchored
row, with two independent divisor constructors.

Artifacts:

```text
04-computation/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.py
05-knowledge/results/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.out
```

Hashes over raw LF bytes:

```text
script cb4f7118398c0b8f2049669742fe18b99688dc679bd72f01b59e444cea77ae79
output ff3df7ff3d6412a9597e1acf63bc4d32e24ffb2e53d6d603c8510d3228917383
```

Replay from the repository root:

```bash
python3 -B 04-computation/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.py \
  | diff -u 05-knowledge/results/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.out -
python3 -B -O 04-computation/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.py \
  | diff -u 05-knowledge/results/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.out -
PYTHONHASHSEED=5209 python3 -B \
  04-computation/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.py \
  | diff -u 05-knowledge/results/lrc14_p5_owner_permutation_apex_obstruction_probe_20260902.out -
```
