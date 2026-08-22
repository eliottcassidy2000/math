---
id: THM-3234
title: "Singer owner compactification and pointed Heisenberg carrier gate"
status: >
  PROVED + VERIFIED-EXACT.
  For every odd prime p, adjoining one point to a cyclic set of size p^2-1
  realizes it as the affine plane F_p^2 with its punctured set one Singer
  orbit.  Combining that Singer cycle with any standard minimal affine
  Heisenberg action generates the full AGL_2(F_p), and the Singer cycle does
  not normalize the Heisenberg subgroup.  A center-faithful Heisenberg
  permutation action with a global fixed point needs at least p^2+1 points,
  sharply.  At p=13 this turns the abstract 168-owner completion into a
  full-affine mixing gate: 169 points work only if the added point is mobile;
  a fixed head requires at least 170.  No physical LRC owner map is supplied.
source: root/creative-reframes/2026-08-03
audit: >
  The assertion-independent exact companion checks the Singer order, radial
  and projective factors, transvection-line orbit, full GL_2 closure, affine
  group order, and nonnormality at p=3,5,7,13.  It independently replays and
  pins THM-3224's explicitly finite 168-owner scout, then exhausts every
  Singer gauge to refute the tempting negative-projective-line reading.  The
  finite scout remains outside the proof graph because it has no branch
  stabilization bound.  An independent theorem and replay audit verified
  the uniform proof, the p=13 specialization, and the exact companion.
depends_on:
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-3224-complete-lrc-orbit-bernoulli-gcd-carry-and-owner-hodge-splitting
  - THM-3228-four-jet-heisenberg-minimal-faithful-permutation-carrier-gate
related:
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2614-punctured-target-root-cosupport-and-principal-deck-no-go
  - THM-2877-semilinear-endpoint-rectangle-classification-and-rank-one-defects
  - THM-3090-affine-projective-prime-power-handshake-and-septimal-counterfeit
script: 04-computation/lrc_singer_heisenberg_pointed_carrier_thm3234.py
output: 05-knowledge/results/lrc_singer_heisenberg_pointed_carrier_thm3234.out
script_sha256: 7d1ed1c5f0632c31964e4441810b61994109fb3a72b23965cac7995cc07cde20
output_sha256: 1e024a5a506febca6bfa22046aca7872098962fad512bfcb7855140f7ca95795
hash_basis: LF-normalized bytes
---

# THM-3234 -- Singer owner compactification and pointed Heisenberg carrier gate

**PROVED + VERIFIED-EXACT.**

THM-3224 leaves a maximal cyclic LRC owner orbit of size `168`; THM-3228
shows that a center-faithful mod-13 four-jet permutation carrier needs at
least `169` points.  The numerical equality

```text
168+1=13^2
```

has one exact group-theoretic realization.  A Singer cycle makes the 168
owners the punctured affine plane, but combining that cycle with the minimal
Heisenberg action forces the entire affine general linear group.  The added
point cannot remain a distinguished fixed head.

The theorem is a carrier gate, not an LRC owner construction.

## 1. Uniform Singer--Heisenberg statement

Let `p` be an odd prime and write

```text
V=F_p^2,                         X=V.                    (1)
```

Identify `V` additively with `F_(p^2)`.  Choose a primitive element
`alpha in F_(p^2)^*` and let

```text
sigma(v)=alpha v.                                         (2)
```

Let `H_p` act on `X` through THM-2779's standard affine action

```text
rho(x,y,z)(r,w)=(r+x,w+z-yr).                            (3)
```

Then:

1. `sigma` fixes zero and is one cycle of length `p^2-1` on `X\{0}`;
2. its projective action is one cycle on the `p+1` lines of `P^1(F_p)`,
   while `sigma^(p+1)` is a scalar of exact order `p-1`;
3. the two permutation groups generate

   ```text
   <rho(H_p),sigma>=AGL_2(F_p);                           (4)
   ```

4. `sigma` does not normalize the displayed `H_p` subgroup; and
5. every center-faithful `H_p`-set with a global fixed point has size at
   least `p^2+1`.  This bound is sharp.

Thus the unpointed center-faithful minimum is `p^2`, while the pointed
minimum is

```text
mu_center(H_p; one global fixed point)=p^2+1.             (5)
```

## 2. The punctured Singer orbit

The multiplicative group `F_(p^2)^*` is cyclic of order `p^2-1`, so `(2)`
fixes zero and acts regularly on every nonzero vector.  Moreover

```text
alpha^(p+1) in F_p^*,
ord(alpha^(p+1))=p-1.                                    (6)
```

Indeed, `alpha^(p+1)` is fixed by Frobenius, and its order is

```text
(p^2-1)/gcd(p^2-1,p+1)=p-1.
```

The quotient

```text
F_(p^2)^*/F_p^*
```

is the projective line and has order `p+1`.  Hence the projective class of
`sigma` is regular on `P^1(F_p)`.  Every nonzero vector has the unique
Singer-coordinate form

```text
alpha^(r+(p+1)k),       0<=r<=p,       0<=k<=p-2,        (7)
```

where `r` is a projective direction and `k` is a nonzero scalar gauge.
This proves assertions 1 and 2.

## 3. The joint group is the full affine group

The action `(3)` contains the full translation group `V`: take `y=0` and
vary `x,z`.  Its linear image contains the nontrivial transvection

```text
U:(r,w) |-> (r,w+r),                                     (8)
```

obtained from `rho(0,-1,0)`.  Thus it is enough to prove

```text
<U,sigma>=GL_2(F_p).                                     (9)
```

Let `L` be the fixed line of `U`.  The conjugate

```text
U_k=sigma^k U sigma^(-k)                                 (10)
```

is a transvection with fixed line `sigma^kL`.  Section 2 shows that these
fixed lines run through all `p+1` projective lines.  The powers of any one
nontrivial transvection give its entire root subgroup

```text
{I+tN:t in F_p}.                                         (11)
```

Choose two distinct fixed lines from `(10)` and use them as a basis.  Their
two root subgroups become the upper and lower elementary unipotent groups.
Those groups generate `SL_2(F_p)`.  Finally

```text
det(sigma)=Norm_(F_(p^2)/F_p)(alpha)=alpha^(p+1)          (12)
```

has order `p-1` by `(6)`.  Therefore `(9)` contains every determinant class
and is `GL_2(F_p)`.  Adding the translations already in `(3)` proves `(4)`.

The linear part of the displayed `H_p` is the root subgroup fixing `L`.
Conjugation by `sigma` changes it to the root subgroup fixing `sigma L`,
which is different from `L`.  Hence `sigma` does not normalize this
Heisenberg subgroup.  The joint carrier is full-affine, not a small
semidirect extension of one fixed four-jet chart.

## 4. A fixed head costs one additional point

Let `H_p` act center-faithfully on a finite set `Y`, and suppose `h in Y` is
fixed by all of `H_p`.  If `|Y|<=p^2`, every orbit in `Y\{h}` has size
strictly below `p^2`.  Since `H_p` is a `p`-group, those orbit sizes are only
`1` or `p`.

On a `p`-point orbit, the stabilizer has index `p`, hence is normal and the
image has order `p`.  That image is abelian and kills

```text
Z(H_p)=[H_p,H_p].                                        (13)
```

The center also fixes every one-point orbit.  It would therefore fix all of
`Y`, contradicting center faithfulness.  This proves the lower bound in
`(5)`.  It is sharp: take the disjoint union of the faithful `p^2`-point
action `(3)` and one fixed point.

Equivalently, every center-faithful action of the minimum unpointed degree
`p^2` is transitive and has no distinguished global fixed point.  This is
the exact extra invoice hidden by the cardinal identity `p^2-1+1=p^2`.

## 5. The p=13 owner compactification

For the maximal orbit in THM-3224,

```text
T=168=13^2-1.                                             (14)
```

As a bare cyclic set it therefore has an equivariant compactification

```text
j in Z/168Z |-> alpha^j in F_169^*,
head_*               |-> 0 in F_13^2.                   (15)
```

The owner shift becomes the Singer cycle.  Formula `(7)` gives the exact
factorization

```text
168=(13+1)(13-1)=14 projective directions * 12 gauges.   (16)
```

In the deterministic control `F_169=F_13[u]/(u^2-2)`, the primitive element
`alpha=1+2u` acts by

```text
S=[1 4; 2 1],              ord(S)=168,
det(S)=6,                  ord_13(6)=12.                 (17)
```

Together with `(8)`, this generates

```text
GL_2(F_13),       order 26,208,
AGL_2(F_13),      order 4,429,152.                        (18)
```

Thus `(15)` supplies an exact **abstract** common set for the cyclic owner
shift and the minimal Heisenberg carrier, but it has two unavoidable costs:

- the Heisenberg action moves the added zero through all 169 points; and
- the two actions jointly generate the full affine group.

If `head_*` is required to remain fixed, Section 4 raises the minimum from
`169` to

```text
170.                                                       (19)
```

The word “head” here denotes a proposed distinguished carrier point.  It is
not THM-3224's pre-stabilization asymptotic head, and `(15)` does not identify
either one with a physical LRC state.

## 6. FINITE-EXACT hostile: the sign line is counterfeit

THM-3224 records a **FINITE-EXACT scout**, explicitly without a branch
stabilization bound, whose 168 fitted second-owner coefficients have signs

```text
156 positive,       12 negative,       0 zero.            (20)
```

The multiplicities tempt a stronger Singer reading: after giving the added
point the negative sign, the negative class would have size 13, exactly one
affine line.  Exact replay refutes that geometry.  The negative owners are

```text
N={0,1,2,3,4,5,162,163,164,165,166,167},                 (21)
```

all with value `-17/3087000`.

Under every Singer-equivariant gauge

```text
j |-> c alpha^(aj),               gcd(a,168)=1,           (22)
```

the nonzero points of a line through zero pull back to one congruence class
modulo `14`, because `F_13^*=<alpha^14>`.  The twelve consecutive boundary
owners in `(21)` occupy twelve distinct residue classes modulo 14.  Hence

```text
max_(Singer gauges, vector lines through zero) |N intersect line|=1. (23)
```

So `N union {head_*}` is never an affine/vector line through zero.  (Its
punctured part represents one projective direction, not a projective line.)
The sign-count match in `(20)` is a genuine counterfeit: multiplicity alone
forgets the cyclic placement consumed by the Singer action.

Equations `(20)--(23)` are a finite hostile probe only.  They do not promote
the fitted owner word to an eventual all-dilation theorem.

## 7. Connection contract and scope

```text
source:      the maximal cyclic 168-cell owner orbit of THM-3224;
map:         j -> alpha^j, followed by adjoining zero;
target:      the punctured affine plane plus one completion point;
preserved:   cyclic owner shift and the exact 168+1 cardinality;
destroyed:   interval geometry, phase endpoints, q-values, physical ancestry;
new action:  the minimal center-faithful H_13 affine action;
joint cost:  full AGL_2(F_13), with the completion point mobile;
sidecar:     a lawful owner-to-endpoint-origin map carrying the LRC predicate;
hostiles:    fixed-head lower bound and the refuted sign-line scout.
```

The theorem does not identify THM-3224's Bernoulli gcd-carry curvature with
the four-jet center, construct the mixed transition square required by
THM-2591, or select a graph in THM-2614's punctured cospan.  It gives no
cellwise safety, finite-head stabilization, physical entry, row exclusion,
or proof of `LRC(14)`.

Its exact gain is a sharper design rule.  A 169-state proposal must let the
extra state participate in one transitive Heisenberg orbit and must survive
full affine mixing.  A proposal that insists on a fixed basepoint needs at
least 170 states before any semantic test begins.

## 8. Exact companion

Run

```text
python3 04-computation/lrc_singer_heisenberg_pointed_carrier_thm3234.py
python3 -O 04-computation/lrc_singer_heisenberg_pointed_carrier_thm3234.py
```

and compare LF-normalized bytes with the declared output.  The companion
uses explicit `require` gates and exact finite-field arithmetic.  At
`p=3,5,7,13` it verifies the Singer order, scalar and projective factors,
all transvection directions, full `GL_2` closure, affine orders, and failure
of normalization.  It pins and byte-replays THM-3224 before reading the
finite scout, then exhausts every gauge in `(22)` for `(23)`.  Ordinary and
optimized modes must byte-match the stored transcript.

QED.
