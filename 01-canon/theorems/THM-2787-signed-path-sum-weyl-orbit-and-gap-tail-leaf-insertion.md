---
id: THM-2787
title: "Signed path-sum Weyl orbit and gap-tail leaf insertion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A tree with m
  edges is graceful exactly when a signed permutation
  of 1,...,m integrates to root-path potentials for which every nonempty
  simple path has signed sum of absolute value between 1 and m.  Equivalently,
  reduced tree incidence carries the A_m vertex-permutation orbit of
  (0,...,m) into the B_m signed-permutation orbit of (1,...,m).  For a fixed
  graceful labeling, lifting across a vacant label gap and attaching one
  marked leaf is graceful exactly when the gap-crossing column is a suffix in
  the edge-difference order and the new edge fills the missing difference.
  Every marked vertex of every nonisomorphic tree through nine vertices has
  such a gap-tail witness.  The corresponding all-tree marked extension
  statement remains a conjecture; it would imply, but does not prove, the
  Graceful Tree Conjecture.
source: graceful-signed-path-audit/gap-tail-2026-07-28
audit: >
  root/2026-07-28 independently rederived the translation-quotient
  A_m/B_m orbit typing, suffix-column iff, marked-vertex quantifier, and
  conjecture/QED boundary, and replayed normal/-O/stored exactly: ACCEPT.
depends_on:
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
related:
  - THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge
  - THM-2765-rooted-nullstellensatz-linear-range-distinct-edge-labeling
  - THM-2783-weighted-long-wall-binary-null-avoidance-and-ternary-state-reconstruction
  - THM-2786-binary-golomb-universal-edge-difference-separation-and-graceful-compression-boundary
script: 04-computation/signed_path_gap_tail_graceful_thm2787.py
output: 05-knowledge/results/signed_path_gap_tail_graceful_thm2787.out
script_sha256: fa87047f0ac54b1767583e9f099faf99b91a89eada4b1f5938cbd4a5895aba04
output_sha256: d7d3a31af3e2f4c65a7fd94f5fb56178e3cee62cf2ee2f1d3890f2dcda4a8c42
hash_basis: LF-normalized bytes
---

# THM-2787 -- graceful labeling is bounded signed path integration

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2770 proves that reduced tree incidence is a unimodular map from vertex
potentials to edge gradients.  The inverse is a root-path sum.  Applying that
inverse not to arbitrary gradients but to the signed-permutation orbit of
`(1,...,m)` gives an exact finite reformulation of graceful labeling:

```text
every nonempty tree path must have a nonzero signed sum of size at most m.
                                                                    (1)
```

This separates the real graceful debt into two independent constraints.
Binary radix weights solve the nonzero-sum half universally but at
exponential span; the compact alphabet `1,...,m` must solve nonvanishing and
the upper cap simultaneously.

The same language yields an exact one-leaf extension rule.  Insert one vacant
label into a graceful vertex order.  Every old edge difference either stays
fixed or increases by one.  The new differences remain collision-free
exactly when the edges crossing the gap form an upper tail in the old
edge-difference order.  This produces a concrete marked induction target,
but its all-tree existence is still open.

## 1. Signed path sums and the `A_m/B_m` orbit intersection

Let `T` be a finite tree with `m>=1` edges and `m+1` vertices.  Choose a root
`r`, orient every edge away from it, and let

```text
partial_T: Z^V/Z1 -> Z^E                               (2)
```

be reduced incidence.  By THM-2770, `(2)` is a unimodular lattice
isomorphism and

```text
(partial_T^-1 d)_v
 =sum_(e in path(r,v)) d_e.                             (3)
```

Choose a bijection

```text
pi:E(T)->{1,...,m}                                      (4)
```

and signs `epsilon_e in {+1,-1}`.  Put

```text
d_e=epsilon_e pi(e),
x_r=0,
x_v=sum_(e in path(r,v)) d_e.                           (5)
```

When a path is traversed against the root orientation, reverse the
corresponding sign.  Thus for every two vertices

```text
x_v-x_u
 =sum_(e in path(u,v)) eta_(u,v;e) epsilon_e pi(e),      (6)
```

where every displayed `eta` is `+1` or `-1`.

The following are equivalent.

1. After a translation, `x` is a graceful labeling of `T`.
2. The potentials are distinct and

   ```text
   max_v x_v-min_v x_v <=m.                             (7)
   ```

3. Every nonempty simple path `P` satisfies

   ```text
   1<=|sum_(e in P) eta_(P,e)epsilon_e pi(e)|<=m.        (8)
   ```

Indeed, `(6)` makes distinctness exactly the nonzero half of `(8)` and makes
`(7)` exactly its upper half.  The edge carrying magnitude `m` forces the
span in `(7)` to be at least `m`, hence exactly `m`.  There are `m+1`
distinct integer potentials in an interval of length `m`, so they are all
the consecutive integers in that interval.  Translating its minimum to zero
gives vertex labels `{0,...,m}`, while `(4)--(5)` already give edge
differences `{1,...,m}`.  The converse is immediate from the difference of
two graceful vertex labels.

There is also an exact Weyl-orbit statement.  Write

```text
lambda=(0,1,...,m) in Z^(m+1),
rho=(1,2,...,m) in Z^m.                                 (9)
```

Let `W(A_m)=S_(m+1)` permute vertex coordinates and let
`W(B_m)=C_2^m semidirect S_m` signed-permute edge coordinates.  In the
translation quotient,

```text
T is graceful
iff partial_T(W(A_m)lambda) intersects W(B_m)rho.        (10)
```

This is a lattice-orbit intersection, not the arrangement nonvanishing of
THM-2770.  The latter uses the `A_m/D_m` discriminants to avoid collisions;
the finite target orbit in `(10)` is of type `B_m` because every individual
edge sign is allowed.

Changing the root adds a common constant before taking the root-zero section.
Reversing every sign negates every potential and, after translation, is the
usual complementary graceful labeling.  Thus `(8)--(10)` do not depend on
the auxiliary root or orientation.

## 2. Both path conditions, and all vertex pairs, are load-bearing

Neither half of `(8)` can be removed.

- On `P3`, use successive signed edge values `+1,+2`.  The potentials are
  `(0,1,3)`.  Every nonempty path sum is nonzero, but the full path has size
  `3>m=2`.
- On `P4`, use successive signed edge values `+1,-3,+2`.  The potentials are
  `(0,1,-2,0)`.  Every path sum has absolute value at most `m=3`, but the
  full path sum is zero.
- Checking only root paths also fails already on the three-vertex star
  rooted at its centre.  Signed leaf values `+1,-2` are individually
  nonzero and bounded by two, while their leaf-to-leaf difference is three.

So the upper cap reduces to the two global extrema, but nonvanishing is not
an extrema statistic.  It is the exact avoidance of all vertex-pair path
zeros.

THM-2783 explains what binary and ternary do here.  Assigning distinct powers
of two to the edges makes every nonempty signed path sum nonzero, regardless
of `T`, `pi`, or the signs: the largest power dominates all smaller powers.
Ternary powers reconstruct the full signed path state.  Both solve a
strictly stronger coding problem and pay exponential span.  Graceful
labeling instead seeks path-restricted null avoidance and discrepancy at
most `m` inside the dense alphabet `{1,...,m}`.  This is why the radix
theorems do not close Graceful Tree.

## 3. The consecutive-ones gap matrix

Now fix a graceful labeling

```text
f:V(T)->{0,...,m}.                                      (11)
```

For `1<=d<=m`, let `e_d={u_d,v_d}` be the unique edge whose difference is
`d`.  For a prospective vacant label

```text
t in {0,...,m+1},                                       (12)
```

define

```text
C_f(d,t)=1
 iff min(f(u_d),f(v_d))<t<=max(f(u_d),f(v_d)).           (13)
```

Each row has consecutive ones and exact row sum

```text
sum_t C_f(d,t)=d.                                       (14)
```

Let

```text
iota_t(a)=a                 if a<t,
           a+1              if a>=t.                    (15)
```

Then the lifted difference of `e_d` is exactly

```text
|iota_t(f(u_d))-iota_t(f(v_d))|=d+C_f(d,t).             (16)
```

The column in `(13)` has a second exact interpretation.  Since `f` is a
permutation, put

```text
S_t=f^-1({0,...,t-1}).                                  (17)
```

Then

```text
{e_d:C_f(d,t)=1}=delta_T(S_t),                           (18)
```

the tree cut of the first `t` vertices in label order.  Thus the new leaf
problem is a compatibility question between an integral vertex order and a
`Z/2` cut/coboundary.  No claim is made that an arbitrary consecutive-ones
matrix with row sums `(1,...,m)` is realizable by a tree labeling.

## 4. Exact gap-tail leaf insertion

Mark a vertex `v` of `T`, add a new leaf `z` adjacent to it, lift every old
label by `(15)`, and give `z` the vacant label `t`.  Let

```text
k=|t-iota_t(f(v))|.                                     (19)
```

Then the enlarged labeling is graceful if and only if

```text
C_f(d,t)=0 for d<k,
C_f(d,t)=1 for d>=k.                                    (20)
```

Equivalently,

```text
delta_T(S_t)={e_k,e_(k+1),...,e_m}.                      (21)
```

To prove this, write `c_d=C_f(d,t)`.  The old lifted differences are
`d+c_d`.  If `c_d=1` and `c_(d+1)=0`, then two old edges both have difference
`d+1`.  Therefore they are distinct only if the zero-one word
`c_1...c_m` is a suffix `0^(k-1)1^(m-k+1)`, including the all-zero case
`k=m+1`.  In that case the old differences are exactly

```text
{1,...,m+1}\{k}.                                        (22)
```

The new edge fills the missing difference precisely when `(19)` has the
same `k`.  This proves both directions.

There is a useful rank form.  If `v in S_t`, condition `(19)` is

```text
f(v)=t-k;                                               (23)
```

if `v notin S_t`, it is

```text
f(v)=t+k-1.                                             (24)
```

Thus a valid insertion is a Ferrers/suffix column together with one marked
rank on the appropriate side of its threshold cut.

The largest-new-edge boundary `k=m+1` is exact.  Equation `(20)` makes the
cut empty.  Since `T` is connected, `t` must be `0` or `m+1`; equation `(19)`
then says that the marked parent has label `m` or `0`, respectively.  Hence

```text
new edge gets m+1
iff the marked parent is an old extreme vertex.          (25)
```

This is the familiar extreme-parent induction, isolated as only the two
outer columns of the more general gap matrix.

## 5. Sharp extension and descent boundaries

The gap columns add genuine flexibility beyond `(25)`.

Take the six-vertex tree

```text
E={(1,0),(1,2),(1,3),(0,4),(4,5)}                       (26)
```

and mark vertex `0`.  It has exactly twelve graceful labelings, and in none
of them does vertex `0` receive label zero or five.  Thus the extreme-parent
route is impossible.  Nevertheless

```text
f=(1,4,2,3,5,0),       t=3,       k=2,
C_f(:,3)=01111                                           (27)
```

satisfies `(19)--(20)`.  The lifted seven-vertex labels are

```text
(1,5,2,4,6,0,3),                                       (28)
```

and their edge differences are `(4,3,1,5,6,2)`.

The quantifier over the graceful labeling is also essential.  On the star
`K_(1,3)` with centre label zero and leaf labels one, two, and three, the
leaf carrying label two has no valid gap-tail insertion.  Other graceful
labelings of the same marked tree do work.

Finally, graceful deletion does not automatically invert gap insertion.
The path with edges

```text
{(1,0),(1,2),(0,3),(3,4)}                               (29)
```

has graceful labels `(0,2,3,4,1)`.  Delete leaf `2`, whose label is three,
and compress labels above the gap.  The remaining edge differences become

```text
(2,2,3),                                                (30)
```

not `(1,2,3)`.  Therefore one cannot prove the marked extension statement
by taking an arbitrary graceful labeling of the larger tree and deleting
its leaf.

## 6. Finite marked census and the exact open induction target

The exact companion checks every nonisomorphic tree through nine vertices.
For every tree shape and every marked attachment vertex it finds some
graceful labeling and some pair `(t,k)` satisfying `(19)--(20)`:

```text
n:                    2  3  4  5  6   7   8   9
tree shapes:          1  1  2  3  6  11  23  47
marked vertices:      2  3  8 15 36  77 184 423.        (31)
```

This motivates, but does not prove:

> **Marked gap-tail extension conjecture.** For every finite tree `T` and
> every marked vertex `v`, some graceful labeling `f` and some vacant label
> `t` satisfy `(19)--(20)`.

The conjecture would imply the Graceful Tree Conjecture by leaf induction.
Given a larger tree, delete one leaf, apply the marked statement to its
parent in the smaller tree, and reinstall the leaf by Section 4.  The
six-vertex example `(26)--(28)` shows why restricting this induction to
extreme-parent labelings is already too strong.

This is a precise new target, not a proof beyond `(31)`.  It asks for one
suffix column in the consecutive-ones matrix of a suitably chosen graceful
labeling, coupled to the marked rank `(23)` or `(24)`.

## 7. Exact verification

Run

```bash
python 04-computation/signed_path_gap_tail_graceful_thm2787.py
python -O 04-computation/signed_path_gap_tail_graceful_thm2787.py
```

The exact companion uses explicit exceptions and no truth-bearing Python
assertions.  It checks the signed-path equivalence on all `48,794` signed
edge-permutation assignments over all labelled trees through five vertices.
It checks `(13)--(22)` on `2,650` graceful labelings of every nonisomorphic
tree through seven vertices: `15,314` consecutive-ones rows, `20,614`
threshold cuts, and `140,516` marked gap cases.  It verifies all `748`
marked vertices in `(31)`, freezes their first-witness bank by SHA-256, and
checks every hostile in Sections 2 and 5.  Normal and optimized runs
byte-match the stored transcript.

```text
PROVED HERE:              signed path-sum characterization of gracefulness;
                          exact A_m/B_m lattice-orbit intersection;
                          separate nonzero and upper-cap hostiles;
                          consecutive-ones gap matrix and threshold cut;
                          necessary-and-sufficient gap-tail insertion rule;
                          exact extreme-parent boundary;
                          six-vertex non-extreme positive repair;
                          fixed-label and naive-descent hostiles;
                          marked gap-tail census through nine vertices.

NOT PROVED:               the marked gap-tail conjecture for all trees;
                          the Graceful Tree Conjecture;
                          that an arbitrary consecutive-ones matrix is a
                          realizable tree gap matrix;
                          an improved THM-2765 range from Nullstellensatz;
                          a modular-group, tournament, LRC, Keller, JC(2),
                          or DC(2) consequence.                            (32)
```

QED.
