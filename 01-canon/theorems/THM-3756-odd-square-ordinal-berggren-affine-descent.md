---
id: THM-3756
title: "Odd-square ordinal chart and affine descent for the Berggren tree"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a primitive
  Pythagorean triple with ordered even leg B, the odd square roots of B+C and
  C-B define two ordinals (r,s).  These pairs are exactly r>s>=1 with
  gcd(2r-1,2s-1)=1; the outer-rank fibre has size phi(2r-1)/2; and the three
  Berggren children are affine in (r,s), while the unique inverse descent is
  piecewise affine on three disjoint cones.
  Without the gcd filter, the full triangular cone is a forest of
  odd-square-scaled Berggren trees indexed by odd-root content.  The algebraic
  proof was independently audited, including a separate brute enumeration of
  all 792 primitive ordered triples with C<=5000 and zero converse failures.
source: opus / 2026-08-23
depends_on:
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
related:
  - THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses
  - THM-3499-regular-shortlex-languages-have-logarithmic-density
  - THM-3510-binary-shortlex-equal-level-count-log-density-boundary
script: 04-computation/odd_square_ordinal_berggren_affine_descent_thm3756.py
output: 05-knowledge/results/odd_square_ordinal_berggren_affine_descent_thm3756.out
script_sha256: 9fed96e425c77fbc4f263a7980cb8dcb2de25cce00375ecb532c94b3b0d261b0
output_sha256: 690f4020ba06f9d2ce8633bb48f2075d411177f76ba305ebd0f45a272f675e70
semantic_sha256: 5f631c9e700911ce515a7393bb8bc6c5a1551355441834c74544606e0922fc6f
independent_script: 04-computation/odd_square_ordinal_berggren_independent_audit_thm3756.py
independent_output: 05-knowledge/results/odd_square_ordinal_berggren_independent_audit_thm3756.out
independent_script_sha256: b5cffb683649ecf792805c06be42a191b90fc62fd97eaf55939560b1a93dcc7b
independent_output_sha256: a511a5c533afabd7d393b1f2fe79b593a3a9ccba1e6d970921fdfb6a6e1c66bb
hash_basis: raw LF bytes
---

# THM-3756 -- odd-square ordinals on the Berggren tree

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The displayed
integer algebra and strict descent prove the global claims.  The exact
companion passes identically in ordinary and optimized Python, and an
independent hostile audit accepted the statement, branch convention, inverse
cones, equality boundary, address domains, and scope after reproducing the
primitive converse by a separate enumeration through `C<=5000`.

No literature-priority or global-novelty claim is made.  The point is to type
an especially simple natural-number rank suggested by the two legs adjacent
to the hypotenuse, then identify exactly what its overlaps preserve and lose.

## 1. Inheritance pass and conventions

[THM-3333](THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md)
extends

```text
T(z)=z(z+1)/2
```

to every integer and uses it as a quadratic refinement on the Pythagorean
light cone.  [THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
identifies the consecutive-parameter U-spine, while
[THM-3357](THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md)
fixes the positive parameter convention and the three branch labels used
below.  [THM-3382](THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses.md)
is the canonical hostile: the same tree object has different density and
harmonic behaviour under time, depth, and heap injections into the naturals.

The ordered triple convention is load-bearing:

```text
(A,B,C)=(odd leg, even leg, hypotenuse).                 (1)
```

Without it, even the root gives `3+5=8` rather than a square.  The closest
normalization warning is MISTAKE-418: a raw parameter recurrence must be
primitive-normalized before a slope-dependent Berggren claim.  The canonical
hostile here occurs one step above the root:

```text
(5,12,13): B+C=25,          (15,8,17): B+C=25.           (2)
```

Thus the outer odd square supplies a useful rank, but not a node address.
The least-used sidecar is the second identity `C-B=(odd)^2`.

## 2. The integer triangular fold and its centered sidecar

For every `z in Z`, direct expansion gives

```text
T(z)=T(-z-1),
T(z+h)-T(z-h)=h(2z+1)                    (h in Z).       (3)
```

The first formula is an exact two-to-one fold: if `n>=0`, then

```text
T(z)=T(n)  iff  z=n or z=-n-1.                           (4)
```

Indeed, the difference factors as `(z-n)(z+n+1)/2`.  For every nonzero
`h`, the second formula restores the oriented integer `z`; it is an
antisymmetric sidecar for the fold.  At the scale in the motivating identity,

```text
[T((r-1)+2)-T((r-1)-2)]/2=2r-1.                         (5)
```

Consequently the `r`th positive odd square is

```text
(2r-1)^2=[T(r+1)-T(r-3)]^2/4.                           (6)
```

Equations (3)--(6) are identities on integers, not merely on the
nonnegative triangular sequence.  Their role below is coordinate extraction;
shared triangular syntax alone supplies no cross-problem theorem under
MISTAKE-222.

## 3. The lossless two-ordinal chart

Let

```text
Omega={(r,s) in Z^2 : r>s>=1 and gcd(2r-1,2s-1)=1}.     (7)
```

For `(r,s) in Omega`, put

```text
q=2r-1,                         d=2s-1,                  (8)
Xi(r,s)=(qd, (q^2-d^2)/2, (q^2+d^2)/2).                 (9)
```

### Theorem 3.1 (odd-square ordinal bijection)

`Xi` is a bijection from `Omega` to the positive primitive Pythagorean
triples in the ordered convention (1).  Its inverse is

```text
r=(sqrt(B+C)+1)/2,               s=(sqrt(C-B)+1)/2.     (10)
```

In particular,

```text
B+C=(2r-1)^2,                    C-B=(2s-1)^2.           (11)
```

*Proof.*  If `q>d` are coprime positive odd integers, (9) is integral and

```text
(qd)^2+((q^2-d^2)/2)^2=((q^2+d^2)/2)^2.                (12)
```

An odd prime dividing all three coordinates would divide both `q^2` and
`d^2`; two cannot divide the odd first coordinate.  Hence the triple is
primitive.

Conversely, let `(A,B,C)` be primitive with `B` even.  Then `A` is odd, since
otherwise the two legs share a factor two, and `C` is odd by the Pythagorean
identity.  The coordinates are pairwise coprime: a prime common to any two
would divide the third and contradict primitivity.  The two positive odd
integers `C+B` and `C-B` are coprime: a common prime would divide both `2C`
and `2B`, hence (being odd) both `C` and `B`.  Their product is `A^2`, so each
is a square,

```text
C+B=q^2,                         C-B=d^2,                (13)
```

for coprime positive odd `q>d`.  Solving (13) gives (9), while writing
`q=2r-1,d=2s-1` gives the unique pair (10) in `Omega`.  QED

There is a useful ambient content law.  If the coprimality condition in (7)
is omitted and `g=gcd(2r-1,2s-1)`, then `g` is the odd-root content and the
triple in (9) has exact coordinate content `g^2`.
After dividing by `g^2`, (12) is the primitive chart for the coprime odd pair.
Section 6 proves more: every fixed-`g` locus is a complete copy of the
primitive Berggren tree, rooted at `g^2(3,4,5)`.  Thus rejected primitive
slots are structured scaled components, not numerical waste.

## 4. Coarse rank, overlap multiplicity, and three natural addresses

Project `Omega` to its first coordinate:

```text
pi(r,s)=r.                                               (14)
```

This is a surjection onto `{2,3,...}`.  The section `s=1` gives

```text
Xi(r,1)=(2r-1,2r(r-1),2r^2-2r+1),                      (15)
```

the consecutive-parameter U-spine of THM-3334.  Under the THM-3357 branch
labels below it is the repeated `L` branch.

### Theorem 4.1 (exact shell multiplicity)

For every `r>=2`,

```text
|pi^(-1)(r)|=phi(2r-1)/2.                               (16)
```

*Proof.*  Put `q=2r-1`.  The possible `d=2s-1` are precisely the odd reduced
residues with `1<=d<q`.  The involution `d -> q-d` pairs all reduced residues
and, because `q` is odd, places exactly one odd residue in each pair.  QED

Thus prime `2r-1` gives the maximal overlap `r-1`; the coarse rank has
unbounded fibres.  This is useful compression when the consumer asks only
for the `B+C` shell, and a destructive quotient otherwise.

There are three canonical ways to continue from the overlapping rank.

1. **Ambient triangular address.**  On every pair `r>s>=1`, define

   ```text
   a(r,s)=T(r-2)+s.                                     (17)
   ```

   The shell `r` is the consecutive interval
   `[T(r-2)+1,T(r-1)]`, so (17) bijects the full triangular cone with the
   positive naturals.  By integer reflection it can equally be written
   `a(r,s)=T(1-r)+s`.  Restricted to `Omega` it is injective but has holes.
   The first is `a(5,2)=8`: here `(q,d)=(9,3)` and (9) is
   `(27,36,45)=9(3,4,5)`.

2. **Selected-shell address.**  Rank the allowed inner ordinals by

   ```text
   j_r(s)=|{t:1<=t<=s, gcd(2t-1,2r-1)=1}|,              (18)
   A(r)=sum_(k=2)^(r-1) phi(2k-1)/2,
   c(r,s)=A(r)+j_r(s).                                  (19)
   ```

   Equations (16), (18), and (19) make `c:Omega -> N_{>0}` a bijection.
   This is the literal “which admissible choice in this shell?” ordinal.  It
   removes the primitivity holes, but the totative-selection step replaces
   the affine child law by a decoding operation.

3. **Ternary heap address.**  The inverse descent in Section 6 assigns a
   unique word in `{L,M,R}*`; composing it with THM-3382's heap recursion
   gives a bijection with `N_0` whose branch append is affine base three.  It
   preserves ancestry operations but not the order of the odd-square shells.

The addresses solve different problems.  Their loss ledger is:

| carrier | exact gain | information lost or made nonlocal |
|---|---|---|
| outer rank `r` | `B+C` shell; U-spine section | inner square, triple, depth, word, children |
| pair `(r,s)` | full primitive triple and affine branches | no scalar total order chosen |
| ambient `a(r,s)` | closed triangular shell packing | not onto on primitive nodes; decode before branching |
| selected `c(r,s)` | bijection of primitive nodes with `N` | holes/content and affine branch visibility |
| heap `h(word)` | branch append and ancestry | arithmetic shell order; analytic weights change |

## 5. Berggren children become affine ordinal maps

Use THM-3357's parameter convention

```text
0<m<n,
Psi(m,n)=(n^2-m^2,2mn,n^2+m^2),                        (20)
```

and its branches

```text
L(m,n)=(n,2n-m),
M(m,n)=(n,2n+m),
R(m,n)=(m,2m+n).                                        (21)
```

The ordinal and parameter charts are mutually inverse:

```text
m=r-s,                       n=r+s-1,                   (22)
r=(m+n+1)/2,                 s=(n-m+1)/2.               (23)
```

Substitution in (21) gives the three exact affine laws

```text
L(r,s)=(r+2s-1,s),
M(r,s)=(2r+s-1,r),
R(r,s)=(2r-s,r).                                        (24)
```

Equivalently, on the coprime odd roots `(q,d)` they are the linear maps

```text
L(q,d)=(q+2d,d),
M(q,d)=(2q+d,q),
R(q,d)=(2q-d,q).                                        (25)
```

The matrices in (25) have determinants `1,-1,1`, respectively.  They
preserve coprimality, oddness, and the positive chamber.  Therefore (24)
acts internally on `Omega` and agrees with the canonical Berggren children,
not merely with an abstract ternary labelling.

The outer child ranks satisfy

```text
r_L-r=2s-1,               r_M-r=r+s-1,       r_R-r=r-s,
r_L+r_R=r+r_M,
r_M-r_L=r-s>0,
r_M-r_R=2s-1>0,
r_L-r_R=3s-r-1.                                        (26)
```

Thus the middle child always has the largest coarse rank.  The two outer
children tie only when `r=3s-1`, which is primitive only at the root.  The
root therefore retains an honest tie: its `L` and `R` children are exactly
the two rank-three triples in (2).  This rank preorder is not THM-3357's
hypotenuse tournament; the latter distinguishes those two children.

Every child raises `r` by at least one.  Therefore a node of tree depth `h`
satisfies

```text
h<=r-2.                                                  (26a)
```

Equality requires an increment of one at every edge.  In (26), this happens
only for `L` on `s=1` or for `R` on `s=r-1`; hence equality holds exactly on
the two boundary rays `L^h` and `R^h`.  Rank is a well-founded grading, but
not depth.

## 6. Intrinsic cone descent and the full content forest

Let

```text
C={(r,s) in Z^2:r>s>=1},                                (27)
g(r,s)=gcd(2r-1,2s-1).
```

### Theorem 6.1 (affine Euclidean forest)

Every `(r,s) in C` off the line `r=3s-1` lies in exactly one of the following
three cones, with the displayed parent:

```text
r>=3s:             L-parent=(r-2s+1,s),
2s<=r<=3s-2:       M-parent=(s,r-2s+1),
s<r<=2s-1:         R-parent=(s,2s-r).                    (28)
```

Each parent is in `C`, has smaller first coordinate, has the same `g`, and
applying the named map in (24) returns `(r,s)`.  Repetition terminates at

```text
root_g=((3g+1)/2,(g+1)/2),                              (29)
```

where `g` is the original odd content.  Hence the full cone is the disjoint
union, over positive odd `g`, of ternary trees rooted at `root_g`.  Under
`Xi`, the `g`-component is exactly `g^2` times the primitive Berggren tree.
In particular, `Omega` is the `g=1` component rooted at `(2,1)`.

*Proof.*  The images of `L,M,R` in (24) satisfy respectively the three
disjoint inequalities in (28), and solving the affine formulas gives the
displayed parents.  On each cone the inverse is one of the unimodular affine
maps inverse to
(25), so it preserves the gcd of the two odd roots.  The inequalities make the new
coordinates positive with first greater than second, and the parent first
coordinate is strictly smaller.

Only the line `r=3s-1` is not covered by the three integer intervals.  There

```text
2r-1=3(2s-1),                                           (30)
```

so its content is `g=2s-1`; solving for `(r,s)` gives (29), and (9) gives
`Xi(root_g)=g^2(3,4,5)`.  Strict descent must terminate, content fixes its
unique terminal, and disjointness of the cones makes the word unique.  QED

For a primitive pair, `g=1`, so Theorem 6.1 independently gives the usual
unique Berggren word and proves that (24) generates every primitive triple
exactly once.  The line `r=3s-1` meets `Omega` only at `(2,1)`.

The coarse rank alone cannot run this descent.  The minimal collision (2)
already has different labelled child-rank signatures:

```text
(3,1) -> (4,6,5),                  (3,2) -> (6,7,4),     (31)
```

where the entries are `(r_L,r_M,r_R)`.  Rank also fails to determine tree
depth: rank four occurs at depth one at `(4,2)` and at depth two at `(4,1)`
and `(4,3)`.  The inner ordinal `s`, or an equivalent word sidecar, is
load-bearing.

## 7. Connection contract and frontier boundary

The actual connection is

```text
primitive ordered triple
  <-> coprime odd roots (sqrt(B+C),sqrt(C-B))
  <-> ordinal pair (r,s)
  <-> affine Berggren word by (28).                      (32)
```

It preserves the entire primitive triple, branch labels, and ancestry when
both ordinals are kept.  Projecting to `r` preserves exactly the chosen odd-
square shell and destroys the data witnessed in (31).  The native missing
coordinate is `s`; a scalar repair is either (17), (19), or the heap address, chosen by
the operation the consumer needs.

This does **not** supply an LRC owner, root in a relation plane, phase, clock,
endpoint, semantic arrival, or maximum.  Replacing a selected odd square by
its ordinal is sound bookkeeping, not an LRC reduction.  It does not turn
ties into a tournament, make density invariant under reindexing, identify
level masses with nodewise selectors, or transfer the triangular syntax to
JC/FC.  An arbitrary nonprimitive scale is also outside the two-square chart:
for example `2(3,4,5)=(6,8,10)` has `B+C=18`; Section 6 covers the odd-square
scales `g^2`, not every scaled triple.  THM-3333, THM-3357, MISTAKE-209,
MISTAKE-222, and THM-3382 remain the relevant stopping boundaries.

## 8. Reproduction and controls

Run

```bash
python 04-computation/odd_square_ordinal_berggren_affine_descent_thm3756.py
python -O 04-computation/odd_square_ordinal_berggren_affine_descent_thm3756.py
```

Both modes byte-match the LF-normalized stored transcript literally on every
supported platform.  The exact companion checks:

- (3) for `-4000<=z<=4000` and `1<=h<=9`; the factorization and exact fibre
  equivalence (4) for `0<=n<=600` across the same `z` range; and (5)--(6) for
  `1<=r<=600`;
- all `145,861` pairs in `Omega` with `2<=r<=600`, including the triple,
  inverse, primitivity, `phi/2` fibre, ambient and selected addresses;
- all `437,583` labelled child/parent round trips in that universe, both in
  ordinal and THM-3357 parameter coordinates;
- strict descent and word replay for every pair, with maximum audited length
  `598` on the consecutive section;
- the full `179,700`-pair triangular cone through `r=600`, including
  `539,100` forest branch round trips, exact `g^2` content, word replay, and
  all `200` component roots visible in that box;
- the complete ternary tree through depth ten, with `59,049` nodes on the
  last level and no same-level or earlier-level collision;
- the root line, rank-three overlap, differing child signatures, cross-depth
  rank collision, first ambient hole `(r,s)=(5,2)`, prime full fibres, and
  exact nonprimitive content hostile `(27,36,45)`.
- an independent implementation's direct census of all `792` primitive
  ordered triples with `C<=5000`, with zero failures of (10)--(13).

The semantic digest is
`5f631c9e700911ce515a7393bb8bc6c5a1551355441834c74544606e0922fc6f`.
The computation verifies consequence objects and conventions; the global
quantifiers follow from the factorization, totient pairing, affine formulas,
and strict descent above.
