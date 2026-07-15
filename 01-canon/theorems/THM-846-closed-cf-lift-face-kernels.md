---
id: THM-846
title: Closed continued-fraction lift/face composites have distinct apex, seam, theta-sheet, and node-operation kernels
status: PROVED EXACT COORDINATE ALGEBRA + FINITE-EXACT 58-CELL OPERATION-KERNEL CLASSIFICATION
source: codex-2026-07-15-S13f
depends_on: [THM-828, THM-834, THM-838, THM-840, THM-842]
related: [THM-813, THM-829, THM-843, HYP-6880]
verification:
  - 04-computation/n9_closed_cf_lift_face_kernels_codex_S13f.py
  - 05-knowledge/results/n9_closed_cf_lift_face_kernels_codex_S13f.out
  - 05-knowledge/results/n9_closed_cf_lift_face_kernels_codex_S13f.json
---

# THM-846 — the first closed CF lift/face word

Let `Phi:X_9->X_10` be THM-838's centered continued-fraction coordinate
copy, and let `A_10,B_10,C_10:X_10->X_9` be the three marked-path faces of
THM-842.  Define the three closed words

```text
F_R=R_10 Phi,                    R in {A,B,C}.             (1)
```

Their complete coordinate algebra is explicit.  The ordered endpoint word
`(F_A,F_C)` loses exactly the common apex bit.  The middle word `F_B` loses
exactly one different seam bit.  On THM-828's 58 collision cells, forgetting
the order of the endpoint word loses a third bit, the theta sheet.  These are
three distinct equivalence relations.  At the merged-node level the bare
direct and partner marginals have still different operation kernels, whereas
the endpoint-coupled complement-line coordinate `bar P` is operation-safe on
the entire finite bank.

This is the first exact closed continued-fraction lift/delete calculation in
the n=9 defect programme.  It proves no closure for arbitrary CF words.

## 1. Literal coordinate theorem

Write `rho:T_10->T_9` for the source-coordinate map of `Phi`; thus
`Phi(x)_(a,b)=x_(rho(a,b))`.  Then

```text
F_A(x)_(a,b)=x_(rho(a,b)),
F_B(x)_(a,b)=x_(rho(a+1,b)),
F_C(x)_(a,b)=x_(rho(a+1,b+1)).                           (2)
```

The coordinate ranks are

```text
rank(F_A),rank(F_B),rank(F_C) = 21,27,21.                (3)
```

More precisely, `F_A` misses exactly the seven top-row coordinates

```text
(9,1),(9,2),...,(9,7),
```

and `F_C` misses exactly the seven first-column coordinates

```text
(3,1),(4,1),...,(9,1).
```

Their used-coordinate union therefore has size 27 and misses only the apex
`a=(9,1)`.  Choosing the first occurrence of every other source coordinate
gives an explicit left inverse

```text
x <-> (F_A(x),F_C(x),x_a).                              (4)
```

All 116 endpoints in THM-828's canonical representative bank have `x_a=0`.
Hence `(F_A,F_C)` is literally injective on that apex-zero representative
slice without an extra token, but not on the full 28-dimensional cube.  This
corrects the plausible but false
rank-28 endpoint-union guess.

The middle composite is even simpler:

```text
F_B(x)_(6,4)=x_(7,3),
F_B(x)_(a,b)=x_(a,b) for (a,b)!=(6,4).                  (5)
```

Consequently `F_B^2=F_B`.  Its image is the hyperplane
`y_(6,4)=y_(7,3)`, and its unique lost literal coordinate is

```text
s_B(x)=x_(6,4) xor x_(7,3).                             (6)
```

The inverse is explicit:

```text
x_(6,4)=F_B(x)_(7,3) xor s_B(x),                        (7)
```

with every other coordinate read directly from `F_B(x)`.  Thus
`(F_B,s_B)` is injective on the whole cube.  This seam is not the apex token
of (4).

For reference, `Phi` itself has rank 28 inside the 36-dimensional target
cube.  Its image has codimension eight, with equality on the two target
copies of each source coordinate

```text
(3,1),(4,2),(7,3),(5,3),(6,4),(7,5),(8,6),(9,7).        (8)
```

The stored JSON gives the eight target-coordinate pairs exactly.

## 2. Symmetry and the rank-four defect

Every `F_R` commutes with complement.  Reflection exchanges the endpoint
words and fixes the middle word:

```text
F_A sigma = sigma F_C,       F_C sigma = sigma F_A,
F_B sigma = sigma F_B.                                      (9)
```

Moreover `F_B theta=theta F_B` for THM-842's apex-relative antipode
`theta`.  The seam `s_B` is theta-invariant.

The four THM-832 defect basis words and all eleven occupied defect sectors
have seam zero and are fixed pointwise by `F_B`.  Under each of A, B, and C,
the four basis images still have rank four and the eleven sectors still have
eleven distinct images.  Thus none of the three composites loses the linear
defect coordinate.  All losses below occur in the nonlinear decorated
fibres.

## 3. The exact 58-cell decks

Let `u_i,v_i`, `i=1,...,58`, be THM-828's ordered source pairs.  Let `Q` be
the literal reflection-orbit line coordinate and let

```text
bar P(x)={nu(x),nu(kappa x)}                             (10)
```

be the unordered pair of converse-merged endpoint nodes.  On the chosen
representatives, every individual role has 58 singleton values for both Q
and `bar P`:

```text
{Q(F_R u_i):i}                 58 x 1,
{bar P(F_R u_i):i}             58 x 1,       R=A,B,C.   (11)
```

But rolewise descent across each source collision pair is

```text
Q(F_R u_i)=Q(F_R v_i)          for 0,58,0 cells,
bar P(F_R u_i)=bar P(F_R v_i)  for 0,58,0 cells,          (12)
```

in A,B,C order.  Thus B is the only closed word constant on the old
58-cell relation.  Individual injectivity in (11) and descent in (12) answer
different questions.

Endpoint order is precisely the next boundary:

```text
ordered (Q(F_A u),Q(F_C u))              58 x 1,
ordered (bar P(F_A u),bar P(F_C u))      58 x 1,
unordered {Q(F_A u),Q(F_C u)}            29 x 2,
unordered {bar P(F_A u),bar P(F_C u)}    29 x 2.          (13)
```

The two unordered partitions are exactly the 29 free theta orbits.  The seam
is constant on each orbit—thirteen orbits have seam zero and sixteen have
seam one—so it cannot recover endpoint order.  The ordered endpoint word is
the theta-sheet coordinate.  If the sheet is oriented by the repository's
canonical integer order on the two source Q cells, the B-Q images preserve
that orientation in all 29 orbits; A/C role order induces the opposite
convention.  Adjoining `bar P(F_B u)` to the unordered endpoint deck gives
58 singletons:

```text
({bar P(F_A u),bar P(F_C u)},bar P(F_B u))   58 x 1.      (14)
```

This corrects a second tempting false summary: the unordered Q endpoint deck
has 29 doubletons, not 58 singletons.  Only its ordered form is injective.

## 4. Operation kernels at merged nodes

THM-840's criterion must be applied separately to every proposed source and
target observation.  The following is a one-step test on the 58 canonical
representatives.  Here “direct” means the representative chart and “partner”
its complemented chart; this convention can swap the two atlas labels used
elsewhere.  The bare source direct node has fibre profile

```text
54 cells = 51 x 1 + 2 x 2 + 1 x 3,
```

and one of those source cells splits after `F_B`.  The bare source partner
node has

```text
53 cells = 48 x 1 + 5 x 2,
```

and two source cells split.  Hence neither marginal is an `F_B` congruence.

Adjoining the input seam repairs the direct source operation kernel but not
the partner source operation kernel:

```text
source observation                 split source cells after F_B
direct node                                      1
partner node                                     2
(direct node,input seam)                         0
(partner node,input seam)                        1
bar P                                            0
literal Q                                        0.                       (15)
```

The two last zeroes are operation-safe on this finite bank.  In fact source Q
and source `bar P` already have 58 singleton values, so these two kernel tests
pass by source injectivity.  Thus `bar P(u)` determines `bar P(F_Bu)` on this
bank, and Q likewise descends; no global quotient functor follows.

This statement must not be confused with static target injectivity.  The
direct and partner `F_B` node marginals each have

```text
55 cells = 52 x 1 + 3 x 2.                       (16)
```

Even appending the *input* seam to either target marginal leaves exactly the
same profile.  By contrast the coupled target `bar P(F_Bu)` has 58
singletons and no loops.  Thus the seam repairs literal inversion in (7) and
one named source congruence in (15), but it is not a universal node-fibre
separator.  Coupling the complement partner is essential here.

Across both endpoints of all 58 input pairs, the B image contains 232 masks
in 105 locally classified merged nodes.  The exact classifier examines 592
masks and 276 local merged nodes; no hashed node label enters any equality
claim.

## 5. What information must be preserved

The calculation exposes a typed state rather than one scalar address:

```text
literal inversion under A/C        apex token;
literal inversion under B          B seam;
unordered-to-ordered endpoint deck theta sheet;
merged-node continuation under B   coupled complement partner;
linear defect continuation          four-bit sector coordinate.           (17)
```

No one entry in (17) substitutes for every other entry.  This is a concrete
Myhill--Nerode warning for the recursive metagraph: minimize states only
relative to a named operation and a named output predicate.  Static recovery,
descent of a quotient, target injectivity, and future continuation are four
different tests.

The challenged vertex assumption in Tournament Analysis is therefore useful
here.  The seven tournament vertices are information carriers—defect sector,
seam, direct-node marginals, endpoint deck, coupled node pair, and Q deck—not
runners or tournament vertices.  The pairwise observable is separated source
pairs; switching from raw retention to retention per logarithmic cost flips
20 of 21 edges.  Both gauges are transitive, have no directed triangle,
singleton SCCs, and one tie Hamiltonian path.  This ranks carrier economy; it
does not turn the carrier tournament into the proof object.

The theorem preserves literal coordinate maps, complement/reflection action,
the stated theta action, the rank-four finite defect bank, literal Q, and
locally exact merged-node equality.  It destroys LRC gaps, owners, walls,
metric loneliness, and all unaudited CF words.  In particular, it proves no
implication for LRC(14). ∎
