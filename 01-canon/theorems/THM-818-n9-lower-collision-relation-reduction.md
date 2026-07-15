---
id: THM-818
title: The n=9 lower-codec collision search is an exact three-way Cech join over a 5,997,416-row n=8 equality relation
status: PROVED REDUCTION + FINITE-EXACT R8 CENSUS; exact n=9 join completed by THM-828
source: codex-2026-07-15-S13
depends_on: [THM-553, THM-796, THM-801, THM-809]
related: [THM-811, THM-814, THM-825, THM-828, HYP-6880]
verification:
  - 04-computation/mobius_cech_n9_relation_preflight_codex_S13.py
  - 05-knowledge/results/mobius_cech_n9_relation_preflight_codex_S13.out
  - 05-knowledge/results/mobius_cech_n9_relation_preflight_codex_S13.json
  - 04-computation/mobius_cech_join_regression_codex_S13.py
  - 05-knowledge/results/mobius_cech_join_regression_codex_S13.out
---

# THM-818 — the `n=9` frontier is a relation join, not a line scan

Let `X_8` be the `2^21` oriented staircase tilings at size eight.  Write
`kappa` for bitwise complementation, `sigma` for staircase reflection,
`chi_8(x)` for the blue/black colour of the complement line through `x`, and
`nu_8(x)` for the merged isomorphism-class node containing `x`.  Thus `nu_8`
already identifies a tournament class with its `sigma`-converse class.  Define

```text
H_8(x) = (nu_8(x), nu_8(kappa x), chi_8(x)),
R_8    = {(x,y) in X_8^2 : H_8(x)=H_8(y)}.                 (1)
```

Then the search for distinct oriented `n=9` tilings with equal lower face-`H_8`
data is exactly a three-way equijoin of copies of `R_8` along their pairwise
`n=7` overlaps.  It is not necessary to enumerate and sort all `2^27`
apex-zero complement-line representatives (or all `2^28` oriented tilings).

The equality relation in (1) has exactly

```text
|R_8| = 5,997,416,       |R_8 \ Delta| = 3,900,264.       (2)
```

Of its rows, `5,216` are blue and `5,992,200` are black.  This theorem proves
the reduction and relation census.  THM-828 subsequently executes the
complete join and classifies its 58 final doubletons.

## Exact three-face reduction

Use the three lower B3 faces `A,B,C` of an oriented size-nine tiling.
Each pair has a shared size-seven subface.  For a relation row
`r=(x,y)`, let

```text
rho_AB^A(r), rho_AB^B(r), rho_AC^A(r), rho_AC^C(r),
rho_BC^B(r), rho_BC^C(r)
```

denote the ordered pair of literal subfaces obtained by applying each
role-specific intersection map to both `x` and `y`.  The superscript matters:
the two face coordinate systems identify the same global cells through
different local positions.  These are ordered 30-bit pairs (two
literal size-seven masks), never unordered lines or node ranks; forgetting
their order would reintroduce the phase holonomy that the join is meant to
control.  The literal B3 cover has the elementary gluing property:

> A compatible triple of labelled `A,B,C` faces has a unique upper tiling.

Indeed, every upper bit occurs in at least one face and the overlap equalities
make repeated occurrences agree; equality on the triple overlap is then
automatic.  Consequently ordered pairs `(u,v)` of
upper tilings with equal lower `H_8` data are in bijection with triples

```text
r_A,r_B,r_C in R_8

rho_AB^A(r_A)=rho_AB^B(r_B),
rho_AC^A(r_A)=rho_AC^C(r_C),
rho_BC^B(r_B)=rho_BC^C(r_C).                              (3)
```

The inverse sends `(u,v)` to its three ordered face pairs.  This proves both
surjectivity and injectivity of the Cech gluing description; no tournament
canonicalization at size nine is hidden in it.

Thus the raw join is exactly the collision relation for the three face-`H_8`
keys.  To recover the exact lower-codec collision search from (3), perform the
following filters in this order:

1. remove the all-diagonal locus `r_A,r_B,r_C in Delta`, which is precisely
   `u=v`;
2. retain pairs whose reconstructed upper apex bits are **both** zero
   (equivalently the two marked B-face apex coordinates vanish), because the
   codec domain consists of canonical apex-zero representatives.  The exact
   relation does contain `285,244` rows in each mixed-apex orientation, but
   those rows lie outside this canonical comparison domain;
3. require equality of the upper colour bit `chi_9(u)=chi_9(v)`.  Together
   with the three face-colour coordinates already equal in `H_8`, this is
   exactly equality of the full `UABC` word;
4. require equality of the full raw `S2` word (`tau=3,...,8` and fixed
   `tau=9`), or refine it by positional
   moments when studying the repair.

Rows with only one or two diagonal faces must not be discarded: they can glue
to a distinct upper pair.  Likewise, equality of the three face colours does
not replace the upper `UABC` test.  These two cautions are the main orientation
and quotient boundaries of the reduction.

## Exact `R_8` census

The existing exact size-eight classifier gives `6,880` ordinary tournament
isomorphism classes.  Grid reflection realizes tournament converse, so
ranking

```text
(min(class(x),class(sigma x)), max(class(x),class(sigma x)))
```

gives exactly `3,528` merged nodes.  No second converse classifier is needed.
Grouping the `2,097,152` oriented tilings by (1) gives `876,512` nonempty
fibres, with maximum fibre size `26`.  If `a_j` is the number of fibres of
size `j`, the complete histogram is

```text
j : a_j
1 : 3,116       2 : 733,838       3 : 36          4 : 115,633
6 : 17,680      8 : 4,209        10 : 1,034       12 : 457
14: 242        16 : 130          18 : 66          20 : 45
22: 20         24 : 4            26 : 2.
```

The two independent checksum identities

```text
sum_j j a_j   = 2,097,152,
sum_j j^2 a_j = 5,997,416                                      (4)
```

prove the line and relation totals.  In particular the observed relation is
only `2.859790802...` rows per oriented tiling, far smaller than the ambient
square.

Packed 64-bit rows require `47,979,328` bytes.  The conservative layouts
used for planning require `239,896,640` bytes for rows plus two aligned
16-byte projection indexes, and `335,855,296` bytes with three indexes.  The
exact census therefore passes the declared single-machine gate
`|R_8| <= 30,000,000` comfortably.

## One-level-down exact replay

The same ordered, role-specific join was executed one level lower without
using THM-809's collision list as input.  The `n=7` observation kernel has
`113,632` rows, `80,864` off-diagonal.  Its exact three-way join gives

```text
nontrivial ordered triangles                 1,672
after requiring both upper apex bits zero      836
after requiring equal upper colour              836
canonical unordered upper pairs                 418.
```

The final `418` is exactly THM-809's independently obtained base-collision
count.  The fact that upper colour removes nothing at this size is a
regression result, not a reason to omit that logically necessary filter at
`n=9`.  Runtime is about 2.2 seconds with maximum RSS about 101 MB in the
stored Python implementation.

## Why `n=9` is the last count-plus-first-moment size

The exact lower key at size nine can be packed in 98 bits:

```text
six n=8 merged-node ranks        6*12 = 72 bits
UABC colour                                  4 bits
S2 mixed radix                               22 bits
total                                        98 bits.
```

On canonical apex-zero representatives, the raw `S2` radix is

```text
4^2 * 10^2 * 20^2 * 4 = 2,560,000 < 2^22.                 (5)
```

The last factor is four, not five: the fixed layer has four positions but
the apex bit is already fixed to zero, so its one-count ranges only from zero
through three.  The earlier `*5=3,200,000` preflight was a safe unoriented
upper bound; it was not the exact canonical radix.

More importantly, `n=9` is the last size at which adjoining the first
positional moment to every `S2` state is unconditionally literal-exact.
Every nonfixed layer has at most three positions.  The fixed layer has four,
but apex-zero orientation leaves only three free.  A subset of at most three
ordered positions is determined by its cardinality and position sum.  At
`n=10`, the two subsets `{0,3}` and `{1,2}` first have equal count and sum, so
second moments or a literal layer word can become necessary.

This separates two questions which were previously conflated:

- whether raw counts already distinguish every `n=9` collision;
- whether a small, guaranteed-exact positional repair exists if they do not.

THM-828 answers the first: raw counts leave exactly 58 disjoint reflection
doubletons.  It answers the second more sharply: one antisymmetric chirality
bit, smaller than the full first-moment sidecar, repairs all 58.

## Structural interpretation

The crucial object is not merely the set of metagraph nodes.  It is the
kernel pair `R_8` of the observation map `H_8`, equipped with three literal
face projections.  Nodes say which observations agree; relation rows retain
*which tilings* agree; the projections retain whether three such agreements
can coexist as the boundary of one larger tiling.  This is the first carrier
in the present recursion that simultaneously records nodes, their tiling
fibres, and composable edge witnesses.

In database language (3) is a sparse three-way join.  In sheaf language it is
the Cech equalizer, or zero-cocycle compatibility condition, for the B3
cover; the phase one-cocycles are exactly what ordered literal overlaps keep
from reappearing.  In tournament
language its rows are parallel coloured complement edges whose endpoint
classes agree, while the overlap projections remember marked Hamiltonian-path
sections that node isomorphism erased.  These are descriptions of the same
finite object, not analogies.

For the required Tournament Analysis, the verifier treats four storage
carriers as vertices.  The pairwise observable is how many declared join
obligations each representation can execute, with switches for retention and
retention per logarithmic storage cost; the reported fingerprints are a
planning diagnostic, not evidence for (3).  The mathematical tournament
objects remain the marked tilings and complement edges.

## Preservation boundary and next exact target

`R_8` preserves both members of every `H_8` fibre, orientation, line colour,
and the literal face intersections needed by (3).  It deliberately does not
preserve metric LRC gaps, owner labels, wall-crossing chronology, affine
carry, or future lift/deletion behaviour.  Passing from rows to fibre sizes
would also destroy the gluing data and is therefore illegal for the next
step.

The challenged assumption is that tournaments or metagraph nodes are the
only useful vertices.  Here the correct vertices for recursion are equality
witnesses `(x,y)`; nodes are merely their quotient labels.  This reduction is
finite and exact, but it is neither an all-size Markov theorem nor a proof of
the fourteen-runner case.

THM-828 carries out this target by a stronger difference-syndrome
factorization: raw-S2 equality restricts `u xor v` to a 14-dimensional code,
and the `A/C` faces then glue the upper pair uniquely.  The complete join has
9,540 overlap matches, 636 B-compatible candidates, and 58 final pairs.  Its
survivor differences span a punctured rank-four reflection-defect cube. ∎
