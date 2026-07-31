---
id: THM-2847
title: "Q3/Q11 transverse endpoint edge and rank-one E3 realization horn"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The target-side
  endpoint minor, affine address edge, exact 20-cell E3-only horn, and split
  rank-one macro-truth mapping-cone algebra are exact.  The imported 449
  columns remain coefficient/support data, not a physical carry-to-response
  action.  No row or LRC(14) conclusion follows.
source: root/lrc-q3-q11-transverse-endpoint-horn-2026-07-28
depends_on:
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2813-affine-lift-transvection-and-projective-horn-decoder
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
  - THM-2829-q11-semantic-reselection-and-fine-ancestry-phase-obstruction
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
related:
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
script: 04-computation/lrc14_q3_q11_transverse_endpoint_horn_thm2847.py
output: 05-knowledge/results/lrc14_q3_q11_transverse_endpoint_horn_thm2847.out
script_sha256: 258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72
output_sha256: 155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d
hash_basis: LF-normalized bytes
---

# THM-2847 -- q3/q11 transverse endpoint edge and rank-one E3 horn

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
The exact companion passes normal, `-O`, and stored-output replay after the
hostile scope and encoding repairs.  LF-normalized SHA-256:

```text
script  258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72
output  155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d
```

This theorem does not exclude an LRC row or identify the semantic response
basepoint required by THM-2835.

## 1. The new synthesis

THM-2807/2813 give the relative affine-lift torsor

```text
A_t(y)=y+t 13^5(y-7 mod13)                    mod13^6.
```

On the residue-eight sheet this is simply

```text
A_t(y)=y+t 13^5.
```

The physical target translate with allocation label `q` has fine source
address

```text
n(q)=6716+beta(q)13^5,       beta(q)=2q mod13.
```

Since `2*7=1 mod13`, the physical allocation defect `q=7t` and the normal
decoder `beta(q)=t` are exact inverses.  More generally,

```text
n(q_1)=A_(beta(q_1)-beta(q_0)) n(q_0),
q_1-q_0=7(beta(q_1)-beta(q_0)).
```

Thus the target-label torsor and the transverse residue-eight address torsor
are not merely equinumerous: they are canonically semiconjugate through the
already proved normal label.

For `(q_0,q_1)=(3,11)`,

```text
n(3)=2234474,
n(11)=3348353,
n(11)-n(3)=3*13^5,
beta(11)-beta(3)=3,
11-3=8=7*3 mod13.
```

Hence the address edge underlying the new q3/q11 physical cospan is exactly
the transverse lift `A_3`.  This identifies the address labels; it does not
assert that `A_3` acts physically on the endpoint or semantic current.

## 2. A genuine target endpoint boundary lies on that edge

Exact reconstruction of the THM-2820/2829 bank gives:

```text
q3  full-six-factor cells:  56,
q11 full-six-factor cells:  63,
intersection:               42.
```

On each common cell, the two typed right-origin columns `(0,0),(12,0)`
and the two rows `q=3,11` have the same nonzero scalar `c` and matrix

```text
c [[1,0],
   [1,1]].
```

The determinant is `c^2`; after dividing by the common scalar the occupancy
minor is exactly one.  Across the full q-bank,

```text
R_00=c(delta_0+delta_3+delta_11),
R_12=c(delta_0+delta_11),
R_00-R_12=c delta_3.
```

The signed endpoint boundary is therefore the monomial `c z^3` over the
inherited endpoint field.  Put

```text
K_0=Q(zeta_1183).
```

In fact

```text
c=zeta_2366^2203-zeta_2366^65 in Q(zeta_1183).
```

Since `1183` is odd, `Q(zeta_2366)=Q(zeta_1183)`, and in the latter field
the same coefficient is

```text
c=zeta_1183^624-zeta_1183^510.
```

It does **not** descend to `Q(zeta_13)`: the Galois automorphism `u=27`
fixes `zeta_13` but sends these exponents to `(286,757)`.  Their difference
is

```text
x^624-x^510-x^286+x^757,
```

a nonzero polynomial of degree `757<phi(1183)=936`.  After normalizing
by `c`, the unit monomial `z^3` does lie in `Q[C_13]`; its inverse exponent
is ten and every one of its thirteen Fourier characters is nonzero.  This
removes relative two-origin ambiguity, but it does not normalize the overall
response phase.

This is the first physical target-relative endpoint boundary known to sit
on a proved transverse affine-lift edge while retaining all six native
factors.  It resolves the target ambiguity left by the older q3-versus-q11
either/or statement.

There is a further exact triangular clutch.  THM-2835's proved
`QB(source)->QA(target)` whole-cylinder column is zero at `q=0,3` and has
count `449` at `q=11`.  On the same 42 cells, adjoining this column to the
two endpoint-origin columns gives, in row order `q=0,3,11`,

```text
[[1,1,  0],
 [1,0,  0],
 [1,1,449]],
```

with determinant `-449`.  Restoring the common nonzero endpoint scalar
only multiplies this determinant by a nonzero power.  Thus q11's semantic
QA response gives a coefficient/support basepoint datum on the target normal
edge.  It is not merely inferred from the q label, but neither this column
nor its count supplies a physical QB response action or an equivariant
carry-to-response map.

The next obstruction can now be displayed as one rank collapse.  Add the
proved `QB(source)->QAB(target)` filler of count `449` at `q=7`.  Before
the macro gate, rows `q=0,3,11,7` and columns
`(origin00,origin12,QA,QAB)` form

```text
[[1,1,  0,  0],
 [1,0,  0,  0],
 [1,1,449,  0],
 [0,0,  0,449]],
```

with determinant `-449^2`.  This is a complete coefficient-support frame.
The physical four-row statement has an exact 20-cell scope:

```text
s in {0,3,8,9,12},  t in {5,6,9,10},  clock=1.
```

On precisely these cells q7 loses only `E3`; on the other 22 common
q3/q11 cells it also loses at least one of the ordinary q1/q2 factors.
THM-2829 proves that the physical `E3` support is exactly `{0,3,11}`.
Multiplying by the retained macro gate therefore deletes precisely the q7
filler row on this 20-cell horn and drops its determinant to zero.

The q7 interval itself is wholly in the complementary `E3` truth block.
Consequently the Boolean split `E3 direct-sum (not E3)` retains all four
rows of the 20-cell horn, but in two different `E3` macro-truth blocks.  The
outer `QA/QAB` semantic-word grading is separate.  A polarized macro-block
reference retaining and aligning those semantic columns, or a lawful
q11-to-q7 action that transports/recharts the `E3` factor, would glue that
coefficient frame.  Without one, taking its determinant as one physical
current conflates independently typed blocks.  This links the obstruction
to the older polarized cross-word sidecar (THM-2380/2460), while explaining
why another unpolarized support census cannot help.

More precisely, this is a rank-one mapping-cone boundary.  Projecting the
four-frame to the `E3` block leaves rank three and kernel

```text
K_0 * (0,0,0,1),
```

the pure `QAB` column.  Projection to the complementary block sends that
kernel generator to `449` times the q7 basis vector.  Hence

```text
0 -> K_0*QAB -> K_0^4 --P_E3 A--> K_0^{q0,q3,q11} -> 0
```

is exact, and the complementary block identifies its kernel with the
one-dimensional missing row.  The `E3`-graded direct sum already has full
rank four.  No additional support dimension is missing; what is absent is a
lawful physical functor that contracts or intertwines the two macro-truth
blocks while retaining the independent `QA/QAB` response/current semantics.
This is a sharper target than another determinant or support census.

Over the field `K_0` this short exact sequence splits; there is no ordinary
linear `Ext` or homology obstruction.  The nontrivial boundary is entirely
in the realization category: the algebraic splitting has not been lifted to
a typed positive/current morphism preserving both the macro truth and the
outer-word phase.  "Mapping cone" below always refers to that realization
boundary, not to a nonsplit vector-space extension.

## 3. Exact gain and exact stopping boundary

The result supplies:

- a physical target-side base vector `delta_3` in a regular `C_13` bank;
- the exact transverse address lift carrying it (`A_3`);
- the normal decoder identifying the lift index from the allocation label;
- a unit normalized two-origin determinant, hence no relative two-origin
  ambiguity.
- a nonzero three-row endpoint/semantic coefficient clutch selecting the q11
  QA support column.

It does **not** supply the THM-2835 basepoint

```text
physical carry delta_0  <->  canonical QB response.
```

The reason is now sharper.  The raw weighted source atom `I` is nonempty,
but the q3 and q11 target intervals are disjoint and the inherited
target-relative endpoint construction supplies neither a common
endpoint-source allocation coordinate nor a common `H` face.  The attached
QB/QA column is a one-way coefficient/support cospan rather than an action
on the response orbit.  The joint
`(target label, transverse lift)` state also lies on the graph

```text
q=7t,
```

so it has only thirteen states, not an independent `13 x 13` phase square.
The normalized two-dimensional Fourier transform of the *unweighted graph
indicator* is supported exactly on

```text
{(a,b):7a+b=0},
```

one thirteen-point annihilator line.  For an arbitrary weighted function on
the graph, the transform is constant on the thirteen cosets indexed by
`7a+b`; it can be nonzero at all 169 frequency pairs (a point mass is the
sharp hostile), but it still contains only thirteen distinct restricted
characters.  Full q-character support of the signed endpoint edge therefore
does not create a second independent normal character.
It chooses the target torsor leg but gives no second, semantically typed
source/response leg.

On the 20-cell horn, the cheapest decisive sidecar is a lawful contraction
of the rank-one `E3` mapping cone: a polarized `E3`/complement reference
retaining the three endpoint/QA/QAB columns, or a q11-to-q7 action that
transports/recharts `E3`.  Globally, a second debt remains: turn the imported
QB-to-QA
coefficient column into a physical carry-to-response
basepoint/intertwiner.  Neither is supplied by the present count.

## 4. Connection contract

```text
source:
  THM-2820 q3 full-factor target and THM-2829 q11 two-origin target;

target:
  THM-2807/2813 residue-eight transverse affine-lift torsor;

map:
  q |-> beta(q)=2q |-> n(q)=6716+beta(q)13^5;

preserved:
  clock one, 42 common semantic cells, all six native factors,
  weighted target value, right origin (0,0), and the q11 (12,0) edge;

destroyed:
  common endpoint-source allocation/common H face, independent target phase,
  physical QB response action, owner/current provenance, and global
  equivariant action;

needed sidecar:
  on the 20-cell horn, one polarized E3/complement macro-truth reference
  retaining/aligning the QA/QAB columns, or a q11-to-q7 action that
  transports/recharts E3; globally, one lawful carry-to-response
  basepoint/intertwiner for the imported QB-to-QA coefficient column;

cheapest falsification:
  exact fibre product of the 20 horn cells with the existing polarized
  cross-word sidecars, followed separately by the 42-cell QB
  carry-to-response basepoint test.  Empty support records the respective
  horn or global stopping boundary.
```

## 5. Exact evidence and audit boundary

The companion reconstructs the full exact integer interval geometry rather
than reading the displayed counts as input.  Its finite universe is

```text
13 allocation labels x 13 s-labels x 13 t-labels x 7 clocks,
two endpoint origins,
all 13 affine-address labels,
and the six literal native-factor containment predicates.
```

It pins the LF hashes of the THM-2806 allocation constructor,
THM-2782 endpoint constructor, and THM-2835 semantic-horn companion.
The value `449` is explicitly a **CITED/pinned THM-2835 theorem value**;
the present script does not pretend to recompute that ancestry census.
It independently reconstructs the 56/63/42/20/22 cell split, both endpoint
supports and their common nonzero coefficient in two exact finite-field
embeddings, the constructor-linked conductor/exponent descent, the
two- and three-row determinants, the four-frame and rank-one kernel, all
169 affine-torsor identities, and the exact graph-character quotient.

Run

```text
python 04-computation/lrc14_q3_q11_transverse_endpoint_horn_thm2847.py
python -O 04-computation/lrc14_q3_q11_transverse_endpoint_horn_thm2847.py
```

Both modes byte-match the stored output after LF normalization.  The proof
uses the exact interval reconstruction for the finite geometry and the
displayed algebra for the torsor, Fourier, determinant, field-descent, and
split mapping-cone conclusions.

The independent hostile audit reconstructed the endpoint field from the
literal q3 interval and caught the false provisional descent to
`Q(zeta_13)`; the repaired conductor-1183 statement above is exact.  It also
forced the 42-cell bank to split into the exact 20-cell E3-only horn and 22
cells with extra q1/q2 losses, distinguished arbitrary weighted graph
profiles from the unweighted graph indicator, separated the raw source atom
from the absent endpoint-source/common-`H` datum, and separated the `E3`
macro truth from the `QA/QAB` semantic-word grading.  Finally it replayed
normal, optimized, and UTF-8/LF stored output and accepted the theorem after
the output-encoding hash repair.  No remaining mathematical or evidence
defect was found.
