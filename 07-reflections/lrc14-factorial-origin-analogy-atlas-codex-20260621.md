---
source: codex-2026-06-21
status: exploratory analogy atlas for HYP-2718
tags:
  - lrc14
  - factorial-origin-atom
  - relation-height
  - analogy-atlas
  - tournament-analysis
---

# LRC14 Factorial-Origin Analogy Atlas

## Target

The current sharp proof obligation is the HYP-2718 origin-atom estimate.  For a
moderate-span balanced split row `E`, after the HYP-2717 carrier-relation
filter, prove

```text
|Q_0(E)| <= cap_k - ProductCover(E),
```

where THM-561 identifies

```text
Q_0 = ProductCover(E) - p0(E) = M_6.
```

So the task is not to dominate every miss-zeta coordinate, every missed-count
atom, or every product-vs-actual discrepancy.  It is to bound one protected
boundary atom after the right quotient.

## Search Principle

The repo-wide analogy search suggests a simple filter:

```text
keep analogies that preserve the hidden carrier until Q_0 is forced;
discard analogies that scalarize before the signed transform.
```

This is the same discipline that kept reappearing in older work under different
names: side-channel ledgers, packet addresses, endpoint ownership, finite
difference handoffs, transfer taxes, root packets, and operator rather than
scalar shadows.

## Closest Analogies

### Relation-Lattice Signed Wall

THM-532, THM-537, and THM-538 already warned that literal majorants and
coordinatewise `L1` bounds are too expensive in the seven-sector route.  The
right object was a signed relation-lattice correction with a finite low-height
ledger and a high-height tail.

S68 sharpened this for the present target.  In
`lrc14_top_character_relation_height_codex_s68.py`, offsets `(30,60)` have raw
gap `27` but primitive relation `(2,-1)` and relation height `3`.  The pressure
leader is arithmetic-phase sensitive, not raw-gap monotone.

Transfer to HYP-2718: the tail variable should be a carrier relation address,
not a distance or span scalar.  A useful ledger key is

```text
(relation vector, support, height, n.M, residue/phase, sector mask, atom t,
 x-cell state word, sign chamber).
```

### Bonferroni Transfer Tax

THM-558 proves the exact one-speed insertion identity

```text
Delta U4 = mass(1 -> 0) - mass(5 -> 4) - 4 mass(6 -> 5).
```

The only positive source is a genuine one-missed-sector closure; high-tail
transitions pay the tax.

Transfer to HYP-2718: search for a `Q_0` tax identity under an interpolation
between the shared-`x` actual carrier law and the independent product carrier
law.  Positive origin-atom creation should have to cross a named relation
packet, and those packets should be charged to tail or non-origin terms before
absolute values are taken.

### Generated Miss-Zeta Words

HYP-2698 says residual profiles must be represented as pointwise generated
miss-zeta words

```text
Z_context,x(A) = product_i Z_i,x(A)
```

before the final `x`-average.  The full residual cone is too large; generated
words are the object.

Transfer to HYP-2718: do not average away carrier provenance before the
factorial-origin transform.  Prove the singleton-product context first using
the explicit hit-count kernel

```text
g_r(t) = 7^-r sum_j (-1)^j binom(t,j)(7-j)^r,
```

then prove that coherent merging of singleton carriers does not increase the
origin-atom risk.

### Finite-Difference Handoff

The THM-438 handoff reflection moved from an opaque signed transform to genuine
counts `t(k,m)`, where a finite alternating sum becomes a root at zero and
falling-factorial divisibility.  THM-561 is the LRC version: miss-zeta moments
invert to missed-count atoms, and `Q_0` is the alternating boundary finite
difference.

Transfer to HYP-2718: look for a cofactor polynomial or divisibility statement
for the generated `Q_t` profiles.  The promising computation is not
`max |W_j|`; it is the factor left after the origin finite difference is
evaluated on admissible generated words.

### Reverse-Pair Cancellation

THM-560 proves a degree drop by pairing every odd cycle with its reverse before
scalarization.  The top-degree terms cancel because the paired leading
coefficients differ by `(-1)^(odd)`.

Transfer to HYP-2718: pair carrier Fourier modes `n` and `-n`, or pair sector
mirror/complement state words, before absolute values.  A plausible small lemma
is that the leading boundary contribution of `Q_0` is odd under this paired
mode involution and therefore cancels inside each generated relation packet.

### Cut/Cycle Seam

THM-559 makes `c3` an exact two-body Ising energy on the cut-space face, while
the higher OCF data live in cycle space.  That is the same structural warning
seen in the LRC cover route.

Transfer to HYP-2718: `ProductCover` is the cut-like main term, but `Q_0=M_6`
is a cycle-like top character.  A two-body scalar statistic can rank rows, but
should not be expected to prove the atom estimate without the many-body packet
labels.

### Operator Ledger Versus Scalar Shadow

HYP-2446's p-curvature toy atlas shows scalar mod-`p` shadows can lie in both
directions: a naive scalar rank can be full when the operator curvature is zero,
and zero when the operator curvature is full.

Transfer to HYP-2718: raw relation height, raw gap, and raw atom size are only
shadows.  The proof object must be an operator-style carrier ledger that records
which exact relation, phase, support, and cell word produced the scalar `Q_0`
contribution.

### Root Packets

THM-513 shows the FKN neighborhood is root-packeted: one-flip atoms are interval
roots, and two-flip `c3` is controlled by packet incidence.

Transfer to HYP-2718: low-height carrier relations should be packet-atlased,
not pooled.  The finite ledger should group by relation support, mod-7 phase,
sign chamber, and state-word address, with a separate exact charge for each
packet family.

### Endpoint-Protected Boundary

The S450 LRC analogy atlas identified the native LRC hard object as a protected
boundary of the bad cover, not a visual circle-cover resemblance.

Transfer to HYP-2718: `Q_0` is the protected boundary atom of the transformed
missed-sector law.  Useful external analogies must preserve that boundary atom
and its endpoint debt, or they are only metaphors.

### Tiling Recursion Warning

THM-442 and the half-tiling addendum say the full and half tiling recurrences
are real cell-affine recurrences, but they do not compute non-additive
cycle-space invariants after scalarization.

Transfer to HYP-2718: half-tilings are useful as an address quotient and for
organizing state words.  They are not, by themselves, a scalar recurrence proof
for `Q_0` unless the carrier labels survive the quotient.

### Divisor-Copy / Squarefree Profile

The Euler-copy reframe is Mobius inversion on divisor zeta data: the copy
function satisfying `sum_{d|n} c(d)=n` is `phi(n)`.  The useful part is not the
totient number alone but the fact that the correct basis is found by inverting
the incidence algebra.

Transfer to HYP-2718: miss-zeta coordinates are the incidence-zeta side, and
`Q_0` is obtained only after the correct finite Mobius/binomial inversion.
Squarefree profiles are a warning to keep support masks until after inversion.

## Proof Moves To Try Next

1. **Origin-atom transfer tax.**  Interpolate from actual shared-`x` carriers
   to independent product carriers one carrier at a time, track transitions in
   `T=#missed sectors`, and search for an identity whose positive `Q_0` terms
   are paid by named tail transitions.

2. **Generated singleton base.**  Prove the HYP-2698 singleton-product kernel
   inequality in the factorial-origin basis, then show that merging singleton
   carriers into coherent blocks cannot increase the admissible `Q_0` pressure.

3. **Mode-pair cancellation.**  Expand `Q_0` in carrier Fourier modes and pair
   `n` with `-n` before taking absolute values.  Test whether the leading packet
   terms cancel by mirror/complement parity.

4. **Low-height packet atlas.**  Extend the S68 relation-height scout to print
   the whole `Q_t` atom profile by relation-height bucket, first for two
   carriers and then triples.  The finite ledger should store exact fractions.

5. **Cofactor search.**  For generated rows, fit the atom profile to a small
   polynomial in the falling basis and look for a cofactor whose endpoint
   values mirror the THM-438 wild/tame bridge.

6. **Operator ledger.**  Record each tail term by
   `(n, support, height, n.M, phase, sector mask, t, x-cell word, sign chamber)`
   and only then ask whether the signed sum fits under `cap_k-ProductCover`.

## False Friends

- A product upper envelope is false; S68 and HYP-2718 both see mixed signs.
- Raw gap monotonicity is false; relation height and arithmetic phase survive.
- Full residual-coordinate dominance is too strong; HYP-2698's generated cone is
  smaller than the full orthant.
- Literal full-torus Weyl is false for multiple integer carriers because exact
  carrier relations always remain.
- Full or half tiling recurrences are cell-affine tools, not scalar proofs for
  non-additive `Q_0` after labels have been erased.

## Tournament Analysis

Vertices are proof-route carriers, not runners:

```text
generated_miss_zeta_word
relation_height_top_character_tail
bonferroni_origin_tax
finite_difference_cofactor
low_height_packet_atlas
operator_ledger
endpoint_boundary
tiling_address_quotient
raw_scalar_analogies
```

Pairwise observable: route `A` beats route `B` if it preserves more of
`Q_0`, generated labels, finite-ledger address, and tail-bound compatibility,
with penalties for known false friends.  The resulting pressure path is

```text
generated_miss_zeta_word
> relation_height_top_character_tail
> bonferroni_origin_tax
> finite_difference_cofactor
> low_height_packet_atlas
> operator_ledger
> endpoint_boundary
> tiling_address_quotient
> raw_scalar_analogies.
```

This ranking is intentionally transitive: score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}`, directed `3`-cycles `0`.

## Assumption Challenge

This pass considered runners, carrier blocks, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
relation vectors, residual masks, missed-count atoms, low-height packets, and
proof obligations as possible vertices.

The productive vertex set is proof-route carriers that preserve the predicate
`Q_0=ProductCover-p0` and retain enough hidden address data to make a signed
estimate.  This quotient destroys sector ownership, block provenance, and
exact cell history unless they are explicitly reattached in the ledger.

Challenged assumption: the final analytic step should look like a scalar
equidistribution, span, or product-envelope estimate.  The repo analogies point
instead to a generated-word, relation-filtered, finite-difference boundary
atom estimate.

