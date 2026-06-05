---
id: HYP-2248
status: OPEN method hypothesis with finite antidiagonal and regressive toys
source: codex-2026-06-05-S673
related:
  - HYP-2247
  - HYP-2246
  - HYP-2245
  - HYP-2243
  - HYP-2242
  - HYP-2241
  - HYP-2240
  - HYP-2189
---

# HYP-2248: Borel Antidiagonalization Needs Embedded Constructive Addresses

## Claim

Borel antidiagonalization, constructive mathematics, tangible incompleteness,
and embedded maximality are one proof-design pattern:

```text
diagonal escape exists
  -> constructive witness needs a named address
  -> embedded maximality names the ambient extension/cut where the witness lives
  -> tangible incompleteness appears when the finite extension law requires
     a uniform bound stronger than the local quotient can prove
```

In the HYP-2247 four-face language:

```text
fraction  = diagonal / cut / owner / Borel code address
recursion = extension law that constructs, preserves, or diagonalizes past it
```

The slogan:

```text
Do not merely prove the antidiagonal exists; retain the address that computes it.
```

## Source Anchors

Cantor diagonalization supplies the primitive pattern: from a proposed list,
construct a diagonal set/word disagreeing with each row at its own coordinate.
MathWorld records the general `S -> P(S)` form and the contradiction produced
by the diagonal subset.

Constructive mathematics supplies the witness-extraction discipline.  The
Stanford Encyclopedia of Philosophy article distinguishes constructive
practice by its algorithmic/intuitionistic methodology and records Bishop/CZF/
type-theory foundations as modern carriers.

Friedman's finite-functions paper is the tangible-incompleteness anchor: it
advertises finite mathematical theorems that require large-cardinal strength
and presents that issue in concrete natural-number terms.

Friedman's publication list gives nearby source territory: proof theory and
intuitionism, set-theoretic foundations for constructive analysis, reverse
mathematics of homeomorphic embeddings, internal finite tree embeddings,
elementary descent recursion, maximality, and Borel reducibility.

The opened Friedman subtle-cardinals/linear-ordering manuscript gives a useful
technical miniature: a `k`-critical linear ordering has no endpoints and forces
shifted equalities for every regressive map on increasing `k`-tuples.  HYP-2248
uses only a finite toy shadow of that idea, not the theorem.

Links used:

- <https://mathworld.wolfram.com/CantorDiagonalMethod.html>
- <https://plato.stanford.edu/entries/mathematics-constructive/>
- <https://arxiv.org/abs/math/9811187>
- <https://u.osu.edu/friedman.8/foundational-adventures/publications/>
- <https://bpb-us-w2.wpmucdn.com/u.osu.edu/dist/1/1952/files/2014/01/subtlecardinals-1tod0i8.pdf>

## S673 Antidiagonal Channel Audit

Given an `m x m` binary matrix `M`, define the constructive antidiagonal witness:

```text
a_i = 1 - M_ii.
```

S673 exhausts all binary matrices for `m<=4` and audits which quotient channels
determine `a`.

For `m=4`:

| Channel | Groups | Mixed antiword buckets | Max bucket | Pure? |
|---|---:|---:|---:|---|
| row weight multiset | 70 | 68 | 6912 | no |
| row weight sequence | 625 | 609 | 1296 | no |
| column weight sequence | 625 | 609 | 1296 | no |
| row and column weights | 16145 | 8769 | 90 | no |
| diagonal vector | 16 | 0 | 4096 | yes |
| diagonal plus row weights | 4096 | 0 | 81 | yes |
| full matrix | 65536 | 0 | 1 | yes |

So even exact row/column scalar summaries fail to construct the witness.  The
diagonal vector is the minimal embedded address.

The mixed examples are tiny and vivid.  For `m=4`, two matrices can have the
same row-weight sequence `(1,0,0,0)` while requiring different antiwords:

```text
diag 1000 -> antiword 0111
diag 0000 -> antiword 1111
```

This is the finite version of the constructive warning:

```text
existence by diagonalization is programmatic only if the diagonal address
survives the quotient.
```

## Regressive Shift-Collision Toy

For increasing `k`-tuples `t` in `[N]` with `min(t)>0`, a finite regressive map
has:

```text
f(t) < min(t).
```

The shifted-collision target for a `(k+1)`-block
`b_1 < ... < b_{k+1}` is:

```text
f(b_1,...,b_k) = f(b_2,...,b_{k+1}).
```

This is only a finite toy inspired by Friedman's `k`-critical linear-ordering
definition.  It deliberately shows why endpoints/embedding matter.

On finite intervals through the checked windows, there is exactly one avoiding
assignment:

```text
f(t) = min(t) - 1.
```

For `k=1,2,3`, S673 finds one avoiding assignment at every checked `N`.  The
predecessor map avoids all shifted collisions.  When the predecessor boundary
value is forbidden, domains become empty immediately and the toy collapses.

The moral is not that the finite toy proves anything about subtle cardinals.
It says:

```text
the boundary/predecessor address is the whole avoiding branch.
```

A no-endpoint/large-cardinal/critical-order theorem is exactly the kind of
embedded maximality statement where that cheap predecessor escape cannot be
used globally.

## Integration With Embedded Maximality

HYP-2242 says a maximum is meaningful only as:

```text
maximal(object, ambient embedding, allowed extensions).
```

HYP-2248 says a diagonal witness is meaningful only as:

```text
witness(object, diagonal address, allowed quotient/extension).
```

The two are the same shape.  A local maximum can be destroyed by a new cut; a
diagonal witness can be destroyed by quotienting away the coordinate that
constructs it.  In both cases, the proof-relevant object is not the scalar
status but the embedded address.

## Constructive Reading

Constructive mathematics is not just a philosophical constraint here.  It is a
test for whether a repo invariant is usable:

```text
Can the invariant produce the next witness/child/extension without consulting
hidden labels?
```

This recasts HYP-2246:

```text
half-filter purity through n=8
```

is a strong classification result, but a production enumerator still needs a
constructive child-selection rule:

```text
parent class + incident orbit + half-filter trace -> selected child class.
```

Purity without selection is existence; purity with selection is constructive
enumeration.

## Tournament Analysis

Vertices are proof-design lanes:

1. `Cantor antidiagonal`
2. `Borel diagonalization`
3. `Constructive mathematics`
4. `Tangible incompleteness`
5. `Embedded maximality`
6. `LRC14 owner/carry rank`
7. `A000568 endpoint recursion`
8. `Unit-distance spine recursion`
9. `Raw quotient/cardinality`

Pairwise observable:

```text
constructive witness,
embedded address need,
recursion strength,
repo transfer,
risk control,
classical existence strength.
```

Switch: majority.  Tie Hamiltonian path: listed priority order.

Fingerprints:

- `score_histogram={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}`
- `directed_3cycles=0`
- `scc_sizes=[1,1,1,1,1,1,1,1,1]`
- `hamiltonian_paths=1`

Top order:

```text
A000568 endpoint recursion
Cantor antidiagonal
Embedded maximality
Tangible incompleteness
Borel diagonalization
Constructive mathematics
LRC14 owner/carry rank
Unit-distance spine recursion
Raw quotient/cardinality
```

The current leader is A000568 because HYP-2246 already has exact finite
constructive machinery.  The deeper target is to move LRC14 and unit distance
up the same ladder by giving them constructive embedded addresses.

## Transfers

### LRC14

Treat the `C=27` owner/carry fiber as a diagonal matrix:

```text
rows/columns = visible residues, cover owners, carry lifts
diagonal address = owner-private deletion / cut where a witness is forced
antiword = explicit lonely-time or strict-tax witness
```

Next proof target:

```text
define diagonal/cut addresses for every AP/Vstar/2AP local carry extension,
then prove the rank drop constructively rather than by post-hoc classification.
```

### A000568

HYP-2246's half-filter trace is the address.  HYP-2248 asks for a constructive
selection theorem:

```text
If a raw endpoint-orbit candidate has a pure half-filter bucket,
then the bucket has a canonical selected child representative
computable from parent deck data.
```

### Unit Distance

Unit-spine ownership is the diagonal address.  Bulk edge counts and direction
counts are row/column weights: useful, but insufficient to construct the
mandatory witness spine.  A unit-distance proof should retain:

```text
point-deletion frontier owner,
unit-spine owner,
bulk tail extension rank.
```

### Borel / CH / Generic Threads

The repo should be careful with set-theoretic analogies.  The productive
finite question is not "is this literally a Borel theorem?" but:

```text
which code/address must remain Borel/constructible/natural under the quotient?
```

That is the bridge from diagonalization to embedded maximality.

## Next Tests

1. For HYP-2246, implement a child-selection audit: does the half-filter trace
   identify a canonical representative before full child canonicalization?
2. For HYP-2241, define the LRC14 diagonal/cut address table over AP/Vstar/2AP
   carry extensions and test whether it constructs the strict witness.
3. For unit distance, build a small Moser-spine antiword lab: same edge-count
   and direction-count rows, different spine-owner witnesses.
4. Replace the finite regressive interval toy with a no-endpoint finite
   surrogate, such as cyclic or bi-infinite-window regressive constraints, and
   measure when the predecessor escape disappears.
