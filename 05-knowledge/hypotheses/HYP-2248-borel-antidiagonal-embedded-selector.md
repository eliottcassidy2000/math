---
id: HYP-2248
status: OPEN companion/addendum method hypothesis with finite invariant-selector evidence
source: codex-2026-06-05-S673b
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

# HYP-2248 Addendum: Borel Anti-Diagonalization as Embedded Selector Obstruction

Companion to `HYP-2248-borel-constructive-embedded-maximality.md`.  The main
S673 lane proves that diagonal witnesses require retained diagonal addresses;
this addendum measures a complementary address tax for invariant outside
selectors under finite symmetry quotients.

## Claim

The repo should distinguish two opposite-looking diagonal moves:

```text
raw diagonalization:
  given a named list, construct an object outside it

Borel / invariant anti-diagonalization:
  under the right invariance requirement, a uniform outside selector is forced
  to hit the named object somewhere
```

The useful bridge to the repo is:

```text
outside selector + invariance => address tax
outside selector + outer extension => recursive capture
embedded maximality = the statement that the named address survives every
                     allowed capture/extension move
```

So "constructive mathematics" is not the opposite of diagonalization here.  It
is the discipline of recording the witness procedure, its named address, and
the stage/rank at which outer extension can absorb the witness.

This sharpens HYP-2247.  Recursion is the fourth representation face, but
HYP-2248 says why a recursive face needs an address before it can become a
proof: otherwise the selected outside point is either not invariantly
selectable, or it becomes the next named point of the ambient extension.

## Source Anchors

Friedman's FOM note on Borel/DST records the core shape:

```text
There is a Borel function F:R^infinity -> R such that F(x) is off x.
But if F is invariant under the listed sequence-equivalences, then there is an
x such that F(x) is a coordinate of x.
```

The same source explicitly frames this as revisiting his old Borel
diagonalization work in `On the Necessary Use of Abstract Set Theory`.

The tangible-incompleteness notes put this near invariant maximality,
emulation theory, and sequential constructions: finite rational-vector
statements can remain transparent while the uniform construction principle
carries large proof strength.  The embedded-maximality note then reintroduces
embeddings and maximal squares over rational boxes, giving the local model for
this repo:

```text
maximality is meaningful only relative to the embedding and allowed maps.
```

External source URLs used:

- `https://fomarchive.ugent.be/2003-May/006637.html`
- `https://fomarchive.ugent.be/2019-September/021671.html`
- `https://fomarchive.ugent.be/2020-May/022140.html`
- `https://u.osu.edu/friedman.8/foundational-adventures/publications/`

## S673 Finite Selector Evidence

S673b adds `04-computation/borel_antidiagonal_embedded_maximality_s673.py` and
stores output in
`05-knowledge/results/borel_antidiagonal_embedded_maximality_s673.out`.

The finite model is deliberately modest:

- `U={0,...,n-1}`;
- a state is a named subset `A`;
- an outside selector must choose `y in U\A`;
- a symmetry group `G` acts on `U`;
- an invariant outside selector can choose at `A` only if some outside point is
  fixed by every symmetry stabilizing `A`;
- naming anchors shrinks `G`.

The minimum number of anchors making every non-full state selectable is the
finite "address tax."

For `n=6`:

| Symmetry group | Selectable states | Blocked states | Address tax |
|---|---:|---:|---:|
| ordered/trivial | `63/63` | `0` | `0` |
| path reflection | `56/63` | `7` | `1` |
| cyclic rotations | `54/63` | `9` | `1` |
| dihedral cycle | `36/63` | `27` | `2` |
| full symmetric | `6/63` | `57` | `5` |

The full-symmetric line is the finite warning: if the quotient forgets almost
all address, only the co-singleton states can make a canonical outside choice.
Most apparent anti-diagonal choices are not invariantly meaningful.

After the needed anchors are named, the constructive least-outside selector is
still captured by outer extension:

```text
A=() -> select 0; outer extension names 0
A=(0,) -> select 1; outer extension names 1
...
A=(0,1,2,3,4) -> select 5; outer extension names 5
```

This is the finite skeleton of the Borel anti-diagonal intuition: a selector
does not become a proof merely by choosing outside the current stage.  The
extension law can name what was just selected.

## Tournament Analysis

Vertices are proof/selector lanes:

```text
embedded_address_tax
borel_invariant_antiselector
ph_bad_branch_rank
lrc_owner_carry_rank
a000568_half_filter_rank
constructive_named_order
raw_cantor_diagonal
raw_count_shadow
```

The pairwise observable is:

```text
invariance_respect,
constructive_witness,
needs_address,
recursion_rank,
embedded_maximality_fit,
LRC_transfer,
overclaim_safety
```

Fingerprints:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`
- `directed_3cycles=0`
- one Hamiltonian path:

```text
embedded_address_tax
> lrc_owner_carry_rank
> borel_invariant_antiselector
> ph_bad_branch_rank
> a000568_half_filter_rank
> constructive_named_order
> raw_cantor_diagonal
> raw_count_shadow
```

This ranking says the repo should now treat "add the address tax" as more
primary than either raw diagonal construction or raw recursion.

## Transfers

### LRC14

HYP-2241's owner-private deletion flag is an address tax already paid.  The
next theorem should not say only that this coordinate separates sampled local
fibers.  It should say:

```text
after the owner-private address is attached,
every coherent +27 carry extension either
  (a) remains in the global floor family,
  (b) pays positive maximin tax,
  (c) or drops a bad-child rank.
```

This is Borel anti-diagonalization in LRC clothing: a would-be counterexample
cannot keep selecting its loneliness witness into the unnamed tail if the
owner/carry address is invariantly retained.

### A000568

HYP-2246 found that unpaired deletion-card data almost reconstructs tournament
classes, but loses an `L/U` owner side in two regular strong collision buckets.
HYP-2248 reads the half-filter trace as an address tax: without the side, the
quotient has no invariant selector for the owner; with it, the endpoint
recursion becomes a finite candidate for a constructive enumerator.

The next refinement is:

```text
half-filter trace + child-count / extension-rank profile.
```

### Paris-Harrington

HYP-2247 already supplies the recursion-rank face.  HYP-2248 supplies the
selector language:

```text
bad coloring = a selector that keeps choosing outside the homogeneous trace
relative largeness = a named initial-segment gate
PH witness growth = the uniform bound for recursive capture of all bad selectors
```

The exact pair miniature is small, but the formulation is correct: the raw
count of bad colorings is less important than the profile of bad children.

### Unit Distance

For unit-distance constructions, symmetric point sets can make the "next
outside point" non-canonical.  A spine-owner or point-deletion frontier is the
address tax.  The proof target should ask whether every slab/ear extension
either names a unit-spine witness or drops a bad-construction rank.

## Assumption Challenge

For LRC/Tournament Analysis, vertices need not be runners, arcs, or even
colorings.  S673 uses proof lanes and selector states.  Candidate vertex sets
considered included named subsets, symmetry groups, anchors, bad coloring
branches, endpoint owner states, carry fibers, proof obligations, and extension
laws.

Preserved predicate:

```text
whether an outside/witness selector remains meaningful and pure under the
allowed extension relation.
```

Destroyed information:

```text
raw labels and scalar counts when they do not survive the invariance quotient.
```

## Next Tests

1. Define the LRC14 `bad-child rank` over the HYP-2241 owner-private fibers and
   test it on coherent `+27` carry subspaces rather than all local Hamming
   perturbations.
2. Extend the A000568 endpoint enumerator with
   `half-filter trace + child-count profile` and audit whether it reduces
   `8->9` fallback canonicalization.
3. Build a finite "Borel selector" toy where the group is not a value-symmetry
   group but a coordinate finite-permutation group acting on bounded lists;
   compare address tax with the subset model here.
4. For unit distance, test whether point-deletion spine-owner profiles separate
   symmetric extension choices in the `n=21/22` Moser carriers.
