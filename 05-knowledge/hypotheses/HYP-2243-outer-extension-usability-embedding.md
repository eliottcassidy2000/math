---
id: HYP-2243
status: OPEN method hypothesis with bounded finite-tree evidence
source: codex-2026-06-05-S668
related:
  - HYP-2242
  - HYP-2241
  - HYP-2240
  - HYP-2171
  - HYP-2167
  - HYP-2165
  - HYP-2164
---

# HYP-2243: Outer Extension Usability and Embedding

## Claim

Extend HYP-2242 from "what is an embedded maximum?" to "when is an embedded
maximum usable under ambient growth?"

Local terminology:

```text
outer extension = an allowed growth of the ambient object that preserves the
                  internal bonds of the current object;

usable quotient = a quotient/address channel whose fibers remain pure for the
                  target predicate after every allowed outer extension;

embedding address = the retained finite downset/profile that tells which probes
                    embed into the extended object, plus the extension address.
```

The resulting repo-form theorem targets are:

```text
Outer Extension Usability Theorem:
  a quotient is usable for a predicate only if allowed outer extensions have
  retained addresses that keep the predicate pure.

Outer Extension Embedding Theorem:
  for tree-like or proof-obligation-like carriers, the natural retained address
  is a bounded embedding/downset profile together with the insertion/extension
  address.
```

This is a finite, operational cousin of the Friedman/Kruskal/reverse-math
phenomenon: simple-looking finite tree and embedding statements can hide very
high proof strength when the bad-sequence/extension bound is made uniform.

## Source Anchors

The user pointed to Harvey Friedman's publication page:

- `https://u.osu.edu/friedman.8/foundational-adventures/publications/`
- `https://u.osu.edu/friedman.8/foundational-adventures/downloadable-manuscripts/`

Relevant official anchors from the page include recursion theory, `Reverse
Mathematics of Homeomorphic Embeddings`, `Internal finite tree embeddings`,
`Finite Functions and the Necessary Use of Large Cardinals`, `Invariant
Maximality and Incompleteness`, and `Concrete Mathematical Incompleteness:
Basic Emulation Theory`.

The downloadable manuscript list includes:

- `Finite Trees and the Necessary Use of Large Cardinals`, where insertion
  domains/rules place new higher vertices into finite trees.
- 2025 `Invariant Maximality Derivations` and `Invariant Maximality Reversals`,
  where order constraints, usability, invariant maximality, and reversals to
  strong consistency statements are developed.

Important caution: S668 did not find the exact phrase "outer extension" in the
opened Friedman PDFs.  The phrase is used here as repo terminology for the
extension operation suggested by the prompt and by the finite-tree insertion
machinery.

## S668 Evidence

S668 adds `04-computation/outer_extension_usability_friedman_s668.py` with
stored output in
`05-knowledge/results/outer_extension_usability_friedman_s668.out`, plus
reflection `07-reflections/outer-extension-usability-friedman-s668.md`.

The script builds a bounded toy universe of colored rooted finite trees:

- `colors=3`
- `max_nodes=5`
- `tree_count=1788`
- `outer_extension_rows=3204`
- target predicate: the extension contains a rooted homeomorphic color chain
  `0 -> 1 -> 2` somewhere.

It treats a one-step outer extension as adding one colored leaf at a chosen
parent, while recording the insertion address:

```text
(parent_path_colors, new_leaf_color, parent_degree_before)
```

It then audits which quotient channels keep the target predicate pure across
all extension rows.

Projection audit:

| Channel | Groups | Mixed fibers | Max bucket |
|---|---:|---:|---:|
| `size_height` | `10` | `6` | `1215` |
| `color_hist` | `52` | `10` | `332` |
| `frontier` | `633` | `46` | `21` |
| `outer_address` | `1899` | `55` | `9` |
| `small_embedding_profile` | `1500` | `0` | `14` |
| `full_embedding_profile` | `1785` | `0` | `3` |
| `address_plus_small_embedding` | `3120` | `0` | `2` |

The key signal is sharp: coarse size/color/frontier data leaks under outer
extension, and even the raw insertion address leaks.  A small homeomorphic
embedding downset profile repairs the bounded toy completely.

N-color sweep at `max_nodes=4`:

| Colors | Trees | Extension rows | `0->1->2` rows | `size_height` mixed | `small_embed` mixed |
|---:|---:|---:|---:|---:|---:|
| `1` | `8` | `8` | `0` | `0` | `0` |
| `2` | `72` | `96` | `0` | `0` | `0` |
| `3` | `303` | `441` | `21` | `3` | `0` |
| `4` | `876` | `1328` | `29` | `3` | `0` |

So the "tree of 3" is the first nontrivial color threshold for this predicate;
once it exists, coarse quotients leak, while embedding profiles remain stable in
the bounded window.

## Tournament Analysis

Vertices are quotient/address channels, not tree vertices.  The observable is:

```text
(exactness,
 max_bucket score,
 compression,
 embedding naturality,
 LRC transfer,
 actionability)
```

Fingerprints:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`
- `directed_3cycles=0`
- `scc_sizes=[1,1,1,1,1,1,1]`
- `hamiltonian_paths=1`

Top order:

1. `address_plus_small_embedding`
2. `small_embedding_profile`
3. `full_embedding_profile`
4. `outer_address`
5. `frontier`
6. `color_hist`
7. `size_height`

The LRC-shaped channel wins: keep the extension address, but pair it with an
embedding/downset profile rather than treating the address as sufficient by
itself.

## LRC14 Transfer

The dictionary is:

| Finite-tree toy | LRC14 |
|---|---|
| internal rooted tree | visible `Res_27` shell |
| one-leaf outer extension | `+27` carry / owner-route lift |
| homeomorphic embedding downset | owner/deletion/carry proof-obligation profile |
| color chain `0->1->2` | floor-vs-strict obstruction predicate |
| small embedding profile pure | private-owner deletion bit separates S666 fibers |

HYP-2241 already found the first LRC instance: visible `Res_27` data leaks, but
the paired owner-private deletion bit gives zero mixed local carry fibers.
HYP-2243 says to generalize this from one bit to an embedding profile over proof
obligations:

```text
For each coherent carry/owner extension, record which bounded proof-obligation
probes embed into the extended certificate.  The extension is usable only if
this profile separates floor from strict rows, except for globally coherent
scalar floor lifts.
```

## Micro-Incompleteness Reading

The script includes a deliberately tiny greedy controlled bad-sequence toy.  It
is not a theorem and does not model Friedman's full strength.  It records the
right warning: finite bad sequences can be built locally, and the hard content
is a uniform no-infinite-bad-sequence or controlled-miniature bound.

This is the repo's "micro incompleteness" moral: a quotient may look stable on a
small finite window, but a uniform outer-extension theorem asks for a bound over
all finite stages.  That is exactly where homeomorphic embeddings, recursion
theory, and finite combinatorics begin to carry proof strength.

## Next Tests

1. Build the LRC14 proof-obligation embedding profile: probes should be small
   D/U/N owner-route certificates, endpoint protector fragments, and carry
   cocycle fragments.
2. Audit coherent carry subspaces with keys:
   `visible Res_27`, `outer carry address`, `owner-private deletion profile`,
   and `bounded proof-obligation embedding profile`.
3. Repeat the finite-tree toy with target predicates other than `0->1->2`:
   antichain of three leaves, three-color path with gaps, and `n`-color chains.
4. For tournament decks, replace "tree embeds into extension" with "rooted card
   probe embeds into the deleted-card deck"; test whether paired deleted-score
   is the first small embedding profile.
