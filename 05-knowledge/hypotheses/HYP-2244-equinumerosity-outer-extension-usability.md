---
id: HYP-2244
status: OPEN method hypothesis with finite cardinal-collision evidence
source: codex-2026-06-05-S669
related:
  - HYP-2243
  - HYP-2242
  - HYP-2187
  - HYP-2186
  - HYP-2236
---

# HYP-2244: Equinumerosity Is Not Outer-Extension Usability

## Claim

Equinumerosity is a cardinal/count shadow.  It becomes a usable proof transfer
only after attaching a retained profile that keeps the target predicate pure
under the allowed outer extensions.

```text
same cardinality or same count
  -> attach deletion/link/embedding/fiber profile
  -> audit target-predicate purity under outer extension
```

This extends HYP-2243 from finite-tree embedding profiles to arbitrary
count-equivalent worlds.

## Source Anchors

The prompt pointed again to Harvey Friedman's official pages:

- `https://u.osu.edu/friedman.8/foundational-adventures/publications/`
- `https://u.osu.edu/friedman.8/foundational-adventures/downloadable-manuscripts/`

Relevant official publication anchors include `Reverse Mathematics of
Homeomorphic Embeddings`, `Internal finite tree embeddings`, `Finite Functions
and the Necessary Use of Large Cardinals`, and `Invariant Maximality and
Incompleteness`.

The two outside count/dimension anchors used this session:

- Royle et al., `Tournaments and even graphs are equinumerous`
  (`https://arxiv.org/abs/2204.01947` and the St Andrews repository page).
- The standard topological distinction that finite-dimensional Euclidean
  spaces have the same continuum cardinality but are not homeomorphic when
  dimensions differ; dimension is retained by homeomorphism/invariance-of-domain
  style data, not cardinality alone.

## S669 Evidence

S669 adds `04-computation/equinumerosity_outer_extension_s669.py` and stores
`05-knowledge/results/equinumerosity_outer_extension_s669.out`.

### Equal-Cardinality Euclidean Toy

The finite toy uses three carriers with identical cardinal shadow `|V|=27`:

| Carrier | Dimension label | Graph |
|---|---:|---|
| `line_P27` | `1` | path on 27 vertices |
| `plane_3x9` | `2` | rectangular grid on 27 vertices |
| `space_3x3x3` | `3` | cubic grid on 27 vertices |

Target predicate: the dimension label.

Projection audit:

| Channel | Groups | Mixed fibers | Max bucket |
|---|---:|---:|---:|
| `cardinality` | `1` | `1` | `3` |
| `cardinality+deletion` | `2` | `1` | `2` |
| `cardinality+degree` | `3` | `0` | `1` |
| `growth_profile` | `3` | `0` | `1` |
| `small_embedding_profile` | `3` | `0` | `1` |

Cardinality alone mixes all dimensions.  Single-point deletion splits the
`1D` carrier from the higher-dimensional carriers, because path punctures can
disconnect while the grid and cube remain connected.  But deletion alone still
mixes `2D` and `3D`.  Degree, local ball-growth, and small embedding profiles
separate all three.

### Outer-Extension Check

Allowed extension: add an outer layer while preserving the internal carrier
type:

```text
line_P27     -> line_P29
plane_3x9    -> plane_5x9
space_3x3x3  -> space_5x3x3
```

The full profile can change because boundary and growth counts change.  The
dimension signature

```text
("branch3", "branch5", "C4", "puncture_splits")
```

remains stable in all three extensions.  This is the finite version of the
claim: raw size is extension-unstable, while the right embedding/link signature
is usable.

### Royle/Tournament Count Shadow

S669 imports S617's count audit:

- degree-even graph classes agree with tournament counts only at `n=3`;
- naive Burnside-even is `graphs - tournaments`, not the Royle-even property;
- even when the true Royle-even equinumerosity theorem is used, the repo still
  lacks a predicate-preserving fiber functor.

S617's retained-fiber splits remain the key warning:

- at `n=5`, `H=9` splits into `beta1=0` and `beta1=1`;
- at `n=6`, `H=17,23,33,37,45` split by `beta1`;
- adding disjoint directed-3-cycle packet data refines the quotient further.

So the tournament/even-graph bridge is not "same number, therefore same proof
object."  It is:

```text
same number + functor preserving H/beta1/packet/embedding side channels.
```

## Tournament Analysis

Vertices are proof-transfer channels:

1. `address_plus_embedding`
2. `small_embedding_profile`
3. `local_growth_profile`
4. `(H,beta1,packet)_fiber`
5. `deletion_profile`
6. `Royle_even_cardinal_theorem`
7. `bijection_without_profile`
8. `raw_cardinality`

Pairwise observable:

```text
(cardinal fidelity,
 extension purity,
 dimension retention,
 tournament fiber retention,
 Friedman embedding fit,
 actionability)
```

Switch: majority.  Tie Hamiltonian path: the listed priority order.

Fingerprints:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`
- `directed_3cycles=0`
- `scc_sizes=[1,1,1,1,1,1,1,1]`
- `hamiltonian_paths=1`

The transitive order puts `address_plus_embedding` first and
`raw_cardinality` last.

## LRC14 Transfer

For LRC14, equinumerous/cardinal shadows correspond to visible scalar quotients
such as `Res_27` shell size, raw support, or raw count of covers.  These are
useful guides but not proof objects.

HYP-2244 says the proof object should be a retained profile:

- owner-private deletion bits;
- carry cocycle labels;
- endpoint protector fragments;
- pair-pinch and D/U/N owner fragments;
- bounded proof-obligation embedding profiles.

The theorem target becomes:

```text
Fixed visible scalar/cardinality shadow plus fixed embedding/fiber profile
forces AP, Vstar, nonprimitive 2AP, or strict looseness,
except for globally coherent scalar floor lifts.
```

## Next Tests

1. Replace the finite grid toy by actual point-set/unit-distance deletion
   profiles: equal edge counts should be split by point-deletion gain and
   direction-support signatures.
2. Build the Royle-even fiber functor target explicitly: search whether any
   graph-side invariant preserves `H`, `beta1`, and disjoint `c3` packets.
3. Build the LRC14 proof-obligation embedding profile over coherent carry
   subspaces rather than local Hamming perturbations.
