# HYP-2467 - Small-factor resonance capacity separates the 27/28 unit-distance gate

**Status:** OPEN proof program with exact connected-factor atlas.

**Source:** codex-2026-06-13.  Extends OPEN-Q-057, THM-493, THM-494,
THM-495, HYP-2461, HYP-2462, and HYP-2466.

## Claim

The resonant-product obstruction at `N=27=3*9` is stronger than the curated
THM-493 search suggested.  Exhausting every connected triangular-lattice factor
patch through size `9`, modulo translation and `D6` symmetry, shows:

```text
N=27: best connected two-factor resonant product = 75 < 81 = 3N.
N=28: best connected two-factor resonant product = 85 > 84 = 3N.
```

The 27 lane fails because the edge-dense size-3 factor is `K3`, and `K3` has no
non-degenerate norm-`t` displacement (`t>1`).  Size-3 factors that do carry a
transverse displacement lose an internal edge; their best resonance bonus
against all connected 9-patches does not repay that generic edge loss.

This complements THM-495's chord-spectrum theorem.  THM-495 proves the bonus can
survive only on shared chord norms and that `t=3` is forced at the `4*7`
crossing.  HYP-2467 adds the finite capacity side: among all connected factors
of the relevant sizes, buying a non-unit chord at size `3` costs too many unit
edges, while size `4` is the first place where edge density and a `sqrt(3)`
chord coexist.

Thus, inside the exact connected two-factor resonant carrier class, `28=4*7`
is the first `3N` crossing.  This does not prove `u(27)<=81`, because an
arbitrary 27-point planar configuration need not compress to two connected
triangular factors.  It gives a sharper proof target: any 82-edge counterexample
at 27 must be genuinely irreducible with respect to this small-factor
resonance-capacity quotient.

## Exact Evidence

Script:

```text
04-computation/unit_distance_resonance_capacity_atlas_codex.py
```

Stored output:

```text
05-knowledge/results/unit_distance_resonance_capacity_atlas_codex.out
```

Connected triangular patch counts, up to translation and `D6`, are:

```text
size 1: 1
size 2: 1
size 3: 3
size 4: 7
size 5: 22
size 6: 82
size 7: 333
size 8: 1448
size 9: 6572
```

For each factorization `N=a*b` with `a,b<=9`, the script computes

```text
U(A (+)_t B) = e(A)|B| + |A|e(B) + Delta_t(A,B),
Delta_t = 1/2 * sum_v m_A(v) m_B(v),
```

where the sum runs over non-degenerate shared norm-`t` displacement vectors
with `t>1`; every relative `D6` orientation of the second factor is tested.

Capacity scan:

```text
N=24: best 68 < 72   [4x6, generic 66, bonus 2, t=3]
N=25: best 72 < 75   [5x5, generic 70, bonus 2, t=3]
N=27: best 75 < 81   [3x9, generic 75, bonus 0]
N=28: best 85 > 84   [4x7, generic 83, bonus 2, t=3]
N=30: best 90 = 90   [5x6, generic 87, bonus 3, t=3]
N=32: best 99 > 96   [4x8, generic 96, bonus 3, t=3]
```

Size-3 stress against all connected 9-patches:

```text
K3:                         total 75, generic 75, bonus 0
collinear 3-path:           total 69, generic 66, bonus 3, t=4
bent path with norm-3 pair: total 70, generic 66, bonus 4, t=3
```

This is the finite resource inequality the next proof should try to lift:
resonance-bearing size-3 factors are too edge-poor, while the edge-rich size-3
factor is resonance-free.

## Interpretation

THM-493 says the crossing bonus is a displacement-spectrum correlation, not a
mystery edge.  HYP-2467 says that correlation has a small-factor capacity law:

```text
edge density + transverse displacement support
```

is a joint resource.  At size `3`, the two resources are mutually exclusive.
At size `4`, the rhombus carries both resources, and pairing it with the
7-point rosette gives exactly the known `85`-edge crosser.

This makes the `27 -> 28` transition look less like a lucky construction and
more like the first point where the carrier has enough room to hold both an
edge-dense factor and a norm-3 displacement packet.

## Assumption Challenge / Tournament Analysis

This session did not take unit-distance points or unit edges as the default
tournament vertices.

Candidate vertex sets considered:

- points,
- unit edges,
- connected patches,
- factor sizes,
- displacement vectors,
- Loeschian rungs,
- carrier lattices,
- proof obligations.

The selected Tournament Analysis uses proof obligations as vertices:
`small_factor_capacity`, `generic_edge_budget`, `resonance_bonus_spectrum`,
`free_patch_annealing`, `carrier_ladder_gate`, and `global_upper_bound`.

It preserves route leverage toward `u(27)<=81`, while destroying individual
coordinates.  The challenged assumption is that the proof should be organized
around points or directions; here the resource capacity itself is the vertex.

Stored tournament fingerprints:

```text
score histogram: {0:1, 1:1, 2:1, 3:1, 4:1, 5:1}
directed 3-cycles: 0
SCC sizes: [1,1,1,1,1,1]
Hamiltonian path count: 1
leader: small_factor_capacity
```

The transitivity is itself useful: in this local proof quotient, the atlas
capacity lemma should be attacked before the global upper bound.

## Next Proof Route

1. Prove the size-3 capacity lemma without enumeration:
   any 3-point triangular factor is either `K3` and has no `t>1`
   displacement, or has at most two unit edges and cannot recover the lost
   generic budget by any norm-`t` correlation against nine points.
2. Lift the lemma from exact triangular factors to Moser-lattice patch
   decompositions: find a compression theorem saying every putative 27-point
   `82`-edge Moser patch either descends to a forbidden connected-factor
   capacity lane or exposes a genuinely irreducible obstruction.
3. Add a branch-and-bound side channel to free-patch search: track
   displacement-spectrum capacity, not only edge count, so candidate
   27-point patches with impossible bonus budgets can be pruned early.
4. Transfer the resource-coordinate viewpoint back to LRC14:
   small factor size here plays the same role as shell-27/13-clock debt there;
   both proofs should count incompatible resources rather than only scalar
   witnesses.
