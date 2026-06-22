---
id: HYP-2847
title: LRC lower-threshold origin-neighborhood width checksum -- HYP-2842's exact-Farey failure is repaired by the b=1 collapse neighborhood for bounded dense branches N<=14
status: FORMALIZATION CHECKSUM; exact-rational scout positive for bounded dense N=6,8,10,12,14; weaker than HYP-2844/HYP-2846 but Lean-facing
source: codex-2026-06-22-S88
related:
  - HYP-2842
  - HYP-2848
  - HYP-2849
  - HYP-2844
  - HYP-2846
  - HYP-2827
  - HYP-2835
  - HYP-2838
  - THM-527
---

# HYP-2847 -- origin-neighborhood width certificate

## Claim

The exact-center version of the Farey witness-floor idea is the wrong quotient
(HYP-2842), and KPS S32's HYP-2844/HYP-2846 incoming work shows the stronger
q-uniform resonant-neighborhood witness route.  This note records an
independent exact-width checksum for the bounded dense lower-threshold branch:
even the proper Farey-neighborhood quotient by itself is too small in the
deepest bounded dense branch, and the missing piece is the
denominator-one collapse: an origin neighborhood of radius `(q-1)/(q*S)` for
`N=2q` and cluster span `S`.

With this `b=1` origin neighborhood included, the exact-rational scout
`04-computation/lrc_lower_threshold_nbhd_width_codex_s88.py` gives positive
bounded dense residual for every even lower threshold checked:

| N=2q | deepest bounded k | proper-only minimum | origin+proper minimum | worst P for origin+proper |
|---:|---:|---:|---:|---|
| 6 | 4 | 0 | `1/15` | `(1,)` |
| 8 | 5 | 0 | `5/168` | `(1, 6)` |
| 10 | 6 | 0 | `1/60` | `(1, 3, 4)` |
| 12 | 7 | 0 | `5/264` | `(1, 3, 4, 5)` |
| 14 | 8 | 0 | `2/273` | `(1, 2, 3, 4, 5)` |

The N=14 row is the tightest reported origin+proper bounded floor.  A Lean
checksum for these rational readouts is in
`TournamentH7.LRCLowerThresholdNeighborhood`.

S89 addendum (HYP-2849): the proper-neighborhood part reduces exactly to a
divisor-depth ledger at bounded span `S=2q-1`; non-divisor `G_P` holes have
positive separation margin from every proper `a/b` neighborhood, so only the
smallest selected multiple of each denominator `b` matters.  The origin
neighborhood is a decoy: if `p=1` is absent, origin mass survives; if `p=1`
is present, one `P` slot is spent and some proper denominator obligation
survives.  For the q=7 stress row this gives the compact split
`1∉P => 11/182` and `1∈P => 2/273`.

## Method

For `N=2q`, a proper Farey center `a/b` with `2 <= b < q` collapses all cluster
phases to at most `b` values, hence gives a cyclic gap at least `1/b`.  If
`S=max(E)-min(E)`, the same gap remains `>1/q` throughout the conservative
neighborhood

```text
|x-a/b| < (1/b - 1/q)/S = (q-b)/(b*q*S).
```

The origin collapse is the same calculation with `b=1`:

```text
|x| < (q-1)/(q*S).
```

A small speed `p in P` removes holes around `j/p` of radius `1/(2*q*p)` from
the safe set `G_P = {x : ||p*x|| >= 1/(2q) for every p in P}`.  The script
unions all guaranteed lonely neighborhoods, subtracts these exact rational
holes modulo one, and scans all `P subset {1,...,2q-1}` with
`|P|=2q-1-k` for `k=q+1,...,2q-1`.

## Interpretation For LRC14

This is weaker than the incoming HYP-2844/HYP-2846 resonant-neighborhood floor:
those scripts use low-denominator centers `{0,1/2,1/3,2/3}` and report much
larger q-uniform `G2` floors.  The value here is different: it isolates an
exact interval-subtraction checksum and a compact Lean-facing bounded-width
formalization target.  After HYP-2842 killed exact Farey-center survival, the
right bounded-dense object is neighborhood width, and the missing neighborhood
is the `b=1` origin collapse.  The N=14 bounded dense stress point
`k=8, |P|=5` has exact residual `2/273` in this conservative model.

The wide-span rows in the same scout can still have zero residual, so the
route should split:

1. bounded dense branch: HYP-2844/HYP-2846 give the stronger resonant route;
   this HYP-2847 checksum gives a smaller exact origin+proper interval target
   that is suitable for Lean formalization;
2. wide/unbounded branch: keep using the existing spreading/decorrelation,
   `L_y`, Part-A, and period-bounded arc machinery.

## Values Smaller Than 14

The same certificate is already positive for N=6,8,10,12.  HYP-2846 gives the
stronger q-uniform witness-route validation for N=6,10,14; this checksum gives
a conservative, exact-rational interval target for formalizing the bounded
dense part.  In these smaller cases the obstruction is not the bounded dense
witness floor; it is only whether the proof system is organized to expose the
right neighborhoods before switching to the wide-span tools.

## Lean Formalization Contract

Current Lean state:

- `TournamentH7.LRCLowerThresholdNeighborhood` records the exact finite
  readouts, proves proper-only failure, origin+proper positivity for
  N=6,8,10,12,14, and proves that `2/273` is the smallest reported
  origin+proper floor.
- `TournamentH7.Verify` imports and audits this checksum.

Next Lean atom:

1. Define a finite list of rational intervals modulo one for the
   origin+proper neighborhoods and the `G_P` holes.
2. Prove interval subtraction lower-bounds `slowμ (goodSet E ∩ safeSet P)`.
3. Prove the containment lemma: if `x` lies in the neighborhood of `a/b` with
   radius `(q-b)/(b*q*S)`, then the cluster phases have a gap `>1/q`.
4. Instantiate the exact bounded dense table for N<=14.

This keeps the formal boundary honest: the new Lean file certifies the exact
arithmetic ledger, while the remaining theorem is the real interval-geometry
bridge into `GOOD(E) ∩ G_P`.

## Tournament Analysis And Assumption Challenge

Pairwise observable: certified residual mass after subtracting `G_P` holes.

Switch/gauge: larger residual wins; ties use proof-route specificity.

Score histogram: `{6:1, 5:1, 4:1, 2:1, 1:1, 0:1}`.

Hamiltonian path:

```text
origin_neighborhood_width > bounded_width > consec_width > wide_width > raw_runner_vertices > exact_farey_points
```

Challenged assumption: tournament vertices do not need to be runners or arcs.
For this branch, better vertices are exact Farey centers, origin neighborhoods,
proper Farey neighborhoods, bounded-width certificates, wide-width
certificates, raw runners/arcs, and proof obligations.  The chosen quotient
preserves interval ownership and residual lower bounds; it destroys exact
center identity, speed ownership, and full component geometry.
