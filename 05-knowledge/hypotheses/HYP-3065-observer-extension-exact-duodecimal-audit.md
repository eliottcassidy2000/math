---
id: HYP-3065
title: Observer-extension exact duodecimal deletion-sector audit
status: SYNTHESIS / exact small tournament audit and LRC carrier abstraction; not a proof
source: codex-2026-06-26-S230
tangent: T1147
script: 04-computation/observer_extension_exact_duodecimal_audit_codex_s230.py
result: 05-knowledge/results/observer_extension_exact_duodecimal_audit_codex_s230.out
related:
  - HYP-3059
  - HYP-3058
  - HYP-3057
  - HYP-3056
  - HYP-3055
  - HYP-3054
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-3045
  - HYP-3043
  - HYP-3041
  - HYP-3040
  - HYP-3039
  - HYP-3038
  - HYP-3037
  - HYP-3022
  - HYP-3021
  - HYP-2991
  - HYP-2989
  - HYP-2120
  - HYP-2121
  - THM-381
  - THM-385
  - OPEN-Q-108
---

# HYP-3065: Observer-extension exact duodecimal deletion-sector audit

## Claim

Layered over HYP-3054's observer-extension cut-payload calculus, HYP-3055's
duodecimal first-defect ledger, HYP-3056's observer-cut orbit ledger,
HYP-3057's value-origin ledger, and the sibling HYP-3058 warning that
attractive sidecar numbers still need proof payloads, with HYP-3059 as the
incoming twelve-layer cut-payload stress test, the `48`, `56`, and `12`
surface around the first rooted-perspective failure should be read as an exact
observer-extension overlap audit, not as raw addition.

The literal count is:

```text
48 + 12 = 60, not 56.
```

The exact small-tournament identities are:

```text
P(4)  = 12  rooted/node perspectives on 4-tournament classes
U(5)  = 12  unrooted 5-tournament classes
SC(6) = 12  self-converse 6-tournament classes
K_{4,5} rectangle redundancy = 12
```

But at the first failure,

```text
P(5) = 48
U(6) = 56
U(6) - P(5) = 8.
```

So the theorem-shaped ledger is:

```text
U(6) = P(5) + 8
     = P(5) + (2/3)*12
     = P(5) + SC(6) - U(4)
     = 48 + 12 - 4.
```

The hidden term is the overlap kernel `U(4)=4`, one third of a dozen.  Thus
the user's "48 plus a dozen gives 56" clue is structurally right only after
controlled forgetting exposes the overlap: a dozen-sized self-converse /
rectangle-residue reservoir contributes net eight because four base classes
are already counted by the rooted cache.

This is currently a numerical identity plus a carrier target, not a proved
bijection.  The next theorem target is to construct the actual inclusion-
exclusion object whose overlap is `U(4)`.

## Exact Audit

The S230 script recomputes tournament isomorphism classes through `n=6` by
canonical labelling and Burnside:

```text
n  U(n)  rooted P(n)  extension E(n->n+1)  self-converse
1     1            1                     2              1
2     1            2                     4              1
3     2            4                    12              2
4     4           12                    48              2
5    12           48                   296              8
6    56          296                                   12
```

The fractions at the first failure are:

```text
P(5)/U(6)      = 6/7
defect/U(6)    = 1/7
defect/dozen   = 2/3
overlap/dozen  = U(4)/12 = 1/3
U(6)/dozen     = 14/3
P(5)/dozen     = 4
```

So the count `56` is not four dozen plus one dozen.  It is four dozen plus a
two-thirds-dozen correction, or equivalently four dozen plus one dozen minus a
one-third-dozen overlap.

## Burnside Source

The exact Burnside terms explain where the numbers enter.

For `n=5`:

```text
type=(5,)             perms=24  fixed=4
type=(3,1,1)          perms=20  fixed=16
type=(1,1,1,1,1)      perms=1   fixed=1024
total/120 = 12
```

For `n=6`:

```text
type=(5,1)            perms=144 fixed=8
type=(3,3)            perms=40  fixed=32
type=(3,1,1,1)        perms=40  fixed=128
type=(1,1,1,1,1,1)    perms=1   fixed=32768
total/720 = 56
```

The `[3,3]` term is important: it has no fixed vertex, yet contributes
`40*32` labelled fixed shadows to `U(6)`.  This matches the older S211/S212
reading that the first `P(5)->U(6)` defect is rootless/cyclic, not deeper
node memory.

## Deletion And Sector Payloads

The self-converse branch is not isolated from the ordinary rooted/deletion
ledger.  Among the `56` six-classes, S230 finds:

```text
self_converse: 12 classes, rooted_mult_sum=60
chiral:        44 classes, rooted_mult_sum=236
total rooted 6-count = 296
```

Deletion-parent counts split differently on the two sides:

```text
self_converse deletion_parent_count_hist = {1:3, 2:2, 3:1, 4:5, 5:1}
chiral        deletion_parent_count_hist = {2:4, 3:12, 4:8, 5:16, 6:4}
```

HYP-3049's ordered-pair sector fact is reverified:

```text
sector_size_deck          separates 55/56
sector_size_internal_deck separates 55/56
sector_full_cross_deck    separates 56/56
```

The collision before cross-sector orientation is exactly the converse pair
`344/345`.  Thus the compact observer-extension payload after rooted memory is
again `cross_sector_orientation_word`.

## Rectangle / Hourglass Reading

HYP-3053 supplies a second source of the dozen.  In the diagonal flow bridge
between layers of sizes `k` and `k+1`,

```text
lines = k(k+1)
rank = 2k
rectangle redundancy = k(k-1).
```

For `k=4`, the redundancy is exactly:

```text
K_{4,5}: lines=20, rank=8, rectangle_redundancy=12.
```

This dozen is not class mass.  It is a rectangle-cycle reservoir.  For full
adjacent-layer flow,

```text
n  lines  rank  redundancy  local_rectangles  hourglass_residues
4      8     5           3                 2                   1
5     20     9          11                 8                   3
6     40    14          26                20                   6
7     70    20          50                40                  10
```

The fixed-path diagonal-flow rule says to keep rectangle and hourglass cycle
residues as proof payloads.  If they vanish, the quotient descends to layer
potentials.  If they do not vanish after LRC labels are attached, the residue
is a named hidden coordinate.

## Observer-Extension / Cut Payload Abstract

The common abstraction is:

```text
observer-extension creates a cut payload;
controlled forgetting may remove the observer only when that payload is
fiber-constant, reconstructible, dual-annihilated, descended familywise,
or emitted as named residual debt.
```

Concrete repo applications:

```text
HYP-3057 value-origin ledger:
  quotient = small integer reused across adjacent tournament levels
  payload  = value_origin_type plus the coordinate it came from

HYP-3056 observer-cut orbit ledger:
  quotient = visible-fiber automorphism quotient
  payload  = orbit of boundary slice / incidence word / extended shadow

HYP-3039 controlled forgetting:
  quotient = coarse quotient fiber
  payload  = hidden coordinate stage / residual debt

HYP-3038 q=23 drop/add square:
  quotient = fixed-margin rectangle
  payload  = exact-M zeta and endpoint-owner strip

HYP-3037 residual capacitors:
  quotient = two-plate route cut
  payload  = first-cut stage / capacitor id

HYP-3041 AP-tail repair:
  quotient = AP-core one-tail observer
  payload  = q=13 puncture or reciprocal fixed-point clock

HYP-3045 endpoint-owner transfer:
  quotient = B18Z6 endpoint shadow
  payload  = external owner strip and transfer delta

HYP-3021/HYP-3022 pair-good decoys:
  quotient = blocker generator quotient
  payload  = barcode / normal-fan active owner

HYP-3049 ordered-pair sectors:
  quotient = old-root/new-observer cut
  payload  = cross-sector orientation word

HYP-3053 diagonal flow:
  quotient = K_{k,k+1} line quotient
  payload  = rectangle/hourglass cycle residue
```

This is the duodecimal synthesis: the dozen keeps recurring when the quotient
has just enough symmetry to look closed, but the actual proof object is the
overlap/residue/cut payload that says what the quotient forgot.

## Tournament Analysis

S230 uses proof-carrier vertices, not runners.

Pairwise observable:

```text
retained LRC predicate payload minus quotient loss.
```

Switch:

```text
orient toward the carrier that explains an actual lost coordinate; exact ties
follow the sidecar order below.
```

Carrier tournament:

```text
exact_canonical_audit
> burnside_odd_cycle_sidecar
> cross_sector_orientation_word
> deletion_parent_fiber_profile
> rectangle_hourglass_cycle_residue
> self_converse_branch_locus
> duodecimal_overlap_ledger
> rooted_perspective_cache
> raw_orbit_count_coincidence
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

## Assumption Challenge

Vertices considered:

```text
runners, arcs, rooted vertices, ordered pairs, edge sectors, diagonal tiles,
K_{k,k+1} line sectors, rectangle cycles, hourglass cycles, deletion fibers,
self-converse branch points, fixed circle sections, section boundaries,
wall-crossing events, cover arcs, barcode bars, endpoint owners, residues,
Fourier modes, matroid cocircuits, and proof obligations.
```

Chosen vertices are observer-extension proof carriers and their cut payloads.

Preserved predicate:

```text
whether the quotient keeps enough information to separate or reconstruct LRC
route/status-changing packet pairs.
```

Destroyed information:

```text
labelled runner identity, old-root/new-observer role, cross-sector
orientation, deletion-parent multiplicity, endpoint-owner names, and
rectangle/hourglass residues when they are not sidecars.
```

Challenged assumption:

```text
the dozen should not be added as independent proof mass.  At this scale the
live statement is the duodecimal overlap law 56 = 48 + 12 - 4, with net
correction 8 = (2/3)*12.
```

## Next Pulls

1. Construct an actual inclusion-exclusion object for
   `U(6)=P(5)+SC(6)-U(4)` or refute it as only a numerical identity.
2. Emit rectangle and hourglass cycle bases whose residues can be compared
   against the ordered-pair sector collision `344/345`.
3. Attach `observer_extension_cut_signature`, `duodecimal_overlap_kernel`,
   `value_origin_type`, `observer_cut_payload_orbit`,
   `self_converse_branch_locus`, `cross_sector_orientation_word`,
   `deletion_parent_profile`, `rectangle_cycle_residue`, and
   `hourglass_cycle_residue` to LRC packet experiments.
4. Test the observer-extension abstraction on residual capacitors,
   endpoint-owner transfer, q=23 drop/add squares, AP-tail q=13 repair, and
   pair-good decoy generator decks before promoting another scalar quotient.
5. Decide whether the overlap kernel is literally `U(4)` as a contraction/
   deletion boundary, a self-converse fixed-set overlap, or a rectangle-cycle
   cohomology kernel.
