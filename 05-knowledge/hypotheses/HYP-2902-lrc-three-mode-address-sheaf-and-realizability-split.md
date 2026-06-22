---
id: HYP-2902
status: STRUCTURAL PROOF-TARGET / exact correction and guardrail
source: codex-2026-06-22-S114
tags: [lrc14, tournament-recursion, half-tiling, legendre, eisenstein, mobius, exact-period, finite-atlas, weyl, OPEN-Q-108]
depends_on:
  - HYP-2901
  - HYP-2899
  - HYP-2900
  - HYP-2886
  - HYP-2889
  - HYP-2876
  - THM-565
related:
  - THM-549
  - THM-550
  - THM-566
  - HYP-2856
  - HYP-2870
  - HYP-2875
  - HYP-2890
  - HYP-2895
  - OPEN-Q-108
results:
  - 04-computation/lrc_legendre_venn_correction_kps.py
  - 05-knowledge/results/lrc_legendre_venn_correction_kps.out
  - 04-computation/lrc_three_mode_address_sheaf_codex_s114.py
  - 05-knowledge/results/lrc_three_mode_address_sheaf_codex_s114.out
---

# HYP-2902: the three recursion modes form a parity-local address sheaf, not a scalar recurrence

This hypothesis refines HYP-2901's denominator-wall guardrail with the address
sheaf interpretation of the owner's correction and KPS S31s.  The corrected
recurrence data should be read as a local address system:

```text
Mobius full mode      : all sizes, full inclusion-exclusion skeleton
Eisenstein half mode  : even sizes, pronic/complement fold
Legendre half mode    : odd sizes, square three-corner geometry
```

The letters `A..G` are therefore not global objects.  They are local slots in
the current mode, and their subtournament sizes differ by mode.

## Exact slot atlas

The S114 audit records the scalar size projections:

```text
Mobius full       : {N-1: +3, N-2: -3, N-3: +1}
Eisenstein even   : {N-1: +2, N-2: -1}
Legendre odd      : {N-1: +2, N-2:  0, N-3: -2, N-4: +1}
```

The corrected Legendre odd slots are:

```text
A,B : N-1
C,D : N-2
E,F : N-3
G   : N-4
```

with geometry:

```text
corners: A, D, B
edges:   A+D-E,  B+D-F,  A+B-C
center:  A+B+D-C-E-F+G
```

This is the crucial correction.  `C` and `D` cancel in scalar cardinality
because both are size `N-2`, but they are different addresses.  `C` belongs to
the `AB` edge; `D` is a corner and also feeds the `AD`, `BD`, and center
regions.  Any proof that collapses `C-D` too early loses the geometry the
owner's half-tiling correction is trying to preserve.

## LRC14 composition

For `N=14`, the half-tiling is even/pronic:

```text
h(14)=42=7*6.
```

The top half-mode is Eisenstein:

```text
A,B size 13; C size 12.
```

The apex coordinate `7` is not a child size in that recurrence.  It is the
pronic fold parameter `N/2`, where the odd Legendre character coordinate
appears.  Thus the composition has two size coordinates:

```text
recurrence children of N=14: 13 and 12
fold/apex coordinate:        7
```

This prevents a common false shortcut.  The odd Legendre chart at apex size
`7` has its own local slots:

```text
A,B size 6; C,D size 5; E,F size 4; G size 3.
```

The odd half chart at `N=15` has:

```text
A,B size 14; C,D size 13; E,F size 12; G size 11.
```

So `14=2*7` should be read as a parity-stratified address sheaf:

```text
Mobius always
  + Eisenstein even/pronic quotient at N=14
  + Legendre odd/square quotient on the exposed apex coordinate
  + exact-period packet labels before scalar cap/floor.
```

## Finite-atlas guardrail

S114 also tests the prompt's divisor-loaded family

```text
S_X = {1,...,11,13,lcm(2,...,X)}.
```

For every `D<=X`, the last speed is `0 mod D`, hence

```text
N(S_X,D)=0.
```

This gives the theorem-level obstruction needed to refute any fixed finite
denominator basis.  The stronger wording "the first witness is nextprime(X)"
is false in this raw family.  Exact scans found:

```text
X=14: first witness D=41, not 17
X=24: first witness D=53, not 29
X=41: first witness D=53, not 43
X=90: first witness D=97
```

The correct statement is:

```text
q_min(S_X) > X,
```

with the first opening governed by the AP-core safe-denominator ladder after
the divisor-loaded tail has killed all small denominators.  This is exactly the
finite/analytic split already visible in HYP-2876 and HYP-2886: fixed residue
bases are useful atlases, not closures.

## Proof-route implication

The current LRC14 architecture should split as follows.

Node 2 is finite/algebraic:

```text
three-gap / AP-hull / anchored Fejer majorization
```

but it must retain sector and address labels.  Scalar additive-energy
monotonicity, local compression, and unlabelled AP comparisons have already
failed in HYP-2889 and related audits.

Node 3 is analytic/exact-period:

```text
delete divisibility-killed denominators
keep unit packets a mod D and CRT/chi7/affine labels
prove a signed low-frequency floor R' >= c > 0
route coherent packets to AP/GAP finite atlases
route incoherent packets to L2/Dirichlet/Weyl tails
then use THM-565 to convert the floor into finite witnesses.
```

The three recurrences connect these branches by addresses.  They are not a
closed scalar recurrence for witness counts, coprime density, or `p0`.

## Tournament Analysis

Vertices are proof carriers:

```text
exact_period_packets
> address_sheaf
> legendre_CD_seam
> eisenstein_pronic_fold
> mobius_scalar_skeleton
> three_gap_AP_hull
> raw_runner_vertices
```

Pairwise observable: preserved labels, common labels, and declaration order.
The S114 tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
Hamiltonian path as above.
```

Assumption challenged: the scalar half-tile recurrence preserves the LRC
predicate.  It does not.  The quotient must retain exact-period packets, local
slot labels, and parity mode until the final cap/floor scalarization.

## What this does not prove

This is not a proof of LRC14.  It improves the proof order:

```text
address-labelled recurrences
  -> exact-period packet floor
  -> AP-facing finite majorization
  -> signed residual leak / L2 tail
  -> THM-565 finite-ruler conversion.
```

The sharp remaining theorem remains the uniform Node-3 floor after coherent
low packets are separated from the analytic tail.
