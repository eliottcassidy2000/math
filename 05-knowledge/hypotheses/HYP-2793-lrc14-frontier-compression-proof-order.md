---
id: HYP-2793
title: LRC14 frontier compression leaves genuine-wide slack plus formal glue after bounded and single-far closure
status: UPDATED proof order; bounded and single-far legs computationally closed
source: codex-2026-06-21
depends_on:
  - THM-563
  - HYP-2788
  - HYP-2790
  - HYP-2792
  - THM-557
  - THM-548
  - THM-546
related:
  - HYP-2791
  - HYP-2780
  - HYP-2781
  - HYP-2775
  - HYP-2701
  - HYP-2684
  - OPEN-Q-108
---

# HYP-2793: LRC14 Frontier Compression Proof Order

## Claim

The current LRC14 proof frontier should no longer be treated as one large
search space.  After the bounded-leg certificate, the THM-563 period-max
identity, and the HYP-2788 single-perturbation dichotomy, the proof is now a
three-leg DAG with two computationally closed legs:

```text
bounded span <= 14
  -> exact cap certificate and formal split completeness

single-perturbation / single-far
  -> CLOSED by full bounded-base periodmax(B) < 15*(cap_k - Plat(B)) audit

genuine-wide / at least two far from every bounded scaffold
  -> scale-free slack floor or survival-currency signed deviation bound
```

The useful next step is not another scalar search.  It is to formalize the two
closed legs and close the remaining genuine-wide/slack leg in its own
relation-lattice or survival-currency coordinates.

## Incoming Evidence Integrated

### 1. Bounded Leg Is Computationally Closed

Incoming result:

```text
05-knowledge/results/lrc14_THREAD3_bounded_leg_certificate_macmini.out
```

gives three independent engines agreeing on every `0 in E` shape, zero
scale-invariance failures, and an exhaustive exact bounded-span certificate for
primitive sets with `span <= 14`.  The cap check covers the LRC14 sizes
`k=8..12`; every row is below cap.  The binding bounded row is k=8 consecutive:

```text
k=8 max = 481/1470 at (0,1,2,3,4,5,6,7)
cap_8 = 2243/5880
margin = 319/5880
```

The same file also records split completeness: the reduction routes exactly to
`BOUNDED` iff primitive reduced span is at most `14`, and to `WIDE` iff it is
greater than `14`; no gap or overlap was found.  This is now a formalization
target, not a search target.

Thread 5 independently rechecked the beginning of this bounded leg in

```text
05-knowledge/results/lrc_bounded_full_recheck_thread5.out
```

with the same k=8 binding row.

### 2. Single-Far Is Closed By A Complete Period-Max Ledger

`THM-563` rewrites the one-far correction as an exact periodic endpoint sum:

```text
w*Delta_w(B) = sum_j sum_{endpoint t of A_j(B)} +/- S_j({w*t})
period = 7*lcm(B)
```

so the branch closes for a bounded base `B` once

```text
periodmax(B) < 15 * (cap_k - Plat(B)).
```

`HYP-2790` tested the tempting Boolean/type bridge and found it weak: among
`135` high-plateau bounded bases k=8..12, there were zero periodmax failures,
but the Boolean/type and containment deficits do not stably rank the pressure.
The persistent coordinate is the endpoint-period numerator itself.

The latest incoming complete ledger is:

```text
05-knowledge/results/lrc_periodmax_THM563_general_check_COMPLETE_macmini_0621s7.out
```

It checks every primitive bounded base `B subset [0,14]` with `0 in B`,
`|B|=k-1`, for every `k=8..13`:

```text
total bases = 12805
failures = 0
skipped = 0
periodmax(B) < 15*(cap_k - Plat(B)) everywhere
```

Per-k summary:

```text
k=8:  3003/3003 PASS, worst exact ratio 10.8188 at (0,2,4,6,8,10,12)
k=9:  3432/3432 PASS, worst exact ratio 13.2805 at (0,2,4,6,8,10,12,14)
k=10: 3003/3003 PASS, worst exact ratio 10.8140 at (0,2,4,6,8,10,11,12,14)
k=11: 2002/2002 PASS, worst exact ratio  9.7616 at (0,3,7,8,9,10,11,12,13,14)
k=12: 1001/1001 PASS, worst exact ratio  8.2059 at (0,5,6,7,8,9,10,11,12,13,14)
k=13:  364/364 PASS, worst exact ratio  5.9838 at (0,3,5,6,7,8,9,10,11,12,13,14)
```

The global binding base is the k=9 even AP:

```text
B = (0,2,4,6,8,10,12,14) = 2*consec_7
Plat = 621/1715
margin = 129643/980980
periodmax = 86/49
15*margin - periodmax = +0.22725...
ratio = 13.2805 < 15
```

Therefore the single-far bounded-base leg is closed by exact finite check.
The HYP-2792 signed endpoint/Dedekind reciprocity theorem remains valuable as
proof compression and conceptual explanation, but it is no longer required to
finish this branch computationally.

### 3. Genuine-Wide Should Be Proved In Its Own Currency

The HYP-2775/HYP-2788 lesson is that independent suprema are false currency:

```text
decorr_sup + err_sup <= cap
```

is too strong, because large decorrelated value and large resonance error occur
in different row families.  The correct comparison is pointwise:

```text
err(E) < cap_k - p0_decorr(E).
```

Equivalently, the proof should split scale-reducible rows from genuine-wide
rows.  The live genuine-wide carrier is not a runner tournament; it is a
relation-lattice / scale-cluster object.  The strongest current currencies are:

```text
survival middle mass:
  C = p1+p2+p3+p4-4*p6 = 1-U4

two-far boundary death-chain:
  only depths 1,2,5,6 are live for the signed deviation

Freiman scale reducibility:
  low relation-height far pairs route to finite atlases
  high relation-height pairs should obey signed Abel/BV cancellation
```

This is the same structural lesson as HYP-2790: keep the signed phase or
relation coordinate until the final cap comparison.

## Sharp Proof Obligations

### Obligation A: bounded split and cap table

Formalize:

```text
primitive reduced span <= 14  <=>  BOUNDED route
primitive reduced span > 14   <=>  WIDE route
```

and import or rederive the exact bounded cap ledger for k=8..12.  The
formalization should include the canonical cap table:

```text
cap_8..13 = 2243/5880, 1979/4004, 55/91, 66/91, 6/7, 1.
```

### Obligation B: bounded-base periodmax theorem

The finite exact audit is complete.  Formalize or import:

```text
for all primitive B subset [0,14], 0 in B, |B|=k-1:
periodmax(B) < 15*(cap_k - Plat(B))
```

The signed endpoint-packet theorem is now an optional compression target.  It
should be written in terms of the endpoint orbit under `w`:

```text
F_B(w) = sum_r eps_r S_{j_r}({w*a_r/q_r}).
```

External calibration: generalized Dedekind/Hardy sums defined by sawtooth or
Bernoulli functions admit finite trigonometric representations via the
discrete Fourier transform; see Rassias-Toth
`https://arxiv.org/pdf/1512.01466`.  Shifted Dedekind-Rademacher sums and
multiple Dedekind sums also have reciprocity laws; see
`https://people.mpim-bonn.mpg.de/zagier/files/acta/73-4/fulltext.pdf` and
`https://www.numdam.org/item/10.5802/aif.2261.pdf`.

The LRC-specific theorem is not quoted from these papers.  The inference is
that the endpoint packet can probably be converted to a signed cotangent or
Dedekind-Rademacher packet and bounded by reciprocity/arc pairing.  This would
explain the finite certificate and may simplify a Lean or paper proof.  A
generic Erdos-Turan or Koksma absolute discrepancy estimate destroys the
cancellation and returns to the 125x loss.

### Obligation C: genuine-wide slack theorem

Prove a true-wide branch statement in one of two equivalent forms:

```text
p0(E) <= Q(k-1) + signed_error(E) < cap_k
```

with a pointwise room-vs-error comparison, or

```text
C(E) = p1+p2+p3+p4-4*p6 >= 1 - cap_k
```

with signed two-far and r>=3 margins handled separately.  The two-far proof
should keep the live-depth kernel `{1,2,5,6}` and relation-height data visible.

### Obligation D: Lean nucleus

The Lean work should start with the pieces that are pure finite arithmetic:

```text
cap constants and order relations
death-chain live-depth identities
bounded split definitions
periodmax inequality statement shape
```

The full row atlas can remain Python-generated until the theorem boundary is
stable.

### S77 addendum: periodmax arithmetic kernel formalized

Codex S77 added a focused certificate kernel for Obligation B/D:

```text
04-computation/lrc_periodmax_worstrow_certificate_codex_s77.py
04-computation/lean/TournamentH7/TournamentH7/LRCPeriodmaxCertificate.lean
```

The Python kernel recomputes the exact `Plat(B)` margins for the six per-k
worst rows of the completed THM-563 audit, brute-verifies the k=8 and k=9
endpoint-period maxima, and records exact headrooms:

```text
k=8   303/392
k=9   44585/196196
k=10  464839/535080
k=11  5417609/2942940
k=12  2664689/1177176
k=13  2207057/588588
```

The Lean module is intentionally not the 12805-row exhaustive proof.  It is the
small arithmetic nucleus surrounding the completed finite audit: all six
worst rows have positive headroom, the k=9 even AP is the largest ratio among
the six per-k worst rows, the count checksum is
`12805` bases / `3995` trivial / `8810` scanned / `0` skipped / `12805` passed /
`0` failed, and the k=8 dilated AP normalization is `periodmax=2`.
Focused build:

```text
lake build TournamentH7.LRCPeriodmaxCertificate
```

succeeds and is stored at
`05-knowledge/results/lrc_periodmax_certificate_lean_codex_s77.out`.
This moves the integer single-far bounded-base leg from "computed result" to a
Lean-facing arithmetic import boundary.  It does not touch the remaining
genuine-wide theorem or the continuous dilation handoff.

## Tournament Analysis

Vertices tested as proof carriers:

```text
bounded_span_certificate
singlefar_endpoint_periodmax
genuinewide_relation_lattice
survival_middle_mass
finite_lowheight_atlas
boolean_type_cut
raw_erdos_turan_koksma
runner_vertices
```

Pairwise observable:

```text
preserves exact LRC predicate,
keeps signed phase or relation data,
has a closed certificate path,
is Lean-local enough to formalize,
does not merge known opposite-risk families.
```

The resulting risk ranking is transitive:

```text
bounded_span_certificate
> singlefar_endpoint_periodmax
> genuinewide_relation_lattice
> survival_middle_mass
> finite_lowheight_atlas
> boolean_type_cut
> raw_erdos_turan_koksma
> runner_vertices.
```

This ranking is not saying the bounded leg is "hardest"; it is saying it is the
cleanest closed proof object.  After the complete THM-563 check,
`singlefar_endpoint_periodmax` is also closed computationally; its remaining
role is formalization or proof compression.  The live mathematical work is now
concentrated in `genuinewide_relation_lattice` / `survival_middle_mass`, plus
the formal glue connecting the closed finite legs to the dichotomy.

## Assumption Challenge

The old default assumption that vertices should be runners is now actively
harmful.  Runner vertices preserve speed ownership but destroy the quotient
that the proof actually uses.  Bounded bases preserve the single-far predicate
but forget genuine-wide interactions.  Endpoint arcs preserve the periodmax
numerator but forget the cap plateau.  Scale clusters and relation equations
preserve true-wide resonance but forget bounded finite cap data.

The current proof should therefore use proof obligations themselves as the
top-level vertices, with local vertex choices inside each obligation:

```text
bounded leg      -> bounded bases / cap rows
single-far       -> endpoint arcs and w-orbits
genuine-wide     -> scale clusters, relation equations, live-depth packets
formal glue      -> exact finite arithmetic objects
```

That quotient preserves the information needed for each branch and avoids
forcing a single lossy scalar across incompatible regimes.
