---
id: THM-2720
title: "Unshifted deepest-source and canonical present-packet global disjointness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  typed row, every source-one present packet
  F_(ell,s) contains the unshifted c3-safe factor, while the exclusive
  deepest source E3 contains the complementary unshifted c3-danger factor.
  Hence E3 intersect F_(ell,s) is empty for every one of the 7*13 labels,
  and remains empty after arbitrary rail, phase, clock, carry, root, or unit
  restrictions.  This is uniform set-theoretic disjointness for the
  unchanged present grammar, not a scalar-row exclusion or LRC(14) result.
source: root/unshifted-deepest-source-present-wall-2026-07-28
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
related:
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2742-full-two-target-present-sheet-deepest-source-semantic-current
script: 04-computation/lrc14_unshifted_deepest_source_present_wall_thm2720.py
output: 05-knowledge/results/lrc14_unshifted_deepest_source_present_wall_thm2720.out
script_sha256: b8e7619f2ca56759794b90794e936935c0f41012206cfeacf3f869e25ed3c0db
output_sha256: fd3b7f8d47efbe0faa43cb8db168eccb2d4ed7671b667593faee5a3020123456
hash_basis: working-tree bytes (LF)
audit: thm2705-2709-audit-2026-07-28 (independent proof audit; endpoint hostile check; independent 91-cell reconstruction; normal/-O/hash/docs replay)
---

# THM-2720 -- the canonical present grammar cannot host a deepest-owner source

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The full half/`C_221` census found positive scalar `E_3` cospans but no
attachment to a canonical present packet.  The reason is not a rare phase or
root accident.  It is one complementary pair of literal factors already
built into the definitions.

## 1. The two incompatible factors

Work on the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5).                 (1)
```

Write `D_v(theta)` for the danger comb centred at `theta` and
`G_v(theta)` for its safe complement.  Let `A0` be the common guard/unit
safe set.  The exclusive deepest-owner source of THM-2305 is

```text
E3=A0 intersect D_c3(0) intersect G_c1(0) intersect G_c2(0).
                                                                  (2)
```

In particular,

```text
E3 subset D_c3(0).                                      (3)
```

For `ell in Z/7` and `s in F_13`, let `F_(ell,s)` be the canonical
source-one present packet used by the lawful response tables and all of the
later rail carriers.  Its literal interval definition has the form

```text
F_(ell,s)=B_(ell,s) intersect G_c3(0),                 (4)
```

where `B_(ell,s)` contains the shifted `c1` source, guard and five unit
safeties, and the target-shifted `q1/c2` factors.  The last factor in `(4)`
is unshifted because the present target slice has deepest coordinate `t=0`.
It is exactly the final `c3`-safe cut in the canonical `build_F(ell,s)`
definition.

Equations `(3)` and `(4)` give, for every label pair,

```text
E3 intersect F_(ell,s)
 subset D_c3(0) intersect G_c3(0)
 =empty.                                                 (5)
```

This proves the asserted disjointness.  No counting, genericity, or sampled
phase is used.

## 2. Every unchanged sidecar re-root remains empty

Suppose a refined carrier vertex is obtained by adding any combination of
rail, phase, clock, carry, private-root, primitive-unit, or affine-arrow
conditions to the same present packet.  Its physical support `C` satisfies

```text
C subset F_(ell,s).                                     (6)
```

Therefore `(5)` implies `C intersect E3=empty`.  Changing a half phase,
choosing another carry, replacing root `6` by roots `1` or `12`, or changing
a rail cannot repair the source vertex while the factor `G_c3(0)` is kept.

This upgrades the finite observations in the half/`C_221` source-fibre
census to a uniform stopping theorem for the **unchanged** present grammar.
It does not make the raw semantic cospan empty: strict open instances of

```text
E3 intersect D^(-6) Q_(3,{1,2})                        (7)
```

exist in both full fibres.  It proves only that such an instance cannot also
be a source vertex of a carrier whose support is contained in `(4)`.

## 3. The wall is sharp under deletion

The exact companion reconstructs all `7*13=91` packets and factors each one
as in `(4)`.  It verifies all `91` intersections in `(5)` are empty.  As a
hostile control, it deletes **only** the final `c3`-safe factor and recomputes
the intersection with `E3`.  Exactly

```text
78 of 91 label pairs become positive,                       (8)
```

with exact unnormalized grid masses between

```text
528530047500 and 1406405080080.                            (9)
```

The first positive deletion control is `(ell,s)=(1,0)`, of grid mass
`1060449843180`.  Thus the wall is load-bearing rather than a redundant
factor.  The remaining thirteen empty deletion rows also show that removing
this factor is necessary here, not by itself a universal construction.

## 4. Consequence and scope

For the deepest-owner lane the next operation must alter the object, not the
remaining sidecars.  One must delete, shift, or relativize the unshifted
`c3`-safe present factor and then rebuild the rail/unit/current data on the
resulting support.  A present-free or relative-present mapping cone is the
minimal honest formulation.  The old packet cannot simply be relabelled as
an `E3` packet.

Nothing here supplies the rebuilt primitive unit, a nonnegative endpoint
current, a common middle fibre, a scalar-row exclusion, or LRC(14).  The
result is confined to row `(1)` and the canonical `F_(ell,s)` grammar.

## 5. Exact reproduction

Run

```bash
python3 04-computation/lrc14_unshifted_deepest_source_present_wall_thm2720.py
python3 -O 04-computation/lrc14_unshifted_deepest_source_present_wall_thm2720.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_unshifted_deepest_source_present_wall_thm2720.out
```

The exact interval logic uses no truth-bearing Python `assert`.  It pins the
audited canonical present-packet implementation by SHA-256,
reconstructs `E3`, checks the exact danger/safe partition, factors all
packets before and after the final `c3`-safe cut, and supplies the deletion
controls `(8)--(9)`.

An independent hostile audit rederived the pointwise inclusion
`E3 subset D_c3(0)` and the literal carrier factor
`F_(ell,s) subset G_c3(0)`, checked that both half-open and strict-open
endpoint conventions leave the disjointness unchanged, and independently
rebuilt all `91` labelled cells.  It reproduced the exact `78/91` deletion
control, its mass range, and both extrema in `(9)`.  A second immutable audit
then found that the promoted file's pinned carrier hash and its two declared
artifact hashes were stale despite the first audit's hash claim.  The carrier
pin was repaired to the current dependency hash before replay; normal and
optimized executions now byte-match the stored transcript, the LF-normalized
artifact hashes are those declared in the frontmatter, and the documentation
check passes.  This was an evidence repair, not a change to the set-theoretic
theorem.

QED.
