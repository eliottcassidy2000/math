---
id: HYP-3146
title: The n=4 tournament tables expose a filler/canary shift-package quotient useful for LRC packet design
status: SYNTHESIS / exact n=4 scout and proof-packet proposal; not a proof
source: codex-2026-06-27-S274
tangent: T1211
technique: LTI-272
tournament_technique: LTT-170
script: 04-computation/tournament_n4_shift_package_erdos870_codex_s274.py
result: 05-knowledge/results/tournament_n4_shift_package_erdos870_codex_s274.out
reflection: 07-reflections/n4-shift-package-erdos870-filler-canary-codex-s274.md
external:
  - https://github.com/davidturturean/erdos-870
related:
  - HYP-3145
  - HYP-3143
  - HYP-3142
  - HYP-3141
  - HYP-3140
  - HYP-3139
  - HYP-3137
  - HYP-3134
  - HYP-3133
  - HYP-3124
  - HYP-3054
  - HYP-3053
  - HYP-3049
  - OPEN-Q-108
---

# HYP-3146: n=4 Shift Packages and Erdos-870 Filler/Canary Logic

## Claim

The two user-proposed `n=4` tournament tables are not merely two notations for
the same object.  After the S276/S274 rebases, HYP-3143 carries the exact-order
subbasis audit, HYP-3145 names the filler-core interface, and HYP-3146 carries
the complementary canary/scaffold policy.  The packet chain exposes a
proof-design distinction that matches the recently formalized negative answer
to Erdős Problem #870:

```text
raw fixed-path cover with redundant representations
  versus
finite-filler scaffold that turns the quotient into a shift package.
```

In `davidturturean/erdos-870`, the transferable architecture is:

```text
sparse order-2 core
+ finite filler gadget
+ finite shift packages
+ clustered canaries
+ deletion-stable nonminimality.
```

The LRC/tournament analogue is:

```text
fixed-path half-tiling cube
+ finite arc scaffold
+ two-bit class shift package
+ S-fiber canary cluster
+ controlled-forgetting/deletion-stability guardrail.
```

This does not make additive bases directly into tournament classes.  It gives
a disciplined way to decide when to keep redundant fiber mass and when to add
finite scaffold data so a quotient becomes proof-usable.

## Exact n=4 Scout

With fixed Hamiltonian path `0->1->2->3`, the free chords can be named

```text
a=(0,2), b=(1,3), c=(0,3).
```

The restricted generator table is exactly the user's tiling table:

```text
      E   a   b   c
E     T   +   -   S
a     +   T   S   S
b     -   S   T   S
c     S   S   S   T
```

But the full fixed-path cube has eight states, with class fibers:

```text
T: 1
+: 1
-: 1
S: 5
```

More precisely:

```text
S = {c, ab, ac, bc, abc}.
```

So the fixed-path table is a cover, not a group quotient.  Its class
generating polynomials by flip weight are:

```text
T(z)=1
+(z)=z
-(z)=z
S(z)=z + 3z^2 + z^3.
```

The `S` fiber is the canary cluster.  The representative `abc` is
delete-one-stable inside `S` because deleting any one of `a,b,c` leaves a
two-flip `S` representative.  The singleton `c` and two-flip representatives
are not delete-one-stable.  Thus redundancy is not noise; it is the local
nonminimality resource.

The second model fixes four arcs with partial outscore vector `(2,1,1,0)`,
sorted `(0,1,1,2)`, and leaves endpoint variables

```text
x=(0,1), y=(2,3).
```

One exact scaffold is:

```text
0->2, 0->3, 2->1, 1->3 fixed.
```

It gives the user's two-bit table:

```text
      E   x   y
E     T   +   -
x     +   T   S
y     -   S   T
```

After naming `xy=S`, this is the Klein four table:

```text
      E   x   y  xy
E     T   +   -   S
x     +   T   S   -
y     -   S   T   +
xy    S   -   +   T
```

So the finite scaffold turns the four class labels into an exact shift package.
It sacrifices the fixed-path `S` redundancy but gains quotient legality.

## The Hidden Compression

The class-preserving map from the fixed-path cover to the two-bit scaffold is
not linear.  It is the monotone Boolean compression

```text
x = a OR c,
y = b OR c.
```

The chord `c` is the clustered canary: activating it forces both endpoint
variables live.  The pair flip `ab` also forces both variables live.  This is
why all mixed or canary states collapse to `S`.

In controlled-forgetting language:

```text
fixed-path quotient:
  keeps representation multiplicity and deletion robustness,
  but is not a congruence;

finite-filler scaffold:
  removes the hidden fiber mass,
  but turns the class labels into a true shift package.
```

The proof question becomes: for the next operation, do we need deletion
stability or quotient legality?

## LRC14 Transfer

HYP-3141 says directed tournament edges must carry tail payload, tip payload,
observer-cut orbit, commutator defect, coordinate-resurrection cover, and
terminal exit.  HYP-3142 says the bounded-core `k=8` node now has an exact
`U4` exit candidate.  HYP-3140 says the multi-far `Rprime` scalar is really a
14-sheet conditional first moment.  HYP-3053/HYP-3054 say fixed-path tilings
and observer cuts are legal only after rectangle/hourglass and cut payloads
are named.

HYP-3146 adds the finite quotient recipe on top of HYP-3143's
exact-order/no-lower-order audit and HYP-3145's filler-core boundary:

```text
if a quotient has a large unsafe fiber,
  either keep it as a canary cluster when deletion stability is needed,
  or add finite filler/scaffold data until the quotient is a shift package.
```

For LRC packet manifests this suggests new fields:

```text
shift_package_scaffold_id
fixed_path_cover_fiber_pgf
canary_cluster_fiber
delete_one_stable_representative
monotone_or_compression_word
finite_filler_arc_set
quotient_congruence_status
deletion_stability_status
terminal_exit_or_named_debt
```

The immediate use is not a new scalar bound.  It is a better checklist before
using any A000568/fixed-path/edge-witness quotient inside the HYP-3141 ->
HYP-3142 packet theorem.

## Tournament Analysis

The executable scout uses proof carriers as vertices, not runners or raw arcs.

Pairwise observable:

```text
majority over predicate retention, exactness, quotient legality, fiber memory,
deletion stability, and LRC transfer.
```

Switch:

```text
orient A -> B when A wins more axes;
ties prefer lower proof cost, then the listed guardrail path.
```

Fingerprint:

```text
score_hist={0:1, 1:1, 3:2, 4:2, 6:1, 7:1}
directed_3cycles=2
edge_flips_against_listed_path=6
selected path =
  edge_tip_tail_witness_packet
  -> finite_filler_scaffold_shift_package
  -> clustered_canary_S_fiber
  -> monotone_OR_shift_package
  -> fiber_PGF_moment_packet
  -> fixed_path_half_tiling_cube
  -> raw_score_sequence_class
  -> raw_fixed_path_class_count.
```

The directed cycles are useful: the scaffold, canary, and monotone compression
are not globally ordered.  They answer different next-operation questions.

## Guardrail

This is not evidence that Erdős #870 resolves an LRC subproblem.  The external
theorem is used as a proof-architecture source: finite fillers make shift
packages, canary clusters preserve deletion-stability, and representation
multiplicity can be load-bearing.  Every LRC transfer still needs a named
predicate, destroyed coordinate, sidecar, and terminal exit.
