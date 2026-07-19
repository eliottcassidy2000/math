---
id: THM-1262
title: A COHERENT BLOCKER TWO-CYCLE CANNOT BE A MARKED INVERSION — ascent protection forces an aligned corridor through a third owner
status: PROVED (two-cycle speed orientation; protected ascent tooth disjoint from the reverse target tooth; nonconsecutivity; binary mismatch-to-adjacency contradiction; aligned phase corridor; forced third-owner protected handoff). This eliminates the inverted two-cycle cell but does not exclude aligned two-cycles, six-comb coverage, or LRC(14)
source: codex-2026-07-19-S78 coherent-word continuation
depends_on: [THM-1240, THM-1248, THM-1252, THM-1253, THM-1256]
related: [THM-1156, THM-1238, THM-1254, THM-1260, HYP-7870]
script: 04-computation/lrc14_blocker_two_cycle_alignment_thm1262.py
output: 05-knowledge/results/lrc14_blocker_two_cycle_alignment_thm1262.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCBlockerTwoCycleAlignment.lean
script_sha256: a0e0163a193e577e0f71724713b2e5b91f53d78b57e1bda1ee830416e17d3dc6
output_sha256: 6729240fdb64cc5e09a7dfc2e7dcaac82a3633fd1b8559d08f75317d64b9dae5
formalization_sha256: af0ede0412ac7dbe7cf18c2edcbde4d6cb5e0a5fe21288e638da50536ba19daa
---

# THM-1262 — a coherent blocker two-cycle is aligned and has a third-owner bridge

## 1. Setup

Let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)]
```

be a complete `c`-safe gap covered by the six strict danger combs of

```text
c<d_1<...<d_6.
```

Choose THM-1253's deletion-minimal chronological tooth cover and choose all
THM-1240 centered-spoke blockers from this one word, as in THM-1254.  Suppose
two labels `l,h`, with

```text
d_l<d_h,                                                (1)
```

form a blocker two-cycle

```text
l -> h -> l.                                           (2)
```

Write

```text
H = the selected d_h-tooth containing t_l,
L = the selected d_l-tooth containing t_h.             (3)
```

Thus `H` is the marked target tooth of the speed ascent `l->h`, while `L` is
the marked target tooth of the speed descent `h->l`.  They are distinct
selected teeth because they have distinct owners.

## 2. Ascent protection separates the two marked teeth

Let `S_l` be the closed `d_l`-safe component through `t_l`.  The ascent half
of (2) satisfies the protected-containment clause of THM-1252:

```text
closure(H) subset int(G intersect S_l).                (4)
```

But `L` is a strict `d_l`-danger tooth.  Therefore

```text
H intersect L = empty.                                 (5)
```

Consecutive teeth in the irredundant chronological word have a positive raw
overlap.  Equation (5) consequently implies that `H` and `L` are not
consecutive.  If their word positions are `A_H,A_L`, then

```text
|A_H-A_L|>=2.                                          (6)
```

This is the point that is unavailable for a general cycle edge: the reverse
edge in a two-cycle identifies the other marked tooth with a tooth of the
very source comb that (4) keeps safe.

## 3. The binary descent cannot be inverted

For the descent `h->l`, the target satisfies

```text
d_l<d_h<c+d_h,
```

so THM-1248 gives a binary relative digit.  THM-1256 then gives the exact
landing dichotomy for the pair `(L,H)`:

1. their tooth order agrees with the order of the centered phases
   `(t_h,t_l)`; or
2. the orders disagree, in which case `L` and `H` are consecutive selected
   teeth and overlap in a directly charged seam.

The second alternative contradicts (5)--(6).  Hence only the first is
possible:

```text
sign(A_L-A_H)=sign(t_h-t_l).                           (7)
```

The phase separation is quantitative:

```text
5/[28(c+d_l)]<|t_h-t_l|<23/[28(c+d_l)].               (8)
```

Thus a coherent blocker two-cycle can never occupy THM-1256's marked
adjacent-inversion branch.  Reflection reverses both orders and preserves
(7), so this is not a choice of lower-gap gauge.

## 4. Alignment forces a third-owner bridge

By (7), the chronological subword from `H` to `L` is oriented in the same
direction as the actual phase segment from `t_l` to `t_h`.  Consecutive
selected teeth overlap, so this subword literally covers that closed phase
segment.  By (6) it contains at least one intermediate tooth.

Let `J` be the immediate neighbour of `H` on the side of `L`.  Its owner is
not `h`, because distinct teeth of one danger comb are disjoint.  Its owner
is not `l`, because `J` overlaps `H` while (4) places all of `H` inside the
`d_l`-safe component.  Therefore

```text
owner(J) notin {l,h}.                                  (9)
```

The overlap `J intersect H` is the corridor-facing wall handoff of `H`.
It is one of THM-1253's pairwise-disjoint raw seams and satisfies

```text
|J intersect H| >= 1/[14 lcm(owner(J),d_h)].          (10)
```

At its wall the carrier `c` and source `d_l` are strictly safe, while `d_h`
is at equality.  Hence (10) is precisely THM-1252's protected third-support
handoff, now forced to point into the aligned two-cycle corridor.

> **Two-cycle alignment theorem.**  Every coherent blocker two-cycle has an
> aligned centered-phase corridor of at least three selected teeth.  The
> corridor leaves the ascent target through a literal protected seam owned
> by a third fast label.  The marked-inversion alternative is impossible.

The seam in (10) already belongs to the full chronological invoice; it is
not additional measure.  The gain is operational: a two-cycle cannot close
inside a two-label adjacent cell and must export its continuation through a
third owner.

## 5. Tournament and carrier audit

The speed tournament on `{l,h}` is transitive and contains only `l->h`; it
forgets the blocker two-cycle.  The blocker relation retains (2) but forgets
why its two target teeth cannot touch.  Completing the remaining four labels
by speed order produces gauge-dependent tournament edges and no proof.

The proof-bearing vertices are instead the two marked target teeth and the
corridor-facing wall event:

```text
(H,t_l,source-safe label l),
(L,t_h,source-safe label h),
(J intersect H, owner(J), chronological side).        (11)
```

The pairwise observable is positive tooth overlap; the switch is circle
reflection, which reverses both phase and tooth order; the chronological word
is the tie Hamiltonian path.  Projecting (11) to runners destroys (4), while
projecting to wall seams destroys the identification of `L` as the reverse
target.  The challenged assumption is that the shortest blocker cycle is the
easiest terminal cell: after placement on the common word, it is forced to
open into a third-owner corridor.

## 6. Verification and scope

The dependency-free referee checks the finite order logic, reflection
covariance, every two-label speed/phase/order orientation, and the
nonconsecutive-to-intermediate-word implication.  The Lean module checks the
abstract disjointness/nonconsecutivity and binary mismatch contradiction,
the alignment sign law, and the third-owner exclusion from the protected
overlap.  Irredundant-cover extraction, strict danger/safe disjointness, and
the geometric identification of the corridor-facing wall remain the named
paper topology providers.

Frozen artifact hashes are

```text
source         a0e0163a193e577e0f71724713b2e5b91f53d78b57e1bda1ee830416e17d3dc6
output         6729240fdb64cc5e09a7dfc2e7dcaac82a3633fd1b8559d08f75317d64b9dae5
formalization  af0ede0412ac7dbe7cf18c2edcbde4d6cb5e0a5fe21288e638da50536ba19daa
```

THM-1262 does not rule out an aligned two-cycle.  It removes one complete
local branch and sharpens its survivor to a placed three-owner transport
problem.  Together with THM-1260's single-fork `chi_7` universality no-go, it
also shows what a useful Fano consumer would have to see: incidence among
several such exported bridges, not the sign of one isolated fork.  ∎
