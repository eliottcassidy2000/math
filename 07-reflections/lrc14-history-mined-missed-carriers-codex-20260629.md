---
date: 2026-06-29
tags: [lrc14, random031, history-mining, seam-complement, fiber-pgf, sidecars, tournament-analysis]
hypothesis: HYP-3494
tangent: T1454
technique: LTI-454
tournament_technique: LTT-354
---

# LRC14 History-Mined Missed Carriers

This pass started as a search for overlooked niche topics in the repo.  The
surprise was that the search did not merely recover more analogies.  It turned
HYP-3486's fiber graph into a smaller coefficient object:

```text
P_escape(y) = 24 + 226 y + 8 y^2.
```

That polynomial is the random031-local version of HYP-3140's sheet-count PGF.
It says the seam-complement witness set is mostly one rank-2 exit per occupied
fiber, with eight two-exit fibers and twenty-four zero-exit fibers.  The
zero-exit fibers are not failed routes; their signatures are free-hole packets.

The clearest reframe is seam cobordism.  HYP-3490 is the blocked-current side:
every E/branch-touched label on random031 is private, so dead-projection
current deletion cannot move.  HYP-3486 is the phase-flow side: after deleting
the hard seam, phase flow routes through rank-2 exits, free-hole packets, and a
single pure bypass component.  The missing proof object is the legal bridge
between these sides, not another scalar wall.

The incoming HYP-3490 Lean firewall ledger sharpens that bridge.  The module
`TournamentH7.LRCPrivateLabelFirewall` can serve as the formal left boundary of
the cobordism; HYP-3494 is then the right-boundary problem of making the
random031 phase-flow fiber graph, PGF buckets, and sidecars meet it without an
illegal quotient.

The concurrently landed HYP-3493 seam-sheaf scaffold is exactly the kind of
right-boundary partial execution HYP-3494 was asking for.  Its no-mixed-stalk
readout means the next obstruction is sharper than "build topology": prove
owner-boundary persistence for the pure bypass, or emit named owner/two-adic
debt.

The newer HYP-3510/HYP-3511 seam packet makes that even sharper.  Coarse phase
transport is connected after deleting the hard seam, and the `40` free-hole
cells are bracketed by ordinary rank-2 boundaries except for two controlled
same-branch doublets.  So HYP-3494's quotient-price table should treat
free-hole packets as a local lemma and reserve true debt for the pure bypass
owner-boundary.

The older carriers now fall into roles:

```text
HYP-3422/HYP-3428: do not identify n*2 with legal topology; carry branch loss.
HYP-3034/HYP-3243: build relative H1/topes for seam complement and boundary.
HYP-3402/HYP-3451/HYP-3437: price owner-current and Menger/Green cut debt.
HYP-3023/HYP-3300: make illegal quotienting visible as a missing sidecar.
HYP-3140/HYP-3511: convert raw witness/free-hole counts into coefficient and
  bracket lemmas.
HYP-3480: give the Lean-facing dispatch split, six singleton rows plus random031.
```

The most actionable next experiment is not a huge new search.  It is a small
observability table for random031 itself.  Columns should be:

```text
u_index
branch
mirror mate
cell class
endpoint rank
owner word
private-label firewall bit
vertical half-turn sidecar
relative-H1 class
PGF coefficient bucket
```

If those columns separate rank-2 exits, HYP-3511 free-hole brackets, and the
pure bypass owner-boundary with no illegal quotient, then HYP-3494 becomes a
formal interface theorem: a random031 terminal proof can forget raw runner
order only after this sidecar tuple is attached.

Assumption challenge: I deliberately did not let runners or gates be the only
vertices.  I considered witness cells, u-fibers, branch sheets, owner labels,
punctures, relative cycles, component cuts, PGF coefficients, observability
columns, Lean dispatch packets, and proof obligations.  The chosen carrier is
the proof-interface bundle because that is what preserves the terminal
predicate after the seam is deleted.
