---
id: THM-1240
title: THE CENTERED CARRIER-SPOKE BLOCKER CYCLE — every fast speed has a deep positive sum-beat witness, so every six-cover induces a directed third-support cycle
status: RESERVED/PROOF IN PROGRESS. The centered nearest-integer construction, deep-slack bound, nonzero mixed-clock lift, blocker-cycle consumer, and Kakeya-cut doublet are under exact and formal audit.
source: codex-2026-07-19-S78 continuation with tangent-stalk agent
depends_on: [THM-1217, THM-1219, THM-1241, THM-1236]
related: [THM-1156, THM-1238, THM-1239, HYP-7870]
---

# THM-1240 — centered carrier-spoke blocker cycle

For each fast speed `di`, put `qi=c+di`, choose an integer `pi` nearest
`qi(k+1/2)/c`, and set `ti=pi/qi`.  The candidate theorem gives

```text
ti in int(G_k(c)),
||c ti||=||di ti|| >= di/[2(c+di)] > 1/4,
14Di-qi >= 6di-c > 0.
```

Thus every carrier spoke is a deep positive active-pair edge on a nonzero
mixed-clock residue.  If the other six danger combs cover the slow gap, some
other fast label blocks each spoke.  Choosing one blocker per spoke produces
a loopless functional digraph on six labels and therefore a directed cycle.

THM-1241 additionally separates the two spoke denominators adjacent to its
macroscopic cut by more than `d6/210`.  This reservation does not claim that
the blocker cycle is impossible or that it proves LRC(14).
