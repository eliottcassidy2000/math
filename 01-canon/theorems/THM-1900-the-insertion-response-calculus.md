---
id: THM-1900
title: "THE INSERTION-RESPONSE CALCULUS -- generalizing THM-1865's source/sink diagnostic to an arbitrary insertion, and unifying it with boxeph THM-1855 (order-join) + kind-pasteur THM-1860 (H=prod H(SCC)). OPERATION: add a vertex u to T that beats EXACTLY the subset P of V (u->j for j in P, j->u otherwise); source = P=V, sink = P=empty. TWO EXACT LAWS (verified exhaustively n=3,4,5). (A) THE c3-VELOCITY IS THE FORWARD CUT: Delta c3(T,P) = e(P -> V\\P) = the number of arcs from P to its complement. PROOF: a 3-cycle through u is exactly u->a->b->u with a in P, b notin P, a->b. So c3 is insertion-NEUTRAL iff P is a CLOSED set (no forward cross-arcs), i.e. a union of TERMINAL strong components = a down-set of the condensation. (B) H-NEUTRALITY = CONDENSATION DOWN-SETS: adding u beating P leaves H unchanged iff P is a down-set of the condensation poset (the linear order of strong components), and the number of H-neutral patterns equals #SCC(T)+1. VERIFIED n=3,4,5: the count of tournaments with exactly 2 neutral patterns (only source+sink, i.e. strongly connected) is 2, 24, 544 = A051337 (strongly-connected labelled tournaments) exactly; the fully-transitive tournaments (n! of them, #SCC=n) have the maximal n+1 neutral patterns. Both H and c3 are neutral on the SAME family (the condensation down-sets). CONSEQUENCE (the inflation-response frame, THM-1865): H and c3 are inflation-NEUTRAL except on the rigid skeleton of down-sets (which collapse to just source/sink when T is strongly connected), so their extremals are inflation-hard; and this EXPLAINS the forgotten HYP-260 (the delta-inequality fails only at source/sink -- those are exactly the two trivial down-sets that survive on a strongly-connected core)"
status: PROVED. (A) c3-velocity = forward cut is elementary. (B) the all-order neutrality proof was repaired on 2026-10-20 using the exact insertion-position count, nonnegative orphan contributions, and Camion's Hamilton-cycle theorem on a split strong component. The old deletion-of-HP premise was false. FINITE-EXACT controls now include all 4054 attachments of classes through order six; see zero_heavy_ports_20261020g/profile_growth.md.
author: opus-2026-07-20-S439
generalizes: THM-1865 (source/sink = the two trivial down-sets P=V, P=empty)
unifies: THM-1855 (boxeph: order-join T1|>T2), THM-1860 (kind-pasteur: H=prod H(SCC), c3=sum c3(SCC))
explains: HYP-260 (delta-inequality fails only at source/sink)
depends_on: [THM-1865 (inflation-response diagnostic), THM-1855 (order-join algebra), THM-1860 (SCC decomposition), A051337 (strongly-connected labelled tournaments 2,24,544)]
---

# THM-1900 — The insertion-response calculus

**Correction, 2026-10-20:** the original mechanism for B assumed that deleting
the inserted vertex from a Hamiltonian path always leaves a Hamiltonian
path. A directed triangle already refutes that premise. The theorem is
unchanged; the proof below counts insertions into old paths without
discarding the additional orphan paths. See
[the full repair and endpoint-cover response](../../05-knowledge/results/zero_heavy_ports_20261020g/profile_growth.md).

The generative-engine seed for [THE ZOO](../../00-navigation/THE-ZOO.md) (S439). Push the
inflation-response frame (THM-1865, §3.19 of the zoo) on the **insertion operation**: add a new
vertex `u` to `T` on `V = {0,…,n−1}` that **beats exactly the subset `P`** (`u→j` for `j∈P`,
`j→u` otherwise). Source = `P=V`, sink = `P=∅` (THM-1865's two cases). Two exact laws.

## A. The c₃-velocity is the forward cut

> **`Δc₃(T,P) = e(P → V∖P)`** — the number of arcs from `P` to its complement.

**Proof.** A directed 3-cycle through the new vertex `u` is `u→a→b→u`, which forces `a∈P`
(so `u→a`), `b∉P` (so `b→u`), and `a→b`. Conversely every such `(a,b)` gives one 3-cycle. No
3-cycle avoiding `u` is created. Hence `Δc₃ = #\{(a,b): a∈P,\ b∉P,\ a→b\} = e(P→V∖P)`. ∎

So **c₃ is insertion-neutral ⟺ `e(P→V∖P)=0` ⟺ `P` is a *closed set*** (all cross-arcs point
`V∖P → P`), i.e. `P` is a union of **terminal** strong components — a **down-set of the
condensation**. Verified exactly for all `(T,P)`, `n=3,4,5`.

## B. H-neutrality = condensation down-sets

The condensation of a tournament is a **linear order** of strong components
`C_1 ▷ C_2 ▷ ⋯ ▷ C_k` (`C_1` beats all later). Its **down-sets** are the suffixes
`C_j ∪ ⋯ ∪ C_k` (the terminal components) — there are `k+1` of them (`∅` up to `V`).

> **Adding `u` beating `P` is H-neutral ⟺ `P` is a down-set of the condensation; and the number
> of H-neutral patterns is `#SCC(T) + 1`.**

**Repaired proof.** On an old Hamiltonian path put bit 1 at a vertex in P
and bit 0 elsewhere. The allowed insertion positions number
`[first=1] + #(0->1 transitions) + [last=0] = 1 + #(1->0 transitions)`.
Each old path therefore has at least one distinct extension. Deleting u
from any of these recovers its particular old path, so there is no overlap.
Other paths, whose deletion leaves an invalid joining arc, contribute
nonnegative additional counts: they are the orphans of
[THM-012b](THM-012b-insertion-decomposition.md).

When P is a suffix of the SCC order, the new tournament is the same ordered
join with u placed at that cut. Every old path has exactly one extension
and every new path is one of these, so H is unchanged. If P splits an SCC,
take a Hamilton cycle of that strong component (Camion's theorem). It has
a P-to-complement transition; deleting a different cycle edge retains that
transition in a spanning path. Concatenate arbitrary Hamiltonian paths of
the other components. If P is a union of SCCs but is not a suffix, arbitrary
concatenated component paths already have a P-to-complement transition.
In either case an old path has at least two insertions, proving strict
increase. Camion's theorem is the same classical dependency used in
[THM-4526](THM-4526-shaved-tournaments-path-plus-chord-and-the-largest-unavoidable-graphs.md).
Thus both directions hold at every order; for a strong parent only the
empty and full beaten sets are neutral.

**Verification (n=3,4,5, exhaustive).** The number of tournaments with exactly **2** neutral
patterns (only source+sink) is `2, 24, 544` = **A051337** (strongly-connected labelled
tournaments) — exactly. The **`n!`** transitive tournaments (`#SCC=n`) attain the maximal `n+1`.
The full histograms (`#neutral : #T`): `n=5 → {6:120, 4:120, 3:240, 2:544}`, matching
`#neutral = #SCC+1` with the SCC-count distribution. The families `{H-neutral}`, `{c₃-neutral}`,
and `{condensation down-sets}` **coincide** at every tested `(T)`.

## C. Why this unifies four threads

| thread | relation to the calculus |
|---|---|
| **THM-1865** (source/sink inflation-response) | the two *trivial* down-sets `P=∅, P=V`; the general down-set family is the refinement |
| **boxeph THM-1855** (order-join `T₁▷T₂`) | a join is an insertion at a component boundary; the velocity table is the boundary case of `Δc₃=e(P→V∖P)` |
| **kind-pasteur THM-1860** (`H=∏H(SCC)`, `c₃=Σc₃(SCC)`) | down-sets are exactly the cuts of the SCC product; neutrality = "the insertion lands between factors" |
| **forgotten HYP-260** ("δ-inequality fails only at source/sink") | source/sink are the *only* down-sets on a strongly-connected core — so the boundary failures are exactly the trivial down-sets |

**The frame it feeds (zoo §6, generator 4).** For any operation and invariant, classify
NEUTRAL / PUMPED. Here: `H, c₃` are neutral on the down-set skeleton (rigid) while `Δc₃` off it is
the exact forward-cut velocity. Next targets named in the zoo: is **min-feedback-arc-set**
insertion-neutral under `+source`? (a source adds `0` back-arcs ⟹ predicted neutral, which would
put min-FAS in the rigid-extremal class and fill one of the two empty invariant columns).

## Verification

`04-computation/insertion_response_calculus_opus_S439.py` (+ `.out`) — `Δc₃=e(P→V∖P)` for all
`(T,P)`, `n=3,4,5`; the H-neutral histograms matching `#SCC+1` and A051337; `c₃`-neutral = closed
sets.
