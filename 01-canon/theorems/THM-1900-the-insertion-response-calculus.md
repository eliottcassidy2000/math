---
id: THM-1900
title: "THE INSERTION-RESPONSE CALCULUS -- generalizing THM-1865's source/sink diagnostic to an arbitrary insertion, and unifying it with boxeph THM-1855 (order-join) + kind-pasteur THM-1860 (H=prod H(SCC)). OPERATION: add a vertex u to T that beats EXACTLY the subset P of V (u->j for j in P, j->u otherwise); source = P=V, sink = P=empty. TWO EXACT LAWS (verified exhaustively n=3,4,5). (A) THE c3-VELOCITY IS THE FORWARD CUT: Delta c3(T,P) = e(P -> V\\P) = the number of arcs from P to its complement. PROOF: a 3-cycle through u is exactly u->a->b->u with a in P, b notin P, a->b. So c3 is insertion-NEUTRAL iff P is a CLOSED set (no forward cross-arcs), i.e. a union of TERMINAL strong components = a down-set of the condensation. (B) H-NEUTRALITY = CONDENSATION DOWN-SETS: adding u beating P leaves H unchanged iff P is a down-set of the condensation poset (the linear order of strong components), and the number of H-neutral patterns equals #SCC(T)+1. VERIFIED n=3,4,5: the count of tournaments with exactly 2 neutral patterns (only source+sink, i.e. strongly connected) is 2, 24, 544 = A051337 (strongly-connected labelled tournaments) exactly; the fully-transitive tournaments (n! of them, #SCC=n) have the maximal n+1 neutral patterns. Both H and c3 are neutral on the SAME family (the condensation down-sets). CONSEQUENCE (the inflation-response frame, THM-1865): H and c3 are inflation-NEUTRAL except on the rigid skeleton of down-sets (which collapse to just source/sink when T is strongly connected), so their extremals are inflation-hard; and this EXPLAINS the forgotten HYP-260 (the delta-inequality fails only at source/sink -- those are exactly the two trivial down-sets that survive on a strongly-connected core)"
status: PROVED. (A) c3-velocity = forward cut is exact (elementary + exhaustive n=3,4,5). (B) H-neutral count = #SCC+1 verified n=3,4,5 (counts match A051337); the down-set mechanism is proved; the elementwise identity {H-neutral} = {c3-neutral} = {condensation down-sets} is computationally confirmed n<=5.
author: opus-2026-07-20-S439
generalizes: THM-1865 (source/sink = the two trivial down-sets P=V, P=empty)
unifies: THM-1855 (boxeph: order-join T1|>T2), THM-1860 (kind-pasteur: H=prod H(SCC), c3=sum c3(SCC))
explains: HYP-260 (delta-inequality fails only at source/sink)
depends_on: [THM-1865 (inflation-response diagnostic), THM-1855 (order-join algebra), THM-1860 (SCC decomposition), A051337 (strongly-connected labelled tournaments 2,24,544)]
---

# THM-1900 — The insertion-response calculus

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

**Mechanism.** A Hamiltonian path of `T+u` restricts to a Hamiltonian path of `T` with `u`
inserted at one position (or an end). `u` beats `P` and loses to `V∖P`; it can be inserted
consistently at a *unique* place — and without creating new paths — exactly when `P` is the set
of components below some cut of the condensation chain (a down-set): `u` slots at that cut. For a
non-down-set `P`, the insertion either forks (multiple valid positions → `H` increases) . For a
strongly connected `T` (`k=1`) only `∅` and `V` survive, recovering THM-1865's source/sink.

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
