---
source: opus-2026-06-07-S705 (user: work on u(22)∈{60,61}; creative proof/disproof tricks)
status: SHARPENS the reduction + Moser-ring evidence for 60; does NOT resolve u(22). The 61-edge UDG
  (if any) deletes a degree-4-or-5 vertex to a 57- or 56-edge 21-core (Alexeev–Mixon–Parshall:
  u(21)=57, exactly 5 densest-21 graphs). The δ extension vertex needs δ core points concyclic at
  circumradius EXACTLY 1 (a "unit-cocyclic δ-set"). Since u(22)≤61 is PROVEN, extension degree on a
  57-core is ≤4, on a 56-core ≤5. VERIFIED in the Moser ring M_L: the δ=4 route is EMPTY (max
  extension degree 3 on the 57-edge cores; W6⊕Δ family), and a hill-climb caps at U=60 — so within M_L
  u(22)=60. Any 61 lives in the δ=5 route (56-core) or OUTSIDE M_L (like the n=17 densest exception).
  THM-440, HYP-2310. Full trick menu (tested vs open) below.
tags: [unit-distance, erdos, u22, alexeev-mixon-parshall, moser-ring, unit-cocyclic, extension,
  totally-unfaithful, rigidity, cyclotomic-cayley, additive-energy, two-value, frontier, honest]
---

# u(22) ∈ {60,61}: the unit-cocyclic extension, and a menu of two-value tricks

**Prompt (user):** work on the unit distance problem for n=22; since it can only be one of two values,
think creatively about tricks for proof/disproof of one or the other.

**State of the art (Alexeev–Mixon–Parshall 2024 + Engel et al.):** `u(16..21)=41,43,46,50,54,57`
(exact; the **5** densest-21 UDGs enumerated), and `60 ≤ u(22) ≤ 61` — the first open case. The bit:
build an embeddable 61-edge UDG, or kill all 61-edge candidates.

## The reduction made sharp (THM-440)

A 61-edge UDG on 22 vertices has `Σdeg=122<132`, so min degree `δ≤5`; deleting a min-degree vertex
leaves `61−δ≤u(21)=57` edges, so `δ≥4`. Thus `δ∈{4,5}` and the 21-core has **57** (δ=4, a *densest*
21-UDG — one of 5) or **56** (δ=5) edges. The new vertex's `δ` neighbours sit on its **unit circle**:
a **unit-cocyclic δ-set** (δ points concyclic at circumradius *exactly* 1). And because `u(22)≤61` is
*proven*, that extension degree is `≤4` on a 57-core, `≤5` on a 56-core — so:

> **u(22)=61 ⟺ the maximal unit-cocyclic extension degree is attained: 4 on a densest-21 core, or 5
> on a 56-edge 21-core.** A clean, finite, geometric trigger.

## What I tested (the Moser ring) — δ=4 is empty, M_L caps at 60

Exact arithmetic over `ℚ(√3,√11)` in `M_L=ℤ[ζ₆]+w₃ℤ[…]`, `w₃=(5+i√11)/6` (the CM carrier of all
known dense configs): for the 12 `W₆⊕Δ` densest-21 cores, **every** candidate center sees at most **3**
core points (distribution `{1:165,2:48,3:1}`) — **the δ=4 route is empty**, core+1 reaches only `U=60`.
A vertex-swap hill-climb over an `M_L` patch also caps at `U=60`. **Within the Moser ring, u(22)=60.**
(Honest: M_L search can *prove* 61, never *disprove* it; this corroborates 60.)

## The menu of tricks (what's creative, tested, and open)

### To PROVE u(22)=61 (a construction — only positive evidence counts)
- **(C1) δ=5 route [OPEN, the live one].** Find a 56-edge 21-UDG with 5 points unit-cocyclic and a
  faithful center. Not yet tested; the natural next computation (generate 56-edge M_L cores by
  one-vertex perturbation of the 57-cores and run the extension census for degree 5).
- **(C2) Escape the Moser ring [OPEN].** The lone n=17 densest graph already lives *outside* M_L. A
  61-config may need a **new CM direction** `√−d` (a third norm-1 generator beyond `ζ₆, w₃`). Search
  `ℚ(√−3,√−d)` for `d∈{7,15,19,…}` for a 22-point `U=61` set — the "new-quantum" trick (cf. S648
  Moser-fixed-quantum, S699m Heegner tower). Adding a CM direction adds unit vectors → richer
  unit-cocyclic sets.
- **(C3) Glue two cores on a shared unit-cocyclic spine [OPEN].** Instead of +1, overlap two dense
  sub-clusters along 4 concyclic-unit points (a "unit-circle seam"), the additive-energy/Minkowski
  trick that built the 57 = W₆⊕Δ itself.

### To DISPROVE 61 (prove u(22)=60 — must kill *all* candidates)
- **(D1) Totally-unfaithful extension [the experts' weapon, here localized].** Show every degree-4
  unit-cocyclic extension of each of the 5 densest-21 UDGs (and every degree-5 of each 56-core)
  creates a **totally-unfaithful pair**: two non-adjacent vertices forced to unit distance in every
  embedding ⟹ not a faithful 22-UDG. This is *finite* (5 cores × their unit-cocyclic 4-sets) and is
  the most promising disproof route. My M_L census is the geometric shadow: no faithful 4-set exists
  in M_L at all.
- **(D2) Unit-cocyclic non-existence [tested in M_L: holds].** Show no near-densest 21-core has 4
  (resp. 5) points concyclic at circumradius exactly 1 with a faithful center. Verified empty in the
  M_L W₆⊕Δ family (max 3); needs the 5 *actual* enumerated cores (some may be outside M_L) to be a
  proof.
- **(D3) Rigidity/self-stress count [OPEN, creative].** A 61-edge UDG on 22 vertices has self-stress
  dimension `61−(2·22−3)=20` (if generically rigid); a 60-edge one, 19. Totally-unfaithful pairs are
  *self-stresses*. Bounding the realizable self-stress space of a 22-point UDG (via the rigidity
  matroid over `ℚ(√−3,…)`) could cap edges at 60 — the algebraic form of D1.
- **(D4) Hereditary double-deletion [tested: not tight].** Every 20-subset has `≤u(20)=54` edges; two
  adjacent degree-4 vertices force the 20-core to exactly 54 (extremal). Gives constraints on the
  degree sequence but no contradiction alone (slack remains) — useful only combined with D1.

## The bridge (Tournament Analysis directive)

The extension vertex's unit circle is the locus of CM roots of unity around it; `M_L`'s 18 unit
vectors are the `Cay(M_L, U)` generators (S599: unit-distance = additive energy on a cyclotomic
Cayley graph). **A degree-`d` extension ⟺ `d` of a center's 18 unit-translates lie in the core** = a
translation-richness / additive-energy condition. "δ=4 route empty" = no core point has 4 of its 18
translates present — the same `r_min`/additive-energy scarcity that makes the LRC residual hard (S700)
and the same `ζ₆` cyclotomic floor as the worry-set (THM-403). The two-value gap `60↔61` is one unit
of additive energy at the boundary of what the CM ring supports — the unit-distance face of the
Vitali/depth wall (S704).

## Honest status

- **PROVED:** the unit-cocyclic extension reduction + the `≤4/≤5` sharpening (given the cited exact
  values & proven `u(22)≤61`).
- **VERIFIED:** δ=4 route empty in the M_L W₆⊕Δ densest-21 family; M_L hill-climb caps at 60.
- **NOT done:** resolving `u(22)` (open). The M_L evidence corroborates 60 but cannot disprove 61
  (could be outside M_L or in the untested δ=5 route).
- **Most promising next:** D1 (totally-unfaithful extension) on the 5 *actual* densest-21 graphs, and
  C1/C2 (δ=5 route, and a `√−d` ring escape) for a construction.

**Artifacts:** `04-computation/unit_distance_u22_extension_census_s705.py` (+`.out`). Theorem
**THM-440**. New **HYP-2310**. Builds on S614, THM-432/HYP-2301, S599 (Cayley bridge), S648/S699m,
Alexeev–Mixon–Parshall, Engel et al.
