# THM-440 — u(22)∈{60,61} reduces to a faithful unit-cocyclic extension of a near-densest 21-core; the δ=4 route is empty and the Moser ring caps at 60

**Status:** PROVED (the extension reduction, given Alexeev–Mixon–Parshall u(21)=57 & the proven
u(22)≤61) + VERIFIED (the Moser-ring census: δ=4 route empty, M_L caps at 60). Does **NOT** resolve
u(22) (an open problem); it sharpens the reduction and gives Moser-ring evidence for 60.
**Source:** opus-2026-06-07-S705. Builds on Alexeev–Mixon–Parshall (arXiv:2412.11914, u(16..21)=
41,43,46,50,54,57; 60≤u(22)≤61; densest graphs enumerated, #densest-21 = 5), Engel et al.
(arXiv:2406.15317, Moser ring, u(22)≥60), S614 (the deletion reduction), THM-432/HYP-2301 (M_L).

## The unit-cocyclic extension reduction (PROVED)

> Let `G` be a hypothetical 61-edge unit-distance graph on 22 vertices. `Σdeg = 122 < 132 = 6·22`,
> so `δ(G) ≤ 5`; deleting a min-degree vertex `v` leaves a UDG on 21 vertices with `61−δ ≤ u(21)=57`
> edges, so `δ ≥ 4`. Hence `δ∈{4,5}` and the 21-core `C=G−v` has **57** (δ=4) or **56** (δ=5) edges.
> The `δ` neighbours of `v` lie on a common **unit circle centred at `v` (circumradius exactly 1)** —
> a faithful **unit-cocyclic `δ`-set** (no other core vertex at distance 1 from `v`).
>
> **Therefore `u(22)=61` ⟺**
> - **(δ=4)** one of the **5 densest 21-vertex UDGs** admits a faithful 4-point unit-cocyclic
>   extension, **OR**
> - **(δ=5)** some **56-edge** 21-vertex UDG admits a faithful 5-point unit-cocyclic extension.
>
> Both are *finite* checks given the (enumerated) embeddings.

**Sharpening via the proven `u(22)≤61`.** A center `v` added to a 57-edge core can have at most **4**
core-neighbours (5 would give `U≥62 > u(22)`); added to a 56-edge core, at most **5**. So each route
is the single question "is the maximal extension degree attained?" — degree 4 on a 57-core, or 5 on a
56-core.

## The unit-cocyclic obstruction

> A degree-`d` extension vertex requires `d` core points **concyclic with circumradius exactly 1**.
> Since `U_count` (exact unit-distance count) equals the faithful edge count, **any 22-point subset
> of a carrier lattice with `U=61` proves `u(22)=61`**; conversely the routes above must be empty for
> `u(22)=60`.

## Moser-ring census (VERIFIED, this session)

In `M_L = ℤ[ζ₆]` extended by `w₃=(5+i√11)/6` (the biquadratic CM ring `ℚ(√−3,√−11)` carrying all
known dense small configs), exact arithmetic over `ℚ(√3,√11)` (no float decides adjacency):
- the `W₆⊕Δ` family yields **12** distinct 57-edge 21-vertex cores; for **every** one, the candidate
  extension centers have unit-neighbour-count distribution `{1:165, 2:48, 3:1}` — **max extension
  degree = 3**. So the **δ=4 route is EMPTY in M_L**: a densest-21 Moser core + 1 vertex reaches only
  `U=60`.
- a vertex-swap hill-climb over an `M_L` patch caps at **`U=60`** (no 61 found).

> **So within the Moser ring, `u(22)=60`** (61 unreachable by +1 on a densest-21 core, and 60 is the
> hill-climb ceiling). Any realisation of `61` must therefore lie in the **δ=5 route** (a 56-edge
> core) **or outside `M_L`** — exactly the situation of the lone n=17 densest graph whose embedding
> is *not* in the Moser ring (Alexeev–Mixon–Parshall, Table 1 footnote).

## Scope / honesty

- The reduction is rigorous **given** the cited exact values; the `≤4/≤5` sharpening uses the proven
  `u(22)≤61`.
- The Moser census is exact but **carrier-restricted**: an `M_L` search can **prove** `u(22)=61` (by
  exhibiting a 61-set) but **cannot disprove** it (61 could live outside `M_L`). It found only 60, so
  it **corroborates `u(22)=60`** without proving it. The `W₆⊕Δ` family need not be all 5 enumerated
  densest-21 graphs.
- This resolves no open case. The genuine disproof of 61 is graph-side: kill all 61-edge candidates
  via the **totally-unfaithful** obstruction (Globus–Parshall / Alexeev–Mixon–Parshall) applied to
  the unit-cocyclic extension.
- S129 adds a complementary negative lemma for the spine route: the five exact `n=21`, `57`-edge
  graph6 cores are endpoint-universal. Every vertex deletion is traceable, every incident edge can
  be the reattachment edge, and every graph-only degree-4 or degree-5 one-vertex extension preserves
  a Hamiltonian unit spine. Therefore the missing obstruction in the `delta=4` route cannot be plain
  traceability; it must be the faithful unit-cocyclic geometry or a totally-unfaithful subgraph.

**Tournament/LRC bridge (directive):** the extension vertex's unit circle is the locus of 6th/CM
roots of unity around it; the 18 `M_L` unit vectors are the `Cay(M_L,U)` generators. A degree-`d`
extension ⟺ `d` of the center's unit-translates lie in the core = an **additive-energy/translation
richness** condition on the core (the S599 unit-distance ↔ cyclotomic-Cayley bridge). "δ=4 route
empty" = the densest-21 Moser core has no point with 4 of its 18 unit-translates present.

**Artifacts:** `04-computation/unit_distance_u22_extension_census_s705.py` (+`.out`), reusing the exact
M_L arithmetic of `unit_distance_moser_lattice_u21_monad_s4.py`. Reflection
`07-reflections/u22-the-unit-cocyclic-extension-and-the-two-value-tricks-s705.md`. New: **HYP-2310**.
Builds on S614, THM-432/HYP-2301, Alexeev–Mixon–Parshall, Engel et al.
