# Smith diagrams square the staircase: the Tutte nexus, the harmonic aspect ratio, and a comprehensive map of the tournament/metagraph threads

**Session:** kind-pasteur-2026-07-15-S128 (cont.7)
**Prompt (owner):** consider Smith diagrams and squaring the square; search back through ALL previous
tournament and metagraph work; go down niche underexplored threads for connections.
**Computations:** 04-computation/smith_diagrams_staircase_kps_S128c7.py (+.out) — all exact.
**Companion:** HYP-6865.

---

## 1. Why Smith diagrams belong here

Brooks–Smith–Stone–Tutte (1940) turned a squared rectangle into an electrical network — the **Smith
diagram**: nodes = maximal horizontal segments, edges = squares, Kirchhoff currents = square sizes,
the whole thing solved by spanning-tree counts. Two facts make this OUR machinery and not an analogy:

1. **The correspondence is a cut⊕cycle statement.** Potentials live on the cut space, currents on
   the cycle space; solving the network is exactly reconciling the two — and the cut⊕cycle
   decomposition of K_n (base path = cut, tiles = cycle) is this project's master frame (the owner's
   (n−1)² = 13 + 2·T(12); opus-S270's perspective family; the wiggly/waggly layers).
2. **Tutte is already in the repo, unexecuted.** kps-S11's even-graph equinumerosity block
   (INVESTIGATION-BACKLOG ~line 798) recorded "G_n ↔ E_n = TUTTE chromatic↔flow duality (even
   subgraphs = flow support)" and never ran it. Squaring the square and the flow/chromatic duality
   are the same man's two hats; this session connects them through the staircase.

## 2. The staircase HAS a Smith diagram, and it is exactly solvable (new, verified n=4..10)

The staircase Δ_{n−2} dissected into its unit cells is a dissection like any squared rectangle, so it
carries a canonical Smith network **N_n**: nodes = maximal horizontal segments, edges = the m =
C(n−1,2) cells — i.e. **a canonical planar graph structure on the tile set** (each tile = one arc slot
of the tournament). Computed exactly (n = 4..10) and then read off structurally:

- The staircase's rows are left-aligned, so each integer height is ONE segment: **N_n is the series
  chain of parallel bundles** of multiplicities (n−2, n−3, …, 1) — row y is a bundle of n−1−y
  parallel unit resistors.
- **κ(N_n) = (n−2)!** (spanning trees = row transversals: pick one tile per row).
- **Pole-to-pole resistance R_n = H_{n−2} = 1 + 1/2 + … + 1/(n−2), the HARMONIC NUMBER**
  (3/2, 11/6, 25/12, 137/60, 49/20, 363/140, 761/280 — exact Fractions from the Kirchhoff solve).
- **Self-duality is exact:** the vertical-segment (BSST-dual) network has identical κ, R, and degree
  data at every n — the grid anti-diagonal reflection (= the project's grid-symmetry involution) IS
  Smith-diagram duality. Blue tilings are orientations of a SELF-DUAL network.
- The BSST current through row y splits equally over its bundle: the canonical **electrical tile
  weight w(x,y) = 1/(n−1−y)**. The "squared rectangle" of the staircase is the classical harmonic
  dissection: a 1 × H_{n−2} rectangle in rows of k equal squares, k = n−2, …, 1.

**The γ resonance.** R_n = H_{n−2} ~ ln n + γ: the Euler–Mascheroni constant — already the fourth
constant of the triangle foundation (Gamma-function asymptotics in Burnside/Stirling) — reappears as
the staircase's electrical growth rate. The triangle keeps only four constants and keeps reusing them.

## 3. The flow laws are electrical laws (verified n=5,6)

The d=m complement-line layer, which the fleet just mapped (codex THM-785 C3-flux, opus THM-787/790
leg law, my THM-791 H-stratum), reads verbatim as a circuit on the merged metagraph:

- **The leg law is an EMF statement.** ΔE4 = 8(e₁ − e_n) per line: only the two staircase LEGS carry
  axis current ("the interior is silent") — the legs are the battery terminals, opus's apex GF factor
  (z + 1/z) is the unit current generator, in their own words.
- **Blue lines satisfy e₁ + e_n = n − 2 identically** (verified on every grid-sym tiling, n=5,6):
  grid symmetry conjugates vertex-1 wins into vertex-n losses. One line then gives
  Δx = 16·e₁ − 8(n−2) ≡ 8n (mod 16) — the THM-785/787 blue parity dichotomy re-derived instantly.
- **The battery table.** Class flux Φ(K) = Σ_fiber (e₁ − e_n): global ΣΦ = 0 (conservation; the layer
  is a closed circuit); the strongest batteries are MIXED classes (±32 at n=6, vs ±4 for pure-blue);
  PBk fluxes are small and — unexplained — come in EQUAL PAIRS across distinct classes at n=6
  (−28,−28; 37,37; …). The interface (mixed) classes don't just connect the phases; they POWER the flow.
- Boundary marker: a KCL absorption inequality was tried and WITHDRAWN on the LRC side
  (THM-767/MISTAKE-146) — the electrical reading here is bookkeeping-exact, not an inequality import.

## 4. The comprehensive thread map (what exists, what is dormant, where Smith connects)

**Core frames (active):** triangle foundation (staircase = tournaments; Modes A/B); cut⊕cycle
perspective family (opus-S270; THM-790 leg law); the transitivity flow atlas (THM-785 codex C3-flux +
interface law + flow address; THM-787 opus E4 axis; THM-791 H-stratum mod-4 + majorization; THM-793
opus Mode-B line tower "axis current fiber-measurable"). Smith reading: flow = current, proved.

**Dormant/niche (surveyed this session, ranked for connection potential):**
1. **kps-S11 even-graph equinumerosity** (07-01, unexecuted "Next:" items): G_n↔E_n Tutte duality;
   E_n = Curie–Weiss Ising (verified identity!); Seidel/two-graphs/equiangular lines with the
   28 = χ(E_7) = max-equiangular-lines(ℝ⁷) pinch at the special primes {3,7,23}. ⭐ the natural
   sequel: N_n gives the PLANAR substrate on tiles; the Tutte polynomial of N_n specializes to both
   sides of the duality — compute T(N_n; x,y) (series-parallel ⟹ closed form) and test the kps-S11
   "S_n-symmetrized Tutte specialization" guess against V(E_n), V(G_n), bridge rank B.
2. **Hadwiger on G_n/Z₂** (klein-S198): minor theory on the metagraph — now has an electrical cousin
   (effective-resistance embeddings bound minors).
3. **Wiggly-class anisotropy** (kps-S20er: self-loop rates differ by tile position): the electrical
   tile weight w = 1/(n−1−y) is a candidate explanatory covariate — silent-mutation rate vs w.
4. **Circulant homology/THM-125, path homology β₂** — Laplacian eigenspaces; the mod-2 even-graph
   projection T_cycle = (I + L(K_n))T already uses the K_n Laplacian; N_n is its planar shadow.
5. **LRC-side electrical motifs** (codex T1336 Green conductance, Thomson profiles on k=8 rows) —
   the OTHER pillar already speaks circuit language; the staircase now does too, from the geometry.
6. **Square complexes / locally testable codes lead** (kps-S22 external): the tiling as a square
   complex; Smith networks are its 1-skeleton economics.

## 5. Honest boundaries

- N_n's structure theorem is elementary (series-parallel); its value is CANONICITY (the dissection
  functor applied to the project's central object) and the exact constants ((n−2)!, H_{n−2}, γ), not
  difficulty. K_n − P (the tile graph on VERTICES) stops being planar at n=8 — vertex-side Smith
  diagrams are small-n only; the CELL-side network N_n is planar for all n. Perfect-squared-square
  PERFECTNESS (all sizes distinct) fails for the staircase (equal sizes within a row) — the harmonic
  dissection is canonical, not perfect.
- The equal-pair PBk flux phenomenon and the Tutte polynomial computation are logged, not done.

## 6. What transcends

Every time this project touches a classical duality it lands on the same three-way coincidence: the
staircase's geometric reflection = the tournament's converse = a functorial duality (Smith/Tutte/
planar). The blue sector is always the self-dual locus, and the self-dual locus is always where the
laws are cleanest (H ≡ 1 mod 4; e₁+e_n = n−2; zero majorization-incomparables). The mathematics keeps
saying: symmetry classes of the dissection ARE the arithmetic strata of the invariants.
