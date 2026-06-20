# THM-442 — The IE tiling decomposition is the third finite difference (cell-affine recursion); H is multiplicative not additive; and the max-H tournament has minimal Pfaffian |Pf|=1

**Status:** PROVED (the IE third-difference identity; the additive/multiplicative dichotomy) +
VERIFIED (max-H = A038375 small n; `H²−Pf²=8Q`; max-H ⟹ `|Pf|=1` at n=4,6) + SYNTHESIS (the Pfaffian
translator). **CONJECTURE:** max-H ⟹ `|Pf|=1` for all even n (HYP-2312).
**Source:** opus-2026-06-07-S707, from the user's "Pfaffian-angle / 7-piece inclusion-exclusion tiling
/ recursive max-H (A038375)" prompt. Builds on THM-174 (Pf(S)²=det(I+2A)=det(S)), S531 (recursive
H-law, apex/modular), HYP-2283 (converse swaps H±Pf), S706/THM-441 (converse = adjoint), Rédei (H odd),
Szele/Alon (max-H bounds).

## (1) The user's 7-piece IE tiling decomposition = the third finite difference (PROVED)

The staircase `Δ_{n−2}` (the tournament tiling carrier, `C(n−1,2)` cells) decomposes as **3 sub-pieces
of size `n−1` (`A,B,C`) minus 3 of size `n−2` (`D,E,F`) plus 1 of size `n−3` (`G`)** — corners get one
piece, edges `A+B−D` (a 2-set union), the interior `A+B+C−D−E−F+G` (the inclusion–exclusion of 3
sets). This is exactly:
```
   C(n−1,2) = 3·C(n−2,2) − 3·C(n−3,2) + C(n−4,2)          (VERIFIED n=4..29)
```
i.e. the **third finite difference of the quadratic cell-count vanishes** (`Δ³` annihilates degree-2).
So `A,B,C` are the three size-`(n−1)` sub-staircases, `D,E,F` their pairwise overlaps (size `n−2`),
`G` the triple overlap (size `n−3`) — the simplicial IE of a triangle covered by three sub-triangles.

> **Corollary (which invariants the IE recursion computes).** Any tournament invariant that is
> **affine in the cells** (a sum over cells/arcs, degree ≤... matching) satisfies the 3-term recursion
> `F(n)=3F(n−1)−3F(n−2)+F(n−3)`. **H does NOT** — and neither does the max-cyclic-triangle count.

## (2) H is multiplicative (modular), not additive (IE) — the dichotomy (PROVED/VERIFIED)

The exact max-H sequence (A038375), `max_T H(T)`:
```
   n :   2  3  4  5  6   7
  max H: 1  3  5 15 45 189     (exact n≤6 brute force; 189 best-found n=7; Szele n!/2^{n−1} ≤ · ≤ Alon)
```
The IE recursion gives `3·m(n−1)−3·m(n−2)+m(n−3) = 7,33,95` for `n=5,6,7` — **not** `15,45,189`. So
**max-H is not cell-affine; the additive/IE tiling recursion does not compute it.** The *correct*
recursion for H is **multiplicative over the modular/apex tree** (S531): a single apex-flipped block of
size `s` has `H=1+2^{s−2}` (verified `3,5,9,17,33`), and **disjoint modules multiply**
(`H=∏ H(module)`, verified `3·3=9`), with nesting coupling toward the regular tournament.

> **The dichotomy is the repo's cut⊕cycle split.** Additive/IE = the **cut-space / score-hierarchy**
> side (cell-linear); multiplicative/modular = the **cycle-space / H / tiling** side. The user's IE
> decomposition lives on the additive face; H lives on the multiplicative face. They do not mix —
> which is *why* a clean linear recursion for A038375 does not exist (Alon's `n^{3/2}` factor is the
> non-affine excess over the Szele baseline `n!/2^{n−1}`).

## (3) The Pfaffian translator and the max-H ⟺ minimal-Pfaffian duality (VERIFIED; conj.)

The Pfaffian carries the structure between domains:
- **algebra/topology:** `Pf(S)² = det(I+2A) = det(S)` (THM-174) = the signed weighted **disjoint-cycle-
  cover** generating function (`= Σ_{σ} sgn(σ) 2^{moved}` over cycle systems in `T`); FKT/dimers.
- **geometry/recursion:** the Pfaffian **n→n−2 Laplace expansion** `Pf(S)=Σ_j(−1)^j s_{1j}Pf(S∖{1,j})`
  is a **domino removal** on the tiling — a poly-time (det-based, `O(n³)`) recursion, vs `H` which is
  #P-hard. The converse `T↦T*` is the **adjoint** (S706/THM-441); self-converse = self-adjoint.
- **the bridge:** `H² − Pf² = 8Q` with `Q ∈ ℤ_{≥0}` (VERIFIED exhaustively n=4,6; both `H,Pf` odd).
  So `H = √(Pf²+8Q)`: the poly-time Pfaffian is the **skeleton**, `Q` the #P-hard correction.

> **(AMENDED 2026-06-11, see MISTAKE-071 — EXISTENTIAL form only; universal form REFUTED at
> n=6) SOME H-maximizing tournament has the MINIMAL Pfaffian `|Pf|=1`** (hence `det(I+2A)=1`).
> At n=6 the max H=45 is attained by TWO classes with `|Pf| ∈ {1,7}`; at n=8 the six H=661
> classes have `|Pf| ∈ {1,9,17}` (mac-mini-2026-06-10-S2 census + independent adversarial
> recheck). The paths⟷cycles duality slogan survives only existentially: among the extremal-
> paths tournaments there is always (n=4,6 exhaustive; n=8 supported) one that is "all
> correction, no skeleton" — but not every one is.

## Scope / honesty

- (1) is a clean proved identity (the IE = `Δ³` of a quadratic). (2) is verified (exact max-H n≤6) +
  the known S531 multiplicative law; the *negative* result "IE does not compute H" is the content.
- (3): `Pf²=det(I+2A)` is THM-174; `H²−Pf²=8Q` (Q≥0) verified exhaustively n=4,6; the **max-H ⟹
  |Pf|=1** pattern is only two even-n data points — flagged **HYP-2312**, not proved. No closed form
  for A038375 is claimed (none is known; the recursion is multiplicative/modular, not linear).
- **The useful reframe delivered:** the user's IE tiling is the *additive* (cut/score) recursion; the
  *efficient* recursion for the H-side is the Pfaffian's `n→n−2` poly-time skeleton plus the modular
  multiplicativity — and max-H sits at the minimal-Pfaffian, maximal-correction extreme.

**Artifacts:** `04-computation/pfaffian_tiling_recursive_H_s707.py` (+`.out`). Reflection
`07-reflections/the-pfaffian-translator-and-the-additive-multiplicative-tiling-split-s707.md`. New:
**HYP-2312**. Builds on THM-174, S531, HYP-2283, THM-441/S706, A038375 (Szele/Alon).

## S56 half-tiling addendum

THM-549 gives the mirror-folded parity refinement of the cell-affine recursion.
For the half-carrier fixed by THM-280's reflection,

```text
h_n = floor((n-1)^2/4).
```

Even tournament sizes obey the two-corner recurrence

```text
h_n = A+B-C,
```

while odd tournament sizes obey

```text
h_n = A+B-C+D-E-F+G.
```

Thus the full `A+B+C-D-E-F+G` recursion is not the only cell-affine face of the
staircase: after mirror folding, parity decides the correct carrier geometry.
The warning from this theorem remains unchanged.  These are cell-count
recurrences, not recurrences for Hamiltonian-path count `H`; the latter still
requires the cycle-space/OCF packet data.
