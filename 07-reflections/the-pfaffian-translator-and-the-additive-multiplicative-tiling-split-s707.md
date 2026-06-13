---
source: opus-2026-06-07-S707 (user: Pfaffian angle as translator; 7-piece IE tiling; recursive max-H)
status: SYNTHESIS + a new pattern. The user's 7-piece inclusion-exclusion tiling decomposition
  (A,B,C size n-1; −D,−E,−F size n-2; +G size n-3) IS the third finite difference annihilating the
  quadratic staircase cell-count: C(n-1,2)=3C(n-2,2)-3C(n-3,2)+C(n-4,2) (verified n=4..29). It computes
  any CELL-AFFINE invariant by F(n)=3F(n-1)-3F(n-2)+F(n-3) — but H is MULTIPLICATIVE (modular, S531),
  not additive, so the IE does NOT compute max-H (A038375=1,3,5,15,45,189; IE gives 7,33,95). This is
  the repo's cut⊕cycle split: additive/IE = cut/score side, multiplicative/modular = cycle/H side.
  PFAFFIAN as translator: Pf(S)²=det(I+2A)=det(S) (THM-174, cycle-covers); Pf recurses n→n-2 (domino,
  poly-time) vs #P-hard H; H²−Pf²=8Q (Q≥0, verified n=4,6); NEW: the max-H tournament has MINIMAL
  |Pf|=1 (max-paths ⟺ min-cycle-cover-Pfaffian, a paths/cycles duality) — HYP-2312. THM-442.
tags: [pfaffian, hamiltonian-paths, A038375, max-H, tiling, inclusion-exclusion, third-difference,
  finite-difference, cut-cycle, modular-decomposition, additive-multiplicative, szele, alon,
  cycle-cover, paths-cycles-duality, translator, recursion, tournaments]
---

# The Pfaffian translator, and why the tiling inclusion-exclusion computes scores but not paths

**Prompt (user):** use the Pfaffian angle to translate between topology/geometry/graphs/tournaments/
algebras; come to a more efficient recursive understanding of max-H (A038375) via the tiling model;
the 7-piece IE decomposition (A,B,C size n−1 positive; D,E,F size n−2 negative; G size n−3 positive);
find other recursive truths.

Two things came out: the user's decomposition is **exactly** a clean classical object (the third
finite difference), and chasing *why it fails for H* led to a genuine new pattern (max-H ⟺ minimal
Pfaffian).

## 1. The 7-piece decomposition IS the third finite difference

Corners get one piece (`A`, `C`, `E`); edges a 2-set union (`A+B−D`); the interior the 3-set
inclusion–exclusion (`A+B+C−D−E−F+G`). That is the simplicial IE of a triangle covered by three
sub-triangles, and on the cell-count it reads
```
   C(n−1,2) = 3·C(n−2,2) − 3·C(n−3,2) + C(n−4,2)      (verified n=4..29)
```
— the **third finite difference `Δ³` of a quadratic vanishes**. So `D,E,F` are the pairwise overlaps,
`G` the triple overlap, and the alternating `+3,−3,+1` are the binomial IE signs. Clean and exact.

> **What it computes:** any invariant **affine in the cells/arcs** satisfies `F(n)=3F(n−1)−3F(n−2)+
> F(n−3)`. The cell-count, arc-count, score-sum — anything cut-space/linear — obeys it.

## 2. Why it does NOT compute max-H — the additive/multiplicative wall

The exact answer (brute force): `max H = 1,3,5,15,45,189` for `n=2..7` (A038375). The IE recursion
gives `3·m(n−1)−3·m(n−2)+m(n−3) = 7,33,95` for `n=5,6,7` — **wrong**. Max-H is **not cell-affine**.

The right recursion for `H` is **multiplicative over the apex/modular tree** (S531): a size-`s` apex
block has `H=1+2^{s−2}`, and disjoint modules **multiply** (`H=∏H(module)`), nesting couples toward
the regular tournament. So:

> **The user's IE tiling is the *additive* (cut-space / score-hierarchy) recursion; H lives on the
> *multiplicative* (cycle-space / tiling) face.** This is the repo's GF(2) **cut ⊕ cycle** split
> (CLAUDE.md): base-path arcs = cut = score (additive, IE-friendly); wiggly arcs = cycle = H
> (multiplicative). They don't mix — which is *why* A038375 has no simple linear recursion (Alon's
> `n^{3/2}` factor over the Szele baseline `n!/2^{n−1}` is precisely the non-affine excess).

So the "more efficient recursion for max-H" is honestly: **there is no additive one** — the tiling IE
is the wrong face. The efficient handle is the multiplicative-modular law (already S531) and the
Pfaffian skeleton (below), not the IE.

## 3. The Pfaffian as translator — and the paths/cycles duality

The Pfaffian carries structure across the user's five domains:
- **algebra ↔ topology:** `Pf(S)² = det(I+2A) = det(S)` (THM-174) `=` the signed weighted
  **disjoint-cycle-cover** generating function (Clifford-algebra top form; FKT dimers).
- **geometry ↔ recursion:** the Pfaffian **n→n−2** expansion `Pf(S)=Σ_j(−1)^j s_{1j}Pf(S∖{1,j})` is a
  **domino removal** on the tiling — and it is **poly-time** (`det`, `O(n³)`), unlike the #P-hard `H`.
- **graphs ↔ tournaments:** the converse `T↦T*` is the **adjoint** of the adjacency-convolution
  (S706/THM-441); self-converse = self-adjoint = the worry-set.
- **the bridge:** `H² − Pf² = 8Q`, `Q∈ℤ_{≥0}` (verified exhaustively n=4,6; both `H,Pf` odd). So
  `H = √(Pf²+8Q)`: the **Pfaffian is the poly-time skeleton, `Q` the #P-hard correction.**

> **New pattern (HYP-2312, verified n=4,6):** the **H-maximizing tournament has minimal `|Pf|=1`**
> (so `det(I+2A)=1`, one net signed cycle cover; `H²=1+8Q` with `Q` maximal). **Maximizing
> Hamiltonian PATHS forces minimizing the cycle-cover Pfaffian** — an open⟷closed (paths⟷cycles)
> duality: the extremal-paths tournament is *all correction, no skeleton.* If this holds for all even
> `n`, the max-H search restricts to the `Pf=±1` ("Pfaffian-trivial") tournaments — a real reduction
> of the A038375 extremal problem to a `det(I+2A)=1` constraint.

## 4. Other recursive truths (the harvest)

- **`Δ³` IE recursion** for every cell-affine invariant (the user's decomposition).
- **Pfaffian `n→n−2` poly-time recursion** (domino removal) — the efficient skeleton.
- **Modular multiplicativity** `H=∏H(module)`, apex block `1+2^{s−2}` (S531, re-confirmed).
- **`H²−Pf²=8Q`** bridge; **max-H ⟺ `|Pf|=1`** (new).
- **The dichotomy is canonical:** cut/score = additive (`Δ³`); cycle/H = multiplicative (modular) —
  the two faces never satisfy each other's recursion.

## 5. Honest status

- **Proved:** the IE = `Δ³` identity (cell count); the additive/multiplicative dichotomy framing.
- **Verified:** max-H `1,3,5,15,45,189` (exact n≤6, best-found 189 at n=7); `Pf²=det(I+2A)`;
  `H²−Pf²=8Q` (Q≥0) exhaustive n=4,6; max-H ⟹ `|Pf|=1` at n=4,6; single-apex `1+2^{s−2}`; modular
  `3·3=9`.
- **Conjecture (HYP-2312):** max-H ⟹ `|Pf|=1` for all even `n` (two data points only).
- **Honest negative:** the user's IE tiling does **not** give a recursion for max-H (it is additive,
  H is multiplicative); no closed form / linear recurrence for A038375 is claimed (none known). The
  deliverable is the *reframe* (which face computes what) + the Pfaffian skeleton + the new
  max-H ⟺ minimal-Pfaffian duality.

**Artifacts:** `04-computation/pfaffian_tiling_recursive_H_s707.py` (+`.out`). Theorem **THM-442**.
New **HYP-2312**. Builds on THM-174 (Pfaffian identity), S531 (recursive H-law / apex / modular),
HYP-2283 (H±Pf), S706/THM-441 (converse=adjoint), Szele/Alon (max-H bounds), CLAUDE.md cut⊕cycle.
