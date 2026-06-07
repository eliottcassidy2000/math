# HYP-2319 — Towards u(21) = 57: the maximizer escapes the Eisenstein lattice (√−3, ≤47) into the Moser-spindle field (√−11, =57); a rigorous KST upper bound ≤71

**Session:** S641
**Status:** PARTIAL (lower bound cited+contextualized; rigorous upper bound proved; exact 57 open for us)
**Provenance forward:** math-lean `Math/UnitDistance/EisensteinNeighbors.lean` (sorry-free)
**Goal:** work towards a proof that `u(21) = 57`, the maximum number of unit distances among 21 points.

---

## 0. The two sides, honestly

`u(21) = 57` is an extremal value with a wide hardness gap between its two halves:
- **Lower** `u(21) ≥ 57`: a *construction* (the Moser slab `P_2^-`, HYP-2224/codex S648). Established.
- **Upper** `u(21) ≤ 57`: requires ruling out every 21-point config with ≥58 unit distances — heavy case
  analysis. The clean *general* method (KST) gives only `≤ 71` (proved below); closing `71 → 57` is the
  genuinely open part for us. I do **not** claim to prove the exact upper bound.

This session contributes: (i) a rigorous upper bound `u(21) ≤ 71`; (ii) the formalized **lattice
rigidity** that explains the structure; (iii) the new framing that the **maximizer escapes the
Eisenstein lattice into the √−11 spindle field** — a small-`n` shadow of the grid-optimality disproof.

---

## 1. The lattice is rigid — and far from optimal (formalized)

Two triangular-lattice points are at unit distance iff their difference has Eisenstein norm
`N(a,b) = a² − ab + b² = 1`, which has **exactly the 6 solutions** `±(1,0), ±(0,1), ±(1,1)` — the
hexagon, `6 = 2·3`, the cube-root units (S637/S638/S640). So **every lattice point has degree exactly
6**, capping any `n`-point lattice subset at `≤ 3n` and, by the boundary deficit, at Harborth's
`⌊3n − √(12n−3)⌋`:

```
  n = 21:  best compact triangular-lattice subset = 47 unit distances  (= Harborth floor).
```

Verified (`u21_unit_distance_bounds_s641.py`): the densest 21-point lattice cluster realizes exactly
**47**, matching `⌊3·21 − √249⌋ = 47`. So the Eisenstein lattice is **10 short** of `u(21) = 57`.

**Formalized (math-lean, sorry-free):** `eisenstein_unit_neighbours` —
`a² − ab + b² = 1 ↔ (a,b) ∈ {±(1,0), ±(0,1), ±(1,1)}` (via `nlinarith` bounding `3a² ≤ 4` ⟹ `|a| ≤ 1`,
then `interval_cases`/`decide`); `eisenstein_neighbour_count` — the set has card `6`. This is the rigid
degree-6 cap the optimum must leave behind.

---

## 2. A rigorous upper bound: u(21) ≤ 71 (K_{2,3}-free + cherry counting)

The unit-distance graph is **K_{2,3}-free**: two distinct points have **at most 2** common unit-distance
neighbours, because the two unit circles around them meet in ≤ 2 points. Counting cherries (paths of
length 2):
```
  Σ_v C(deg v, 2)  =  Σ_{pairs {u,w}} |common nbrs(u,w)|  ≤  2·C(21,2)  =  420.
```
With `Σ deg v = 2e` and Cauchy–Schwarz `Σ deg² ≥ (2e)²/21`, this gives
`4e²/21 − 2e ≤ 2·420`, i.e. `4e² − 42e − 17640 ≤ 0`, so

> **`u(21) ≤ 71`** (rigorous, general method).

This is the honest rigorous ceiling. The gap `71 → 57` is the case-analysis content the exact value
needs (Schade-style; out of scope here). The geometric core (`≤ 2` common neighbours) and the cherry
bound are the natural next *formal* targets (handoff).

---

## 3. The escape: √−3 (≤47) ≪ √−11 (=57) ≤ KST (≤71)

The `u(21) = 57` optimum is the **Moser slab `P_2^-`** (HYP-2224) — a **non-lattice** configuration. Its
geometry is the Moser **spindle**: two rhombi sharing a vertex, rotated by `θ` with `cos θ = 5/6`, so
`e^{iθ}` is a root of `3z² − 5z + 3`, **discriminant `−11`**. The optimum therefore lives in `ℚ(√−11)`,
**not** the Eisenstein lattice `ℚ(√−3)`.

```
  triangular lattice (ℚ(√−3), rigid degree 6) :  47
  u(21) optimum = Moser slab (ℚ(√−11) spindle) :  57      ← escapes the lattice
  rigorous KST upper bound (K_{2,3}-free)      : ≤ 71
```

> **The unit-distance maximizer escapes the Eisenstein lattice `√−3` into the spindle field `√−11`.**
> This is the small-`n` shadow of the grid-optimality *disproof* (the lattice is not optimal even at
> `n = 21`), and it is the **same `√−11`** as the Heegner chromatic tower: `χ ≥ 4` requires leaving
> `ℚ(√−3)` for `ℚ(√−11)` (the Moser spindle, HYP-2277). So maximizing unit distances and forcing the
> 4th colour leave the Eisenstein lattice through the **same door** — the Moser spindle, `disc = −11`.
> The lattice's formalized rigidity (§1, degree exactly 6) is precisely the cap both must escape.

---

## 4. A repo inconsistency to flag

`HYP-2170` records "`n = 22` max unit distances = 49 (Harborth triangular-lattice optimum)". But `49` is
the **lattice** optimum (`⌊3·22 − √261⌋ = 49`), **not** `u(22)`: the Moser slab `P_2^+` gives `60`
(HYP-2224), and `u(22) ∈ {60, 61}`. Calling `49` "the max unit distances" conflates the lattice optimum
with the global maximum — the very conflation this session corrects (§3: lattice ≠ global max for these
`n`). Recommend annotating HYP-2170. (Not opening a court case: HYP-2170's deeper LRC/Cayley bridge is
unaffected; only its "max = 49" phrasing is wrong.)

## 5. Connections & handoffs
- **HYP-2224 (S648):** the `P_2^-` = 57 lower bound and the Moser-slab carrier; this session adds the
  rigorous upper bound, the lattice foil, and the field identification.
- **HYP-2277 (S688):** the `√−11` Moser-spindle rotation field; §3 shows the unit-distance *maximizer*
  uses the same field as the chromatic `χ ≥ 4` escape — a maximization↔colouring unification through
  `disc = −11`.
- **S637/S638/S640:** the `6 = 2·3` Eisenstein units / cube root are the lattice's rigid neighbourhood.
- **Handoff (the real upper-bound content):** (a) formalize the geometric `≤ 2 common unit-neighbours`
  (two unit circles meet in ≤ 2 points) and the cherry/KST bound; (b) push `71 → 57` via the degree /
  forbidden-subgraph case analysis (the Schade method) — the actual proof of `u(21) ≤ 57`.
