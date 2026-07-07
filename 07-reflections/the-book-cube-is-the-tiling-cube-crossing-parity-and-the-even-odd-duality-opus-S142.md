---
source: opus-2026-07-07-S142
status: excavation + fusion; F1-F4 verified exhaustively/randomized n<=14; T-join repair theorem
  (one line, upgrades HYP-2712a); affine parity law proved; mirror-breaking counts exact
tags:
  - crossing-number
  - two-page
  - tiling-model
  - even-graphs
  - parity
  - kleitman
  - unit-distance
  - even-odd-duality
---

# The book cube is the tiling cube: crossing parity and the even/odd duality

**opus-2026-07-07-S142.** Owner: *consider our past work on the unit distance problem, the
crossing number of K_n vs the bipartite graph, and the even/odd duality.* Excavated:
HYP-2712a (mac-mini-S7 era: the 2-page crossing optimum of K_n is attained with the page
indicator an **even graph** — "cr(K_n) is a cycle-space quadratic form whose optimum lives on
E_n", then flagged *no LRC leverage*), the cr(K₁₄) = 315 squarefree numerology (codex
HYP-2627/2629), and the Zarankiewicz/annulus tangent (HYP-3954). With this week's machinery
the thread fuses into theorems.

## 1. The identification (F1): the book cube IS the tiling cube

A 2-page book drawing of K_n fixes a spine order — **a Hamiltonian path** — and assigns each
chord a page bit. Chords of length 1 (spine edges) interleave nothing, so the crossing form
`Q(x) = #{interleaved chord pairs on the same page}` depends only on the `l ≥ 2` chords:
**exactly the m = C(n−1,2) tiles of the tournament tiling model.** The two GF(2)^m structures
the repo has been working — orientation-tilings and page-assignments — are the same cube
viewed from the same marked Hamiltonian path; the crossing landscape is a canonical quadratic
weight ON the tiling cube. (Verified n = 6..14.)

## 2. The affine parity law (F2) — Kleitman in book form, and the even/odd duality

Over GF(2), `[x_e = x_f] = 1 + x_e + x_f`, so the crossing **parity is affine**:

> `Q(x) ≡ C(n,4) + Σ_{deg_I(e) odd} x_e (mod 2)`, with `deg_I({a,b}) = (l−1)(n−1−l)`,
> `l = b−a`.

Consequences (verified exhaustively):
- **Odd n: no chord has odd interleave degree** — every 2-page drawing of K_n (n odd) has
  crossing parity **≡ C(n,4) (mod 2), drawing-invariant** (the Kleitman-flavored rigidity;
  e.g. every 2-page K₇ drawing has oddly many crossings).
- **Even n: parity is carried exactly by the EVEN-length chords** — at n = 14, the crossing
  parity of any book drawing is `1001 + #(even-length chords on page 1) (mod 2)`. The
  even/odd duality lands in the crossing world precisely as it does in the runner world
  (the ½-anchor/2-adic structure lives on even speeds; here parity lives on even chords).
- The Z(n)-vs-C(n,4) parity table shows one anomaly, **n = 12** (the only n where the
  optimum needs an ODD count of parity-carrying chords on page 1) — the repo's recurring
  k = 12 anomaly appears in the crossing column too.
- Structural placement of the codex numerology: `cr(K₁₄) = 315 ≡ 1001 ≡ C(14,4)` — the
  optimal drawing's even-chord page count is even; the 315-oddness is the parity law, not
  numerology.

## 3. The T-join repair theorem (upgrades HYP-2712a from census to one line)

HYP-2712a verified by sweeps that the 2-page optimum is attained with page-1 an even
subgraph. Mechanism, now a theorem: **spine edges cross nothing and the spine is a tree, so
the unique T-join on the spine repairs any page-1 degree-parity defect at zero crossing
cost.** Any optimal tile assignment completes to an even page set for free (verified
constructively at n = 6, 7, 8: all optima completable). *Cut space repairs; cycle space
optimizes* — the old "cut vs cycle is the universal seam" reflection landing exactly: the
crossing problem's parity lives on the cut side (spine tree), its optimization on the cycle
side (E_n), which is why the E_n lens partitions the landscape.

## 4. Mirror breaking at n = 8 (F4 + census)

`Q` is invariant under the grid involution σ (spine reversal) and under the page flip (the
LINE pairing of the blue/black metagraph — so **the crossing landscape descends to lines**).
But invariance of the form ≠ symmetry of the optimizers: σ-fixed (blue) optima exist at
n = 6 (4 of 12) and n = 7 (12 of 140), and **vanish at n = 8 (0 of 16: every optimal book
drawing of K₈ comes in a σ-pair; the strict-blue minimum is 20 > Z(8) = 18)**. The same
spontaneous mirror-breaking the fleet found for E[U] (mac-mini-S43) and the observer's
personal space (S138) — the SC/NS phenomenon — now in the crossing landscape, with an
even/odd-n alternation.

## 5. K_n vs K_{m,n}, and where unit distance sits (honest scope)

Guy's optimal cylindrical drawings put ⌈n/2⌉ and ⌊n/2⌋ vertices on two concentric circles:
`cr(K_n)` decomposes as two cliques + a bipartite (Zarankiewicz) part — and the runner
world's 2-adic split is the same decomposition: the half-shift `x ↦ x+½` fixes even-speed
runners and rotates odd ones, so the movie's crossings split into within-parity (clique)
and cross-parity (bipartite) braids. The bipartite part is where the ½-anchor machinery
lives; the crossing world's parity-carrying even chords and the runner world's
parity-interlacing families are two shadows of the same ℤ/2. The unit-distance connection
in the repo is currently only the HYP-3954 tangent (annulus/T-crossing ↔ equidistribution
extremality); beyond that tangent the unit-distance problem remains unmined here — flagged
honestly rather than forced.

## 6. What this opens

- The crossing weight `Q` as a canonical potential on the tiling cube: its interaction with
  the orientation semantics (H, μ-floor) — is the Q-landscape correlated with the H-gradient
  on the metagraph? (The line-invariance makes Q a line-metagraph functional — mac-mini-S47's
  THM-646 objects acquire a weight.)
- Prove the σ-pair structure at even n ≥ 8 (why does blue fail exactly there?) — the
  candidate mechanism is the parity law: at even n the parity-carrying chords are the even
  ones, and σ-fixed assignments constrain them in pairs.
- The bipartite/Zarankiewicz version of F1–F3 (spine = the two-circle order; the tiling
  model's analog for K_{m,n}) — natural next census.
