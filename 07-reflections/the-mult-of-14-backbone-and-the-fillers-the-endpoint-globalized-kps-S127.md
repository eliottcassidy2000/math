# The mult-of-14 backbone and the fillers: the endpoint structure globalized

*kind-pasteur-2026-07-11-S127 cont.65. Owner: "push the endpoint mult-of-14 structure further." The cont.64
endpoint fact — near `1/14` only the mult-of-14 runner covers — was local. Globalizing it turns the mult-of-14
runner into a **backbone** that grids the whole interval, and the covering-min into "the deepest backbone gap the
eleven fillers cannot fill." This lens cleanly separates the proved single-killer case (fillers idle → the
2-runner balance) from the open multi-killer case (a filler contests the gap).*

---

## From endpoint to backbone

A covering family has a multiple of 14, say `14m`. It is bad at **every** grid point `j/14` (`14m·(j/14) = jm ∈
ℤ`), and between them it has intermediate arcs at `k/(14m)`. So `14m` tiles `[1/14, 13/14]` on a **grid of arcs of
spacing `1/(14m)`**, each of width `1/(98m)`, leaving inter-arc **gaps of width `13/(196m)`**. The mult-of-14
runner is a *backbone*; the other eleven runners are *fillers*. Since the covering-min lonely point `t*` is where
the whole family is lonely, in particular `14m` is lonely there —

> **the lonely point `t*` always lies in a backbone gap** (verified: `‖14m·t*‖ ≥ M` for every case tested —
> single-killer, ladder, multi-killer).

## The grid-point modular law — why gaps live between grid points

Near a grid point `j/14`, a runner `b` covers iff `14 ∣ bj`, i.e. `b ≡ 0 (mod 14/gcd(j,14))`:

| `j` | `gcd(j,14)` | covered by |
|---|---|---|
| 1,3,5,9,11,13 | 1 | mult of 14 |
| 2,4,6,8,10,12 | 2 | mult of 7 |
| 7 | 7 | even |

A covering family has a multiple of 14, which — since `14 = 2·7` — is simultaneously a multiple of 7 and even, so
it **covers all 13 grid points by itself**. Hence no gap sits at a grid point; every gap is strictly *between*
grid points, inside a backbone slot. The backbone owns the grid; the fight is in the slots.

## Two runners bind — and who they are is the whole story

At `t*` exactly **two** runners bind (achieve the minimum `M` — the Chebyshev equioscillation, opus-S252). Which
two is the entire distinction between proved and open:

| family | `M` | `t*` | binding pair |
|---|---|---|---|
| single-killer `{1..12,182}` | 14/183 | 14/183 | **{1, 182}** = runner 1 + backbone |
| ladder `{1..12,364}` | 28/365 | 28/365 | **{1, 364}** = runner 1 + backbone |
| multi-killer `{1..11,13,84}` | 7/89 | 37/89 | **{5, 84}** = *filler* + backbone |
| multi `{1..10,13,22,84}` | 2/23 | 2/23 | **{1, 22}** = runner 1 + *filler* |

- **Single-killer: the fillers are idle.** Only runner 1 and the backbone bind; the AP fillers `{2,…,12}` are
  strictly lonelier than `M` at `t*`. So the single-killer covering-min is literally a **2-runner problem** —
  runner 1 against the backbone `14m` — which is exactly the slow-fast balance `M = [1/13]·(182c)/(182c+1)`
  **proved and machine-checked** (cont.60/61, `LRCSingleKillerLadder`). The backbone lens *is* the geometric
  content of that proof: in the backbone's first gap `(1/14, 1/13)`, only runner 1 (rising) and the backbone (its
  arc edges) are active.
- **Multi-killer: a filler contests the gap.** A small filler binds instead of runner 1 (runner `5` for
  `{1..11,13,84}`; runner `22` in place of the backbone for `{1..10,13,22,84}`), and `t*` moves to mid-interval.
  Now the backbone gap is *not* empty of fillers — the depth is set by a filler-vs-backbone (or filler-vs-filler)
  balance, and that is precisely the analytic case the Fourier route must close.

## Net — a clean lens, and where it stops

Globalizing the endpoint fact gives the **backbone + fillers** decomposition: the mult-of-14 runner grids
`[1/14,13/14]`, owns every grid point, and the covering-min is the deepest of its `~12m` slots that the eleven
fillers fail to fill. This *proves nothing new*, but it reorganizes the crux cleanly and pins the boundary
exactly: **the single-killer case is where the fillers are idle and the gap is a 2-runner (runner-1 vs backbone)
balance — done, in Lean; the multi-killer case is where a filler enters the backbone slot and contests it — the
open analytic residual.** It also re-derives, geometrically, why exactly two runners bind (equioscillation) and
why the backbone `14m` is the natural "observer grid" for the whole problem. A concrete next question the lens
poses: bound how deep a filler-contested backbone slot can be, uniformly in `m` — the multi-killer analogue of
the single-killer `1/(196m)` slot depth.

*Files: lrc14_backbone_structure_kps_S127.py (+.out). Extends cont.64 (endpoint fact), cont.60/61 (single-killer
balance, Lean), opus-S252 (equioscillation); complements the Fourier route (opus-S259–263) by isolating the
filler-contested slot as its target. HYP-6236.*
