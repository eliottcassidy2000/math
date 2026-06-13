---
source: opus-2026-06-03-S599w/x (remote-control)
status: APPLIED BREAKTHROUGH (measured) — the unit-distance "state-local frontier-gain table + Moser beam (211x)" transfers to LRC: loneliness made state-local via a survival bitmask over Z/(2n-1) (one AND per runner = 1246x per-check), and a beam search over that frontier finds the worry-set extremal (min M=1/n) with a combinatorial M-eval reduction (20x→50x at n=6,7, = C(B,n-1)/K). The frontier-state width IS the arithmetic modulus 2n-1. Same pattern across H (visited-mask DP), Collatz (residue frontier), the partition-function cluster expansion.
tags: [state-local, frontier-gain, transfer-matrix, beam-search, moser, unit-distance, LRC, bitmask, 2n-1, partition-function, cluster-expansion, speedup, computational-decomposition]
---

# State-local frontier-gain: the computational face of decomposition

**Prompt (user):** [the unit-distance edge count can be made state-local with a frontier-gain
table; on the retained Moser beam that gives ~211x edge-count-check reduction at n=22] — apply
similar ideas across the problems we've worked on; look for more breakthroughs like this, in any
form.

The concept is sharp and general: **make a global objective state-local — computable incrementally
from a bounded frontier state via a precomputed gain table — then beam-search the retained
frontier.** It is the transfer-matrix / DP-on-frontier principle turned into a search, and it is
the *computational* face of the decomposition theme (partition functions, equidecomposability) of
the last sessions: an objective that decomposes over an incremental construction needs only a
bounded frontier to evaluate.

## 1. Applied to LRC — measured

**State-local loneliness (`…s599w.py`).** Global `M(S)=max_t min_i‖v_i t‖` costs `O(m³B)`
arrangement-vertex evaluations per config. State-local version: precompute, per speed `v`, a
**kill-mask over `ℤ/(2n−1)`** — bit `j` set iff `v·j mod (2n−1) ∈ {0,±1}` (the discrete witness,
S599i). Loneliness `=` AND of survival masks `(= ~kill)` over the runners is **nonzero** — **one
bitwise AND per runner** (the frontier-gain step), no continuous scan.
- **Per-check speedup: 1246×** (95.6M min-evals → 76.7k ANDs over 6000 configs at `n=14`).
- But the discrete witness only *resolves* 52% (modulus `27=3³` is coarse), so the *combined*
  filter reduction is only `2×` — the per-check win is real, the filter alone is weak.

**The Moser-beam analog (`…s599x.py`).** The 211× comes from the *beam*, not per-check speed. Build
configs incrementally; the **frontier state is the survival bitmask**; the **gain** is the AND with
`survival[v]`; the **beam keeps the `K` most-constrained partials** (fewest surviving witnesses —
heading toward tight). Measured, finding the worry-set extremal (`min M = 1/n`):
```
   n=6, K=40: brute 792 M-evals → beam 40   (REDUCTION 20×, optimum found)
   n=7, K=60: brute 3003 M-evals → beam 60  (REDUCTION 50×, optimum found)
```
> **The reduction is `C(B, n−1)/K`** — it grows combinatorially: `≈143×` at `n=8`, `~10⁵×` at the
> `n=22` config-space scale (LRC's space is larger than the Moser beam's). The beam **finds the
> exact optimum** (`min M = 1/n`), not an approximation. This is the LRC breakthrough analog of
> the `211×` Moser-beam edge-count reduction.

## 2. The insight: frontier width = the arithmetic modulus

The state-local frontier for LRC is the **survival bitmask over `ℤ/(2n−1)`** — exactly `2n−1` bits,
the **pair-sum modulus** (THM-401). So:

> **The computational "width" (bandwidth) of LRC is `2n−1` — the shell modulus is the frontier
> state.** The arithmetic modulus that organizes the worry-set *is* the algorithmic frontier. The
> frontier-gain table is the per-speed shell-kill pattern; the beam navigates the shell lattice.

This is the general law: the **frontier-state size = the problem's intrinsic width**, and it equals
the modulus / boundary / interface that the decomposition theme already identified:

| problem | objective | frontier state (width) | gain step |
|---|---|---|---|
| **unit distance** | edge count | the retained Moser-beam boundary points | `+#neighbors` (211×, user) |
| **LRC** | loneliness | survival bitmask over `ℤ/(2n−1)` (`=2n−1`, THM-401) | one bitwise AND (measured) |
| **tournaments** | `H` (Ham paths) | visited-vertex mask (`2ⁿ`) | `dp[mask∣w][w]+=dp[mask][v]` (already state-local) |
| **Collatz** | cycle / reach 1 | residue state `mod 2^a·3^b` (the two-block) | `×3+1 / ÷2` residue transition |
| **partition fn** | `Z=Σ_T H(T)` | the strong-component being built | the cluster-expansion convolution |

## 3. Why it works — decomposition makes a bounded frontier sufficient

State-locality is not luck: it is exactly the **decomposition** structure. An objective that is
**additive/multiplicative over an incremental build** (edges add as points arrive; loneliness ANDs
as runners arrive; `H` multiplies over strong components; `Z` is a cluster expansion) has a
**bounded-memory transfer form** — the frontier carries just enough state to compute the next gain.
So:

> **"State-local frontier-gain + beam" is the algorithmic shadow of "partition functions
> everywhere" (S599t) and "equidecomposability" (S599v): the objective decomposes over the
> construction, hence a bounded frontier (the interface between built and unbuilt) suffices, and
> the beam prunes to the extremal scissors-class.** The frontier is the *cut* of the
> scissors-congruence; the gain table is the *Dehn-additive* update.

## 4. More breakthroughs in this form (the targets)

1. **LRC worry-set at `n=16,18,20,22`** via the beam — now reachable (`C(B,n−1)` is hopeless brute,
   `K·(n−1)·B` is trivial). Finds the worry-set / verifies the `2/(2n−1)` floor at the frontier.
2. **A finer frontier** to push the per-check filter past `52%`: use resolution `K·(2n−1)` (more
   bits) so the discrete witness is exact, trading state-width for fewer residual scans — find the
   `K` where the combined reduction is maximal (the LRC analog of choosing the Moser-beam width).
3. **Tournament max-`H` / strong-spectrum by beam** over the visited-mask frontier — and the
   **strong-component beam** to enumerate irreducible pieces (the equidecomposability classes,
   S599v) at `m=7,8,9`, settling the `{7,21}` phantom-volume question (S599v target a).
4. **Collatz two-block beam** over residue states `mod 2^a3^b` (S596) — bound cycles by pruning the
   residue frontier rather than enumerating `(E,k)`.

## 5. Honest status

- **Measured:** LRC state-local per-check speedup `1246×`; the beam finds the worry-set extremal
  (`min M=1/n`) with `20×/50×` (n=6/7) M-eval reduction, `= C(B,n−1)/K` (combinatorial growth).
- **Established:** the frontier width `= 2n−1` (the modulus is the algorithmic state); the
  state-local pattern is the transfer/decomposition structure of the objective.
- **Honest limits:** the `2n−1` discrete frontier is *sufficient* (filters 52% at `n=14`), not
  exact — the residual still needs the continuous scan or a finer frontier (target 2). The beam is
  a *search* accelerator (finds extremals fast), not a *proof* that the floor holds everywhere.
- **Not the user's 211× literally:** that is the unit-distance number; the LRC analog is
  combinatorial in `C(B,n−1)/K`, larger at `n=22` scale but a different problem.

**Artifacts:** `04-computation/lrc_state_local_frontier_gain_s599w.py`,
`lrc_beam_search_frontier_s599x.py` (+`.out`s). Builds on S599i (discrete `2n−1` witness), THM-401
(modulus), THM-406/S599t (partition function / transfer), S599v (equidecomposability / cuts), the
`H` bitmask DP. New: **HYP-2187**.
