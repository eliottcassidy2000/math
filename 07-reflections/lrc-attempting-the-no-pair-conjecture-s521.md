# Attempting the no-pair conjecture and LRC-given-the-conjecture: both reduce to the blocker core (S521)

*claudebox-2026-06-01-S521. An honest attempt at (1) proving the no-pair reduction
conjecture and (2) proving LRC assuming it. Result: the no-blocker case is fully
proved; both remaining tasks bottom out at the same blocker / fully-covered core.
Builds on the balanced apex-pair theorem (Thm A), the tight-locus theorem (Thm B),
and the base reduction.*

## The decisive partition

Every primitive speed set falls into one of three cells by two `mod n` predicates
(`n = m+1`): **blocker** = some speed `≡ 0 (mod n)`; **pair** = some pair of
speeds sums to `0 (mod n)`.

| cell | status |
|------|--------|
| no-blocker | **LRC PROVED**: `t = 1/n` is lonely (`||v_i/n|| ≥ 1/n` since no `v_i ≡ 0`) |
| blocker + no-pair | the no-pair conjecture's real content (residual) |
| blocker + pair | not addressed by the conjecture |

Verified (n=6): no-blocker sets never fail (`M ≥ 1/n`, some tight like AP); the
tight/sub-tight risk lives only in blocker cells.  So **a counterexample, if any,
must have a speed divisible by n** — the base reduction, now seen as the frame for
both tasks.

## Task 1 — prove the conjecture (no pair `≡ 0 mod n` ⇒ `M > 1/n`)

**Proved sub-case (rigorous).** If a no-pair set has **no blocker**, then `t = 1/n`
gives `M ≥ 1/n`, and by Theorem B (`M = 1/n ⇒ n | v_i+v_j`) the no-pair condition
forbids `M = 1/n`; hence `M > 1/n`.  So the conjecture holds for all no-pair sets
without a speed divisible by `n`.

**Residual: blocker + no-pair.**  Here the unique blocker `v_b = nw` sits at the
observer at `t = 1/n`.  Key structural fact from no-pair: the other runners'
residues mod `n` cannot contain both `r` and `-r`; in particular **not both `1` and
`n-1`**.  Perturb `t = 1/n + ε`:
- the blocker moves to distance `||nwε|| ≈ nw|ε|` — safe once `|ε| ≥ 1/(n^2 w)`;
- a runner at residue `1` (distance `1/n`) moves *away* for `ε > 0`; one at `n-1`
  moves away for `ε < 0`.  Since no-pair forbids both `1` and `n-1`, one sign of
  `ε` moves all the boundary (`d=1`) runners *outward*;
- a runner at level `d_i = min(r_i, n-r_i) ≥ 2` (distance `d_i/n`) stays safe while
  `|ε| ≤ (d_i - 1)/(n v_i)`.

So a good `ε` exists **iff** `1/(n^2 w) ≤ (d_i-1)/(n v_i)` for every level-`≥2`
runner, i.e. `v_i ≤ (d_i - 1) v_b` for all `i` (`v_b = nw`).  **This proves the
blocker+no-pair case whenever the blocker dominates** the "tight" (`d_i = 2`)
runners (`v_i ≤ v_b`).  It **fails** when a fast non-blocker sits at level 2
(`v_i > v_b`, `d_i = 2`) — a runner close to the observer and faster than the
blocker — which needs a time away from `1/n` (a different binding pair via
Theorem A).  That residual is open.

So: conjecture **proved** for (no-blocker) ∪ (blocker-dominant no-pair); **open**
for no-pair sets with a fast level-2 non-blocker.

## Task 2 — prove LRC assuming the conjecture

Assume the conjecture (so all no-pair sets have `M > 1/n`).  Combined with the
proved no-blocker case, the only sets not yet handled are **blocker + pair**:
a speed `≡ 0 (mod n)` *and* a pair summing to `0 (mod n)`.  The conjecture says
nothing about these (it is about no-pair sets).  So

> **LRC(n)  (given the conjecture)  ⟺  `M ≥ 1/n` for all blocker+pair sets.**

This is a strictly smaller class, but it is exactly where the extremal AP-type sets
live (AP has a blocker only when `n | v_i` for some `i`; e.g. `{1,...,n-1}` has no
blocker — it is no-blocker+pair and handled by `t=1/n`!).  In fact the tight
extremiser `{1,...,n-1}` is **no-blocker**, so it is already covered.  The
blocker+pair residual is the genuinely hard, non-extremal-looking core: speed sets
forced off `t=1/n` by a blocker yet carrying the pair-resonance.  The conjecture
does not crack it.

**Honest verdict.** Assuming the conjecture does **not** yield LRC: it removes the
blocker+no-pair slice, leaving blocker+pair. Both the conjecture's residual and the
LRC-given-conjecture residual are **blocker** cases — the same fully-covered /
divisible-by-`n` core every S521 thread reached. The base reduction already
isolates this core; the pair/no-pair refinement sorts it further but the divisible-
by-`n` blocker is the irreducible obstruction.

## What is genuinely proved this turn

1. **LRC for all no-blocker sets** (`t = 1/n`), strict (`M > 1/n`) exactly when no
   pair sums to `0 mod n` (Theorem B). In particular the AP extremiser and all
   sets without a speed divisible by `n` are done.
2. **The conjecture for blocker-dominant no-pair sets** (the `t = 1/n ± ε`
   perturbation), using the no-pair fact that residues `1` and `n-1` cannot
   coexist.
3. A clean statement of the irreducible residual: **blocker + (pair, or a fast
   level-2 non-blocker)** — speed sets with a speed `≡ 0 (mod n)` and a runner
   pinned near the observer that cannot be freed by a local perturbation.

## Seed

The residual is a **two-gap race** (THM-386) localized at the blocker: free the
blocker (push `|ε|` up) without driving a fast near-observer runner below `1/n`.
The no-pair structure already broke the symmetry on one side (no `1` and `n-1`
together); the open need is a *global* time (some binding pair `(i,j)` via Theorem
A with `v_i + v_j` not a multiple of `n`) that simultaneously clears the blocker
and the fast level-2 runner. Proving that such a pair always exists for blocker
sets would close LRC — and is exactly the fully-covered core in gap-function
clothing.
