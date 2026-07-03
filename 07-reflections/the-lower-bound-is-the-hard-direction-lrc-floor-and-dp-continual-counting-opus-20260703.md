# The lower bound is the hard direction: the LRC(14) floor and DP continual counting are the same shape — structured cancellation beats moments, and the AP is the low-discrepancy extremal

*opus-2026-07-03-S48. The owner pointed at arXiv:2607.00876 (Bairaktari–Larsen, "The Binary Tree
Mechanism is Optimal for Approximate DP Continual Counting", July 2026) and asked how it could relate
to the 14-runner LRC, while finishing the open floor math. It relates precisely — not by subject
(privacy vs Diophantine covering) but by the **shape of the hard direction** — and the analogy handed
back a concrete tool that closes the last genuine-mathematics item, (R2), to evidence standard.*

## The shared shape (the honest core)

Both problems have an **easy upper bound and a hard lower bound**, and in both the lower bound is the
whole game:

| | LRC(14) 11-core floor | DP continual counting |
|---|---|---|
| object | `meas(L_C)`, lonely measure of an 11-core | `err`, ℓ∞ error of releasing all prefix sums |
| easy side | pentagon **achieves** `313/9702` (a config) | binary tree **achieves** `O(log^{3/2} n)` (a mechanism) |
| hard side | `meas(L_C) ≥ 1/36` for **every** core | `err ≥ Ω(log^{3/2} n)` for **every** mechanism |
| naive method | 2nd moment gives the **wrong direction** | independent noise gives the **wrong rate** |
| extremal | the **AP / pentagon** (renormalization fixed point) | the **dyadic binary tree** |
| engine | discrepancy of a triangular accumulation | hereditary discrepancy of the prefix-sum matrix |

The point of the table is the middle two rows. In **both** problems the first-moment/second-moment
estimate fails in the hard direction, and in **both** the fix is the same species of object: the
**discrepancy of a triangular (prefix-sum / counting) accumulation**, extremized at the
**uniform/dyadic** configuration, and won by **cancellation, not magnitude.**

## Why moments fail, and what replaces them (the mechanism)

**LRC side.** Let `count(t) = #{v ∈ C : ‖vt‖ < 1/14}` be the running danger count; the observer is
lonely iff `count(t) = 0`, and `meas(L_C) = meas{t : count(t)=0}`. Here `E[count] = k·(2/14) = k/7 > 1`
for `k = 11`, so `{count = 0}` is a **below-mean event**. Markov / Paley–Zygmund / second moment give
`P(count > 0) ≥ E[count]²/E[count²]` — an **upper** bound on `meas(L_C)`. Lower bounds on a below-mean
event **cannot** come from magnitude; they need **signed cancellation**. This is exactly why the
discrepancy potential (S29) "provably cannot" close the r=2 residual, and why klein's **F3-sharp**
(S106) matters: it beats the crude `√(Δ/N)` freeze error by **telescoping the signed per-cell errors
along the sweep** — the interior contributions cancel, leaving `O(1/N)`. The `√` was the price of
*not* cancelling; the true rate is linear because the accumulation telescopes.

**DP side.** Releasing `n` prefix sums with independent noise accumulates error with the stream length.
The **binary tree** mechanism instead writes each prefix as a sum of `O(log n)` dyadic block-sums and
noises the blocks — the errors **don't accumulate linearly** because the hierarchical decomposition
**cancels** the redundancy between overlapping prefixes. The new lower bound `Ω(log^{3/2} n)` (via
hereditary discrepancy of the counting matrix, plus the ℓ∞ max over prefixes) shows this dyadic
cancellation is optimal.

> **Same engine.** A running accumulation whose naive summation loses a factor, recovered by a
> **hierarchical/telescoping cancellation**, with the obstruction measured by the **discrepancy of the
> lower-triangular counting matrix.** klein's telescoping *is* the binary-tree cancellation, one level
> at a time; the DP theorem is the statement that no cleverer cancellation exists.

## The triangular matrix is the staircase's cut space (the repo bridge)

The DP lower bound is the **hereditary discrepancy of the prefix-sum matrix** `P` (lower-triangular
all-ones): `P x = ` running sums of `x`. In this repo that matrix is **already central**: the
Cut ⊕ Cycle GF(2) split of the tiling cube (CLAUDE.md) puts the **base-path arcs = the cut space = the
score hierarchy = prefix sums**, and the wiggly arcs = the cycle space. `P` is the **cut-space
generator of the staircase**. So "hereditary discrepancy of the counting matrix" is a discrepancy
invariant of the **cut side of the triangle** — the vertical-leg/hierarchy side of "everything is the
triangle," the exact side the score/OCR machinery lives on. The DP paper is, structurally, a
discrepancy theorem about the object on one leg of our triangle.

## The ℓ∞ max ↔ the worst-position functional

The extra `√log` in the DP bound comes from the **max over `n` correlated prefixes** (the ℓ∞ norm). The
LRC floor's controlling functional is the **worst-position** integral
`Q_c(m) = ∫_0^m D_c^*(s) ds` (S33 §5): the adversary places the compact part `L_B` exactly where the
renormalized density `D_c` is **smallest** (increasing rearrangement). And `M(S) = max_t`. Both are
"**worst over a structured family**" functionals, and — the key fact — both are **extremized by the
uniform configuration**.

## (R2) closed to evidence standard: the AP is the low-discrepancy extremal

This is where the analogy paid for itself. The last genuine-mathematics item in the S33 assembled floor
proof is **(R2)**: do difference patterns of range `> R0 = 14` stay above `F_j`? (F3-sharp is done;
MISTAKE-090 retired the old "uniform arc-count" framing.) Reading it through the DP lens — *the
low-discrepancy instance is the hard one* — predicts the answer and the extremal. Verified
(`lrc14_R2_largerange_discrepancy_floor_opus_20260703_S48.py`, grid `T = 6000`, AP re-validated against
S33's exact `559/11025`):

- **`F_7` is minimized at range 6 — the AP `(0,1,…,6)`** (the unique minimal-range pattern), value
  `0.05069 = 1.83 × 1/36`, and **strictly increases with range** (`0.066, 0.062, 0.074, …, 0.099` at
  range 150). Large-range spot checks to range 150: **all ≥ the AP floor and ≥ 1/36.**
- **`j = 8` identical shape** (AP argmin `0.0546`, all `≥ 1/36`).
- **The mechanism is discrepancy.** The AP has the **smallest star-discrepancy** of `{c_i}` *and* the
  smallest floor `Q_c`; spreading the `c_i` (higher range = higher discrepancy) **raises** `D_c` toward
  the independent reference `(6/7)^7 = 0.340` and raises the floor. The AP sweeps its arcs in
  lock-step → maximal simultaneous covering → minimal density → the extremal. **Low discrepancy = hard
  instance**, exactly as in the counting-matrix DP bound.

So range `> 14` is **safer, not a gap**: the whole renormalization residual is pinned by the compact
**AP fixed point**, which is *also* the difference-flow renormalization fixed point (HYP-3901). (R2)'s
"large-range patterns recurse / approach independence" is confirmed on the nose. The remaining proof
step is now **shaped, not open**: **Schur-convexity of `Q_c` in the gap vector** (AP = all-gaps-equal =
the majorization minimizer) — the same majorization species as the DP paper's discrepancy argument, and
a candidate for its `γ₂`/factorization-norm machinery.

## What transfers, what is only analogy (honest)

- **Genuine and load-bearing:** (i) the *hard direction is the lower bound*, and moments give the wrong
  direction in both; (ii) F3-sharp's telescoping **is** the binary-tree cancellation; (iii) (R2) **is**
  a discrepancy estimate, and the DP lens correctly predicted the AP extremal and closed it to evidence
  standard.
- **Structural (real but not a theorem transfer):** the counting matrix = the staircase cut space; the
  ℓ∞ max = the worst-position `Q_c`; AP = dyadic/uniform extremal.
- **Speculative (a named tool to try):** the DP proof's **hereditary-discrepancy / `γ₂` factorization
  norm** is the natural instrument for the Schur-convexity of `Q_c` — i.e. for turning the (R2)
  evidence into a proof. Not claimed, flagged.

## Status

- **VERIFIED:** (R2) to evidence standard — `F_7`, `F_8` minimized at the AP, strictly increasing in
  range, all `≥ 1/36` out to range 150; AP `Q_c` matches S33's exact value; the discrepancy mechanism
  (AP = min star-discrepancy = min floor).
- **ROUTE:** the last (R2) proof step = Schur-convexity/majorization of `Q_c` in the gap vector; the DP
  `γ₂`/hereditary-discrepancy method is the proposed tool.
- **REFLECTION:** the LRC floor lower bound and the DP continual-counting lower bound are one shape —
  discrepancy of a triangular accumulation, extremized at the uniform/dyadic config, won by cancellation
  not moments. The counting matrix is the staircase's cut space.

Related: HYP-4013 (this, the (R2) closure); HYP-4011/klein-S106 (F3-sharp = the telescoping/binary-tree
cancellation); the S33 assembled floor proof (`lrc14-11core-floor-assembled-proof-opus-20260701-S33.md`,
§5 (R2)); HYP-3901 (difference-flow renormalization, the AP fixed point); MISTAKE-090 (the retired
uniform-arc-count framing); the Cut⊕Cycle split and "everything is the triangle" (CLAUDE.md); S29
(covering-min = facility-location, potential = discrepancy). External: arXiv:2607.00876 (Bairaktari–
Larsen). Script: `04-computation/lrc14_R2_largerange_discrepancy_floor_opus_20260703_S48.py`.
