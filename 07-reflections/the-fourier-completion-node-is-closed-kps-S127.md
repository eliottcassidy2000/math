# The LEM-022 Fourier-completion node is closed

*kind-pasteur-2026-07-11-S127. Owner: "prove the Parseval completion identity, build the ZMod q Fourier
infra." opus-S214 had just proved the identity concurrently — so I closed the node the other way: the final
assembly. This note records the closure of what the mapping pass named the single biggest open measure-side
sub-lemma.*

---

## A concurrent landing, and the right pivot

I was sent to prove the Parseval completion identity `C_w = b²/q + (1/q)·Σ_{k≠0} B̂(k)·conj(B̂(w·k))`. The
first thing I did was pull — and found opus-S214 had proved *exactly* that (`completion_identity`) plus the
`k = 0` split (`completion_diff_bound`: `‖C_w − b²/q‖ ≤ (1/q)·Σ_{k≠0} ‖B̂(k)‖·‖B̂(w·k)‖`), concurrent with my
own cont.18 hand-off of the derivation. The owner dispatches the same prompt to several machines; this is
that, and the discipline is: **do not re-prove what just landed.** opus built the `ZMod`-adjacent `eInt`
character infra and the identity. So the identity was done — and the node's *remaining* open piece was the
other half: the final assembly opus stopped short of.

## What was left, and closing it

opus's `completion_diff_bound` ends at `‖C_w − b²/q‖ ≤ (1/q)·Σ_{k≠0} ‖B̂(k)‖·‖B̂(w·k)‖`. Three things turn
that into LEM-022's target `5q(log₂q+1)²/P`:

1. **The termwise coefficient bound** — opus's B.2 `norm_bandCoeff_le`: `‖B̂(k)‖ ≤ q/(2·cdist k)`.
2. **The `range q → ZMod q` transport** — the diff-bound sum is over integer residues; death-star's harmonic
   sum lives over `ZMod q`.
3. **The harmonic aggregation** — death-star's `harmonic_ratio_sum_mul_le`: `S·P ≤ 20(log₂q+1)²`.

I proved the assembly `completion_closed_of_coeffBound` doing all three, kernel-pure. Two bridges carried it:

- **`P > 0` supplies the unit.** B.2 needs `w·k ≠ 0` for the second factor. Rather than assume `w` is a
  unit, the ratio floor does it: if `w·z = 0` for some `z ≠ 0`, then `P ≤ cdist z · cdist 0 = 0` — so
  `P > 0` *forces* `w·z ≠ 0`. The hypothesis you already have hands you the one you need.
- **`natCast`/`val` is the transport.** `Finset.sum_bij'` between `(range q).filter(≠0)` and
  `univ.filter(≠0)`, `k ↦ (k : ZMod q)`, `z ↦ z.val`. The reciprocal summand is invariant, and the
  constants close on the nose: `(q²/4)·(20/P) = 5q²/P`, then `/q` from the identity gives `5q/P`.

Then I discharged the coefficient hypothesis for the actual band shape. The LRC safe set is a *contiguous
arc*, and for `B = {lo, …, lo+len−1}`, `bandDFT q B j` is precisely opus's B.2 interval coefficient at
frequency `−j` — so `norm_bandDFT_interval` gives `‖B̂(j)‖ ≤ q/(2·cdist j)` via `norm_bandCoeff_le` and
`cdist(−j) = cdist j`. Composed:

```lean
completion_closed_interval : ‖C_w − b²/q‖ ≤ 5·q·(log₂q+1)² / P      -- unconditional, interval band
```

kernel-pure, no hypotheses beyond `P > 0` and the ratio floor. **The node is closed.**

## What closing it means

The mapping pass was explicit: this was *the* sole unproven link on the Fourier / OffLine≤f(E3) route, and
"proving it lights up the whole chain." Every downstream consumer was already in Lean — the per-cell `M` it
produces feeds `MultCorrelation.offdiag_mcorr_sq_le`, which with `E3Budget.E3_lt_choose` gives the
OffLine ≤ f(E3) energy bound and the density floor. So a full analytic route to the residual measure floor
now has no gap at this node.

It took four hands and one clean division of labor: death-star's combinatorial `hyperbola_box_count` and
dyadic `harmonic_ratio_sum_mul_le` (S9–S13), opus's orthogonality, coefficient bound, and Parseval identity
(S213–S214), my off-diagonal aggregation and the closing assembly (cont.18–19). Nobody built the whole
thing; each brick was designed to compose, and they did. That is the shape of a frontier falling — not one
heroic proof, but a stack of kernel-pure pieces that meet at seams someone thought about in advance.

## The lesson

When you're sent to prove `X` and pull to find `X` just landed, the move is not to duplicate and not to
stop — it is to ask *what `X` was for*, and prove the next thing. The identity was a means to the node; the
node was still open at its other half; the honest continuation was the assembly, and it composed with the
concurrent work instead of colliding with it. Pull often, and when you collide, pivot forward.

*Files: `LRCFourierClosed.lean` (`completion_closed_of_coeffBound`, `norm_bandDFT_interval`,
`completion_closed_interval`), building on opus's `LRCFourierCompletionC` (identity + diff bound),
`LRCFourierCompletionB` (B.2), and death-star's `LRCHyperbolaBox` (harmonic aggregation). Continues
[[the-fourier-completion-reduces-to-one-parseval-identity-kps-S127]] — the "one Parseval identity" it named
as remaining is now proved (opus) and the node assembled shut (this).*
