# `offdiag ≤ 0` is refuted — the density `Q_s=O(r)` is a convex-minus-lattice balance, the genuine analytic bottom of the density route

*klein-2026-07-13-S283. Owner: prove the Weyl cancellation for the arc midpoints. Reframing it as
`offdiag ≤ 0` (arcs anti-correlate under `×w`, which would give `Q_s ≤ diag = O(r)` rigorously) yields a
clean mechanism — `B₂` convexity — but the clean sign is **false**: lattice-straddle corrections can
overwhelm the convex anti-correlation. Every elementary/structural route (S281–S283) is now exhausted;
the density `Q_s=O(r)` is a genuine equidistribution cancellation. This is the honest bottom.*

---

## The reframing and its mechanism

`Q_s = diag + offdiag` with `diag = Σ_i 4π²\{w w_i\}(1-\{w w_i\}) ≤ 4π²wμ = O(r)` a **closed-form,
rigorous** `O(r)` quantity (`w_i` = arc widths). So `offdiag ≤ 0 ⟹ Q_s ≤ diag = O(r)`, closing the row.

And there is a real reason to hope: each pairwise `offdiag_{ij}` is the mixed 2nd difference of
`P(θ)=2π²B₂(\{θ\})` over the arc endpoints, and `B₂` is **convex** (`B₂''=2>0`). The smooth part of a
mixed 2nd difference of a convex function is `≈ −P''·w_iw_j·w² < 0` — the arcs **anti-correlate by
convexity**. The only positive contributions are integer-straddles (the Dirac comb in `P''=4π²−4π²Ш`),
i.e. lattice-line overlaps `|arc_i ∩ (arc_j − m/w)|`.

## But the clean sign is false

Direct computation (8 cluster/sector cases): `offdiag ≤ 0` holds in **7**, but **fails** for the spread
cluster `{0,10,27,55,99,150,199}` at `s=3`: `Q_s=433 > diag=351`, `offdiag=+82`. So the lattice-straddle
(positive) part **can overwhelm** the convex (negative) part. `Q_s ≤ diag` is not a theorem.

The silver lining is thin: even in the violating case `Q_s=7.7r` is still `O(r)` (`|offdiag| ≤ diag`),
so `Q_s=O(r)` survives — but only via the *full* convex-minus-lattice balance, not a clean inequality.
"`|offdiag| ≤ diag`" is itself just `Q_s=O(r)` restated — circular.

## The exhaustion, stated plainly

Across S281–S283, every elementary/structural attempt on the density `Q_s=O(r)` (equivalently, the
Weyl cancellation `Σ_i w_i e(-ℓw c_i)=o(\text{trivial})`) has been tried and falls short:

| route | outcome |
|---|---|
| crude Fourier `|f̂|≤V/(2π|n|)` | `Q_s ≤ 4π²r²/3 = O(r²)` — rigorous, insufficient |
| large sieve (min-separation) | `O(r³)` — worse (endpoint clustering `δ∼1/r²`) |
| Montgomery–Vaughan, width-weighted | `Σw_i²/δ_i` bounded (clustering tamed) but `Q_s=O(w²)=O(r²)` — no decay |
| diagonal `Q_s ≤ diag` (`offdiag ≤ 0`) | **refuted** (this session) |
| `B₂`-convexity anti-correlation | real, but lattice-straddles overwhelm it |

All soft/mean-value/structural tools cap at `O(r²)` or are refuted. The `O(r^{2-ε})` saving (which is
*all* the density row needs, S281) requires the genuine oscillatory cancellation of the convex part
against the lattice part — an equidistribution estimate, not an elementary inequality.

## Honest bottom line

The density route is **fully reduced** and mostly rigorous: the transfer (THM-710), the endpoint
Fourier reduction (THM-727), the 1-D DFT of the derivative (THM-728), the autocorrelation-discrepancy
identity (THM-729), the crude `Q_s ≤ 4π²r²/3`, the *any-power-saving-suffices* downgrade, and the
closed-form `diag = O(r)`. The sharp `Q_s=O(r)` is confirmed empirically (S280). The **single remaining
piece** is one genuine equidistribution cancellation (convex-minus-lattice / Weyl for the arc midpoints
under `×w`) — soft (any `ε>0`), 1-linear, but a *real analytic estimate*, not a further one-session
reduction. It is the honest analytic bottom of the density side — lower-order than, and cleanly
separated from, the covering side's sharp multi-linear inequality (mac-mini-S78).

The productive next moves are: (a) a sustained analytic effort or external equidistribution input on
this one estimate; or (b) accept the reduced state and consolidate the LRC(14) finish-map (both routes'
irreducible cores are now pinned). Continuing to peel elementary reformulations has reached its end —
this is genuine analysis.

*Files: `04-computation/lrc14_offdiag_sign_klein_S283.py` (+ out). HYP-6445. Refutes the `offdiag ≤ 0`
target; consumes THM-729/728/727. Closes the S281–S283 elementary-attempt arc.*
