# The wall was a decorrelation — and both sides have zero mean

*kind-pasteur-2026-07-11-S127. Owner: "attack the S_c reciprocal-sum bound, the open wall." I attacked it,
and it turned out the wall is not a lattice-geometry problem at all — it is a one-line decorrelation, and the
piece I proved last turn (THM-699) is its mirror image. This note is that attack and the symmetry it
uncovered.*

---

## Two faces of one wall

The support-6 reciprocal-sum bound — the project's single named open lemma for the LRC(14) wide-spread half —
had been carrying a reputation for hardness earned honestly: the lattice sum `corr(E) = Σ_c D7(c) S_c(E)` is
only conditionally convergent (MISTAKE-078), the absolute Minkowski envelope diverges harmonically, and every
successive-minima attempt overshoots 5×–32×. That is all true of the *lattice-sum* face of the wall.

But the wall has a second face, and HYP-2644 had already turned toward it: the **plateau recursion**. The
whole unbounded direction reduces to a single estimate — *the largest offset decorrelates* — and that
estimate has comfortable margin (0.13–0.18, an order looser than the tight finite check). HYP-2644 verified
it to `O(1/w)` and named it the open piece. I proved it.

## The one-line decomposition, and the Fourier bound

`p_0(E) = meas(S7(E))`. For `E = E' ∪ {w}`, `w = max E`, the cover is *exactly* decomposable:
`E ∪ {w}` covers all seven sectors iff `E'` covers all seven, or `E'` misses *exactly one* sector `s` and the
single point `frac(wx)` lands in it. Two more sectors missed and one point cannot fill them; so, pointwise,

`1_cover(E,x) = 1_cover(E',x) + Σ_s f_s(x)·1{frac(wx) ∈ sector s}`,  `f_s = 1{E' misses exactly s}`.

Integrate, subtract the plateau value `Plat(E') = p_0(E') + (1/7)p_1(E')`, and the error is *centered*:

`p_0(E) − Plat(E') = Σ_s ∫_0^1 f_s(x)·[1{frac(wx)∈s} − 1/7] dx`.

Now it is elementary. Each `f_s` is bounded-variation (`|\hat f_s(m)| ≤ V/2π|m|`); each bracket `g_s` is
mean-zero with `|\hat g_s(ℓ)| = |\sin(πℓ/7)|/π|ℓ|`. Orthogonality forces `m + ℓw = 0`, the mean-zero fact
kills `ℓ=0`, and what survives is `Σ_{ℓ≠0} \hat f_s(−ℓw)\hat g_s(ℓ)`, bounded termwise by
`V(f_s)/(2π|ℓw|)·1/(π|ℓ|)`, summing to `V(f_s)/(6w)`. Hence

> **`|p_0(E) − Plat(E')| ≤ V(E')/(6w)`** (THM-700).

No lattice, no conditional convergence, no successive minima — a bounded-variation function correlated
against a single fast Weyl frequency, which by Fourier decay is `O(1/w)`. The wall's fearsome face was the
lattice sum; its tractable face is the same object read as an integral, where the divergence was never real —
it was an artifact of resumming a convergent integral into a conditionally-convergent series.

## The symmetry: both sides have zero mean

The thing that made THM-700 fall is exactly the thing I proved last turn on the *other* side of the same
correction. THM-699 showed the weight has **zero mean**: `Σ_c D7(c) = 0`, so `corr = Σ_c D7(c)(S_c − S̄)`
sees only the coset-*variation* of the reciprocal sums. THM-700's oscillation `g_s = 1{·} − 1/7` also has
**zero mean**, which is exactly what deletes the `ℓ=0` term and leaves an `O(1/w)` tail. The correction
`Σ_c D7(c) S_c(E)` is a pairing of two mean-zero objects — a centered weight against a centered oscillation —
and the far element decorrelates them. Neither zero-mean fact is a coincidence: they are the same cancellation
(the `(−1)^{|T|}` seven-sector alternation) seen once on the residue-character side and once on the
frequency-oscillation side. The wall is a *correlation of two centered signals*, and a correlation of centered
signals across a spectral gap is small. That is the whole content.

## What is closed, and what is not

Closed: the inductive step — the decorrelation `|p_0 − Plat| ≤ C/w`, rigorous, explicit constant, the single
open analytic lemma HYP-2644 named. Not closed: the *recursion*. Three mechanical-but-unwritten pieces remain
— the identical decorrelation for `p_1` (the `Q(m) = max[p_0 + p_1/7]` object), the accumulation of per-peel
errors down to a bounded core (a summable-geometric bookkeeping, since each peeled `w` dominates its residual
frequency), and the sharp constant (the crude `V ≤ 14Σe'` overshoots the true `≈0.2` by ~10², recoverable
from the `|\sin(πℓ/7)|` numerator and cancellation among the seven `f_s`). And the tight margin was never in
this direction at all — it lives in the finite `consec_k` check, where it has always lived.

The lesson: a wall that resists a hard tool is sometimes waiting for a soft one. The reciprocal-sum bound
did not need geometry of numbers; it needed the observation that a convergent integral had been forced into
a divergent series, and that once written as the integral it is an undergraduate Fourier estimate. Attack the
representation, not the reputation.

*Files: `04-computation/lrc14_plateau_decorrelation_kps_S127.py` (+`.out`), canon `THM-700`. Continues the
support-6 thread from [[the-support6-kernel-is-a-permanent-and-it-vanishes-on-average-kps-S127]] — THM-699 is
the weight side, THM-700 the oscillation side, both zero-mean. The remaining recursion closure is downstream
of THM-700; HYP-2644's `Q(k−1) < cap_k` finite margins (0.13–0.18) are the room it has to work in.*
