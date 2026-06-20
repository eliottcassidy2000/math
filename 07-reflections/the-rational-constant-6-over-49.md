# The rational constant 6/49: when the apex prime's reality gives back the sign

**Source:** mac-mini-2026-06-20-S2. Dispatch: think about a signed Erdős–Turán packet estimate,
leverage Ruzsa modeling / Plünnecke–Ruzsa, and the apex-prime / Bertrand / Chebyshev structure;
aim for an LRC(14) proof. Canon: THM-546 (sharpened), THM-547, HYP-2676. Built on HYP-2657 (the
QR reality of the coset kernel) and codex's HYP-2675 (the collar/true-wide split).

## A bound that threw away its own sign

THM-546 bounded the far-element deviation by `|Δ_w| ≤ κV(E')/(π²w)`, `κ = 1.857…`, an
irrational constant gotten by putting absolute values on every Fourier mode. It was rigorous
but `5–76×` loose — and it had to be, because absolute values discard exactly the cancellation
that makes the quantity small. The question the loose constant raised: where did the sign go,
and can the arithmetic give it back?

## The apex prime, twice

The sign was hiding in a fact already proved on a different route. HYP-2657 had shown the
correction is *real*: because `6 ≡ −1 (mod 7)` is a **non-residue**, the coset kernel satisfies
`D7(−c) = conj(D7(c))`, so the `n ↔ −n` pairing cancels all imaginary parts exactly. Read back
into the deviation, this says `Δ_w = Σ_j Σ_{n≥1} 2 Re[ŝ_j(n) 1̂_{B_j}(−nw)]` — manifestly real,
the phases organized rather than discarded.

Summing the `n`-series first (Abel summation) collapses it to a sum over the *arcs* of `B_j`:
```
Δ_w = −(1/w) Σ_j Σ_{arcs (a,b)} [ F_j(wa) − F_j(wb) ],
```
where `F_j` is the centered antiderivative of `1_{sector j} − 1/7` — a sawtooth that climbs at
slope `+6/7` across the one sector and drifts down at `−1/7` across the other six. Its extreme
value is not some transcendental number; it is `sup|F_j| = 3/49`, exactly, the same for every
sector. So
```
|Δ_w| ≤ (6/49) · V(E') / w,      6/49 = (p−1)/p²  at  p = 7.
```
The constant is **rational and arithmetic**. The apex prime that had appeared as the
`7`-vanishing of the kernel `ŝ_j(7m) = 0` now reappears in the *size* of the bound: `p² = 49`
in the denominator, `p − 1 = 6` non-apex sectors in the numerator. The same prime governs both
which frequencies survive and how large the survivor can be. The irrational `κ/π² ≈ 0.188`
collapses to `6/49 ≈ 0.122` — sharper, and finally *legible*.

## What a finite maximum buys

The sharpening pays off immediately where the base is bounded. For a cluster whose
second-largest element is `≤ 14` (codex's "boundary collar"), the part below the single far
runner lives in `{0,…,14}`, so its coverage `Plat(E')` and its arc-complexity `V(E')` are each
bounded by a *finite maximum* — a number you compute once. The deviation is then
`≤ (6/49) V_max / w`, which drops below the recursion margin the moment `w` exceeds an explicit
cutoff `w* = 54, 90, 103` for the three dangerous rows. Above `w*`, THM-546 certifies; below it,
a finite check. An unbounded direction became a finite computation, purely because the constant
was concrete enough to divide by. Two of the three regions of the crux — the finite half
(`max ≤ 14`, kps) and the collar (`2nd-largest ≤ 14`, here) — are now closed.

## The lesson

Looseness in a bound is often a missing symmetry, not a missing technique. The factor we were
losing was a *sign*, and the sign was guaranteed by an arithmetic fact — that `−1` is a
non-residue modulo the apex prime — that had been proved for an entirely different purpose.
When the estimate finally carried its sign, the constant it produced was not analytic debris
but `(p−1)/p²`: the prime itself, stating how much a single far runner can distort the
seven-sector picture. The remaining region, the genuinely two-scale "true-wide" clusters, is
where Ruzsa modeling and Plünnecke–Ruzsa enter — small doubling forces a low-dimensional
generalized progression, where scale-invariance or the Freiman-dimension penalty takes over.
That is the last room in the house, and its door is additive-combinatorial, not Fourier.
