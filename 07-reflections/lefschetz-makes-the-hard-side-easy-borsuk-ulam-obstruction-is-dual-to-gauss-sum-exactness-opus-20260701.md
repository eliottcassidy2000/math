# Lefschetz makes the hard side easy: the Borsuk–Ulam obstruction is dual to Gauss-sum (Weil/Lefschetz) exactness — certify the free-ℤ₂ regime with a trace formula, not a fixed point

*opus-2026-07-01-S23. The owner asked for tricks to make the hard (free-ℤ₂, p≡3 mod 4, Borsuk–Ulam) side easier,
thinking Lefschetz fixed-point counts. The trick is real and clean: where Brouwer (existence of a symmetric
fixed point) fails, the **Lefschetz count** `L(f)=Σ(−1)ⁱ Tr(f_*|Hᵢ)` still works for the right map — and for the
p≡3 mod 4 / Paley structure that count is the **Gauss sum √p**. The topological hardness is *dual* to arithmetic
exactness.*

## The move: Brouwer counts existence, Lefschetz counts fixed points
Brouwer says a continuous self-map of a convex body has *a* fixed point — but a **free** involution has none, so
Brouwer/symmetric-SOS gives nothing (the Borsuk–Ulam wall, HYP-3814). Lefschetz is stronger: for *any*
endomorphism `f` of a compact space, `L(f)=Σ(−1)ⁱ Tr(f_*|Hᵢ)` equals the **signed count of fixed points**, and it
can be nonzero even when the involution itself is free. So the trick is: **don't certify the hard side by finding
a symmetric fixed point; certify it by a Lefschetz trace of the Frobenius / dilation.** Two computations make
this concrete.

### (A) The Frobenius trace = the Gauss sum √p (verified)
The Paley tournament's skew matrix `S_ij = χ(j−i)` (χ = Legendre) is skew for p≡3 mod 4, and its spectrum is
**exactly `{0, ±i√p}`** — verified p=3,7,11,19,23. Those `±√p` are the **Gauss sums**, i.e. the Frobenius
eigenvalues on `H¹` of the associated variety (Grothendieck–Lefschetz: `#Fix(Frobᵏ)=Σ(−1)ⁱ Tr(Frobᵏ|Hⁱ)`, the
Weil point count). So the Paley/QR structure — the flip-rank obstruction (HYP-3805) and the free-ℤ₂ hard core —
has a spectrum that is an **exact Lefschetz/Weil trace**, not a topological unknown. Its Cayley spectrum is
**concentrated** at just `{1, e^{±2i·arctan√p}}` (the fixed vertex + the Gauss-sum pair), *not* spread like the
transitive tournament's roots of unity — so the "hard" object is arithmetically *simple*, just non-symmetric.

### (B) The dilation trace = the resonance count 1−v (verified)
The LRC runners are the dilations `φ_v : t ↦ vt` on the circle `S¹`. As a degree-`v` map, `L(φ_v) = Tr(H₀) −
Tr(H₁) = 1 − v = −(v−1)` = the signed count of its fixed points `t = k/(v−1)` — the **resonances** (times when
runner `v` returns to the origin). So the entire coincidence/loneliness arrangement is a **Lefschetz count**:
the lonely set is the complement of `⋃_v Fix(φ_v)`, and the three-distance / Farey gap structure that governs
`M(S)` for the AP is exactly the combinatorics of these fixed points.

## The duality: Borsuk–Ulam hardness ⟺ Weil/Lefschetz exactness
The punchline is that the *same* arithmetic switch controls both sides:
- `p ≡ 3 (mod 4) ⟺ −1 ∈ QNR ⟺ the QR-negation ℤ₂ is **free**` — this is the **Borsuk–Ulam** obstruction
  (no symmetric SOS certificate; HYP-3814).
- `p ≡ 3 (mod 4)` is *also* exactly the condition for the **Paley tournament to exist** and for its point
  count / spectrum to be the **exact Gauss sum √p** (Weil/Lefschetz).

So the topological wall and the arithmetic key are the *same door*. **The free-ℤ₂ obstruction that kills the
symmetric certificate is precisely what makes the trace formula exact.** The trick to make the hard side easy:
stop looking for a symmetric fixed point and **read off the Lefschetz/Weil trace** — the hardness is the
solvability, seen from the other side.

## The concrete payoff for the LRC endgame
This is not just an analogy; the LRC crux already *is* a trace formula:
- The AP lonely-measure moments are the **Ramanujan sums** `c_N(k)` (klein HYP-3793) — and Ramanujan sums are
  finite character-sum **Lefschetz traces** (fixed points of the dilation/Galois action on `(ℤ/N)*`).
- The **singular series** — the "achievability" crux where my MOSEK/SOS campaign stalled (HYP-3791) — is a sum
  of these Ramanujan/Gauss traces. It has *no* symmetric-SOS certificate (Borsuk–Ulam) but *is* an exact
  Lefschetz trace formula. That is why the far-element/equidistribution route (klein/mac-mini) works where SOS
  doesn't: it is computing the trace, not searching for a symmetric optimum.
- The covering-min binds at the composite `Φ₆ = n²−n+1` (metric-irreducible, klein HYP-3812). "Metric-irreducible
  at a composite" is the statement that the trace does **not** factor over the CRT — i.e. it is a genuine
  Lefschetz trace on the full `Φ₆`, not a product of prime-level traces (THM-503: `L` is not an Euler product).

So the endgame recipe on the hard side: **replace the (obstructed) symmetric SOS with the exact Lefschetz/Weil
trace** — Gauss sums for the Paley obstruction, Ramanujan sums for the LRC singular series, resonance counts
`1−v` for the loneliness arrangement. The three pillars (POCS / flat-extension / Blaschke, HYP-3796) are the
*constructive* face of the same trace: flat-extension **is** the atomic (fixed-point) decomposition the
Lefschetz count enumerates; the Blaschke circle-map fixed points **are** the `Fix(φ_v)`; POCS converges to the
trace's support.

## Status
- **Verified:** Paley skew spectrum `= {0, ±i√p}` = the Gauss sum (p=3,7,11,19,23); Paley Cayley spectrum
  concentrated at the Gauss-sum pair + fixed vertex; dilation `L(φ_v)=1−v` = resonance count.
- **Trick (organizing, with a verified backbone):** on the free-ℤ₂ / p≡3 mod 4 hard side, certify via the
  **Lefschetz/Weil trace** (Gauss/Ramanujan sums) not a Brouwer fixed point; **Borsuk–Ulam obstruction is dual
  to Gauss-sum exactness** (same p mod 4 switch). The three pillars are the constructive face of the trace.
- **Payoff:** the LRC singular series (SOS-obstructed, HYP-3791) and the covering-min at Φ₆ (metric-irreducible,
  HYP-3812) are exact **trace formulas** — this is *why* the arithmetic (far-element/Ramanujan) route succeeds
  where SOS stalls, and it says the endgame is a trace computation, not a symmetric certificate.
- **Honest:** the Gauss-sum/dilation traces are verified and classical (conference matrices, Weil); "the LRC
  singular series = an exact Lefschetz trace that certifies the bound" is the organizing conjecture — the trace
  is exact, but turning it into the `M(S)≥1/n` bound is still the far-element analysis, not yet closed.

Related: HYP-3814 (Brouwer/Borsuk–Ulam — this supplies the hard-side tool), HYP-3805 (Paley flip-rank
obstruction = the √p Gauss-sum spectrum), HYP-3791 (MOSEK: no SOS shortcut = the Borsuk–Ulam wall), HYP-3793/klein
(moments = Ramanujan sums = Lefschetz traces), HYP-3796/mac-mini (three pillars = constructive trace),
HYP-3812/klein (metric-irreducible at Φ₆ = non-Euler-product trace), THM-503. HYP-3815 (this). Script:
04-computation/lefschetz_paley_gauss_dilation_opus_20260701.py.
