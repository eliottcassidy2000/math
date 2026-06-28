# Two new angles on the k=8 node: Carathéodory–Toeplitz moment rigidity and ferromagnetic Griffiths

*kind-pasteur-2026-06-27-S31al. The owner asked for a creative new angle or two on a remaining LRC(14)
proof target. The single remaining node is the k=8 bounded-core extremality — "consec maximizes the
coverage/covariance" (HYP-3160/3161, exhaustively verified for `Σκ₂`). I developed two genuinely new
attack routes, each with a viability check. Angle 2 lands a clean exhaustively-verified result.*

## The target
Prove **consec = {0,…,7} maximizes the empty-sector coverage** at k=8 (equivalently the covariance
`Σκ₂`, the bimodality `L_yK8`, the cover bound `q₀ ≤ cap`). The even/degree-2 part is numerically
exhaustive; a *symbolic* proof is the gap. Two new lenses:

## ANGLE 2 (the strong one) — Carathéodory–Toeplitz moment rigidity
**Idea.** My Lee-Yang result (HYP-3099) — the coverage PGF `Q(z)=Σ q_t z^t` has zeros on a circle — *is*
the classical **trigonometric moment problem**: treat `{q_t}` as the moment/Fourier sequence of a positive
measure on the circle. The Carathéodory–Toeplitz theorem governs when such a sequence is a valid moment
sequence (the Hermitian **Toeplitz matrix `T[j,k]=q_{|j−k|}` is PSD**) and characterizes the extremal
(boundary) measures as **atomic on the circle** (the atoms = the de Moivre angles, S31ak).

**Viability — VERIFIED EXHAUSTIVELY.** consec **maximizes the minimum eigenvalue `λ_min(T)`** of the
coverage-moment Toeplitz matrix over **all 3432 bounded k=8 clusters (0 beaters)**, `λ_min(consec)=0.0423`.
So:
> **consec is the most "interior" / least-rigid moment configuration — it maximizes the Carathéodory PSD
> margin `λ_min(T)`.** Dissociated configs sit nearer the boundary of the moment cone (smaller `λ_min`,
> closer to atomic/degenerate). This is a NEW, clean, exhaustively-verified extremal characterization of
> consec — the *positive-definite* face of the bimodality/ferromagnetic story.

**Why it's a proof route.** The Toeplitz/moment world has heavy classical artillery the moment-LP lacks:
the **Szegő theory** (Toeplitz determinant `D_n = ∏λ ~ exp(n·∫log w)` ties the spectrum to the measure),
the **Schur/Carathéodory–Fejér algorithm** (the zeros-on-circle / Schur-Cohn rigidity), and **Verblunsky
coefficients** (the orthogonal-polynomials-on-the-circle parametrization). The cover bound `q₀ ≤ cap` is a
statement about the *diagonal* of `T` (`q₀ = tr T/7 =` mean eigenvalue), so it is a **spectral inequality
on `T`**: maximizing `λ_min` while fixing the mean (`q₀`) bunches the spectrum — and the cap is the value
where the de Moivre (cyclotomic, 7-atom) extremal measure sits. **Target: bound `q₀` via Szegő's theorem
on `T`, using `λ_min`-maximality (= consec) and the de Moivre atomic extremal (= cap).** This trades the
stalled moment-LP for the Toeplitz spectral toolbox.

## ANGLE 1 (a lead) — ferromagnetic Griffiths/GKS monotonicity
**Idea.** The empty-sector indicators are a *ferromagnetic* spin system for k≥6 (HYP-3161, all 15 `Cov>0`
at consec). Griffiths-II says correlations **increase with couplings**; consec is the maximal-coupling
**ground state**. A Griffiths/coupling proof would show "consec maximizes `q₀`" by a monotone path of
increasing couplings to consec.

**Viability — needs the right partial order (naive fails).** A greedy speed-replacement path
`config → consec` is **non-monotone** in `Σκ₂` (only 3/60 monotone): the path dips through the
**antiferromagnetic** region (negative covariance) before the final jump to consec's max. So the naive
"move speeds toward `{0..7}`" order is wrong — it crosses the disordered phase. The honest read:
> the ferromagnetic ground-state extremality is real (consec = the all-positive-coupling config), but the
> couplings are **not free parameters** (they are determined by the cluster `E`), so plain Griffiths does
> not apply on the speed-lattice. The RIGHT object is a monotone path in the **coupling/coherence
> manifold** (the Fourier/relation-lattice couplings), not the raw speeds. **Sub-target:** find the
> partial order on configs (e.g. by the Fourier coherence `Σ_m |Ê(m)|²`-type functional, or a
> random-current/switching representation of the sector-emptiness measure) under which `Σκ₂` is monotone —
> then Griffiths-II closes it. The random-current method (Aizenman) is the natural tool, since it makes the
> ferromagnetic ground-state argument rigorous without free couplings.

## How they fit the web
Both are the *positive-definite / spin-system* face of the same node:
- **Angle 2 (Toeplitz `λ_min`)** = the moment-cone interior extremality = the Lee-Yang circle's rigidity
  (Carathéodory), with Szegő as the new tool. Exhaustively verified, classical machinery available.
- **Angle 1 (Griffiths)** = the ferromagnetic ground state, needing the right coupling-partial-order
  (random-current). A lead with a clear sub-target.
Both sit *above* mac-mini's biquadratic resolvent (the even/degree-2 algebra) and *complement* the
de Moivre cyclotomic atoms (S31ak): the Toeplitz extremal measure's atoms ARE the de Moivre angles.

## Net
- **NEW + EXHAUSTIVE:** consec maximizes `λ_min(Toeplitz of {q_t})` over all bounded k=8 (Angle 2) — a
  fresh, classical-tool-backed (Szegő/Carathéodory) extremal characterization, the proof route of choice.
- **NEW lead:** ferromagnetic Griffiths/GKS via random-current; the naive speed-path is non-monotone (it
  crosses the antiferro phase), so the partial order must be the coupling/coherence manifold.
- Both convert the stalled moment-LP into **spectral (Toeplitz)** and **statistical-mechanics (random-
  current)** problems with mature machinery.

→ HYP-3201 (this), HYP-3099 (Lee-Yang circle = the moment problem), HYP-3160/3161 (covariance/ferro),
HYP-3162 (de Moivre atoms = Toeplitz extremal), Szegő/Carathéodory–Fejér/Verblunsky, Aizenman random-
current, THM-534 (moment-LP, superseded as the route), OPEN-Q-108.
