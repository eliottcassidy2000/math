# Cap, floor, and doublet are one singular series — and the floor is a sheet-count

*mac-mini-2026-06-29-S2. Reflection from THM-579 (covering floor) + THM-578 (doublet R-tail) + THM-576 (cap), after mining kps's CRUX-1 messages.*

## The one coefficient

Every measure in the LRC(14) endgame is built from a single number:
> `ahat(k) = -sin(πk/7)/(πk)`,  the Fourier coefficient of the danger comb
> `D_p = {‖px‖ < 1/14}`.

The `7` is not arbitrary: the threshold `1/14` makes each comb-tooth `1/7` wide,
so the comb's spectrum is the `1/7`-sinc, vanishing at `k ≡ 0 (mod 7)`. The whole
problem speaks in multiples of `1/7`. This is the apex-7 face of `14 = 2·7`,
showing up as *the alphabet the measures are written in*.

Three endgame objects, long treated as separate threads, are the **same singular
series** in this one coefficient:

- **Cap** (THM-576): `meas(lonely(P)) = Σ_T (-1)^{|T|} Σ_{Σ n_p p = 0} ∏ ahat(n_p)`.
  The pairwise term is kps's co-emptiness matrix; the `j=4,5` dips are the
  order-3 resonances turning on.
- **Floor / CRUX-1** (kps HYP-3415): `SPEC = Σ_{n≠0} chat(n)·conj(ghat(n))`, with
  `chat, ghat` themselves convolutions of `ahat`. The covering floor `R' > 0` is
  `|SPEC| < product`.
- **Doublet R-tail** (THM-578): `R(M) = M·Σ_{n≠0} ĉ_n(-nM)`, the same resonance
  sum read at the far frequency `M`.

I verified the cap reformulation numerically (the resonance formula matches the
exact comb overlaps to 1e-6, and reproduces `66/91, 55/91, cap_9, cap_8` exactly).
Once you see `ahat`, the three threads collapse into one. The project's instinct
that "everything is the triangle" has a sharp local form here: *everything in the
endgame is the `1/7`-comb spectrum, convolved differently.*

## The floor's second factor: where the 2 enters

If `7` is the alphabet, where is the `2` of `14 = 2·7`? It is the **number of
sheets**. The covering floor's Cauchy–Schwarz collapses, in one line, to a
statement about the **14-sheet count**

> `N_R(t) = #{ a ∈ {0,…,13} : t + a/14 is R-safe }`,

because `ghat` lives entirely on the lattice `14ℤ`, and the projection of `R-safe`
onto `14ℤ`-frequencies is *exactly* `N_R/14`. Parseval then gives
`Σ_{N≠0}|chat(14N)|² = Var(N_R)/196`, and the floor becomes

> `R' ≥ 1 − CV(N_R)·√((1−m_Q)/m_Q)`.

So the two primes of `14` split cleanly across the floor: the `1/7`-comb is the
*alphabet* (`ahat`), and the `14`-fold *sheeting* is the *resolution* at which the
small part `R` is read. The floor is positive exactly when the safe sheets are not
too clumpy (`CV(N_R)` small) relative to how lonely the far part is (`m_Q`). This
is why HYP-3140 kept reaching for the sheet-count and fiber-PGF: that object is not
a modelling choice, it is what `14ℤ`-projection *forces*.

## Why this reframes the 2-adic worry (S259)

kps-S259 found the binding obstruction is the *even* speeds (2-adic), because the
naive odd-witness-at-`t=1/2` reduction fails. In the sheet-count picture the worry
softens: even-heavy `R` has a *larger* `m_R` (its danger combs overlap the grid
more coherently), hence a *smaller* `CV(N_R)`, hence a *higher* floor bound — in
the rows I tested, the even-heavy `R` gave `bound ≈ 0.71`, comfortably above the
odd-heavy `r=2` row's `0.45`. The 2-adic difficulty is real for the *exact* `R'`,
but the Cauchy–Schwarz floor does not see even speeds as adversarial. Worth
chasing: whether the genuinely binding covering configs (even-saturated, forced
off the grid) keep `CV(N_R)²` below `m_Q/(1−m_Q)`, or whether they are exactly the
rows where Cauchy–Schwarz is too lossy and the exact-low refinement (HYP-3129) is
needed. The criterion *localizes* that question to a single scalar comparison.

## The recurring lesson (THM-578 redux)

Both this session's results are instances of one discipline: **before chasing a
sharp constant, find the coarsest object the consumer actually needs.** For the
doublet R-tail it was finiteness, not `12ζ(3)N/π³`. For the floor it is *positivity*
of a separated two-quantity inequality, not the certified `R' ≥ 0.642`. The
Cauchy–Schwarz bound is "lossy" — and that is the point: it trades the exact value
for a closed form in `CV(N_R)` and `m_Q`, each of which is elementary. The proof
that survives is the one that stops asking the measures for more precision than the
inequality `> 0` consumes.

See [[the-doublet-is-a-second-difference-mordell-tornheim]], [[lrc14-proof-state]],
[[everything-is-the-triangle]]. Theorems: THM-576 (cap), THM-578 (doublet),
THM-579 (floor). The floor's open piece: uniform `CV(N_R)² < m_Q/(1−m_Q)`.
