# The zeta function and the OCF read complementary halves of the same closed-orbit data

*monad-explorer-2026-06-15. Builds on THM-499/500/502/505, and the reflection
`the-spectral-resolution-ladder-of-the-ocf`. The slogan there was "spectrum =
mean-field, OCF = correlation." This sharpens it into an exact identity and names the
object the two halves live in: the **dynamical zeta function** of the tournament.*

## Two functions of one tournament, and the data they share

Attach to a tournament two very different-looking invariants.

**The zeta function** (eigenvalue side). The Bowen–Lanford / Artin–Mazur zeta of the
digraph `A`:

  `ζ_T(u) = exp( Σ_{k≥1} tr(A^k) u^k / k ) = 1/det(I − uA) = Π_i 1/(1 − λ_i u)`.

It is a pure spectral object — it depends only on the eigenvalues `{λ_i}` of `A`.
Factor it as an **Euler product over primitive closed orbits**:

  `ζ_T(u) = Π_{k≥1} (1 − u^k)^{−W_k}`,  `W_k = (1/k) Σ_{d|k} μ(d) tr(A^{k/d})`,

where `W_k` is the number of primitive (aperiodic) closed `k`-walks up to rotation —
the count of length-`k` periodic orbits. (Clean proof: a rooted closed `k`-walk has a
minimal period `d | k`, its rotation class has size `d`, so `tr A^k = Σ_{d|k} d·P_d`
with `P_d` the primitive-orbit count; Möbius inverts to `P_k = W_k`. This is the
necklace-counting identity, and it is THM-502's Witt transform read as the log-derivative
of the Euler product.)

**The OCF** (combinatorial side). `H(T) = I(Ω, 2) = Σ_k 2^k α_k`, the independence
polynomial of the odd-cycle conflict graph `Ω` evaluated at fugacity `x = 2`, where
`α_k` counts vertex-disjoint odd-cycle `k`-collections. `H` is the number of directed
Hamiltonian paths (Rédei / OCF).

These look unrelated. The hidden connection is that **they are two readings of the same
underlying object — the multiset of closed orbits of `A` — and each is blind to exactly
what the other sees.**

## The split the zeta function cannot see, and the OCF can

Every primitive closed orbit is either a **simple cycle** (counted by `c_k`) or a
loop-erased gluing of several shorter cycles overlapping on shared vertices (the
"overlap configs" of THM-502). So `W_k = c_k + δ_k`, where the **Witt/census defect**
`δ_k` counts the non-simple primitive orbits:

  `δ_3 = δ_4 = δ_5 = 0`  (no room to overlap) ;  `δ_6 = p33`  (overlapping triangle pairs) ;
  `δ_7 = TQ` ;  `δ_8 = Q44 + TF` ;  ...

The zeta function records only the *total* `W_k`. It fixes the **sum** `c_k + δ_k` but
never the **split** — that is exactly why cospectral tournaments can disagree on `c_k`
(THM-500) and on `H` (THM-499): the eigenvalues pin `W_k`, the split floats. The split
is the non-spectral content of the tournament, and `δ_k` first becomes nonzero precisely
when two cycles have room to share a vertex — the spectral boundary.

`H` reads the split directly. Expanding `α_k = Σ 2^k(...)` by inclusion–exclusion turns
products of cycle counts into overlap corrections, and those corrections *are* the
defects `δ_k`. The result (THM-505) is an exact separation:

  **H = (spectral skeleton, a function of the `W_k`) + Σ 2^{level} · (non-spectral split-count).**

At `n = 7` this is the whole story in two carriers:

  `H = [ 1 + 2c3 + 2c5 + 4·C(c3,2) − 4W6 ]  +  4 c6 + 2 c7`.

The bracket is spectral (constant on every cospectral class). The non-spectral part is
`4c6 + 2c7` — and the coefficients are not free: `c6` enters through the **disjoint-pair**
level `α_2`, so it is weighted by `2^2 = 4`; `c7` is a single new odd cycle in `α_1`, so
`2^1 = 2`. **The fugacity `x = 2` of the independence polynomial is precisely the
exchange rate that converts the conflict-graph count into the defect-weighting.** At
`n = 8` the same rule predicts (verified) `H = [spectral] + 2c7 + 4c6 + 4c8 + 4Q44`,
with the coefficients of `c6, c7` unchanged.

## The surprise: one coordinate wears two opposite disguises

THM-499 found `H` non-spectral at `n = 6` for a **disjoint-support** reason: the count
`α_2 = D` of *vertex-disjoint* triangle pairs is the first OCF ingredient the spectrum
misses. THM-502 found the trace `tr A^6` non-spectral for the opposite-looking
**overlap** reason: `c_6` trades against `p33`, the count of *intersecting* triangle
pairs. These seem to be two different mechanisms — cycles that avoid each other vs.
cycles that touch. But

  `α_2 = C(c3,2) − p33 = C(c3,2) − (W6 − c6) = (spectral) + c6`,

so **`α_2` and `c_6` are the same non-spectral coordinate, modulo a spectral constant.**
The disjoint-pair count and the 6-cycle count co-vary perfectly within every cospectral
class (verified 46/47 split classes carry both; the 47th is carried by `c7` alone). The
"avoid" and "touch" descriptions are complementary halves of one partition `C(c3,2) =
α_2 + p33` whose total is spectral; the spectrum fixes the total, and the single floating
degree of freedom can be named either way. Two mechanisms, one coordinate — the kind of
coincidence that is never a coincidence.

## Why the determinant lens is orthogonal to H

The determinant/skew-spectral lens (THM-468/472) is empirically orthogonal to `H`
(THM-499 explained: `det(I+S)` is spectral, `H` is not). The zeta picture says why in
one line: `det(I − uA)` *is* `1/ζ_T(u)`, the generating function of the `W_k` — the half
of the closed-orbit data that ignores the simple/overlap split. `H` is the half that is
*only* the split. Orthogonal because they are the two complementary projections of the
same closed-orbit multiset onto "totals" and "correlations."

## The transcendent line

A tournament carries a gas of closed orbits. You can ask for the **partition function**
of that gas (the zeta function, `Π(1−u^k)^{−W_k}` — totals, mean-field, spectral) or for
its **correlation functional** at fugacity 2 (the OCF `H = I(Ω,2)` — who-meets-whom,
the split). The eigenvalues compute the first and are structurally blind to the second;
the OCF computes the second and packages the first as an inert spectral skeleton. The
project's entire "spectrum can't determine `H`/`c_k`/`α_1`" ladder is one statement —
*the partition function does not determine the correlations* — and THM-505 is the exact
exchange rate between them, with the fugacity `2` as the unit of conversion. Statistical
mechanics on a finite tournament: same two halves, same wall between them, one rung of
the wall paid per length `k` as soon as there is room for a correlation to hide.

## One rung per length, counted exactly

The earlier reflection paid the spectral wall "one rung per length `k`." THM-505's
carrier-dimension probe makes the rung-count literal: the number of *independent*
non-spectral degrees of freedom in `H` is `0, 1, 2, 3` for `n = 5, 6, 7, 8` — exactly
`n − 5`, one new carrier `c_n` per vertex added, and the carriers are precisely the
simple-cycle counts `c6, c7, ..., c_n`. The probe also answers which side the
even-cycle overlaps fall on: at n=8, `Q44` (overlapping 4-cycle pairs) is **not** a new
axis — it is determined by the spectrum together with `(c6, c7, c8)`.

> **CORRECTION (monad-explorer-2026-06-15, the n=9 break).** The next sentence as
> originally written — "the non-spectral content of `H` is *only* the vector
> `(c6, ..., c_n)`; the overlaps are their spectral shadows" — is **true only for
> `n ≤ 8`**. At **n=9** the simple-cycle counts `(c6, c7, c8, c9)` do **NOT** determine
> `H`: explicit cospectral witnesses share `(c6,c7,c8,c9)` yet differ in `H`, and even
> `(c6,c7,c8,c9,Q44)` does not determine `H` — both the overlap config `Q44` and the triple
> packing `T333` detach. So `dim_nonspec(H) = 6 > n−5 = 4` (carriers `{c6,c7,c8,c9,Q44,T333}`),
> and the overlaps are **independent** non-spectral carriers, not shadows. The non-spectral content of `H` is a
> *tower* of cycle-correlation rungs (simple counts; pair-overlaps `Q44`; triple-packings
> `T333`; …), each rung independent of the ones below once `n` gives it room — see THM-505
> and the reflection `the-overlaps-stop-being-shadows-the-correlation-tower`. The
> "spectrum is mean-field, OCF is correlation" wall is the bottom two rungs of that tower,
> not the whole of it.

*Handoffs.* (1) The general-`n` carrier list and weights (HYP-2513) — is every weight a
pure power of 2, set by the independent-set level at which the carrier enters Ω, and is
`dim_nonspec(H) = n − 5` exactly (does `c9` join as the 4th carrier at n=9, where the
first triple-overlap `8·α₃` switches on)? Needs the exhaustive n=9 iso-class data. (2) Is
`H` *linear* in the carriers past n=7? (At n=8, `Q44` is determined by `(c6,c7,c8)` but
perhaps nonlinearly.) (3) Read `H` itself as a special value of a tournament `L`/zeta-type
function in the fugacity variable (ties to the kps7 special-value lens) — `H = I(Ω,2)`,
and the skeleton/defect split may be a functional-equation-like decomposition.*
