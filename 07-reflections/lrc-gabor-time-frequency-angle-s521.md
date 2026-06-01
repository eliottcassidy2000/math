# The Gabor / time-frequency angle, bridged to the tournament work (S521)

*claudebox-2026-06-01-S521. Investigating a "Gabor" (time-frequency / Gabor-frame)
perspective on LRC and connecting it to the tournament / danger-graph framework.
NOTE: the prompt was garbled ("Gabor"; "tournaments"); I read "Gabor" as Gabor/
time-frequency harmonic analysis and flag the parallel "Galois/cyclotomic" reading
(our A000016 / multiplicative-walk side). This reflection is a framing/recasting —
it re-expresses our results in time-frequency language and exhibits one concrete
bridge — rather than new theorems.*

## Runners as frequencies; the danger count as a Gabor frame operator

Read the runners as **frequencies** `v_1,...,v_m` (observer = the reference/DC).
The safe band `chi = 1_{[1/n,1-1/n]}` is a window on `R/Z`; runner `i` sees the
window dilated by `v_i`. The lonely indicator is `L(t) = prod_i chi(v_i t)`, and
the **danger count** `N(t) = sum_i 1_{B_i}(t)` (`B_i = chi(v_i .)^c`) is the
**diagonal of the time-frequency frame operator** of the dilated-window system. Its
extremes are the **Gabor frame bounds**:

> `A = min_t N(t)` (lower), `B = max_t N(t)` (upper); and **LRC <=> A = 0** — the
> time-frequency system has a *hole* (fails to frame from below).

The total time-frequency density is `sum_i measure(B_i) = 2(n-1)/n -> 2 > 1`:
the system is **overcomplete**, yet LRC is the assertion that it still has a hole.
That is a Balian-Low-flavored phenomenon — overcompleteness does not preclude a
covering gap because the time-frequency lattice has rational **resonances**.

## The one concrete bridge: frame bounds = danger extremes

Computed (exact): the upper frame bound `B = max_t N` equals the **congestion** =
max runners inside a `1/n`-arc = (essentially) the **clique number `omega` of the
danger graph** (the perfect unit-circular-arc graph from the coloring thread). So
the Gabor frame bounds of the runner-window system ARE the danger-graph extremes:

> `(A, B) = (min N, max N) = (loneliness, congestion) = (0 iff LRC, omega-of-danger-graph)`.

This ties the time-frequency picture to the tournament/danger-graph picture by a
genuine identity, not an analogy.

## Resonances = the over-correlation we already found

The **ambiguity / time-frequency overlap** of two runner-windows is exactly the
pairwise overlap `measure(B_i ∩ B_j)`; its maximum is the **doubling pair**
`v_j = 2 v_i` with ratio `n/4` (a time-frequency **commensurability** — a rational
point of the lattice). This is precisely the over-correlation the moment-method
investigation isolated as the obstruction. So "doubling resonance" = "time-frequency
commensurability" = the Gabor lattice's degeneracy — three names for one thing.

The **binding pair** (S521 Thm A) lives at the *sum* frequency `v_i + v_j`, and the
half-turn tournament's edges flip at the *beat* frequencies `v_i - v_j` — the
tournament IS the beat structure of the time-frequency system, and the optimal
lonely time `k/(v_i+v_j)` is a sum-frequency time-frequency event.

## Zak transform / quasi-periodization

`v_i t mod 1` is the quasi-periodization that the **Zak transform** formalizes; the
zeros of a Zak-transform-type object are where the Gabor frame fails — here, the
**lonely set**. So LRC = the Zak-zeros of the runner system are nonempty. The
Diophantine location of those zeros (continued fractions of the `v_i`) is the same
data the cutting-sequence / discrepancy thread exposed.

## The parallel Galois / cyclotomic angle (if "Gabor" was "Galois")

Our companion structure is genuinely **cyclotomic**: the menu count `A000016(m) =
(1/2m) sum_{d|m odd} phi(d) 2^{m/d}` is a Galois/Burnside count over roots of unity
(`phi`, odd divisors = the cyclotomic Galois group's cycle structure); the regular
polygon = the `n`-th roots of unity (Galois-symmetric); the discretized covering of
`Z/N_*` is by **multiplicative dilates** `v_i . (band)` = a Galois/Frobenius action
on residues. So:

- **Gabor (additive/time-frequency):** windows TRANSLATED in time, frequencies `v_i`,
  frame bounds, ambiguity/beats — the *additive* walk side.
- **Galois (multiplicative/cyclotomic):** windows DILATED by `v_i` on `Z/N`, roots
  of unity, `phi`/odd-divisor counts, Frobenius — the *multiplicative* sieve side.

LRC sits where the two meet — the program's recurring "addition meets
multiplication / `A000568`-`A000016`" thesis, now wearing the names **Gabor x
Galois** (time-frequency x cyclotomic). The covering is additive (translates) in
time and multiplicative (dilates) in the residue dual; that duality is exactly the
Fourier/Zak transform.

## Assessment

The Gabor/time-frequency angle is a faithful *recasting*: the danger count is a
frame-operator diagonal, LRC is a vanishing lower frame bound (a hole in an
overcomplete time-frequency system), the upper bound is the danger-graph clique
number, and the resonances are time-frequency commensurabilities (doubling pairs)
— all identities with our prior results, no new theorem. Its value is unifying
language and the Balian-Low-style framing ("overcomplete yet holed, due to rational
resonance"), plus the clean **Gabor x Galois = additive x multiplicative**
synthesis. If "Gabor" meant "Galois," the multiplicative/cyclotomic half above is
the place to dig (the Frobenius action on the `Z/N_*` covering, the cyclotomic
structure of `A000016`); say the word and I'll develop that side in depth.

## Lead

Make the Zak-transform object precise: define `Z(t)` whose zero set is the lonely
set, and study its zeros via Balian-Low / Heil-Ramanathan-Topiwala-type results on
time-frequency translates — asking whether an overcomplete (density `~2`) system at
frequencies `{v_i}` can have a hole only on a Diophantine-thin set, recovering the
covering/discrepancy obstruction in time-frequency form.
