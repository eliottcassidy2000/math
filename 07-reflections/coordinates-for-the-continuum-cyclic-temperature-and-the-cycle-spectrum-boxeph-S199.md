# Coordinates for the continuum: cyclic temperature, iso-cyclic shells, and the cycle spectrum

*boxeph-2026-07-21-S199. Object: THM-2013. Owner: find reframes and invent terms/lenses for the
continuum so we can stop enumerating tournaments and focus on the interesting near-regular behavior.
Builds on THM-1979 (the spectrum), THM-1926 (the zeta / cycle moments), and death-star-S84
(H≥disc's binding case = quasirandom = the continuum center).*

## The problem with enumeration

The continuum (the near-regular interior, THM-1979) has 47 iso classes in one n=7 score shell, 6880
tournaments at n=8, ~10⁶ at n=9. Enumerating it is hopeless and, worse, uninformative — the point is
not the list but the *shape*. We need **intrinsic coordinates**: a handful of numbers that locate a
tournament in the continuum and describe near-regular behavior without ever listing the members. Here
are the coordinates, and the two lenses (thermodynamic and harmonic) they come with.

## The coordinates (new terms)

**Cyclic temperature** `τ = c₃/c₃_max = 1 − σ²/σ²_max ∈ [0,1]`. The one macroscopic coordinate,
score-spread rescaled: `τ=0` at the transitive **ground state** (frozen, ordered), `τ=1` at the
regular **hot** center. The continuum is the high-`τ` region. (Because `c₃ = n(n²−1)/24 − (n/2)σ²`
exactly, THM-1979, temperature is literally cyclicity — it needs *no enumeration*, just the scores.)

**Iso-cyclic shell** `𝒮_τ` = the set of classes at fixed `τ` (fixed `c₃`, fixed score-spread).
Tournament space is a stack of shells: the transitive is the `τ=0` **singleton** shell; the shells
**swell** toward `τ=1`. The shell is the right unit — not the individual tournament.

**Structural entropy** `S(τ) = log₂ |𝒮_τ|`. Zero at the ground state, maximal at the hot center. This
is the continuum's "size" as a function of temperature — a smooth macroscopic curve, not a list.

**Cycle spectrum** (a.k.a. **cyclic harmonics**) `(N₄, N₅, …, N_n)`, `N_k = tr(Aᵏ)` = #closed
k-walks (the zeta moments, THM-1926). Structure: `N₁=N₂=0` (no loops/digons), `N₃=3c₃` is the
**fundamental**, *frozen* by `τ`; the **overtones** `N₄,…,N_n` are the free coordinates that resolve
a tournament *within* its shell. A near-regular tournament ≈ its cyclic harmonics.

**Frozen vs free.** *Frozen* invariants are score-determined (scores, `c₃`, `τ`, `σ²`) — they place
the shell. *Free* invariants carry structure (the overtones `N₄⁺`, and the beyond-spectral `|R|`) —
they place the tournament in the shell. Describe the continuum by the **free** coordinates only.

## Lens 1 — thermodynamic

Read `τ` as temperature and `S` as entropy. Then:
- the **transitive tournament is the T=0 ground state** (unique, zero entropy, char_A=xⁿ, ζ=1);
- the **regular/quasirandom tournament is the T=1 hot phase** (maximal entropy, all strong, the
  positive-entropy continuum of tournamentons `W≈½`);
- **score spread is the order parameter** (magnetization); cyclicity is the disorder; the n=7
  perfection-breaking (odd holes, spectral collapse) is a **phase transition** — the temperature at
  which the ordered description (reduction principles) stops covering the phase.

This is why the reduction hierarchy (THM-1862/1926/1960) describes the *cold* rim and the mathematics
lives in the *hot* interior: reductions are the low-temperature expansion; the continuum is beyond its
radius of convergence. death-star-S84's "H≥disc binding case = quasirandom" is the statement that the
hardest inequality is saturated at `τ=1` — the hot center — exactly where enumeration fails.

## Lens 2 — harmonic

Read the tournament as a signal on the cycle basis. `N_k = Σ_j λ_jᵏ` (eigenvalue power sums), so the
cycle spectrum IS the char-poly, and the **overtones resolve structure**: the fundamental `N₃` sets
the temperature, and `N₄, N₅, …` are the timbre that distinguishes near-regular tournaments sharing a
temperature. Where the harmonic lens runs out — cospectral tournaments (spectral collapse, 89% ties at
n=7) — the **beyond-harmonic** coordinate `|R|` (mac-mini THM-1966, first independent at n=7) takes
over. So the coordinate budget of the continuum is layered:
```
   L0  cyclic temperature τ           (1 real; from the scores, no enumeration)
   L1  cycle spectrum N₄…N_n          (= char_A; resolves shells to cospectral classes)
   L2  beyond-spectral |R|            (separates the cospectral twins, from n=7)
```
A near-regular tournament is pinned by **`τ` + a short cycle spectrum + `|R|`** — a low-dimensional
address. The continuum is not `10⁶` objects; it is a low-dimensional coordinate cloud with a
temperature axis and an entropy profile.

## The payoff (how to stop enumerating)

To study near-regular behavior: fix a temperature `τ≈1`, and describe the shell `𝒮_τ` by the
*distribution* of its free coordinates (the cycle-spectrum cloud + `|R|`), and its entropy `S(τ)`.
The interesting questions become continuous: *how does the cycle-spectrum cloud spread as `τ→1`? where
is the diversity peak (THM-1979 saw it slightly off-center)? what is `S(τ)` asymptotically?* — all
answerable per-shell, per-coordinate, without listing tournaments. The transitive point and the
quasirandom horizon are the two boundary conditions; everything between is the temperature flow.

## Two features the coordinates reveal (verified, n≤7)

**The diversity maximum is at intermediate temperature.** Structural entropy `S(τ)=log₂|𝒮_τ|` does
*not* peak at the hot center: at n=7 it peaks at `τ=5/7` (`c₃=10`, 79 classes, `S=6.30`), while the
`τ=1` regular shell holds only 3 classes (`S=1.58`). The continuum is fattest *just inside* the hot
edge — the diversity peak THM-1979 glimpsed, now located: `τ*≈0.7`. (Thermodynamically this is a
specific-heat-like peak between the ordered and quasirandom phases.)

**An all-strong condensation threshold.** The strong-fraction of a shell jumps sharply to 1 at a
critical temperature: every class with `τ ≥ 9/14 ≈ 0.64` is strongly connected at n=7 (`τ≥3/4` at n=6,
`τ≥3/5` at n=5), and reducible classes appear only below it. There is a genuine **condensation
temperature** `τ_c` above which the whole shell is irreducible.

## Verified anchors (n≤7)

- `N₁=N₂=0`, `N₃=3c₃` frozen by the score sequence; the **first free moment is `N₄`** (varies within
  a score sequence from n=5).
- Shell entropy peaks at `τ*≈0.7`; all-strong condensation threshold `τ_c≈0.64` (n=7).
- Coordinate budget in the hot n=7 shells: cycle spectrum alone resolves `21/47` of the biggest shell
  (cospectral collapse), **`+|R|` resolves `36/47`** — and `15/15` of the `τ=13/14` shell. So
  `(τ, char_A, |R|)` pins most near-regular tournaments; a small irreducible residue survives at the
  very center (the deep continuum needs one more coordinate).

## Threads worked (S200 → THM-2016)

- **Local subtournament density is an L3 probe, with score kept as a sidecar:**
  at the `c3=12` shell, 4&5-profiles genuinely improve
  `(char_A,|R|,score)` resolution `36→44` of 47.  At `c3=11`, the apparent
  `41→50` 4-profile gain is already supplied by score, so it is not independent
  local evidence.  **3/47 still survive
  (score, char_A, |R|, 4-profile, 5-profile)**. The deep center is
  *invariant-resistant* even after global, score, and local sidecars.
- **Reducibility ceiling (proved, proof corrected by THM-2000):**
  `max c₃ over reducible = c₃_max(n−1)`.  Cycles live inside SCCs, but their
  counts must be **summed**, not bounded by the largest summand.  The
  nondecreasing increments
  `c₃_max(n+1)-c₃_max(n)=T_floor(n/2)` let us concentrate every SCC-size
  partition at `(n−1,1)`.  Thus
  `τ_c=c₃_max(n−1)/c₃_max(n)` is the attained reducible ceiling
  (1/2, 2/5, 5/8, 4/7 for n=4..7), while the first all-strong discrete shell is
  `(c₃_max(n−1)+1)/c₃_max(n)` (1, 3/5, 3/4, 9/14).
- **H is a thermometer:** mean H rises monotonically with `τ` (1 → 178 at n=7), the spread carried by
  the free coordinates — locating death-star-S84's `H≥disc` binding case at the hot center `τ=1`.

Links: THM-2013, THM-2016, THM-1979, THM-1926, THM-1966, THM-1960,
[[tournament-space-as-a-spectrum-single-point-to-continuum-boxeph-S198]],
[[the-n-ge-7-regime-what-breaks-what-survives-boxeph-S197]].
