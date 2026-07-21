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

**Cyclic temperature**
`τ = c₃/c₃,max = (σ²_tr−σ²)/(σ²_tr−σ²_min) ∈ [0,1]`, where
`σ²_min=0` for odd `n` and `1/4` for even `n`. The one macroscopic
coordinate is score-spread rescaled: `τ=0` at the transitive **ground state**
and `τ=1` at the regular (odd) or near-regular (even) maximum-cyclic edge.

**Iso-cyclic shell** `𝒮_τ` = the set of classes at fixed `τ` (fixed `c₃`, fixed score-spread).
Tournament space is a stack of shells: the transitive is the `τ=0`
**singleton** shell. Shell sizes are not monotone in `τ`; the shell is the
right coarse unit, not the individual tournament.

**Structural entropy** `S(τ) = log₂ |𝒮_τ|`. It is zero at the ground
state and, already at n=7, peaks at an intermediate temperature rather than
the hot edge. This records the shell-size profile instead of pretending it is
monotone.

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
- the **regular/near-regular edge is the T=1 hot phase** (maximum cyclicity and
  all strong, but not necessarily maximum finite shell entropy or quasirandom);
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
n=7) — the **beyond-harmonic** coordinate `|R|` refines the spectrum from n=6
and becomes independent even of `(spectrum,H)` at n=7 (THM-1966). So the
coordinate budget is layered:
```
   L0  cyclic temperature τ           (1 real; from the scores, no enumeration)
   L1  cycle spectrum N₄…N_n          (= char_A; resolves shells to cospectral classes)
   L2  beyond-spectral |R|            (refines cospectral fibers from n=6)
```
A near-regular tournament is partially localized by **`τ` + a short cycle
spectrum + `|R|`**. The unresolved fibers are part of the object, not noise:
the continuum is a layered coordinate cloud with a temperature axis and an
entropy profile, rather than a complete low-dimensional parameterization.

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

**A reducible ceiling and first all-strong shell.** At n=7, reducible classes
attain `τ_red=4/7`; the next shell, `τ_all=9/14`, and every shell above it
are all strong. The corresponding first all-strong values are `3/4` at n=6
and `3/5` at n=5. The one-integer-shell gap matters when naming the threshold.

## Verified anchors (n≤7)

- `N₁=N₂=0`, `N₃=3c₃` frozen by the score sequence; the **first free moment is `N₄`** (varies within
  a score sequence from n=5).
- Shell entropy peaks at `τ*≈0.7`; `τ_red=4/7` and `τ_all=9/14` at n=7.
- Coordinate budget in the hot n=7 shells: cycle spectrum alone resolves `21/47` of the biggest shell
  (cospectral collapse), **`+|R|` resolves `28/47`** — and `13/15` of the `τ=13/14` shell.  The old
  36/47 and 15/15 rows accidentally used signed `R` while labeling it `|R|`.
  Thus `(τ, char_A, |R|)` improves the address but leaves substantial fibers;
  after local profiles, six twin pairs still survive in this shell.

## Threads worked (S200 → THM-2016)

- **Local subtournament density is an L3 probe, with score kept as a sidecar:**
  under the corrected absolute `|R|` carrier, 4&5-profiles improve the final
  address `28→41` of 47 at `c3=12`, `46→50` of 52 after score at `c3=11`, and
  `13→14` of 15 at `c3=13`.  Six, two, and one twin pairs remain respectively.
  The deep center is *invariant-resistant* even after global, score, and local
  sidecars.
- **Reducibility ceiling (proved, proof corrected by THM-2000):**
  `max c₃ over reducible = c₃_max(n−1)`.  Cycles live inside SCCs, but their
  counts must be **summed**, not bounded by the largest summand.  The
  nondecreasing increments
  `c₃_max(n+1)-c₃_max(n)=T_floor(n/2)` let us concentrate every SCC-size
  partition at `(n−1,1)`.  Thus
  `τ_red=c₃_max(n−1)/c₃_max(n)` is the attained reducible ceiling
  (1/2, 2/5, 5/8, 4/7 for n=4..7), while the first all-strong discrete shell is
  `τ_all=(c₃_max(n−1)+1)/c₃_max(n)` (1, 3/5, 3/4, 9/14).
- **H is a thermometer:** mean H rises monotonically with `τ` (1 → 178 at n=7), the spread carried by
  the free coordinates — locating death-star-S84's `H≥disc` binding case at the hot center `τ=1`.

Links: THM-2013, THM-2016, THM-1979, THM-1926, THM-1966, THM-1960,
[[tournament-space-as-a-spectrum-single-point-to-continuum-boxeph-S198]],
[[the-n-ge-7-regime-what-breaks-what-survives-boxeph-S197]].
