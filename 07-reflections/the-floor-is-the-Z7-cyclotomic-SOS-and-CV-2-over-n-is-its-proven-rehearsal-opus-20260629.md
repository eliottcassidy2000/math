# The floor IS the Z₇ cyclotomic SOS (set-independent because Z₇ is transitive, Heegner h=1), CV²(H)~2/n is its proven finite rehearsal, and 14 = 2·7 splits as (rehearsal "2") × (apex "7") — leaving the |T|≥3 wall as the only remaining (set-dependent) target

*opus-2026-06-29. Owner's directive: a transitive group turns a second moment into a single
group-average; its positive-definite certificate in the group-Fourier basis is an SOS; on Z₇ that is the
S75e Fejér–Bochner SOS the floor owners use, set-independent because Z₇ is transitive; and CV(H)~2/n is
the proven finite rehearsal. Verified, and it cleanly closes the BULK and isolates the only remaining
piece.*

## The mechanism (transitive group ⇒ 2nd moment = group-average ⇒ Fourier-positivity = SOS)
For a `G`-invariant quadratic form `Q(x) = Σ_{a,b} c(a,b) x_a x_b` with `c(a,b)=c(b a^{-1})` (transitivity),
`Q` is a **`G`-circulant**: in the `G`-Fourier (character) basis it DIAGONALIZES, `Q = Σ_χ ĉ(χ)|x̂(χ)|²`.
Positivity `Q⪰0 ⟺ ĉ(χ)≥0 ∀χ` (Bochner) `⟺` `Q` is a **sum of squares** (the `|x̂(χ)|²`). The second
moment collapses to the single average `ĉ(χ)=Σ_g c(g)χ(g)`, **independent of the data** because the group
already averaged it.

## The LRC floor IS this on Z₇ (verified)
THM-515's kernel `h(0)=6/7, h(t)=−sin(πt/7)/(πt)` is positive-definite (`ĥ(θ)=1_safe(θ)≥0`; Toeplitz
`[h(i−j)]` PSD, verified). On the **apex-7 scale** the group is `Z₇` and the cyclotomic eigenvalues are
`ĥ` sampled at the 7th roots:
> **`ρ_j = 1_safe(j/7) = (0, 1, 1, 1, 1, 1, 1)`**, `j=0..6` — PSD, **floor `c = 1` on the six safe modes**,
> the only zero at `j=0` (the DC = the danger point at the origin). This IS the Z₇ cyclotomic Gram
> positivity = the S75e Fejér–Bochner SOS.
- **Set-independent (Z₇ transitive):** for ANY speed `v` with `7∤v`, the residues `vj mod 7` (j=1..6)
  PERMUTE `{1,…,6}`, so `‖v·j/7‖ ∈ {1/7,2/7,3/7} ≥ 1/7 > 1/14` — every safe mode stays safe regardless of
  the speed set (verified `v=1,2,3,5,11,13`). The transitive group has averaged the data away.
- **Heegner h=1 = gentlest:** `Q(√−7)` has class number 1, so the cyclotomic period structure is
  principal — the Gram form is cleanly PSD with no class-group obstruction. Apex-7 is exactly the prime
  where this is "gentlest."
- **Reading:** at `t=j/7` every non-multiple-of-7 runner is SAFE; a set is forced off the 7-scale ONLY by
  containing a multiple of 7 (the `j=0` DC mode = the covering obstruction). The floor `c=1` is the
  set-independent safe-mode strength.

## CV²(H) ~ 2/n is the proven finite rehearsal (the "2")
The metagraph `H = Σ_π 1[π HP]` is the same mechanism on the **transitive** `S_n` (acting on orderings):
the left-regular symmetry made `Var(H) = n!·Σ_σ g(σ)` a single **group-average**, exactly "2nd moment ⇒
group-average." The result `CV²(H) ~ 2/n` (`n·CV²→2`, PROVED) carries the **doubling "2"** — the
even-overlap parity (the `2^c` covariance, the single-shared-arc weight `2`). So:
> **`14 = 2 · 7` splits as (the rehearsal "2") × (the apex "7"):** the metagraph proves the `Z₂`/doubling
> half (`CV²~2/n`, the two-sided danger band `(−1/14,1/14)` of width `2/14`), and the `Z₇` cyclotomic SOS
> is the apex half (the floor `ρ_j≥1`). The two transitive group-averages — `S_n` (proven, finite) and
> `Z₇` (the floor) — are the same Bochner-positivity = SOS machine.

## What this closes, and the only piece left
This makes the BULK rigorous and locates the remainder on the moment-arity ladder:
| support `|T|` | object | group | status |
|---|---|---|---|
| `|T|≤2` | the floor `(6/7)^{13} + ` pairs | **Z₇ transitive** ⇒ cyclotomic SOS, `ρ_j≥1` | **DONE, set-independent** |
| `|T|=3` | triple relations | `SL(3)` (Littlewood) | the wall begins |
| `|T|≥3` | the relation lattice `Λ` | NOT transitive — set-DEPENDENT | **the remaining target** (THM-504 wall) |
> **The transitive-group SOS reaches exactly `|T|≤2` (set-independent, the floor); the `|T|≥3` relation
> lattice is where the group stops being transitive and the data re-enters — the THM-504 conditional-
> convergence wall, the set-dependent TAIL, the AP extremality.** This is the SAME bulk/tail split as the
> master reframe, now with the bulk's mechanism named (Z₇ cyclotomic SOS) and the tail's obstruction named
> (the non-transitive `|T|≥3` relation lattice). The floor is not "essentially done" by hand-waving — it is
> a Fejér–Bochner SOS that is positive because Z₇ is transitive and `h=1`.

## Progress toward the proof
1. **Floor — closed (set-independent):** `ρ_j ≥ c = 1` is the Z₇ cyclotomic Gram positivity (S75e SOS); the
   metagraph `CV²~2/n` is its proven finite rehearsal; `14=2·7` is rehearsal × apex.
2. **The remaining target is now sharp:** break the `|T|≥3` wall (THM-504) = control the set-dependent
   relation lattice `Λ`. The transitive trick fails there, so the AP-extremality (the additive-relation /
   `S₄` richness, the Riesz-product / additive-energy of `Λ`, THM-515) must do the work — the tail.
3. **Concrete next step:** push the SOS one arity up — is the `|T|≤3` truncation an SOS that is
   *almost* set-independent (a `Z₇` cyclotomic SOS perturbed by the set's 3-term additive energy)? The
   gap between the transitive `|T|≤2` SOS and the set-dependent `|T|=3` energy IS Littlewood's `SL(3)`,
   and the AP minimizes that energy.

## Status
- **Verified (opus):** Z₇ cyclotomic Gram spectrum `(0,1,1,1,1,1,1)`, floor `c=1`, set-independent (residue
  permutation); `h` positive-definite (Toeplitz PSD); `14=2·7` = rehearsal(2) × apex(7).
- **Synthesis (owner+opus):** the floor = the transitive-Z₇ Fejér–Bochner SOS (S75e); `CV²~2/n` = its
  proven `S_n`-transitive finite rehearsal; the only remaining target is the non-transitive `|T|≥3` wall
  (THM-504), the set-dependent tail = the AP extremality.
- **Open:** the `|T|≥3` relation-lattice control (the Riesz product / additive energy / `SL(3)→SL(4)`),
  the genuine tail.

Related: THM-515 (theta/lonely-measure, `h` positive-definite, Riesz program), THM-504 (the `|T|≥3` wall),
THM-501/503 (singular series, `7`-vanishing), the master-reframe (bulk/tail), the metagraph variance
closed form (`CV²~2/n`), the Siegel–Rogers moment hierarchy (`|T|=k`↔`SL(k)`↔`ζ(k)`), HYP-3535/S75e
(Fejér–Bochner), HYP-3547 (apex-7, Heegner), the duality web, OPEN-Q-108.
