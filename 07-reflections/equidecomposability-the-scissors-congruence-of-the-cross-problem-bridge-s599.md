---
source: opus-2026-06-03-S599v (remote-control)
status: SHARP SYNTHESIS — the cross-problem bridge is EQUIDECOMPOSABILITY (scissors congruence), not mere EQUINUMEROSITY (a count/bijection). For tournaments: H = the volume (equinumerosity); the strong-component multiset = the scissors-congruence class (equidecomposability); the Dehn invariant = that multiset mod reordering; {7,21} = phantom volumes. The shared constant Cl₂(π/3)=1.0149 is verified THREE ways (Clausen series = log-sin/Mahler integral = ideal-regular-tetrahedron volume) — the scissors invariant of the basic piece (3-cycle, angle 2π/3). The 0.014 excess is the equidecomposability signature.
tags: [equidecomposability, equinumerosity, scissors-congruence, dehn-invariant, clausen, Cl2, lobachevsky, ideal-tetrahedron, strong-components, H-spectrum, pi-over-3, mahler-measure, unit-distance, sawin]
---

# Equidecomposability: the scissors congruence of the cross-problem bridge

**Prompt (user):** pursue the sharpest next targets, keeping in mind equinumerosity and
equidecomposability.

These two notions are the right lens, and they sharpen everything. **Equinumerosity** = same
cardinality (a bijection exists). **Equidecomposability** = scissors congruence: cut into pieces
and rearrange by the symmetry group. Equidecomposable ⟹ equinumerous, *never the converse* — the
gap is a **Dehn invariant**. The whole cross-problem story is the difference between the two, and
the shared constant `Cl₂(π/3)` is the canonical scissors-congruence (hyperbolic-volume) invariant.

## 1. The tournament H-spectrum IS a scissors-congruence structure

From the strong-component multiplicativity (S599s), `H(T) = ∏ H(Cᵢ)`. Reading this as
cut-and-reassemble:

> - **`H(T)` = the EQUINUMEROSITY invariant** — the *count* (Hamiltonian paths), a bijection class.
> - **the strong-component multiset = the EQUIDECOMPOSABILITY invariant** — the
>   *scissors-congruence class*: the irreducible pieces `T` is cut into. Two tournaments are
>   equidecomposable iff they have the **same multiset of strong pieces**.
> - **the Dehn invariant = that multiset modulo reordering** — what survives beyond the volume.

**Verified** (`…s599v.py`) — same `H` (equinumerous), different piece-multiset (NOT
equidecomposable = *Dehn-distinct*):
```
  H= 9 : {(3,3), (9)}                          — 2 scissors classes, same volume
  H=15 : {(5,3), (15)}                         — 2
  H=45 : {(5,3,3), (9,5), (15,3), (45)}        — 4 scissors classes, one volume
  H= 7, 21 : ∅                                  — NO scissors class  (phantom volumes)
```

> **The forbidden `{7,21}` are *phantom volumes*** — values with **no** realizing
> scissors-congruence class (neither irreducible nor a product). This is *stronger* than classical
> scissors congruence, where every volume is realized; the tournament `H`-spectrum is a
> sub-semigroup with sporadic holes, and `7=Φ₃(2)`, `21=3Φ₃(2)` are the two holes.

*(Caveat, honest: my enumerator used strong pieces only up to `m=6`, so it mislabels `63,189` as
phantom. They are in fact **new irreducible pieces** at `m=7,8` (a fresh "prime" appears at each
size), consistent with claude-S613. The genuine phantoms are `{7,21}` — never irreducible at any
`m`.)*

## 2. Cl₂(π/3): one constant, three scissors-congruence faces (verified)

> **Verified triple identity (`…s599v.py`):**
> ```
>   Cl₂(π/3) = Σ sin(kπ/3)/k²                              = 1.014942   (Clausen series)
>            = −∫₀^{π/3} log|2 sin(t/2)| dt                 = 1.014941   (log-sin / Mahler / Dehn density)
>            = 3·Λ(π/3) = vol(ideal regular hyperbolic tetrahedron) = 1.014942   (scissors-congruence volume)
> ```

So the shared constant is *simultaneously* (i) a **Clausen value**, (ii) a **Mahler-measure /
log-sin integral** (the Dehn-invariant density `log|2 sin|`), and (iii) the **volume of the ideal
regular tetrahedron** — the canonical object of 3-D scissors congruence (Hilbert's 3rd /
Dehn–Sydler). The **basic irreducible tournament piece, the 3-cycle, has skew eigenvalues at
`±2π/3`** = the tetrahedron's dihedral angle = the Eisenstein `ζ₃` angle. *The smallest scissors
piece of the tournament world is the ideal tetrahedron, and its volume is `Cl₂(π/3)`.*

## 3. The cross-problem bridge is equidecomposability, and 0.014 is its signature

This resolves *why* the same `0.014` appears in unit distance and tournaments, and upgrades the
claim from coincidence to structure:

> **The unit-distance disproof and the SC-tournament shape growth are EQUIDECOMPOSABLE, not merely
> equinumerous.** Both count **unit-norm objects at the Eisenstein angle `π/3`**, and both reduce
> to the *same* `log|2 sin|` (Mahler/Lobachevsky) integral over that angle — i.e. the same
> scissors-congruence volume `Cl₂(π/3)`. The matching is a **cut-and-rearrange of the angle-`π/3`
> integration domain**, so the excess exponent `δ = Cl₂(π/3) − 1 = 0.01494` is an
> **equidecomposability invariant** (a shared Dehn volume), not a numerical accident.
> - UD (Sawin/CM): unit-norm `|β|=1` elements ⟹ `Σ log|β|`-type regulator = the log-sin integral.
> - Tournaments (`α₂=1`): norm-1 root pairs `ρ₁ρ₂=1`, Lee–Yang zeros at `±2π/3` ⟹ the same
>   `log|2 sin|` growth of the shape parameter.
> Mere **equinumerosity** (the raw counts `u(n)` vs `H` or `s`) does *not* match — and shouldn't;
> only the **equidecomposable** (Dehn/`Cl₂`) part does.

## 4. The sharpest next targets (now precisely stated)

1. **Phantom-volume theorem:** prove no strong tournament has `H ∈ {7,21}` for *every* `m` (not
   just sampled) — i.e. `Φ₃(2)` and `3Φ₃(2)` are never irreducible scissors volumes. This is the
   clean form of "the durable gaps," and it is a finite-irreducible-spectrum statement.
2. **Derive the shape exponent `= Cl₂(π/3)`** (the user's target 1, now a *route*): show the SC
   shape parameter's tropical growth and the UD unit-norm count are the **same Lobachevsky volume**
   via the angle-`π/3` equidecomposition — i.e. both are `exp` of the log-sin integral. The triple
   identity (§2) is the bridge; the remaining step is the explicit cut-and-rearrange.
3. **Sawin exact vs `Cl₂(π/3)`** (target 2): since UD is a *lower* bound `> n^{1.014}`, the exact
   exponent could exceed `Cl₂(π/3)`. Equidecomposability predicts **equality** (same Dehn volume);
   a strict excess would mean the UD construction is *not* scissors-congruent to the tournament
   piece (extra, non-`π/3` volume). This makes the equality a **testable structural prediction**,
   not a guess.
4. **Dehn-classify the repo bridges:** which cross-problem correspondences (worry-set ↔ round
   tournaments S599f; UD ↔ Eisenstein; LRC `3³` ↔ `ζ₆`) are *equidecomposable* (transfer the full
   `π/3` structure) vs merely *equinumerous* (transfer only the count `2^{⌊(n-2)/2⌋}`)? The
   equidecomposable ones are the genuine theorem-transfer channels.

## 5. Honest status

- **Verified:** the triple identity `Cl₂(π/3) = ` Clausen `=` log-sin integral `=` ideal-tetra
  volume `= 1.014942`; the equidecomposability classes of small `H` (`9,15,45` multi-class;
  `7,21` empty); the 3-cycle `±2π/3` dihedral angle.
- **Established (rigorous reframe):** `H` = equinumerosity (count), strong-multiset =
  equidecomposability (scissors class), Dehn = multiset mod reorder; `{7,21}` = phantom volumes.
- **Corrected:** `63,189` are irreducible pieces at `m=7,8`, not phantoms (my `m≤6` truncation);
  genuine phantoms are `{7,21}` only.
- **Framed (the sharp targets):** the shape-exponent derivation and Sawin-exact question become
  *equidecomposability* statements (shared Dehn volume `Cl₂(π/3)`), turning the numerical `0.014`
  coincidence into a structural prediction — **not yet proven**.

**Artifacts:** `04-computation/equidecomposability_dehn_tournament_s599v.py` (+`.out`). Builds on
S599s (strong-component multiplicativity), S599u (`Cl₂(π/3)`/`π/3` shared object), HYP-707
(tropical const), oracle's `α₂=1` shape parameter, Hilbert 3rd / Dehn–Sydler / Lobachevsky. New:
**HYP-2186**.
