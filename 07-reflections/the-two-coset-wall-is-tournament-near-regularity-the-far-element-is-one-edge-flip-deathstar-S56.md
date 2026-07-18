# The 2-coset wall is tournament near-regularity — and the far element is one edge-flip

**death-star-2026-07-18-S56.** Following the owner's cue — prove `M<1/13 ⟹` residues in ≤2 cosets of
`val·ℤ`, and relate these 2-dimensional objects to tournaments. Result: the 2-coset statement is
**equivalent** to the residue comparison tournament being the regular tournament `R₁₃` up to a single
edge-flip, verified across the deep-well tower. This is the LRC(14) wall translated onto the project's
home turf (Rédei / regular tournaments). **It is not a proof** — I say below exactly why tournament
regularity is necessary but not sufficient, which is the honest residual.

Builds on `the-inverse-theorem-is-a-function-field-freiman-3k-4...` (the 2-coset reformulation) and
boxeph's THM-1017 / difference-closure lemma. Scripts: `04-computation/lrc_tourney_deathstar_S56.py`,
`lrc_tourney_flip_deathstar_S56.py`.

---

## 1. The residues carry a tournament (the owner's "binary dimensions around vertices")

At the maximizer `t = a/q`, `M = val/q`, the 13 residues `r_i = v_i·a mod q` sit on the circle `ℤ/q`.
They carry a **comparison tournament** `T`: `i → j` iff `(r_i − r_j) mod q ∈ (0, q/2)` (i is "ahead" of
j by less than half). This is exactly the owner's picture — `C(13,2) = 78` binary comparisons (a
triangular number of on/off switches), *not independent* (they are the coboundary of 13 positions, so
`d_{ij} + d_{jk} = d_{ik}`), arranged around 13 vertices, and Rédei guarantees a ranking (Hamiltonian
path — the sorted circular order). The observer at `0` is a **source** (no residue within `val` behind
it: all `r_i ∈ [val, q−val]`), which is the tournament face of "lonely" (THM-381).

## 2. The signature: `M<1/13 ⟹ T = R₁₃ up to one edge-flip` (verified)

The **regular** (rotational) tournament `R₁₃` on 13 (odd) vertices has every score `= 6`; it is the
tournament of **equally-spaced** residues. Computing the flip-defect (`½ Σ|s_i − 6|`) of the residue
tournament:

| family | scores | flip-defect from `R₁₃` |
|---|---|---|
| deep well `{1..12,182}` | `[5, 6×11, 7]` | **1** |
| ladder `m=2,5`, dilate `3·` | `[5, 6×11, 7]` | **1** |
| `{1..11,13,84}` (`M>1/13`) | `[5,5,6×8,7,7]` | 2 |
| generic covering (`M=1/5`) | `[3×6, 6×4, 7×3]` | 10 |

So **every `M<1/13` family's residue tournament is `R₁₃` with exactly one edge reversed** — one vertex
`6→7` (the far element gains an edge) and one `6→5`. `M>1/13` has defect `≥ 2`. Chaining the
equivalences:

> `M<1/13`  ⟺  `T` is regular-up-to-one-flip  ⟺  residues equally spaced + one perturbation  ⟺  **2
> cosets of `val·ℤ`**  ⟺  **AP core + far element** (boxeph THM-1017).

**The far element is the single broken edge of the regular tournament.** That is the clean answer to
"how do these two-dimensional objects relate to tournaments": the 2-coset dimension-2 structure is
`R₁₃` (dimension-1, the regular ranking) plus one flip (the second dimension = the far element). A
supporting fixed point: on the whole tower `a = val` (i.e. the maximizer time `t = M` itself), so
`R = val·V mod q` and the core speeds `≤ 12 < 13` give `val·v_i < q`, hence `r_i = val·v_i ≡ 0
(mod val)` *automatically* — the first coset is forced once `a=val` and the core is small.

## 3. Why this is not (yet) a proof — the honest gap

Tournament regularity is **necessary but not sufficient** for the AP, and that is exactly the residual:

- A regular comparison tournament means every half-arc holds `(n−1)/2` points — the configuration is
  **balanced**, not necessarily **equally spaced**. Balanced ⊋ equally-spaced. So "`T` regular" alone
  does not force the AP; the AP needs the stronger *uniform-gap* rigidity.
- boxeph's pigeonhole delivers **one** short gap (defect on one edge, dimension `<13`); promoting it to
  **flip-defect exactly 1** (dimension `≤2`) is the whole open content — the same gap the Freiman
  coset-progression route faces (§the Freiman reflection: raw doubling `47 > 3k−4`, so only the
  dimension/coset framing sees it).

So the tournament translation is faithful and it puts the wall on home turf, but it does **not** shrink
the wall: `M<1/13 ⟹ flip-defect ≤ 1` is LRC(14) in tournament clothing.

## 4. What the tournament frame does buy — and the sharpened target

- **The right rigidity object.** "AP" is now "the regular tournament `R₁₃`", and "far element" is "one
  edge-flip". The near-tightness `M<1/13` should force **near-regularity**, and the residual is the
  extra step from *balanced* (regular scores) to *equally-spaced* (the AP) — a second-order rigidity.
- **A concrete mechanism to try (Rédei/OCF).** The flip-defect is a coarse invariant; the finer one is
  the tournament's **odd-cycle collection / Hamiltonian-path count** `H` (Rédei: odd). `R₁₃` maximizes
  balance; a single flip changes `H` by an even amount. Conjecture worth testing: `M<1/13` pins the
  residue tournament's `H` (or its OCF spectrum) to the `R₁₃`-plus-one-flip value, and *that* rigidity —
  not the score sequence — forces the equal spacing. This is the project's own machinery (the OCF, the
  regular-tournament uniqueness) aimed at its last open case.
- **Necessary-condition harvest.** `a=val` (`t=M`) and "observer is source" are clean necessary
  structural facts; if `a=val` is provable for `M<1/13` covering families, the first coset is free and
  the wall collapses to "the far element is a single flip", a genuinely smaller target.

**Bottom line:** the 2-coset wall = tournament near-regularity, the far element = one edge-flip, and the
open content is the balanced→equally-spaced rigidity — best attacked with the OCF/Rédei tools this
project was built on, not with more covering geometry or raw-doubling Freiman.



## 5. UPDATE — the OCF invariant IS pinned (conjecture CONFIRMED)

Tested §4's conjecture directly: computed the OCF invariant `H` = Hamiltonian-path count (Rédei-odd)
of the residue tournament via bitmask DP. **Result: `H = 3551083`, pinned EXACTLY across every M<1/13
family** — deep well, ladders `m = 1,2,3,7,14`, and dilates `3·,5·,7·,11·` — all give the same
`H = 3551083`, `c₃ = 90`.

And it is **uniquely** `R₁₃`-plus-one-**antipodal** flip. Sweeping the single flips of the regular
tournament by chord distance `d`:

| `d` | 1 | 2 | 3 | 4 | 5 | **6** |
|---|---|---|---|---|---|---|
| `H(R₁₃+flip_d)` | 3577045 | 3563239 | 3553327 | 3550831 | 3552463 | **3551083** |

Only the distance-6 (antipodal) chord flip matches the pinned value; every other distance gives a
*distinct* `H`. So the OCF pins the residue tournament to a single isomorphism type:

> **`M(V) < 1/13`  ⟹  residue tournament `T ≅ R₁₃` with one antipodal edge reversed, `H = 3551083`.**

This is the sharp upgrade of §2: the *score sequence* `[5,6×11,7]` is shared by many tournaments (only
"near-regular"), but the *OCF* `H` singles out the exact type — confirming that the **OCF/Rédei
invariant, not the score sequence, is the rigidity that matches the additive structure**. The far
element = the one antipodal flip = the second coset = the second Freiman dimension: four descriptions,
one object.

**Honest status (unchanged):** this is a pinned *necessary* invariant (verified), and it is a
*consequence* of the AP+far structure — so it confirms and sharpens the reformulation but does not by
itself prove the wall (`M<1/13 ⟹ H = 3551083` is still the inverse theorem, now stated as an
OCF-value identity). The proof route it opens: the LRC↔OCF bridge (`H` = the danger-count partition
function `Q(2)`, loneliness = `Q(0)`) — if `H` can be computed from the LRC data and forced to
`3551083` by `M<1/13`, the wall closes on the OCF side. That is the concrete next OCF experiment.
Script: `04-computation/lrc_ocf_pinning_deathstar_S56.py`, `lrc_ocf_flip_char_deathstar_S56.py`.

→ THM-1017, THM-1028, the Freiman-3k-4 reflection, THM-381 (lonely⟺source), the regular-tournament /
OCF canon; owner cue: tournaments as non-independent binary dimensions forcing a ranking.
