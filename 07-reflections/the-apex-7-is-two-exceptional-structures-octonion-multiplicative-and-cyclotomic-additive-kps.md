# The apex 7 is TWO exceptional structures: octonion (multiplicative, tournament) and cyclotomic (additive, LRC)

*kind-pasteur-2026-06-27-S31an. Pushing the abstraction: I tested whether the LRC's apex-7 structure is
the octonion/G₂/Fano exceptional object (`14=dim G₂=dim Aut(O)`, `7`=imaginary octonions=Fano=the G₂
7-rep). The test is a productive NEGATIVE that sharpens the whole picture: the octonion (multiplicative)
structure is the TOURNAMENT side's apex; the LRC side's apex is the CYCLOTOMIC (additive) one. The
two-maps split (HYP-3099) is precisely this multiplicative↔additive duality at 7, bridged by Joukowski.*

## The seductive numerology (real, but on which side?)
- `14 = dim G₂ = dim Aut(octonions)`; `7` = the imaginary octonion units = the **G₂ 7-dimensional irrep**
  = the **Fano plane** `PG(2,2)` (`|PSL(2,7)|=168`); the 7-dim cross product exists only in dim 3 and 7.
- The project's **doubling map `z↦2z mod 7`** (orbits `{1,2,4}=QR`, `{3,5,6}=NQR`, order 3) **is** the
  octonion multiplication's base Fano line `e₁e₂=e₄` and its conjugate (the "two Fano planes" of `T₇`,
  mac-mini). So the *tournament* `H=I(Ω,2)` / Burnside `2^orbits` / apex-7 cycle structure lives on the
  octonion multiplication.

## The test (productive negative): the octonion does NOT organize the LRC associator
`lrc_octonion_g2_fano_associator_kps.py`, the associativity-compression failure `κ₃` (3-way joint
empty-sector cumulant, S31ai) over the 20 inner triples of consec_8:
- **`κ₃` is NOT doubling-invariant** (only 2/20 triples invariant under `z↦2z`) — the octonion
  automorphism is **broken** by the LRC.
- The **Fano-line triples are not extremal** (`{1,2,4}` κ₃=0.0215 vs the max 0.030 at `{2,3,4}`); QR vs
  NQR differ but unremarkably.
- The single-sector `p_j` carry the **reflection symmetry** `s↦6−s` (`p₁=p₅`, `p₂=p₄`), **not** the
  multiplicative QR/NQR symmetry.

## The resolution: multiplicative (tournament) vs additive (LRC) at the apex
The LRC orbit `{frac(i x)}` is an **additive** arithmetic progression on the circle; its symmetries are
**translation + reflection** (`s↦6−s` = complement = `T↦Tᵒᵖ`), and its apex structure is the
**cyclotomic** field `ℚ(ζ₇)` (the de Moivre cubic, additive Gaussian periods, S31ak). The tournament's
`Ω`/`H` is **multiplicative** (the `z↦2z` doubling, QR(7), Burnside `2^orbits`), and its apex structure is
the **octonion/G₂/Fano** multiplication. So:

> **The apex prime 7 carries TWO exceptional structures: a MULTIPLICATIVE one (octonion / G₂ / Fano /
> QR-doubling = the TOURNAMENT / cycle / `H` side) and an ADDITIVE one (cyclotomic `ℚ(ζ₇)` / de Moivre /
> Gaussian periods = the LRC / coverage side). The project's TWO MAPS (HYP-3099) — tournament on the real
> axis, LRC on the circle — are exactly this multiplicative↔additive duality of the apex, and JOUKOWSKI
> `w=z+1/z` (S31am) is the conformal bridge between them.**

This is the deep reason the octonion was *seductive but mislocated*: it is genuinely there at 7, but on
the tournament face, where the multiplication lives. The LRC face is the additive (cyclotomic) shadow of
the same prime. Both are "why 7 is exceptional"; they differ by additive vs multiplicative.

## Why this is more than bookkeeping (the synthesis payoff)
- It **explains the two maps' geometry**: real-axis (tournament) = the multiplicative/Galois `ℝ`-locus
  (`I(Ω,x)` real-rooted, the `2cos` = trace form of the multiplicative group); circle (LRC) = the additive
  roots-of-unity locus (`ℚ(ζ₇)`). Joukowski `z+1/z` = the **trace map** `ζ↦ζ+ζ̄=2cos` that sends the
  multiplicative circle to the additive real cubic — literally the Galois trace `ℚ(ζ₇)→ℚ(cos 2π/7)`.
- It **predicts where each tool applies**: octonion/Moufang/alternativity → tournament `Ω` identities
  (the `{7,21}` H-gaps, the BIBD/Paley multiplicative maximizers); cyclotomic/three-gap/Diophantine → LRC
  coverage. Don't import octonion identities into the LRC associator (the test shows they break); do
  import them into the tournament `H`-structure.
- It **refines "why n=14 first-open"** (S31ak): both the additive (cubic CF, non-periodic) AND the
  multiplicative (octonion non-associativity) exceptionalities of 7 switch on at the apex; the LRC feels
  the additive one (the Diophantine cubic wall), the tournament feels the multiplicative one (the H-gap).

## The associativity-failure, relocated
My S31ai "associativity-compression failure = the odd/Worpitzky `κ₃` = the non-associative apex" stands,
but its home is now precise: the **non-associativity is the octonion's** (multiplicative, tournament-side),
and what the LRC `κ₃` measures is the **additive shadow** of that — the cyclotomic 3-cocycle, the
real-cubic Galois obstruction, not the octonion associator itself (which is doubling-symmetric and breaks
here). The even/odd parity split = the additive reflection (`ℤ/2` complement), not the octonion grading.

## Net
- **Productive negative:** the octonion/G₂/Fano structure is real at 7 but lives on the TOURNAMENT
  (multiplicative) face, not the LRC (additive) face — verified by `κ₃`'s broken doubling-invariance.
- **Synthesis:** apex 7 = multiplicative (octonion/G₂, tournament) ⊕ additive (cyclotomic, LRC); the two
  maps are this duality; **Joukowski = the Galois trace `ζ↦ζ+ζ̄`** bridging them.
- **Guidance:** octonion/Moufang tools → tournament `Ω`/`H`; cyclotomic/Diophantine tools → LRC coverage.

→ HYP-3211 (this), HYP-3099 (two maps), HYP-3162 (cyclotomic/Joukowski=trace), HYP-3160/3161 (parity/
associator/ferromagnetic), the-four-faces-of-14, mac-mini two-Fano-planes, G₂/octonions/Fano, OPEN-Q-108.
