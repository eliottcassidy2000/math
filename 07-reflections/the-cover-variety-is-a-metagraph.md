# The cover variety is a metagraph (and the kills are one interval, moved multiplicatively)

**opus-2026-07-05-S83** (HYP-4126). A pondering session: the open l ≥ 7 cases held against
the repo's supposedly-unrelated threads. Four echoes turned out to be load-bearing.

## 1. The reframing: one interval, twelve translates

Everything the 169-grid does is multiplication in disguise. Cell `a` is killed by runner
value `v` iff `a·v ∈ ±[1,12] (mod 169)` — i.e. the kill set is the **multiplicative
translate** `K(v) = v⁻¹·(±[1,12])`. The cover problem ("no strict witness") is: do twelve
translates of ONE fixed set — one translate `v_r⁻¹·B` from each fiber `v ≡ r (mod 13)` —
cover the unit group?

In discrete-log coordinates (ℤ/169)* ≅ ℤ/156 ≅ ℤ/12 × ℤ/13, this is covering ℤ/156 by 12
additive translates of the fixed 24-set `S = log(±[1,12])`, the shift of translate r
confined to a fixed class mod 12. `S` has exactly two elements per column (the per-column
kill pair; verified perfectly equidistributed), and the free part of each shift is one
coordinate of (ℤ/13)¹². **The dilation shadows are one orbit of the unit action** — but
the action is NOT the straight diagonal: dilating by u shifts every translate by `log u`
AND permutes the positions by r ↦ ur mod 13, with the mod-13 representative reduction
inserting a CARRY twist (log of the reduced representative ≠ reduced log). The spine is
a helix, not a line — the same carry that twists it is the carry the tower levels see.
(The naive straight-diagonal test fails, verified; the orbit statement holds by
construction.)

Two immediate dividends. (i) The abundance of covers stops being surprising: per-column
coverage by 24 quasi-random heights has probability ~0.4, twelve nearly-independent
columns give ~1e-5 of pattern space ≈ 10⁸ covers — matching the observed census growth
and the 0/200k random-sampling rarity simultaneously. (ii) The right measure of a cover's
"specialness" is its position relative to the diagonal — see §3.

## 2. The tower has one-way traffic: witnesses lift

A level-ℓ witness (margin ≥ (13^{ℓ−1}+1)/13^ℓ at t = a/13^ℓ) IS a level-(ℓ+1) witness at
t = 13a/13^{ℓ+1}, with margin to spare: (13^ℓ+13)/13^{ℓ+1} > (13^ℓ+1)/13^{ℓ+1}. So
**cover-hood projects down the tower**, and the set of "covers at every level" is a
closed 13-adic object. The corrected dichotomy conjecture is asymptotic:

> the only 13-adic tower-limit covers are the 13-adic dilations λ·{1..12}, λ ∈ ℤ₁₃*.

False at level 2 (S82: 130k+ covers vs 156 shadows); the real question is the FIRST LEVEL
where non-shadow covers die out. If the dichotomy holds at level 3 (probe running as this
is written), then for any family: witness at level ≤ 3 (kernel row at 2197), or level-3
shadow — and a level-3 shadow with values < 2197 must EQUAL λ·(perm) mod 2197, i.e. be
the dilated standard on its bounded part: the alignment-sliver members (values in
[169, ~1605], all < 2197) would be forced to sit at exact positions λ·i mod 2197 — a
FINITE check over λ. **The sliver would die by a λ-sweep.** That is the cleanest available
strategy for the one open l ≥ 7 case, and it came from pondering, not from pushing the
window machinery harder.

## 3. The metagraph echo (the repo's other half was the template)

The cover variety wears the same three structures as G_n/ℤ₂:

* **The principal line / spine**: the shadow family, λ = 1 → 168 — one helical orbit of
  the unit group (shift + carry-twisted position permutation), the "SC backbone". Every analysis should be oriented relative to it (Hamming
  distance to the nearest shadow = distance to the spine).
* **The mirror**: (r, κ) ↔ (13−r, 12−κ) — kill sets coincide in 78 pairs because
  K(−v) = K(v). This is complement symmetry; the honest object is the MERGED cover
  variety (quotient by the mirror), exactly as V_merged = (A000568+SC)/2 taught us.
* **The sea**: the sampled non-shadow covers differ from each other in 1–3 positions —
  a mutation-connected bulk, the wiggly d=1 layer of pattern space. Whether the sea
  CONNECTS to the spine (covers reachable from shadows by single-position mutations) is
  measurable with the S82 tooling and would say whether covers are "shadows + defects".

The repo's standing instruction — compute every invariant along spine/ribs/sea — has a
literal second application here.

## 4. Paley, Erdős–Turán, Monsky: the "noise" topics were arrows

* **Character sums (Paley/QR thread, flip-rank memory)**: the Fourier bias of
  S = log(±[1,12]) in ℤ/156 is exactly an incomplete character sum over the interval
  [1,12] mod 169: Σ_{x≤12} χ(x). Pólya–Vinogradov/Burgess bound these; small bias ⟹
  covering flexibility (many covers), large bias at low frequencies ⟹ rigidity. The
  general-p statement is crisp and publishable-shaped: *tight-family rigidity for the
  Lonely Runner at n = p+1 runners reduces, level by level, to incomplete character sums
  over [1, p−1] mod p^ℓ*. (−1 is a QR mod 13, so B is QR-symmetric — the Paley heptagon
  thread (HYP-3805) and this grid share the quadratic-residue skeleton.)
* **Erdős–Turán**: the covering-by-translates flexibility is a discrepancy statement
  about S. The codex commits waved these words around emptily; the content lands here.
* **Monsky (2-adic ℝ)**: Monsky proves a combinatorial dissection theorem by extending
  a p-adic valuation to ℝ and coloring. Our tower is natively 13-adic; the witness-lift
  monotonicity is the valuation-coloring's "one-way boundary". A Sperner-style argument
  on the 13-adic tree — color each family class by the level of its first witness —
  is a candidate mechanism for proving the tower-limit dichotomy outright rather than
  level-checking it.
* **Observer-lens (memory)**: same-scale cluster teeth form a PENCIL that coheres near
  t = 0 and re-coheres near every t = k/d (d a pairwise difference): the sliver's good
  windows live inside Bohr sets of the difference vector. The base citation point should
  be CHOSEN near the origin's re-coherence points — the lens picks the window, the
  window doesn't pick the lens.

## 5. What this changes about strategy

The window/fee/descent machinery treats runners as geometry (teeth on a line). The grid
treats them as arithmetic (residues in a group). The pondering's net: at scale ≥ 169 the
arithmetic is STRONGER — one multiplicative interval controls everything, the tower has
monotone traffic, and the last geometric holdout (the sliver) is exactly the place where
the level-3 arithmetic should be invoked instead. Geometry for the crowd and the spread;
arithmetic for the aligned. The tools were already partitioned this way (S81 reflection);
what is new is that the arithmetic side has a LADDER (levels), a SPINE (shadows), and a
CLASSICAL LEVER (character sums) — none of which the geometric side can see.
