/-
  TournamentH7.LRCFarPeelGood — THE FAR-PEEL POSITIVITY on the ACTUAL good region
  (kind-pasteur-2026-07-03-S30).

  opus-S49 (LRCFarPeelCore) proved far-peel positivity in the single-`comb`
  formulation; mac-mini/kps carry the good region as `goodRegion2` (the wrapped
  danger, `dangerPair`).  This file lands the positivity DIRECTLY on `goodRegion2`,
  via the already-proved `damped_peel` (the wrapped-`dangerPair` rate) — so no
  comb/dangerPair bridge is needed.  This is opus's steps 2+3 fused in the right
  formulation:

      length (goodRegion2 (E ++ [w]) h) > 0

  once the base has positive length (step 1 — mac-mini's floor) and the far speed `w`
  clears the peel threshold `(#pieces)·4h < (1−2h)·length·w`.  With
  `exists_lonely_of_goodRegion2_pos` (step 4, kps-S14) this is the far-element peel of
  `CoveringFarLonely 22` for large `w`, minus the finite small-`w` window (step 5).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import TournamentH7.LRCPeelAssembly
import TournamentH7.LRC13Citation

namespace LonelyRunner
namespace RatIntervals

/-- **FAR-PEEL POSITIVITY (good-region form).**  If the base good region has positive
length and the far speed `w` clears the peel threshold
`(#pieces)·4h < (1−2h)·length·w`, the peeled good region `goodRegion2 (E ++ [w])`
still has positive length.  Direct consequence of `damped_peel` (the wrapped
`dangerPair` rate). -/
theorem goodRegion2_peel_pos (E : List ℤ) (w : ℤ) (hw : 0 < w) {h : ℚ}
    (hh0 : 0 < h) (hh1 : h ≤ 1 / 2)
    (hbig : ((goodRegion2 E h).length : ℚ) * (4 * h)
      < (1 - 2 * h) * length (goodRegion2 E h) * (w.toNat : ℚ)) :
    0 < length (goodRegion2 (E ++ [w]) h) := by
  have hd := damped_peel E w hw hh0 hh1
  have hwnat : 0 < w.toNat := by omega
  set G : ℚ := length (goodRegion2 E h) with hG
  set c : ℚ := ((goodRegion2 E h).length : ℚ) with hc
  set W : ℚ := (w.toNat : ℚ) with hW
  have hWpos : (0 : ℚ) < W := by rw [hW]; exact_mod_cast hwnat
  -- hd : (1 − 2h)·G − c·(4·(h/W)) ≤ length (good (E ++ [w]))
  -- hbig : c·(4h) < (1 − 2h)·G·W
  have hthr : c * (4 * h / W) < (1 - 2 * h) * G := by
    have he : c * (4 * h / W) = c * (4 * h) / W := by ring
    rw [he, div_lt_iff₀ hWpos]
    linarith [hbig]
  linarith [hd, hthr]

/-- **FAR-PEEL POSITIVITY, division-cleared integer threshold.**  Same, with the
threshold as `(#pieces)·4h < (1−2h)·length·w` — the integer-friendly `w`-bound for
certifying covering-far families. -/
theorem goodRegion2_peel_pos_of_gt (E : List ℤ) (w : ℤ) (hw : 0 < w) {h : ℚ}
    (hh0 : 0 < h) (hh1 : h ≤ 1 / 2)
    (hlen : 0 < length (goodRegion2 E h))
    (hwbig : ((goodRegion2 E h).length : ℚ) * (4 * h)
      < (1 - 2 * h) * length (goodRegion2 E h) * (w.toNat : ℚ)) :
    0 < length (goodRegion2 (E ++ [w]) h) :=
  goodRegion2_peel_pos E w hw hh0 hh1 hwbig

/-- **THE FAR-ELEMENT PEEL CLOSER** (opus steps 2+3+4 composed): a positive-speed
`13`-family whose last runner `w = v (Fin.last 12)` clears the peel threshold against
the base of the first twelve is lonely.  Reduces `CoveringFarLonely` for a large far
runner to the ONE remaining lemma — the base good-region floor (step 1, from
`LRC(≤13)`), which makes the threshold `hbig` satisfiable.  (No separate positivity of
the base is needed: `hbig` forces it, since a zero-length base makes the strict
inequality `(#pieces)·4h < 0` impossible.) -/
theorem far_peel_lonely (v : Fin 13 → ℤ) (hv : ∀ i, 0 < v i)
    (hbig : ((goodRegion2 (List.ofFn (Fin.init v)) (1 / 14)).length : ℚ) * (4 * (1 / 14))
      < (1 - 2 * (1 / 14)) * length (goodRegion2 (List.ofFn (Fin.init v)) (1 / 14))
          * ((v (Fin.last 12)).toNat : ℚ)) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hw : 0 < v (Fin.last 12) := hv _
  have hpeel := goodRegion2_peel_pos (List.ofFn (Fin.init v)) (v (Fin.last 12)) hw
    (by norm_num) (by norm_num) hbig
  have heq : List.ofFn v = List.ofFn (Fin.init v) ++ [v (Fin.last 12)] := by
    rw [List.ofFn_succ', List.concat_eq_append]
    rfl
  rw [← heq] at hpeel
  exact exists_lonely_of_goodRegion2_pos v hv hpeel

/-- **A member forces positive length** (reverse of `exists_mem_of_length_pos`): if a
rational `x` lies in a region `L`, then `L` has positive length — because the piece
containing `x` has `p.1 ≤ x < p.2`, hence is non-degenerate.  This reduces the base
good-region floor (step 1) to producing ONE rational good point. -/
theorem length_pos_of_mem {L : Region} {x : ℚ} (hx : mem x L) : 0 < length L := by
  obtain ⟨p, hpL, hp1, hp2⟩ := hx
  have hpos : (0 : ℚ) < p.2 - p.1 := by linarith
  unfold length
  have hnn : ∀ y ∈ L.map fun q => max 0 (q.2 - q.1), (0 : ℚ) ≤ y := by
    intro y hy
    rw [List.mem_map] at hy
    obtain ⟨q, _, rfl⟩ := hy
    exact le_max_left _ _
  have hmem : (max 0 (p.2 - p.1)) ∈ L.map fun q => max 0 (q.2 - q.1) :=
    List.mem_map_of_mem hpL
  have hle := List.single_le_sum hnn _ hmem
  have heq : max 0 (p.2 - p.1) = p.2 - p.1 := max_eq_right hpos.le
  rw [heq] at hle
  linarith

/-- **dangerPair soundness**: a point of the danger comb is within `h` of some integer
multiple.  (Boundary-inclusive `≤`; the strict good point below beats it.) -/
theorem mem_dangerPair_le {s : ℕ} (hs : 0 < s) {h : ℚ} {x : ℚ}
    (hmem : mem x (dangerPair s h)) : ∃ m : ℤ, |(s : ℚ) * x - m| ≤ h := by
  rw [dangerPair, mem_union] at hmem
  rcases hmem with h1 | h2
  · rw [TournamentH7.CombPatterns.mem_comb hs] at h1
    obtain ⟨k, _, hlo, hhi⟩ := h1
    refine ⟨(k : ℤ), ?_⟩
    rw [abs_le]; push_cast
    constructor <;> linarith
  · rw [TournamentH7.CombPatterns.mem_comb hs] at h2
    obtain ⟨k, _, hlo, hhi⟩ := h2
    refine ⟨(k : ℤ) + 1, ?_⟩
    rw [abs_le]; push_cast
    constructor <;> linarith

/-- **Reverse membership** (good ⟹ in the good region): a rational `x ∈ [0,1)` that is
STRICTLY `h`-far from every integer multiple of every base speed lies in
`goodRegion2 base h`.  With `length_pos_of_mem` this closes the base floor (step 1) from a
single strict-good rational point — which `LRC(≤13)` (margin `1/13 > 1/14`) supplies. -/
theorem good2_mem_of_strict {speeds : List ℤ} {h : ℚ}
    (hpos : ∀ s ∈ speeds, 0 < s) {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1)
    (hgood : ∀ s ∈ speeds, ∀ m : ℤ, h < |(s : ℚ) * x - m|) :
    mem x (goodRegion2 speeds h) := by
  rw [goodRegion2, mem_diffF]
  refine ⟨⟨(0, 1), List.mem_singleton_self _, hx0, hx1⟩, ?_⟩
  intro q hq hxq
  rw [List.mem_flatMap] at hq
  obtain ⟨s, hs, hqs⟩ := hq
  have hspos := hpos s hs
  have hsnat : 0 < s.toNat := by omega
  have hmem : mem x (dangerPair s.toNat h) := ⟨q, hqs, hxq⟩
  obtain ⟨m, hm⟩ := mem_dangerPair_le hsnat hmem
  have hcast : ((s.toNat : ℕ) : ℚ) = ((s : ℤ) : ℚ) := by
    exact_mod_cast Int.toNat_of_nonneg hspos.le
  rw [hcast] at hm
  exact absurd hm (not_le.mpr (hgood s hs m))

/-- **THE BASE GOOD-REGION FLOOR from one strict-good rational** (step 1, member form):
a single rational `x ∈ [0,1)` strictly `h`-far from every base multiple forces
`0 < length (goodRegion2 base h)`. -/
theorem goodRegion2_length_pos_of_strict {speeds : List ℤ} {h : ℚ}
    (hpos : ∀ s ∈ speeds, 0 < s) {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1)
    (hgood : ∀ s ∈ speeds, ∀ m : ℤ, h < |(s : ℚ) * x - m|) :
    0 < length (goodRegion2 speeds h) :=
  length_pos_of_mem (good2_mem_of_strict hpos hx0 hx1 hgood)

/-- **THE DENSITY BRIDGE** (step 1, closed in Lean): the base good-region floor from the
`LRC(≤13)` citation.  Cite the 12 base speeds → a REAL `t0` with margin `1/13`; a rational
`x` within `1/(182·ΣB)` of `t0` (density) is STRICTLY `1/14`-far from every base multiple
(the `1/13 − 1/14 = 1/182` slack); its fractional part is a strict-good point in `[0,1)`,
so `0 < length (goodRegion2 base (1/14))`.  This is mac-mini's THM-609 continuity step,
and the crisp rational/irrational reading: a REAL lonely point (possibly irrational) forces
a RATIONAL one because the good set is open (opus-S50: bounded magnitude = rational census,
large magnitude = irrational sweep). -/
theorem base_floor_of_cite (cite : LonelyRunner.LRCUpTo13) (base : Fin 12 → ℤ)
    (hpos : ∀ i, 0 < base i) :
    0 < length (goodRegion2 (List.ofFn base) (1 / 14)) := by
  obtain ⟨t0, ht0⟩ := cite 12 (by norm_num) base (fun i => (hpos i).ne')
  set V : ℤ := ∑ i, base i with hV
  have hVpos : 0 < V := Finset.sum_pos (fun i _ => hpos i) ⟨0, Finset.mem_univ 0⟩
  have hVR : (0 : ℝ) < (V : ℝ) := by exact_mod_cast hVpos
  have hle : ∀ i, base i ≤ V := fun i =>
    Finset.single_le_sum (fun j _ => (hpos j).le) (Finset.mem_univ i)
  set δ : ℝ := 1 / (182 * (V : ℝ)) with hδ
  have hδpos : (0 : ℝ) < δ := by rw [hδ]; positivity
  obtain ⟨x, hx1, hx2⟩ := exists_rat_btwn (show t0 - δ < t0 + δ by linarith)
  have hxt : |(x : ℝ) - t0| < δ := by rw [abs_lt]; constructor <;> linarith
  -- strict good at x (real), for every base speed
  have hgoodQ : ∀ i : Fin 12, ∀ m : ℤ, (1 : ℚ) / 14 < |(base i : ℚ) * x - m| := by
    intro i m
    have hbR : (0 : ℝ) < (base i : ℝ) := by exact_mod_cast hpos i
    have hble : (base i : ℝ) ≤ (V : ℝ) := by exact_mod_cast hle i
    have htri : |(base i : ℝ) * t0 - m|
        ≤ |(base i : ℝ) * (x : ℝ) - m| + (base i : ℝ) * |(x : ℝ) - t0| := by
      have h1 : (base i : ℝ) * t0 - m
          = ((base i : ℝ) * (x : ℝ) - m) + (base i : ℝ) * (t0 - (x : ℝ)) := by ring
      rw [h1]
      calc |((base i : ℝ) * (x : ℝ) - m) + (base i : ℝ) * (t0 - (x : ℝ))|
          ≤ |(base i : ℝ) * (x : ℝ) - m| + |(base i : ℝ) * (t0 - (x : ℝ))| := abs_add_le _ _
        _ = |(base i : ℝ) * (x : ℝ) - m| + (base i : ℝ) * |(x : ℝ) - t0| := by
            rw [abs_mul, abs_of_pos hbR, abs_sub_comm t0 (x : ℝ)]
    have hm13 := ht0 i m
    have hcast : (1 : ℝ) / (13 : ℕ) = 1 / 13 := by norm_num
    rw [hcast] at hm13
    have hbd : (base i : ℝ) * |(x : ℝ) - t0| < 1 / 182 := by
      calc (base i : ℝ) * |(x : ℝ) - t0| < (base i : ℝ) * δ :=
            mul_lt_mul_of_pos_left hxt hbR
        _ ≤ (V : ℝ) * δ := mul_le_mul_of_nonneg_right hble hδpos.le
        _ = 1 / 182 := by rw [hδ]; field_simp
    have hgR : (1 : ℝ) / 14 < |(base i : ℝ) * (x : ℝ) - m| := by
      have : (1 : ℝ) / 13 - 1 / 182 = 1 / 14 := by norm_num
      linarith
    have hcastabs : |(base i : ℝ) * (x : ℝ) - (m : ℝ)| = ((|(base i : ℚ) * x - m| : ℚ) : ℝ) := by
      push_cast; ring_nf
    rw [hcastabs] at hgR
    have : ((1 : ℚ) / 14 : ℝ) < ((|(base i : ℚ) * x - m| : ℚ) : ℝ) := by push_cast; linarith
    exact_mod_cast this
  -- fractional part: strict good preserved, lands in [0,1)
  set x' : ℚ := Int.fract x with hx'
  have hx'0 : 0 ≤ x' := Int.fract_nonneg x
  have hx'1 : x' < 1 := Int.fract_lt_one x
  refine goodRegion2_length_pos_of_strict ?_ hx'0 hx'1 ?_
  · intro s hs
    obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
    exact hpos i
  · intro s hs m
    obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
    have hcase := hgoodQ i (m + base i * ⌊x⌋)
    have heq : (base i : ℚ) * x - ((m + base i * ⌊x⌋ : ℤ) : ℚ)
        = (base i : ℚ) * x' - (m : ℚ) := by
      rw [hx', Int.fract]; push_cast; ring
    rwa [heq] at hcase

/-- **THE DENSITY BRIDGE, general arity** (`j ≤ 12`): the base good-region floor from the
`LRC(≤13)` citation for ANY base of at most 12 positive speeds.  Cite the `j` base speeds → a
REAL `t0` with margin `1/(j+1) ≥ 1/13`; a rational `x` within `1/(182·ΣB)` of `t0` is STRICTLY
`1/14`-far from every base multiple (the `1/(j+1) − 1/14 ≥ 1/182` slack), so
`0 < length (goodRegion2 base (1/14))`.  This is the arity the peels actually consume: removing
`k` far runners leaves a `(13−k)`-speed base — exactly `j = 13−k ∈ [7,12]` for the `≤ 6`-far
simultaneous peel.  (`j = 0` is the empty base, good region `[0,1)`, length `1`.) -/
theorem base_floor_of_cite_gen (cite : LonelyRunner.LRCUpTo13) {j : ℕ} (hj : j ≤ 12)
    (base : Fin j → ℤ) (hpos : ∀ i, 0 < base i) :
    0 < length (goodRegion2 (List.ofFn base) (1 / 14)) := by
  rcases Nat.eq_zero_or_pos j with hj0 | hjpos
  · -- empty base: goodRegion2 [] (1/14) = [(0,1)], length 1
    subst hj0
    have hnil : (List.ofFn base) = ([] : List ℤ) := by
      apply List.ofFn_zero
    rw [hnil]
    show (0 : ℚ) < length (goodRegion2 [] (1 / 14))
    rw [goodRegion2]
    simp only [List.flatMap_nil]
    show (0 : ℚ) < length [((0 : ℚ), 1)]
    simp only [length, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
    norm_num
  · -- nonempty base: cite → t0 at margin 1/(j+1), then the density argument
    obtain ⟨t0, ht0⟩ := cite j hj base (fun i => (hpos i).ne')
    set V : ℤ := ∑ i, base i with hV
    have hVpos : 0 < V :=
      Finset.sum_pos (fun i _ => hpos i) ⟨⟨0, hjpos⟩, Finset.mem_univ _⟩
    have hVR : (0 : ℝ) < (V : ℝ) := by exact_mod_cast hVpos
    have hle : ∀ i, base i ≤ V := fun i =>
      Finset.single_le_sum (fun j _ => (hpos j).le) (Finset.mem_univ i)
    set δ : ℝ := 1 / (182 * (V : ℝ)) with hδ
    have hδpos : (0 : ℝ) < δ := by rw [hδ]; positivity
    obtain ⟨x, hx1, hx2⟩ := exists_rat_btwn (show t0 - δ < t0 + δ by linarith)
    have hxt : |(x : ℝ) - t0| < δ := by rw [abs_lt]; constructor <;> linarith
    -- margin 1/(j+1) ≥ 1/13
    have hjR : ((j : ℝ) + 1) ≤ 13 := by
      have : (j : ℝ) ≤ 12 := by exact_mod_cast hj
      linarith
    have hjRpos : (0 : ℝ) < (j : ℝ) + 1 := by positivity
    have hmargin : (1 : ℝ) / 13 ≤ 1 / ((j : ℝ) + 1) := one_div_le_one_div_of_le hjRpos hjR
    have hgoodQ : ∀ i : Fin j, ∀ m : ℤ, (1 : ℚ) / 14 < |(base i : ℚ) * x - m| := by
      intro i m
      have hbR : (0 : ℝ) < (base i : ℝ) := by exact_mod_cast hpos i
      have hble : (base i : ℝ) ≤ (V : ℝ) := by exact_mod_cast hle i
      have htri : |(base i : ℝ) * t0 - m|
          ≤ |(base i : ℝ) * (x : ℝ) - m| + (base i : ℝ) * |(x : ℝ) - t0| := by
        have h1 : (base i : ℝ) * t0 - m
            = ((base i : ℝ) * (x : ℝ) - m) + (base i : ℝ) * (t0 - (x : ℝ)) := by ring
        rw [h1]
        calc |((base i : ℝ) * (x : ℝ) - m) + (base i : ℝ) * (t0 - (x : ℝ))|
            ≤ |(base i : ℝ) * (x : ℝ) - m| + |(base i : ℝ) * (t0 - (x : ℝ))| := abs_add_le _ _
          _ = |(base i : ℝ) * (x : ℝ) - m| + (base i : ℝ) * |(x : ℝ) - t0| := by
              rw [abs_mul, abs_of_pos hbR, abs_sub_comm t0 (x : ℝ)]
      have hm := ht0 i m
      push_cast at hm
      have hbd : (base i : ℝ) * |(x : ℝ) - t0| < 1 / 182 := by
        calc (base i : ℝ) * |(x : ℝ) - t0| < (base i : ℝ) * δ :=
              mul_lt_mul_of_pos_left hxt hbR
          _ ≤ (V : ℝ) * δ := mul_le_mul_of_nonneg_right hble hδpos.le
          _ = 1 / 182 := by rw [hδ]; field_simp
      have hgR : (1 : ℝ) / 14 < |(base i : ℝ) * (x : ℝ) - m| := by
        have h14 : (1 : ℝ) / 13 - 1 / 182 = 1 / 14 := by norm_num
        linarith
      have hcastabs : |(base i : ℝ) * (x : ℝ) - (m : ℝ)| = ((|(base i : ℚ) * x - m| : ℚ) : ℝ) := by
        push_cast; ring_nf
      rw [hcastabs] at hgR
      have : ((1 : ℚ) / 14 : ℝ) < ((|(base i : ℚ) * x - m| : ℚ) : ℝ) := by push_cast; linarith
      exact_mod_cast this
    set x' : ℚ := Int.fract x with hx'
    have hx'0 : 0 ≤ x' := Int.fract_nonneg x
    have hx'1 : x' < 1 := Int.fract_lt_one x
    refine goodRegion2_length_pos_of_strict ?_ hx'0 hx'1 ?_
    · intro s hs
      obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
      exact hpos i
    · intro s hs m
      obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
      have hcase := hgoodQ i (m + base i * ⌊x⌋)
      have heq : (base i : ℚ) * x - ((m + base i * ⌊x⌋ : ℤ) : ℚ)
          = (base i : ℚ) * x' - (m : ℚ) := by
        rw [hx', Int.fract]; push_cast; ring
      rwa [heq] at hcase

/-- **The density bridge, list form.**  For any list `B` of at most 12 positive speeds, the
`LRC(≤13)` citation gives `0 < length (goodRegion2 B (1/14))`.  This is the shape the peels
consume: the base of a far-peel (`E = init v`, 12 speeds) or of the simultaneous ≤6-far peel
(`B`, the `13−|far|` bounded runners) is exactly such a list, so its good-region FLOOR is
discharged unconditionally by the citation node — no computation of the base region needed. -/
theorem base_floor_list_of_cite (cite : LonelyRunner.LRCUpTo13) (B : List ℤ)
    (hpos : ∀ w ∈ B, 0 < w) (hlen : B.length ≤ 12) :
    0 < length (goodRegion2 B (1 / 14)) := by
  have hget : List.ofFn (B.get) = B := List.ofFn_get B
  rw [← hget]
  refine base_floor_of_cite_gen cite hlen (B.get) ?_
  intro i
  exact hpos (B.get i) (List.get_mem B i)

/-! ## The quantitative floor: an all-good interval survives with its full length

`base_floor_of_cite` gives `0 < length`.  The tower rung (opus HYP-4046) needs the
QUANTITATIVE floor `length ≥ 2δ`.  The bridge is a pure Region-measure fact:
`length_ge_of_mem_cover` — if every point of `[a, b)` is a member of a `Norm` region, that
region has length at least `b − a`.  Combined with `norm_goodRegion2` (below) this upgrades
the base floor to an explicit lower bound. -/

/-- **Cursor covering lower bound.**  If every point of `[a, b)` is a member of a `Norm`
region `L`, then `length L ≥ b − a`.  The interval is swept left to right: the first piece
holding the cursor either reaches `b` (done) or advances the cursor to its right end, and the
`Norm` gap guarantees the tail covers the remainder.  This is the measure content the fleet
requested as `length_ge_of_safe_interval` (Norm length-extensionality, covering form). -/
theorem length_ge_of_mem_cover : ∀ {L : Region}, Norm L → ∀ {a b : ℚ},
    (∀ x : ℚ, a ≤ x → x < b → mem x L) → b - a ≤ length L := by
  intro L
  induction L with
  | nil =>
      intro _ a b hcov
      rcases lt_or_ge a b with hab | hab
      · exact absurd (hcov a le_rfl hab) (by simp [mem])
      · have h0 : length ([] : Region) = 0 := rfl
        rw [h0]; linarith
  | cons p L' ih =>
      intro hN a b hcov
      rcases lt_or_ge a b with hab | hab
      · have hp12 : p.1 < p.2 := norm_head_lt hN
        have hlen : length (p :: L') = (p.2 - p.1) + length L' := by
          simp only [length, List.map_cons, List.sum_cons]
          rw [max_eq_right (by linarith : (0 : ℚ) ≤ p.2 - p.1)]
        by_cases hp2a : p.2 ≤ a
        · -- p is entirely left of the cursor: it covers nothing in [a, b)
          have hcov' : ∀ x : ℚ, a ≤ x → x < b → mem x L' := by
            intro x hx1 hx2
            rcases hcov x hx1 hx2 with ⟨q, hq, hq1, hq2⟩
            rcases List.mem_cons.mp hq with rfl | hqL'
            · exact absurd hq2 (not_lt.mpr (by linarith))
            · exact ⟨q, hqL', hq1, hq2⟩
          have hrec := ih (norm_tail hN) hcov'
          rw [hlen]; linarith [length_nonneg L']
        · push_neg at hp2a
          have hap : p.1 ≤ a := by
            rcases hcov a le_rfl hab with ⟨q, hq, hq1, hq2⟩
            rcases List.mem_cons.mp hq with rfl | hqL'
            · exact hq1
            · have := norm_head_le hN q hqL'; linarith
          by_cases hp2b : b ≤ p.2
          · rw [hlen]; linarith [length_nonneg L']
          · push_neg at hp2b
            have hcov' : ∀ x : ℚ, p.2 ≤ x → x < b → mem x L' := by
              intro x hx1 hx2
              rcases hcov x (by linarith) hx2 with ⟨q, hq, hq1, hq2⟩
              rcases List.mem_cons.mp hq with rfl | hqL'
              · exact absurd hq2 (not_lt.mpr hx1)
              · exact ⟨q, hqL', hq1, hq2⟩
            have hrec := ih (norm_tail hN) hcov'
            rw [hlen]; linarith
      · linarith [length_nonneg (p :: L')]

end RatIntervals
end LonelyRunner
