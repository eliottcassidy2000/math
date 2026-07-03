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
import TournamentH7.LRCSimulPeel

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

/-! ### `Norm (goodRegion2 …)`: the good region is ordered, disjoint, nondegenerate

The region difference keeps degenerate pieces (LRCRegionDiff design note), so the
un-`F` `diff` is not `Norm`.  But `cutF`/`diff1F`/`diffF` FILTER degenerate pieces
(`length_cutF = length_filter_live`), and cutting a live interval out of a sorted-disjoint
region preserves order — so `diffF` of a `Norm` region by LIVE intervals is `Norm`.  This
chain establishes `Norm (goodRegion2 …)`, the missing hypothesis for the cursor bound. -/

/-- Prepend one interval to a `Norm` region, given it ends at or before every existing
piece's start. -/
theorem norm_cons {p : ℚ × ℚ} {L : Region} (hp : p.1 < p.2) (hL : Norm L)
    (hjoin : ∀ r ∈ L, p.2 ≤ r.1) : Norm (p :: L) := by
  rcases L with _ | ⟨q, L'⟩
  · exact hp
  · exact ⟨hp, hjoin q (List.mem_cons_self), hL⟩

/-- Concatenate two `Norm` regions when every left piece ends at or before every right
piece's start. -/
theorem norm_append : ∀ {A : Region}, Norm A → ∀ {B : Region}, Norm B →
    (∀ a ∈ A, ∀ b ∈ B, a.2 ≤ b.1) → Norm (A ++ B) := by
  intro A
  induction A with
  | nil => intro _ B hB _; simpa using hB
  | cons p A' ih =>
      intro hA B hB hjoin
      have hp : p.1 < p.2 := norm_head_lt hA
      have hnorm : Norm (A' ++ B) :=
        ih (norm_tail hA) hB (fun a ha b hb => hjoin a (List.mem_cons_of_mem _ ha) b hb)
      apply norm_cons hp hnorm
      intro r hr
      rw [List.mem_append] at hr
      rcases hr with hrA' | hrB
      · exact norm_head_le hA r hrA'
      · exact hjoin p (List.mem_cons_self) r hrB

/-- Every piece of `cutF p q` lies within `[p.1, p.2]`. -/
theorem cutF_bounds {p q r : ℚ × ℚ} (hr : r ∈ cutF p q) : p.1 ≤ r.1 ∧ r.2 ≤ p.2 := by
  unfold cutF cut at hr
  rw [List.mem_filter] at hr
  obtain ⟨hmem, _⟩ := hr
  rcases List.mem_cons.mp hmem with rfl | hmem'
  · exact ⟨le_refl _, min_le_left _ _⟩
  · rcases List.mem_singleton.mp hmem' with rfl
    exact ⟨le_max_left _ _, le_refl _⟩

/-- Cutting a live interval out of one nondegenerate interval yields a `Norm` region. -/
theorem norm_cutF {p q : ℚ × ℚ} (hp : p.1 < p.2) (hq : q.1 ≤ q.2) : Norm (cutF p q) := by
  have hmid : min p.2 q.1 ≤ max p.1 q.2 := by
    have h1 := min_le_right p.2 q.1
    have h2 := le_max_right p.1 q.2
    linarith
  unfold cutF cut
  by_cases h1 : (p.1 : ℚ) < min p.2 q.1
  · by_cases h2 : max p.1 q.2 < p.2
    · rw [List.filter_cons_of_pos (by simpa using h1),
          List.filter_cons_of_pos (by simpa using h2), List.filter_nil]
      exact ⟨h1, hmid, h2⟩
    · rw [List.filter_cons_of_pos (by simpa using h1),
          List.filter_cons_of_neg (by simpa using h2), List.filter_nil]
      exact h1
  · by_cases h2 : max p.1 q.2 < p.2
    · rw [List.filter_cons_of_neg (by simpa using h1),
          List.filter_cons_of_pos (by simpa using h2), List.filter_nil]
      exact h2
    · rw [List.filter_cons_of_neg (by simpa using h1),
          List.filter_cons_of_neg (by simpa using h2), List.filter_nil]
      trivial

/-- Subtracting one live interval preserves `Norm`. -/
theorem norm_diff1F : ∀ {L : Region}, Norm L → ∀ {q : ℚ × ℚ}, q.1 ≤ q.2 →
    Norm (diff1F L q) := by
  intro L
  induction L with
  | nil => intro _ q _; trivial
  | cons p L' ih =>
      intro hN q hq
      have hp : p.1 < p.2 := norm_head_lt hN
      have hstep : diff1F (p :: L') q = cutF p q ++ diff1F L' q := by
        unfold diff1F; rw [List.flatMap_cons]
      rw [hstep]
      apply norm_append (norm_cutF hp hq) (ih (norm_tail hN) hq)
      intro a ha b hb
      have hab := cutF_bounds ha
      unfold diff1F at hb
      rw [List.mem_flatMap] at hb
      obtain ⟨p'', hp''L', hbcut⟩ := hb
      have hbb := cutF_bounds hbcut
      have hpp'' := norm_head_le hN p'' hp''L'
      linarith [hab.2, hbb.1, hpp'']

/-- Subtracting a list of live intervals preserves `Norm`. -/
theorem norm_diffF : ∀ {B : Region}, ∀ {L : Region}, Norm L → (∀ q ∈ B, q.1 ≤ q.2) →
    Norm (diffF L B) := by
  intro B
  induction B with
  | nil => intro L hL _; exact hL
  | cons q B' ih =>
      intro L hL hlive
      have hstep : diffF L (q :: B') = diffF (diff1F L q) B' := by
        unfold diffF; rw [List.foldl_cons]
      rw [hstep]
      exact ih (norm_diff1F hL (hlive q (List.mem_cons_self)))
        (fun q' hq' => hlive q' (List.mem_cons_of_mem _ hq'))

/-- **`Norm (goodRegion2 …)`**: for positive speeds and `0 ≤ h`, the good region is a
`Norm` region — ordered, disjoint, nondegenerate.  (The base window `[0,1)` is `Norm` and
every `dangerPair` arc is live.) -/
theorem norm_goodRegion2 {speeds : List ℤ} (hpos : ∀ s ∈ speeds, 0 < s) {h : ℚ}
    (hh : 0 ≤ h) : Norm (goodRegion2 speeds h) := by
  unfold goodRegion2
  have hbase : Norm [((0 : ℚ), 1)] := by unfold Norm; norm_num
  apply norm_diffF hbase
  intro q hq
  rw [List.mem_flatMap] at hq
  obtain ⟨s, hs, hqs⟩ := hq
  have hsnat : 0 < s.toNat := by have := hpos s hs; omega
  exact dangerPair_live s.toNat hsnat hh q hqs

/-- **THE QUANTITATIVE FLOOR (interval-survival form)**: if every point of `[a, b)` is
`h`-good for the base `speeds` (strictly `h`-far from every multiple, `x ∈ [0,1)`), the good
region has length at least `b − a`.  This is the strengthening opus/mac-mini requested for
the TOWER RUNG (`length_ge_of_safe_interval`, HYP-4046): the base floor is not merely
positive but bounded below by the safe-interval width.  Composes `length_ge_of_mem_cover`
(the covering measure bound) with `norm_goodRegion2` (the good region is `Norm`) and
`good2_mem_of_strict` (a strict-good point is a member). -/
theorem length_ge_of_safe_interval {speeds : List ℤ} {h : ℚ}
    (hpos : ∀ s ∈ speeds, 0 < s) (hh : 0 ≤ h) {a b : ℚ}
    (hsafe : ∀ x : ℚ, a ≤ x → x < b → 0 ≤ x ∧ x < 1 ∧
      ∀ s ∈ speeds, ∀ m : ℤ, h < |(s : ℚ) * x - m|) :
    b - a ≤ length (goodRegion2 speeds h) := by
  apply length_ge_of_mem_cover (norm_goodRegion2 hpos hh)
  intro x hx1 hx2
  obtain ⟨hx0, hx1', hgood⟩ := hsafe x hx1 hx2
  exact good2_mem_of_strict hpos hx0 hx1' hgood

/-- **THE QUANTITATIVE BASE FLOOR from the citation** (explicit bound).  Sharpens
`base_floor_of_cite`'s `0 < length` to an EXPLICIT lower bound `1/(400·ΣB) ≤ length`: the
`LRC(≤13)` citation gives a real `t0` at margin `1/13`; a rational `x` within `1/(400·ΣB)` of
`t0` has fractional part `x'` that is `≥ 1/13 − 1/400`-good, and a one-sided rational interval
of width `1/(400·ΣB)` beyond `x'` (chosen inside `[0,1)` by casing on `x' ≤ 1/2`) stays
`> 1/14`-good (slack `1/13 − 1/200 > 1/14`).  `length_ge_of_safe_interval` then gives the
width as the floor.  This is the tower-rung (opus HYP-4046) quantitative floor, discharged. -/
theorem base_floor_quant_of_cite (cite : LonelyRunner.LRCUpTo13) (base : Fin 12 → ℤ)
    (hpos : ∀ i, 0 < base i) :
    1 / (400 * ((∑ i, base i : ℤ) : ℚ)) ≤ length (goodRegion2 (List.ofFn base) (1 / 14)) := by
  obtain ⟨t0, ht0⟩ := cite 12 (by norm_num) base (fun i => (hpos i).ne')
  set V : ℤ := ∑ i, base i with hV
  have hVpos : 0 < V := Finset.sum_pos (fun i _ => hpos i) ⟨0, Finset.mem_univ 0⟩
  have hVQ : (0 : ℚ) < (V : ℚ) := by exact_mod_cast hVpos
  have hV1 : (1 : ℚ) ≤ (V : ℚ) := by
    have : (1 : ℤ) ≤ V := by omega
    exact_mod_cast this
  have hVR : (0 : ℝ) < (V : ℝ) := by exact_mod_cast hVpos
  have hle : ∀ i, base i ≤ V := fun i =>
    Finset.single_le_sum (fun j _ => (hpos j).le) (Finset.mem_univ i)
  set δ : ℝ := 1 / (400 * (V : ℝ)) with hδ
  have hδpos : (0 : ℝ) < δ := by rw [hδ]; positivity
  obtain ⟨x, hx1, hx2⟩ := exists_rat_btwn (show t0 - δ < t0 + δ by linarith)
  have hxt : |(x : ℝ) - t0| < δ := by rw [abs_lt]; constructor <;> linarith
  set ρ : ℚ := 1 / (400 * (V : ℚ)) with hρ
  have hρpos : (0 : ℚ) < ρ := by rw [hρ]; positivity
  set x' : ℚ := Int.fract x with hx'
  have hx'0 : 0 ≤ x' := Int.fract_nonneg x
  have hx'1 : x' < 1 := Int.fract_lt_one x
  -- every rational `y` within `ρ` of `x'` is strictly `1/14`-good for the base
  have hgoodY : ∀ y : ℚ, |y - x'| ≤ ρ → ∀ i : Fin 12, ∀ m : ℤ,
      (1 : ℚ) / 14 < |(base i : ℚ) * y - m| := by
    intro y hy i m
    have hbR : (0 : ℝ) < (base i : ℝ) := by exact_mod_cast hpos i
    have hbleR : (base i : ℝ) ≤ (V : ℝ) := by exact_mod_cast hle i
    set m'' : ℤ := m + base i * ⌊x⌋ with hm''
    have hcit := ht0 i m''
    have hcast13 : (1 : ℝ) / (13 : ℕ) = 1 / 13 := by norm_num
    rw [hcast13] at hcit
    -- the two-perturbation triangle bound
    have hkey : (base i : ℝ) * t0 - (m'' : ℝ)
        = ((base i : ℝ) * (y : ℝ) - (m : ℝ)) + (base i : ℝ) * (t0 - (x : ℝ))
            + (base i : ℝ) * ((x' : ℝ) - (y : ℝ)) := by
      have hfr : ((x' : ℚ) : ℝ) = (x : ℝ) - ((⌊x⌋ : ℤ) : ℝ) := by
        rw [hx']; push_cast [Int.fract]; ring
      rw [hm'']; push_cast; rw [hfr]; ring
    have htri : |(base i : ℝ) * t0 - (m'' : ℝ)|
        ≤ |(base i : ℝ) * (y : ℝ) - (m : ℝ)|
          + (base i : ℝ) * |t0 - (x : ℝ)| + (base i : ℝ) * |(x' : ℝ) - (y : ℝ)| := by
      rw [hkey]
      calc |((base i : ℝ) * (y : ℝ) - (m : ℝ)) + (base i : ℝ) * (t0 - (x : ℝ))
              + (base i : ℝ) * ((x' : ℝ) - (y : ℝ))|
          ≤ |((base i : ℝ) * (y : ℝ) - (m : ℝ)) + (base i : ℝ) * (t0 - (x : ℝ))|
              + |(base i : ℝ) * ((x' : ℝ) - (y : ℝ))| := abs_add_le _ _
        _ ≤ (|(base i : ℝ) * (y : ℝ) - (m : ℝ)| + |(base i : ℝ) * (t0 - (x : ℝ))|)
              + |(base i : ℝ) * ((x' : ℝ) - (y : ℝ))| := by
            have := abs_add_le ((base i : ℝ) * (y : ℝ) - (m : ℝ)) ((base i : ℝ) * (t0 - (x : ℝ)))
            linarith
        _ = |(base i : ℝ) * (y : ℝ) - (m : ℝ)|
              + (base i : ℝ) * |t0 - (x : ℝ)| + (base i : ℝ) * |(x' : ℝ) - (y : ℝ)| := by
            rw [abs_mul, abs_mul, abs_of_pos hbR]
    -- perturbation bounds: each ≤ 1/400
    have hb1 : (base i : ℝ) * |t0 - (x : ℝ)| < 1 / 400 := by
      calc (base i : ℝ) * |t0 - (x : ℝ)| ≤ (V : ℝ) * |t0 - (x : ℝ)| :=
            mul_le_mul_of_nonneg_right hbleR (abs_nonneg _)
        _ < (V : ℝ) * δ := by
            rw [abs_sub_comm]; exact mul_lt_mul_of_pos_left hxt hVR
        _ = 1 / 400 := by rw [hδ]; field_simp
    have hb2 : (base i : ℝ) * |(x' : ℝ) - (y : ℝ)| ≤ 1 / 400 := by
      have hyR : |((x' : ℚ) : ℝ) - ((y : ℚ) : ℝ)| ≤ ((ρ : ℚ) : ℝ) := by
        rw [abs_sub_comm]
        have : |((y : ℚ) : ℝ) - ((x' : ℚ) : ℝ)| = (|y - x'| : ℚ) := by
          push_cast; ring_nf
        rw [this]; exact_mod_cast hy
      calc (base i : ℝ) * |(x' : ℝ) - (y : ℝ)| ≤ (V : ℝ) * ((ρ : ℚ) : ℝ) :=
            mul_le_mul hbleR hyR (abs_nonneg _) hVR.le
        _ = 1 / 400 := by rw [hρ]; push_cast; field_simp
    have hgR : (1 : ℝ) / 14 < |(base i : ℝ) * (y : ℝ) - (m : ℝ)| := by
      have h14 : (1 : ℝ) / 13 - 1 / 400 - 1 / 400 > 1 / 14 := by norm_num
      linarith
    have hcastabs : |(base i : ℝ) * (y : ℝ) - (m : ℝ)|
        = ((|(base i : ℚ) * y - m| : ℚ) : ℝ) := by push_cast; ring_nf
    rw [hcastabs] at hgR
    have : ((1 : ℚ) / 14 : ℝ) < ((|(base i : ℚ) * y - m| : ℚ) : ℝ) := by push_cast; linarith
    exact_mod_cast this
  -- the base speeds are positive as a list, and the interval fits in [0,1) by casing
  have hposL : ∀ s ∈ List.ofFn base, 0 < s := by
    intro s hs; obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs; exact hpos i
  have hρval : (x' + ρ) - x' = ρ := by ring
  rcases le_or_gt x' (1 / 2) with hhalf | hhalf
  · -- right interval [x', x' + ρ] ⊆ [0,1)
    have hsafe : ∀ z : ℚ, x' ≤ z → z < x' + ρ → 0 ≤ z ∧ z < 1 ∧
        ∀ s ∈ List.ofFn base, ∀ m : ℤ, (1 : ℚ) / 14 < |(s : ℚ) * z - m| := by
      intro z hz1 hz2
      have hρsmall : ρ ≤ 1 / 400 := by
        rw [hρ]
        exact one_div_le_one_div_of_le (by norm_num) (by nlinarith [hV1])
      refine ⟨by linarith, by linarith, ?_⟩
      intro s hs m
      obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
      exact hgoodY z (by rw [abs_le]; constructor <;> linarith) i m
    have hlen := length_ge_of_safe_interval hposL (by norm_num : (0:ℚ) ≤ 1/14) hsafe
    rw [hρval] at hlen
    calc 1 / (400 * ((∑ i, base i : ℤ) : ℚ)) = ρ := by rw [hρ, hV]
      _ ≤ length (goodRegion2 (List.ofFn base) (1 / 14)) := hlen
  · -- left interval [x' - ρ, x'] ⊆ [0,1)
    have hsafe : ∀ z : ℚ, x' - ρ ≤ z → z < x' → 0 ≤ z ∧ z < 1 ∧
        ∀ s ∈ List.ofFn base, ∀ m : ℤ, (1 : ℚ) / 14 < |(s : ℚ) * z - m| := by
      intro z hz1 hz2
      have hρsmall : ρ ≤ 1 / 400 := by
        rw [hρ]
        exact one_div_le_one_div_of_le (by norm_num) (by nlinarith [hV1])
      refine ⟨by linarith, by linarith, ?_⟩
      intro s hs m
      obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
      exact hgoodY z (by rw [abs_le]; constructor <;> linarith) i m
    have hlen := length_ge_of_safe_interval hposL (by norm_num : (0:ℚ) ≤ 1/14) hsafe
    have hρval2 : x' - (x' - ρ) = ρ := by ring
    rw [hρval2] at hlen
    calc 1 / (400 * ((∑ i, base i : ℤ) : ℚ)) = ρ := by rw [hρ, hV]
      _ ≤ length (goodRegion2 (List.ofFn base) (1 / 14)) := hlen

end RatIntervals
end LonelyRunner
