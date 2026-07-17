/- LRCLacunaryNest.lean -- opus-2026-07-17-S336 (HYP-7190), rewrite 2.
   THE NESTED-GAP LACUNARY BRANCH: the first UNIVERSAL (quantifier-level)
   loneliness branch beyond the corpus's per-certificate gates.

     `lonely_of_pos_lacunary` : every positive speed family whose consecutive
     ratios are ≥ 7/3 (as 7·v i ≤ 3·v (i+1)) is lonely at 1/14.

   PROOF (pure ℚ floor arithmetic): maintain a rational interval [a, b]
   inside a SAFE GAP of every processed comb (`arcSafe (1/14) w 0`).  Step:
   if w·(b − a) ≥ 2, the period [⌈wa⌉, ⌈wa⌉ + 1] fits in [wa, wb]; its gap
   pulls back to a subinterval of length (6/7)/w.  7/3 = 2/(6/7) sustains
   the recursion (`NestOK` threads the window lengths).  Any point of the
   final interval is a rational witness; `norm_ge_of_arcSafe` +
   `lonely_of_norm_forall` finish.

   Paper: THM-928(A′) (7/3 sharpens the measure cascade's R ≥ 15; TIGHT —
   boundary chains land distance exactly 1/14, nested_gap_verify_opus_S336).
   Kernel-pure, no native_decide. -/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCWitnessCert
import TournamentH7.LRC14CertRoute

namespace LonelyRunner
namespace LRC14
open WitnessCert

/-- One nested-gap step: a window of `w`-length at least 2 contains a full
safe gap of the comb of modulus `w`, of length `(6/7)/w`. -/
theorem nested_gap_step (w : ℤ) (hw : 0 < w) (a b : ℚ)
    (hlen : 2 ≤ (w : ℚ) * (b - a)) :
    ∃ a' b' : ℚ, a ≤ a' ∧ b' ≤ b ∧ b' - a' = (6/7) / (w : ℚ) ∧
      arcSafe (1/14) w 0 a' b' := by
  have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  set j : ℤ := ⌈(w : ℚ) * a⌉ with hjdef
  have hja : (w : ℚ) * a ≤ (j : ℚ) := Int.le_ceil _
  have hjb : (j : ℚ) < (w : ℚ) * a + 1 := Int.ceil_lt_add_one _
  refine ⟨((j : ℚ) + 1/14) / w, ((j : ℚ) + 13/14) / w, ?_, ?_, ?_, ?_⟩
  · rw [le_div_iff₀ hwQ]
    have hcomm : a * (w : ℚ) = (w : ℚ) * a := mul_comm _ _
    linarith
  · rw [div_le_iff₀ hwQ]
    have hcomm : b * (w : ℚ) = (w : ℚ) * b := mul_comm _ _
    have h3 : (w : ℚ) * a + 2 ≤ (w : ℚ) * b := by linarith
    linarith
  · field_simp
    ring
  · have hwa : (w : ℚ) * (((j : ℚ) + 1/14) / w) = (j : ℚ) + 1/14 := by
      field_simp
    have hwb : (w : ℚ) * (((j : ℚ) + 13/14) / w) = (j : ℚ) + 13/14 := by
      field_simp
    have hfl : ⌊(j : ℚ) + 1/14⌋ = j := by
      rw [Int.floor_intCast_add]
      norm_num
    refine ⟨?_, ?_, ?_⟩
    · simp only [sub_zero]
      rw [hwa, hwb]
      linarith
    · simp only [sub_zero]
      rw [hwa, hfl]
      push_cast
      linarith
    · simp only [sub_zero]
      rw [hwa, hwb, hfl]
      push_cast
      linarith

/-- The recursive window ledger: the head needs `w`-length ≥ 2 of the current
window; the tail runs in the produced gap of length `(6/7)/w`. -/
def NestOK : List ℤ → ℚ → Prop
  | [], _ => True
  | w :: rest, len => 2 ≤ (w : ℚ) * len ∧ NestOK rest ((6/7) / (w : ℚ))

/-- The nested-gap chain: a `NestOK` list of moduli admits a common nested
window; each modulus is safe on an interval containing it. -/
theorem nested_gap_chain : ∀ (ws : List ℤ), (∀ w ∈ ws, 0 < w) →
    ∀ a b : ℚ, a ≤ b → NestOK ws (b - a) →
    ∃ aF bF : ℚ, a ≤ aF ∧ bF ≤ b ∧ aF ≤ bF ∧
      ∀ w ∈ ws, ∃ lo hi : ℚ, lo ≤ aF ∧ bF ≤ hi ∧ arcSafe (1/14) w 0 lo hi := by
  intro ws
  induction ws with
  | nil =>
      intro _ a b hab _
      exact ⟨a, b, le_refl _, le_refl _, hab, fun w hw => absurd hw (by simp)⟩
  | cons w rest ih =>
      intro hpos a b hab hnest
      obtain ⟨hlen, hnest'⟩ := hnest
      have hw : 0 < w := hpos w (by simp)
      have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
      obtain ⟨a', b', ha', hb', hlen', hsafe⟩ := nested_gap_step w hw a b hlen
      have hab' : a' ≤ b' := by
        have hpos67 : (0 : ℚ) < (6/7) / (w : ℚ) := by positivity
        linarith
      have hnest'' : NestOK rest (b' - a') := by
        rw [hlen']
        exact hnest'
      obtain ⟨aF, bF, haF, hbF, habF, hall⟩ :=
        ih (fun x hx => hpos x (by simp [hx])) a' b' hab' hnest''
      refine ⟨aF, bF, le_trans ha' haF, le_trans hbF hb', habF, ?_⟩
      intro x hx
      rcases List.mem_cons.mp hx with rfl | hx'
      · exact ⟨a', b', haF, hbF, hsafe⟩
      · exact hall x hx'

/-- The ratio arithmetic: a 7/3-ratio step sustains the window ledger. -/
theorem nest_step_arith (u v : ℤ) (hu : 0 < u) (hr : 7 * u ≤ 3 * v) :
    2 ≤ (v : ℚ) * ((6/7) / (u : ℚ)) := by
  have huQ : (0 : ℚ) < (u : ℚ) := by exact_mod_cast hu
  have hrQ : (7 : ℚ) * u ≤ 3 * v := by exact_mod_cast hr
  rw [mul_div_assoc', le_div_iff₀ huQ]
  nlinarith

/-- The original whole-period nested proof: positive speeds, first at least
two, and every consecutive ratio at least `7/3`. -/
theorem lonely_of_pos_lacunary_of_two_le (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (h0 : 2 ≤ v 0)
    (hchain : ∀ i : Fin 12, 7 * v i.castSucc ≤ 3 * v i.succ) :
    ∃ t : ℝ, Lonely 14 v t := by
  have harith : ∀ i : Fin 12,
      2 ≤ ((v i.succ : ℤ) : ℚ) * ((6/7) / ((v i.castSucc : ℤ) : ℚ)) :=
    fun i => nest_step_arith _ _ (hpos _) (hchain i)
  set ws : List ℤ := List.ofFn v with hws
  have hpos' : ∀ w ∈ ws, 0 < w := by
    intro w hw
    rw [hws, List.mem_ofFn] at hw
    obtain ⟨i, rfl⟩ := hw
    exact hpos i
  have hnest : NestOK ws (1 - 0) := by
    have h00 : 2 ≤ ((v 0 : ℤ) : ℚ) * ((1 : ℚ) - 0) := by
      have : (2 : ℚ) ≤ ((v 0 : ℤ) : ℚ) := by exact_mod_cast h0
      linarith
    rw [hws]
    simp only [List.ofFn_succ, List.ofFn_zero, Fin.succ]
    exact ⟨h00, harith ⟨0, by omega⟩, harith ⟨1, by omega⟩,
      harith ⟨2, by omega⟩, harith ⟨3, by omega⟩, harith ⟨4, by omega⟩,
      harith ⟨5, by omega⟩, harith ⟨6, by omega⟩, harith ⟨7, by omega⟩,
      harith ⟨8, by omega⟩, harith ⟨9, by omega⟩, harith ⟨10, by omega⟩,
      harith ⟨11, by omega⟩, trivial⟩
  obtain ⟨aF, bF, _, _, habF, hall⟩ :=
    nested_gap_chain ws hpos' 0 1 (by norm_num) hnest
  refine ⟨((aF : ℚ) : ℝ), lonely_of_norm_forall ?_⟩
  intro i
  obtain ⟨lo, hi, hlo, hhi, hsafe⟩ := hall (v i) (by
    rw [hws, List.mem_ofFn]
    exact ⟨i, rfl⟩)
  have hx1 : ((lo : ℚ) : ℝ) ≤ ((aF : ℚ) : ℝ) := by exact_mod_cast hlo
  have hx2 : ((aF : ℚ) : ℝ) ≤ ((hi : ℚ) : ℝ) := by
    have : (aF : ℚ) ≤ hi := le_trans habF hhi
    exact_mod_cast this
  have := norm_ge_of_arcSafe (le_of_lt (hpos i)) (by norm_num) hsafe hx1 hx2
  simpa using this

/-- The unit head already has the full safe gap `[1/14,13/14]`.  Nesting the
remaining twelve speeds inside that gap removes the artificial `v 0 ≥ 2`
condition from the original whole-period proof. -/
theorem lonely_of_pos_lacunary_of_head_one (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (h0 : v 0 = 1)
    (hchain : ∀ i : Fin 12, 7 * v i.castSucc ≤ 3 * v i.succ) :
    ∃ t : ℝ, Lonely 14 v t := by
  have harith : ∀ i : Fin 12,
      2 ≤ ((v i.succ : ℤ) : ℚ) * ((6/7) / ((v i.castSucc : ℤ) : ℚ)) :=
    fun i => nest_step_arith _ _ (hpos _) (hchain i)
  let tail : Fin 12 → ℤ := fun i => v i.succ
  set ws : List ℤ := List.ofFn tail with hws
  have hpos' : ∀ w ∈ ws, 0 < w := by
    intro w hw
    rw [hws, List.mem_ofFn] at hw
    obtain ⟨i, rfl⟩ := hw
    exact hpos i.succ
  have hfirst :
      2 ≤ ((tail 0 : ℤ) : ℚ) * ((13/14 : ℚ) - 1/14) := by
    rw [show ((13/14 : ℚ) - 1/14) = 6/7 by norm_num]
    simpa [tail, h0] using harith (0 : Fin 12)
  have hnest : NestOK ws ((13/14 : ℚ) - 1/14) := by
    rw [hws]
    simp only [List.ofFn_succ, List.ofFn_zero, Fin.succ]
    exact ⟨hfirst,
      by simpa [tail] using harith (1 : Fin 12),
      by simpa [tail] using harith (2 : Fin 12),
      by simpa [tail] using harith (3 : Fin 12),
      by simpa [tail] using harith (4 : Fin 12),
      by simpa [tail] using harith (5 : Fin 12),
      by simpa [tail] using harith (6 : Fin 12),
      by simpa [tail] using harith (7 : Fin 12),
      by simpa [tail] using harith (8 : Fin 12),
      by simpa [tail] using harith (9 : Fin 12),
      by simpa [tail] using harith (10 : Fin 12),
      by simpa [tail] using harith (11 : Fin 12), trivial⟩
  obtain ⟨aF, bF, haF, hbF, habF, hall⟩ :=
    nested_gap_chain ws hpos' (1/14) (13/14) (by norm_num) hnest
  have hsafeHead : arcSafe (1/14) 1 0 (1/14) (13/14) := by
    norm_num [arcSafe]
  refine ⟨((aF : ℚ) : ℝ), lonely_of_norm_forall ?_⟩
  intro i
  by_cases hi : i = 0
  · subst i
    have hx1 : (((1/14 : ℚ) : ℝ)) ≤ ((aF : ℚ) : ℝ) := by
      exact_mod_cast haF
    have hx2 : ((aF : ℚ) : ℝ) ≤ (((13/14 : ℚ) : ℝ)) := by
      exact_mod_cast (le_trans habF hbF)
    have hsafe := norm_ge_of_arcSafe (s := (1 : ℤ)) (cTar := 0)
      (by norm_num) (by norm_num) hsafeHead hx1 hx2
    simpa [h0] using hsafe
  · have hiPos : 1 ≤ (i : ℕ) := by omega
    let source : Fin 12 := ⟨(i : ℕ) - 1, by omega⟩
    have hsource : source.succ = i := by
      apply Fin.ext
      dsimp [source]
      omega
    obtain ⟨lo, hi', hlo, hhi, hsafe⟩ := hall (tail source) (by
      rw [hws, List.mem_ofFn]
      exact ⟨source, rfl⟩)
    have htail : tail source = v i := by simp [tail, hsource]
    rw [htail] at hsafe
    have hx1 : ((lo : ℚ) : ℝ) ≤ ((aF : ℚ) : ℝ) := by exact_mod_cast hlo
    have hx2 : ((aF : ℚ) : ℝ) ≤ ((hi' : ℚ) : ℝ) := by
      exact_mod_cast (le_trans habF hhi)
    have hsafe' := norm_ge_of_arcSafe (le_of_lt (hpos i)) (by norm_num)
      hsafe hx1 hx2
    simpa using hsafe'

/-- **THE UNIVERSAL LACUNARY BRANCH**: every positive thirteen-speed family
with consecutive ratios at least `7/3` is lonely at `1/14`, with a rational
witness.  Positivity leaves only `v 0 = 1` or `2 ≤ v 0`; the two preceding
nested-gap constructions cover those cases exactly. -/
theorem lonely_of_pos_lacunary (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i)
    (hchain : ∀ i : Fin 12, 7 * v i.castSucc ≤ 3 * v i.succ) :
    ∃ t : ℝ, Lonely 14 v t := by
  by_cases h0 : v 0 = 1
  · exact lonely_of_pos_lacunary_of_head_one v hpos h0 hchain
  · apply lonely_of_pos_lacunary_of_two_le v hpos (by
      have hp0 := hpos 0
      omega) hchain

/-- **THE WINDOW+TAIL GLUE** (opus-S337; the formal core of gap gluing):
a base window certified safe for `base` plus a `NestOK` lacunary tail give
a common rational time at which every modulus of `base ++ tail` is safe.
The base certificate is exactly what census rows / ladder packs / the
cluster-gap lemma provide; the tail is any 7/3-chain above the window. -/
theorem window_tail_glue (base tail : List ℤ) (a b : ℚ) (hab : a ≤ b)
    (hbase : ∀ w ∈ base, arcSafe (1/14) w 0 a b)
    (htailPos : ∀ w ∈ tail, 0 < w)
    (hnest : NestOK tail (b - a)) :
    ∃ t : ℚ, ∀ w ∈ base ++ tail, ∃ lo hi : ℚ,
      lo ≤ t ∧ t ≤ hi ∧ arcSafe (1/14) w 0 lo hi := by
  obtain ⟨aF, bF, haF, hbF, habF, hall⟩ :=
    nested_gap_chain tail htailPos a b hab hnest
  refine ⟨aF, fun w hw => ?_⟩
  rcases List.mem_append.mp hw with hwb | hwt
  · exact ⟨a, b, haF, le_trans habF hbF, hbase w hwb⟩
  · obtain ⟨lo, hi, hlo, hhi, hsafe⟩ := hall w hwt
    exact ⟨lo, hi, hlo, le_trans habF hhi, hsafe⟩

/-- The glue at the real level: the common time keeps every positive modulus
of `base ++ tail` at circle distance ≥ 1/14. -/
theorem norm_glue (base tail : List ℤ) (a b : ℚ) (hab : a ≤ b)
    (hpos : ∀ w ∈ base ++ tail, 0 < w)
    (hbase : ∀ w ∈ base, arcSafe (1/14) w 0 a b)
    (hnest : NestOK tail (b - a)) :
    ∃ t : ℝ, ∀ w ∈ base ++ tail,
      (1 : ℝ) / 14 ≤ ‖(((w : ℝ) * t : ℝ) : UnitAddCircle)‖ := by
  obtain ⟨t, ht⟩ := window_tail_glue base tail a b hab hbase
    (fun w hw => hpos w (List.mem_append.mpr (Or.inr hw))) hnest
  refine ⟨((t : ℚ) : ℝ), fun w hw => ?_⟩
  obtain ⟨lo, hi, hlo, hhi, hsafe⟩ := ht w hw
  have := norm_ge_of_arcSafe (le_of_lt (hpos w hw)) (by norm_num) hsafe
    (x := ((t : ℚ) : ℝ)) (by exact_mod_cast hlo) (by exact_mod_cast hhi)
  simpa using this

/-! ## Axiom audit -/
#print axioms nested_gap_step
#print axioms nested_gap_chain
#print axioms lonely_of_pos_lacunary
#print axioms window_tail_glue
#print axioms norm_glue

end LRC14
end LonelyRunner
