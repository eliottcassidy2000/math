/-
  TournamentH7.LRCRungFloor — the gate-tower lemma stack, formalized
  (death-star-2026-07-19-S59h; THM-1271 L1 + THM-1270 closure lemmas).

  For the canonical-D single-far family F_D(N) = {1..N}\{N-1} ∪ {D(N-1)} with
  binder p = 2D-1 and slack-1 denominator Q = (N+1)D - 1:

  * `cert_dvd` / `rung_core` (L1, integer heart): if (p·a) % Q = D then no
    family element has v·a within D-1 of a multiple of Q.  The whole modular
    argument is ONE divisibility certificate — the ring identity
        v + r(N-1) = (N+1)·((2D-1)·(va-r) - v·((2D-1)a-D)) - Q·(v-2r)   [mod Q]
    — so `linear_combination` does it with no modular-inverse API; the unique
    in-range representative is the DELETED element N-1 (at r = -1), and every
    other case dies by sign/size.

  * `rung_band` / `rung_floor_single` / `rung_floor_witness`: the residue band
    [D, Q-D] and its lift to the real witness t = a/Q putting EVERY family
    element at circle-distance ≥ D/Q — hence M(F_D(N)) ≥ D/((N+1)D-1).

  * `exists_binder_multiplier`: the multiplier exists iff gcd(2D-1, Q) = 1
    (Bézout) — exactly the coprimality the gate's clause N ≢ 1 (mod 2D-1)
    guarantees.

  * `Q_dvd_iff_binder_dvd` + `rung_dead_of_common_factor` (THM-1270): under
    N ≡ 1 (mod ℓ), ℓ ∣ Q ↔ ℓ ∣ 2D-1; and a common factor g ∣ p, g ∣ Q with
    g ∤ D makes (p·a) % Q = D unsolvable — the tower-closure/composite-skip.

  Non-vacuity (MISTAKE-186 discipline): concrete instantiation at
  (N, D, a) = (31, 4, 55), the THM-1285 member 4/127, closes the file.
-/
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace RungFloor

/-! ## The integer core -/

/-- **The divisibility certificate.**  If `Q ∣ (2D-1)·a - D` and
`Q ∣ v·a - r` then `Q ∣ v + r·(N-1)`, where `Q = (N+1)·D - 1`.
Pure `linear_combination`; no inverses. -/
theorem cert_dvd (D N a v r : ℤ)
    (h1 : ((N+1)*D - 1) ∣ ((2*D-1) * a - D))
    (h2 : ((N+1)*D - 1) ∣ (v * a - r)) :
    ((N+1)*D - 1) ∣ (v + r*(N-1)) := by
  obtain ⟨k1, hk1⟩ := h1
  obtain ⟨k2, hk2⟩ := h2
  refine ⟨(N+1)*((2*D-1)*k2 - v*k1) - (v - 2*r), ?_⟩
  linear_combination (N+1)*(2*D-1)*hk2 - (N+1)*v*hk1

/-- **L1 integer core.**  With `3 ≤ D`, `3D-2 < N`, and a multiplier
satisfying `((2D-1)·a) % Q = D`, no family element `v` (either `1 ≤ v ≤ N`
with `v ≠ N-1`, or `v = D(N-1)`) has `v·a ≡ r (mod Q)` for any `|r| ≤ D-1`. -/
theorem rung_core (D N a : ℤ) (hD : 3 ≤ D) (hN : 3*D - 2 < N)
    (ha : ((2*D-1) * a) % ((N+1)*D - 1) = D)
    (v : ℤ) (hv : (1 ≤ v ∧ v ≤ N ∧ v ≠ N-1) ∨ v = D*(N-1))
    (r : ℤ) (hr1 : -(D-1) ≤ r) (hr2 : r ≤ D-1)
    (hcong : ((N+1)*D - 1) ∣ (v * a - r)) : False := by
  have hN1 : (0:ℤ) ≤ N - 1 := by omega
  have hQpos : (0:ℤ) < (N+1)*D - 1 := by nlinarith
  -- Q ∣ p·a − D from the emod hypothesis
  have hdef := Int.emod_def ((2*D-1) * a) ((N+1)*D - 1)
  have hpa : ((N+1)*D - 1) ∣ ((2*D-1) * a - D) :=
    ⟨((2*D-1) * a) / ((N+1)*D - 1), by linear_combination ha - hdef⟩
  obtain ⟨k, hk⟩ := cert_dvd D N a v r hpa hcong
  -- bounds on v
  have hNx : N ≤ D*(N-1) := by
    nlinarith [mul_le_mul_of_nonneg_right hD hN1]
  have hxv : v ≤ D*(N-1) := by
    rcases hv with ⟨_, hvN, _⟩ | hvx
    · omega
    · omega
  have hv1 : 1 ≤ v := by
    rcases hv with ⟨h1, _, _⟩ | hvx
    · exact h1
    · nlinarith [mul_le_mul_of_nonneg_right hD hN1]
  -- product bounds for r(N-1)
  have hrlo : -(D-1)*(N-1) ≤ r*(N-1) :=
    mul_le_mul_of_nonneg_right hr1 hN1
  have hrhi : r*(N-1) ≤ (D-1)*(N-1) :=
    mul_le_mul_of_nonneg_right hr2 hN1
  -- k ∈ {0, 1}
  have hk0 : 0 ≤ k := by
    by_contra hneg
    push Not at hneg
    have hk1' : k ≤ -1 := by omega
    have := mul_le_mul_of_nonneg_left hk1' (le_of_lt hQpos)
    nlinarith
  have hk1 : k ≤ 1 := by
    by_contra hbig
    push Not at hbig
    have hk2' : 2 ≤ k := by omega
    have := mul_le_mul_of_nonneg_left hk2' (le_of_lt hQpos)
    nlinarith
  interval_cases k
  · -- k = 0 : v + r(N-1) = 0
    have hw0 : v + r*(N-1) = 0 := by linarith [hk]
    rcases hv with ⟨hv1', hvN, hvne⟩ | hvx
    · rcases le_or_gt r (-2) with hrle | hrgt
      · -- r ≤ -2 : v = (-r)(N-1) ≥ 2(N-1) > N
        have h2r : (2:ℤ) ≤ -r := by omega
        have := mul_le_mul_of_nonneg_right h2r hN1
        nlinarith
      · rcases le_or_gt 0 r with hr0 | hrneg
        · -- r ≥ 0 : v = -r(N-1) ≤ 0
          have := mul_nonneg hr0 hN1
          nlinarith
        · -- r = -1 : v = N-1, excluded
          have hreq : r = -1 := by omega
          apply hvne
          rw [hreq] at hw0
          linarith
    · -- far element : (D+r)(N-1) = 0
      have hzero : (D + r) * (N - 1) = 0 := by
        rw [hvx] at hw0; linear_combination hw0
      rcases mul_eq_zero.mp hzero with h | h
      · omega
      · omega
  · -- k = 1 : v + r(N-1) = Q
    have hw1 : v + r*(N-1) = (N+1)*D - 1 := by linarith [hk]
    rcases le_or_gt r 0 with hrle | hrpos
    · -- r ≤ 0 : v ≥ Q > D(N-1) ≥ v
      have := mul_nonpos_of_nonpos_of_nonneg hrle hN1
      nlinarith
    · rcases hv with ⟨_, hvN, _⟩ | hvx
      · -- base : v = Q - r(N-1) ≥ Q - (D-1)(N-1) = 2D + N - 2 > N
        nlinarith
      · -- far : r(N-1) = 2D - 1, but r ≥ 1 gives r(N-1) ≥ N-1 > 2D-1
        have hrN : r * (N-1) = 2*D - 1 := by
          rw [hvx] at hw1; linear_combination hw1
        have h1r : (1:ℤ) ≤ r := by omega
        have := mul_le_mul_of_nonneg_right h1r hN1
        nlinarith

/-! ## The band and the real lift -/

/-- **The residue band.**  Every family element's residue lies in `[D, Q-D]`. -/
theorem rung_band (D N a : ℤ) (hD : 3 ≤ D) (hN : 3*D - 2 < N)
    (ha : ((2*D-1) * a) % ((N+1)*D - 1) = D)
    (v : ℤ) (hv : (1 ≤ v ∧ v ≤ N ∧ v ≠ N-1) ∨ v = D*(N-1)) :
    D ≤ (v * a) % ((N+1)*D - 1) ∧
      (v * a) % ((N+1)*D - 1) ≤ ((N+1)*D - 1) - D := by
  have hQpos : (0:ℤ) < (N+1)*D - 1 := by nlinarith
  have hdef := Int.emod_def (v * a) ((N+1)*D - 1)
  have hs0 : 0 ≤ (v * a) % ((N+1)*D - 1) := Int.emod_nonneg _ (by omega)
  have hsQ : (v * a) % ((N+1)*D - 1) < (N+1)*D - 1 := Int.emod_lt_of_pos _ hQpos
  have hdvd_s : ((N+1)*D - 1) ∣ (v * a - (v * a) % ((N+1)*D - 1)) :=
    ⟨(v * a) / ((N+1)*D - 1), by linear_combination -hdef⟩
  constructor
  · by_contra hcon
    push Not at hcon
    exact rung_core D N a hD hN ha v hv ((v * a) % ((N+1)*D - 1))
      (by omega) (by omega) hdvd_s
  · by_contra hcon
    push Not at hcon
    obtain ⟨c, hc⟩ := hdvd_s
    exact rung_core D N a hD hN ha v hv ((v * a) % ((N+1)*D - 1) - ((N+1)*D - 1))
      (by omega) (by omega) ⟨c + 1, by linear_combination hc⟩

/-- **Single-element real lift** (variable-modulus `sieve19_single`):
a residue in the band `[D, Q-D]` puts the element at circle-distance `≥ D/Q`
at the witness time `t = a/Q`. -/
theorem rung_floor_single (Q D v a : ℤ) (hQpos : 0 < Q)
    (hband : D ≤ (v * a) % Q ∧ (v * a) % Q ≤ Q - D) :
    ∀ m : ℤ, (D : ℝ) / (Q : ℝ) ≤ |(v : ℝ) * ((a : ℝ) / (Q : ℝ)) - m| := by
  intro m
  have hQR : (0 : ℝ) < (Q : ℝ) := by exact_mod_cast hQpos
  have hQne : (Q : ℝ) ≠ 0 := ne_of_gt hQR
  have hreal : (v : ℝ) * ((a : ℝ) / (Q : ℝ)) - m
      = ((v * a - Q * m : ℤ) : ℝ) / (Q : ℝ) := by
    field_simp
    push_cast
    ring
  have hdef := Int.emod_def (v * a) Q
  have hs0 : 0 ≤ (v * a) % Q := Int.emod_nonneg _ (by omega)
  have hint : D ≤ |v * a - Q * m| := by
    rcases le_or_gt D 0 with hD0 | hDpos
    · exact le_trans hD0 (abs_nonneg _)
    have hval : v * a - Q * m
        = Q * ((v * a) / Q - m) + (v * a) % Q := by linear_combination -hdef
    rcases le_or_gt 0 ((v * a) / Q - m) with hk | hk
    · have hnn : 0 ≤ Q * ((v * a) / Q - m) := mul_nonneg (le_of_lt hQpos) hk
      rw [abs_of_nonneg (by linarith [hband.1, hs0])]
      linarith [hband.1]
    · have hle : (v * a) / Q - m ≤ -1 := by omega
      have hmul := mul_le_mul_of_nonneg_left hle (le_of_lt hQpos)
      have hneg : v * a - Q * m ≤ -D := by nlinarith [hband.2]
      rw [abs_of_nonpos (by linarith)]
      linarith
  have hintR : (D : ℝ) ≤ |((v * a - Q * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast hint
  rw [hreal, abs_div, abs_of_pos hQR]
  gcongr

/-- **L1, the rung-floor witness (family form).**  With a valid binder
multiplier, EVERY element of `F_D(N)` is at circle-distance `≥ D/Q` at
`t = a/Q` — hence `M(F_D(N)) ≥ D/((N+1)D-1)`. -/
theorem rung_floor_witness (D N a : ℤ) (hD : 3 ≤ D) (hN : 3*D - 2 < N)
    (ha : ((2*D-1) * a) % ((N+1)*D - 1) = D)
    (v : ℤ) (hv : (1 ≤ v ∧ v ≤ N ∧ v ≠ N-1) ∨ v = D*(N-1)) :
    ∀ m : ℤ, (D : ℝ) / (((N+1)*D - 1 : ℤ) : ℝ)
      ≤ |(v : ℝ) * ((a : ℝ) / (((N+1)*D - 1 : ℤ) : ℝ)) - m| :=
  rung_floor_single ((N+1)*D - 1) D v a (by nlinarith)
    (rung_band D N a hD hN ha v hv)

/-! ## Existence of the binder multiplier (Bézout) -/

/-- The witness multiplier exists whenever `gcd(2D-1, Q) = 1`. -/
theorem exists_binder_multiplier (D N : ℤ) (hD : 3 ≤ D) (hN : 3*D - 2 < N)
    (hco : IsCoprime (2*D-1) ((N+1)*D - 1)) :
    ∃ a : ℤ, ((2*D-1) * a) % ((N+1)*D - 1) = D := by
  have hQpos : (0:ℤ) < (N+1)*D - 1 := by nlinarith
  have hDQ : D < (N+1)*D - 1 := by nlinarith
  obtain ⟨u, w, huw⟩ := hco
  refine ⟨u * D, ?_⟩
  have hkey : (2*D-1) * (u * D) = D + ((N+1)*D - 1) * (-(w * D)) := by
    linear_combination D * huw
  rw [hkey, Int.add_mul_emod_self_left]
  exact Int.emod_eq_of_lt (by omega) hDQ

/-! ## The tower-closure / composite-skip core (THM-1270) -/

/-- Under `N ≡ 1 (mod ℓ)`: `ℓ ∣ Q ↔ ℓ ∣ 2D-1` — the slack-1 denominator is
congruent to the binder. -/
theorem Q_dvd_iff_binder_dvd (D N ℓ : ℤ) (hN : ℓ ∣ (N - 1)) :
    (ℓ ∣ ((N+1)*D - 1)) ↔ (ℓ ∣ (2*D - 1)) := by
  obtain ⟨c, hc⟩ := hN
  constructor
  · rintro ⟨k, hk⟩
    exact ⟨k - c*D, by linear_combination hk - D * hc⟩
  · rintro ⟨k, hk⟩
    exact ⟨k + c*D, by linear_combination hk + D * hc⟩

/-- **The dead-rung lemma** (tower closure / composite skip).  A common factor
`g` of the binder and `Q` with `g ∤ D` makes the binder congruence
unsolvable: no multiplier `a` has `((2D-1)·a) % Q = D`. -/
theorem rung_dead_of_common_factor (D Q g : ℤ)
    (hgp : g ∣ (2*D-1)) (hgQ : g ∣ Q) (hgD : ¬ g ∣ D) :
    ∀ a : ℤ, ((2*D-1) * a) % Q ≠ D := by
  intro a ha
  apply hgD
  rw [← ha, Int.emod_def]
  exact dvd_sub (hgp.mul_right a) (hgQ.mul_right _)

/-! ## L3: the packing seal (disjoint-balls counting, no sorting) -/

/-- **L3, the packing seal.**  If every index difference `d ∈ [1, N]` has its
residue `(d·a) % Q` in the band `[c, Q-c]` (i.e. circle-distance `≥ c/Q` —
this is `min(c, e) = c` when the ghost error is no smaller than the base
clearance), then `(N+1)·c ≤ Q`: the `N+1` orbit points pack the circle.
Proof: `(u, i) ↦ (u·a + i) % Q` is injective on `[0,N] × [0,c)` — a collision
would put some `(d·a) % Q` within `c-1` of `0` — so `(N+1)·c ≤ Q` by counting. -/
theorem packing_seal (N c Q a : ℤ) (hc : 1 ≤ c) (hN : 1 ≤ N) (hQ : 0 < Q)
    (hband : ∀ d : ℤ, 1 ≤ d → d ≤ N →
      c ≤ (d * a) % Q ∧ (d * a) % Q ≤ Q - c) :
    (N + 1) * c ≤ Q := by
  -- the injection grid, in ℕ
  have hNn : ∃ n : ℕ, (n : ℤ) = N := ⟨N.toNat, Int.toNat_of_nonneg (by omega)⟩
  obtain ⟨n, hn⟩ := hNn
  have hcn : ∃ γ : ℕ, (γ : ℤ) = c := ⟨c.toNat, Int.toNat_of_nonneg (by omega)⟩
  obtain ⟨γ, hγ⟩ := hcn
  -- the map
  set f : ℕ × ℕ → ℤ := fun p => ((p.1 : ℤ) * a + (p.2 : ℤ)) % Q with hf
  have hmaps : ∀ p ∈ Finset.range (n+1) ×ˢ Finset.range γ,
      f p ∈ Finset.Ico (0:ℤ) Q := by
    intro p _
    simp only [Finset.mem_Ico]
    exact ⟨Int.emod_nonneg _ (by omega), Int.emod_lt_of_pos _ hQ⟩
  have hinj : Set.InjOn f ↑(Finset.range (n+1) ×ˢ Finset.range γ) := by
    intro p hp q hq hpq
    simp only [Finset.mem_coe, Finset.mem_product, Finset.mem_range] at hp hq
    obtain ⟨hp1, hp2⟩ := hp
    obtain ⟨hq1, hq2⟩ := hq
    -- from equal emods, Q ∣ difference
    have hdvd : Q ∣ ((p.1 : ℤ) * a + (p.2 : ℤ)) - ((q.1 : ℤ) * a + (q.2 : ℤ)) := by
      have h1 := Int.emod_def ((p.1 : ℤ) * a + (p.2 : ℤ)) Q
      have h2 := Int.emod_def ((q.1 : ℤ) * a + (q.2 : ℤ)) Q
      refine ⟨((p.1 : ℤ) * a + (p.2 : ℤ)) / Q - ((q.1 : ℤ) * a + (q.2 : ℤ)) / Q, ?_⟩
      have hfeq : f p = f q := hpq
      simp only [hf] at hfeq
      linear_combination h2 - h1 + hfeq
    -- WLOG p.1 ≥ q.1; kill p.1 ≠ q.1 by the band
    have key : ∀ u₁ u₂ i₁ i₂ : ℕ, u₁ ≤ n → u₂ ≤ n → i₁ < γ → i₂ < γ →
        u₂ < u₁ → ¬ Q ∣ ((u₁ : ℤ) * a + (i₁ : ℤ)) - ((u₂ : ℤ) * a + (i₂ : ℤ)) := by
      intro u₁ u₂ i₁ i₂ hu₁ hu₂ hi₁ hi₂ hlt hdvd'
      set d : ℤ := (u₁ : ℤ) - (u₂ : ℤ) with hd
      have hd1 : 1 ≤ d := by simp only [hd]; omega
      have hdN : d ≤ N := by
        simp only [hd]
        have : (u₁ : ℤ) ≤ (n : ℤ) := by exact_mod_cast hu₁
        omega
      obtain ⟨hb1, hb2⟩ := hband d hd1 hdN
      -- (d·a) % Q = (i₂ - i₁) % Q
      have hcong : Q ∣ d * a - ((i₂ : ℤ) - (i₁ : ℤ)) := by
        obtain ⟨k, hk⟩ := hdvd'
        exact ⟨k, by linear_combination hk⟩
      have hii : -(c-1) ≤ (i₂ : ℤ) - (i₁ : ℤ) ∧ (i₂ : ℤ) - (i₁ : ℤ) ≤ c-1 := by
        constructor
        · have : (i₁ : ℤ) < c := by rw [← hγ]; exact_mod_cast hi₁
          omega
        · have : (i₂ : ℤ) < c := by rw [← hγ]; exact_mod_cast hi₂
          omega
      -- so (d·a) % Q is within c-1 of 0 or Q — contradicting the band
      obtain ⟨k, hk⟩ := hcong
      have hdaQ : d * a % Q = d * a - Q * (d * a / Q) := Int.emod_def _ _
      -- d·a = Q·k + (i₂ - i₁); compare with the band via cases on the residue
      rcases le_or_gt 0 ((i₂ : ℤ) - (i₁ : ℤ)) with hpos | hneg
      · -- residue is i₂ - i₁ itself: < c
        have : d * a % Q = (i₂ : ℤ) - (i₁ : ℤ) := by
          rw [show d * a = ((i₂ : ℤ) - (i₁ : ℤ)) + Q * k by linarith [hk],
            Int.add_mul_emod_self_left]
          exact Int.emod_eq_of_lt hpos (by omega)
        omega
      · -- residue is Q + (i₂ - i₁): > Q - c
        have : d * a % Q = Q + ((i₂ : ℤ) - (i₁ : ℤ)) := by
          rw [show d * a = (Q + ((i₂ : ℤ) - (i₁ : ℤ))) + Q * (k - 1) by
            linarith [hk], Int.add_mul_emod_self_left]
          exact Int.emod_eq_of_lt (by omega) (by omega)
        omega
    rcases lt_trichotomy p.1 q.1 with hlt | heq | hgt
    · exact absurd (dvd_neg.mpr hdvd) (by
        have := key q.1 p.1 q.2 p.2 (by omega) (by omega) hq2 hp2 hlt
        simpa [neg_sub] using this)
    · -- p.1 = q.1 ⟹ Q ∣ i₁ - i₂ small ⟹ i₁ = i₂
      have hi : Q ∣ ((p.2 : ℤ) - (q.2 : ℤ)) := by
        obtain ⟨k, hk⟩ := hdvd
        refine ⟨k, ?_⟩
        have : (p.1 : ℤ) = (q.1 : ℤ) := by exact_mod_cast heq
        linear_combination hk - a * this
      have hp2' : (p.2 : ℤ) < c := by rw [← hγ]; exact_mod_cast hp2
      have hq2' : (q.2 : ℤ) < c := by rw [← hγ]; exact_mod_cast hq2
      have h2c : 2 * c ≤ Q := by
        obtain ⟨hb1, hb2⟩ := hband 1 le_rfl (by omega)
        omega
      have hii0 : (p.2 : ℤ) = (q.2 : ℤ) := by
        obtain ⟨k, hk⟩ := hi
        have hk0 : k = 0 := by
          by_contra hkne
          rcases lt_or_gt_of_ne hkne with hkneg | hkpos
          · have : Q * k ≤ Q * (-1) :=
              mul_le_mul_of_nonneg_left (by omega) (le_of_lt hQ)
            omega
          · have : Q * 1 ≤ Q * k :=
              mul_le_mul_of_nonneg_left (by omega) (le_of_lt hQ)
            omega
        rw [hk0, mul_zero] at hk
        omega
      have : p.1 = q.1 := heq
      have : p.2 = q.2 := by exact_mod_cast hii0
      exact Prod.ext heq this
    · exact absurd hdvd (key p.1 q.1 p.2 q.2 (by omega) (by omega) hp2 hq2 hgt)
  -- the count
  have hcard := Finset.card_le_card_of_injOn f hmaps hinj
  rw [Finset.card_product, Finset.card_range, Finset.card_range] at hcard
  have hIco : (Finset.Ico (0:ℤ) Q).card = Q.toNat := by
    rw [Int.card_Ico]; simp
  rw [hIco] at hcard
  have : ((n+1) * γ : ℤ) ≤ (Q.toNat : ℤ) := by exact_mod_cast hcard
  rw [Int.toNat_of_nonneg (le_of_lt hQ)] at this
  calc (N + 1) * c = ((n : ℤ) + 1) * (γ : ℤ) := by rw [hn, hγ]
    _ ≤ Q := by exact_mod_cast this

/-! ## L2: the small-moduli seal (Dirichlet pigeonhole + parity reflection) -/

/-- Low moduli `q ≤ N`: some family element is a multiple of `q` outright. -/
theorem small_moduli_seal_low (D N q a : ℤ) (hq1 : 2 ≤ q) (hq2 : q ≤ N) :
    ∃ v : ℤ, ((1 ≤ v ∧ v ≤ N ∧ v ≠ N-1) ∨ v = D*(N-1)) ∧ (v * a) % q = 0 := by
  rcases eq_or_ne q (N-1) with hq | hq
  · exact ⟨D*(N-1), Or.inr rfl, by
      rw [← hq]
      exact Int.emod_eq_zero_of_dvd ⟨D*a, by ring⟩⟩
  · exact ⟨q, Or.inl ⟨by omega, by omega, hq⟩, by
      exact Int.emod_eq_zero_of_dvd ⟨a, rfl⟩⟩

/-- **Dirichlet step.**  For `0 < q`, some `u ∈ [1, q/2]` has `(u·a) % q`
within `1` of `0` or `q` (boxes of width two force a collision). -/
theorem dirichlet_step (q a : ℤ) (hq : 4 ≤ q) :
    ∃ u : ℤ, 1 ≤ u ∧ u ≤ q/2 ∧
      ((u * a) % q ≤ 1 ∨ q - 1 ≤ (u * a) % q) := by
  by_contra hcon
  push Not at hcon
  -- every u in [1, q/2] has residue in the middle [2, q-2]
  have hmid : ∀ u : ℤ, 1 ≤ u → u ≤ q/2 →
      2 ≤ (u * a) % q ∧ (u * a) % q ≤ q - 2 := by
    intro u h1 h2
    obtain ⟨hlo, hhi⟩ := hcon u h1 h2
    have h0 : 0 ≤ (u * a) % q := Int.emod_nonneg _ (by omega)
    have hlt : (u * a) % q < q := Int.emod_lt_of_pos _ (by omega)
    omega
  -- box map into [1, (q-2)/2]; source [1, q/2] is strictly bigger
  set K : ℤ := q/2 with hK
  have hKpos : 1 ≤ K := by omega
  set src : Finset ℤ := Finset.Icc 1 K with hsrc
  set tgt : Finset ℤ := Finset.Icc 1 ((q-2)/2) with htgt
  have hmapsto : ∀ u ∈ src, ((u * a) % q) / 2 ∈ tgt := by
    intro u hu
    simp only [hsrc, Finset.mem_Icc] at hu
    obtain ⟨hm1, hm2⟩ := hmid u hu.1 hu.2
    simp only [htgt, Finset.mem_Icc]
    omega
  have hcard : tgt.card < src.card := by
    simp only [hsrc, htgt, Int.card_Icc]
    omega
  obtain ⟨u₁, hu₁, u₂, hu₂, hne, heq⟩ :=
    Finset.exists_ne_map_eq_of_card_lt_of_maps_to hcard hmapsto
  simp only [hsrc, Finset.mem_Icc] at hu₁ hu₂
  -- order them
  rcases hne.lt_or_gt with hlt | hlt
  case _ =>
    -- u₁ < u₂ : use d = u₂ - u₁
    obtain ⟨hm1, hm2⟩ := hmid u₁ hu₁.1 hu₁.2
    obtain ⟨hm1', hm2'⟩ := hmid u₂ hu₂.1 hu₂.2
    have hd1 : 1 ≤ u₂ - u₁ := by omega
    have hd2 : u₂ - u₁ ≤ K := by omega
    obtain ⟨hdm1, hdm2⟩ := hmid _ hd1 hd2
    -- (d·a) % q = the residue difference (mod q), which is within 1 of 0 or q
    have hdvd : q ∣ (u₂ - u₁) * a - ((u₂ * a) % q - (u₁ * a) % q) := by
      have h1 := Int.emod_def (u₁ * a) q
      have h2 := Int.emod_def (u₂ * a) q
      exact ⟨u₂ * a / q - u₁ * a / q, by linear_combination h1 - h2⟩
    obtain ⟨k, hk⟩ := hdvd
    rcases le_or_gt 0 ((u₂ * a) % q - (u₁ * a) % q) with hpos | hneg
    · have : ((u₂ - u₁) * a) % q = (u₂ * a) % q - (u₁ * a) % q := by
        rw [show (u₂ - u₁) * a
            = ((u₂ * a) % q - (u₁ * a) % q) + q * k by linarith [hk],
          Int.add_mul_emod_self_left]
        exact Int.emod_eq_of_lt hpos (by omega)
      omega
    · have : ((u₂ - u₁) * a) % q = q + ((u₂ * a) % q - (u₁ * a) % q) := by
        rw [show (u₂ - u₁) * a
            = (q + ((u₂ * a) % q - (u₁ * a) % q)) + q * (k - 1) by linarith [hk],
          Int.add_mul_emod_self_left]
        exact Int.emod_eq_of_lt (by omega) (by omega)
      omega
  case _ =>
    -- u₂ < u₁ : symmetric with d = u₁ - u₂
    obtain ⟨hm1, hm2⟩ := hmid u₁ hu₁.1 hu₁.2
    obtain ⟨hm1', hm2'⟩ := hmid u₂ hu₂.1 hu₂.2
    have hd1 : 1 ≤ u₁ - u₂ := by omega
    have hd2 : u₁ - u₂ ≤ K := by omega
    obtain ⟨hdm1, hdm2⟩ := hmid _ hd1 hd2
    have hdvd : q ∣ (u₁ - u₂) * a - ((u₁ * a) % q - (u₂ * a) % q) := by
      have h1 := Int.emod_def (u₁ * a) q
      have h2 := Int.emod_def (u₂ * a) q
      exact ⟨u₁ * a / q - u₂ * a / q, by linear_combination h2 - h1⟩
    obtain ⟨k, hk⟩ := hdvd
    rcases le_or_gt 0 ((u₁ * a) % q - (u₂ * a) % q) with hpos | hneg
    · have : ((u₁ - u₂) * a) % q = (u₁ * a) % q - (u₂ * a) % q := by
        rw [show (u₁ - u₂) * a
            = ((u₁ * a) % q - (u₂ * a) % q) + q * k by linarith [hk],
          Int.add_mul_emod_self_left]
        exact Int.emod_eq_of_lt hpos (by omega)
      omega
    · have : ((u₁ - u₂) * a) % q = q + ((u₁ * a) % q - (u₂ * a) % q) := by
        rw [show (u₁ - u₂) * a
            = (q + ((u₁ * a) % q - (u₂ * a) % q)) + q * (k - 1) by linarith [hk],
          Int.add_mul_emod_self_left]
        exact Int.emod_eq_of_lt (by omega) (by omega)
      omega

/-- **L2, the small-moduli seal** (`N < q ≤ 2N`, `N` odd): some family element
of `F_D(N)` sits within circle-distance `1/q` of the integers at `t = a/q` —
so no candidate at a small modulus reaches past the floor.  The `u₀ = N-1`
escapes die by `x = D(N-1)` (residue `0`), by parity (`q` even is impossible
since `N-1` is even), and by the reflected element `u' = q-(N-1)` (odd, hence
`≠ N-1`, and in range). -/
theorem small_moduli_seal (D N q a : ℤ) (hNodd : Odd N) (hN : 8 ≤ N)
    (hq1 : N < q) (hq2 : q ≤ 2*N) :
    ∃ v : ℤ, ((1 ≤ v ∧ v ≤ N ∧ v ≠ N-1) ∨ v = D*(N-1)) ∧
      ((v * a) % q ≤ 1 ∨ q - 1 ≤ (v * a) % q) := by
  obtain ⟨u₀, hu1, hu2, hu3⟩ := dirichlet_step q a (by omega)
  have huN : u₀ ≤ N := by omega
  rcases eq_or_ne u₀ (N-1) with hu | hu
  · -- the Dirichlet element is the deleted one
    rcases le_or_gt ((u₀ * a) % q) 1 with hlow | hhigh
    · rcases eq_or_ne ((u₀ * a) % q) 0 with h0 | h0
      · -- residue 0: the far element x = D(N-1) also has residue 0
        refine ⟨D*(N-1), Or.inr rfl, Or.inl ?_⟩
        have hdvd : q ∣ (N-1) * a := by
          have hdef := Int.emod_def (u₀ * a) q
          rw [hu] at h0 hdef
          exact ⟨(N-1) * a / q, by linear_combination h0 - hdef⟩
        obtain ⟨k, hk⟩ := hdvd
        have : (D*(N-1)) * a % q = 0 :=
          Int.emod_eq_zero_of_dvd ⟨D * k, by linear_combination D * hk⟩
        omega
      · -- residue 1 (q must be odd; reflect)
        have hres1 : (u₀ * a) % q = 1 := by
          have := Int.emod_nonneg (u₀ * a) (show q ≠ 0 by omega)
          omega
        -- parity: N-1 even, residue odd forces q odd
        obtain ⟨t, ht⟩ := hNodd
        have hqodd : ¬ (2 ∣ q) := by
          rintro ⟨s, hs⟩
          have h2q : (2:ℤ) ∣ q := ⟨s, hs⟩
          have h2res : (2:ℤ) ∣ ((u₀ * a) % q) := by
            rw [hu, Int.emod_def]
            exact dvd_sub ⟨t*a, by linear_combination a*ht⟩ (h2q.mul_right _)
          omega
        -- reflected element u' = q - (N-1)
        refine ⟨q - (N-1), Or.inl ⟨by omega, by omega, by omega⟩, ?_⟩
        rw [hu] at hres1
        have hdef := Int.emod_def ((N-1) * a) q
        rw [hres1] at hdef
        right
        have : (q - (N-1)) * a = (q - 1) + q * (a - ((N-1) * a / q) - 1) := by
          linear_combination hdef
        rw [this, Int.add_mul_emod_self_left]
        have hh : (q - 1) % q = q - 1 := Int.emod_eq_of_lt (by omega) (by omega)
        omega
    · -- residue q-1 (mirror of the residue-1 case)
      have hresq : (u₀ * a) % q = q - 1 := by
        have := Int.emod_lt_of_pos (u₀ * a) (show 0 < q by omega)
        omega
      obtain ⟨t, ht⟩ := hNodd
      have hqodd : ¬ (2 ∣ q) := by
        rintro ⟨s, hs⟩
        have h2q : (2:ℤ) ∣ q := ⟨s, hs⟩
        have h2res : (2:ℤ) ∣ ((u₀ * a) % q) := by
          rw [hu, Int.emod_def]
          exact dvd_sub ⟨t*a, by linear_combination a*ht⟩ (h2q.mul_right _)
        omega
      refine ⟨q - (N-1), Or.inl ⟨by omega, by omega, by omega⟩, ?_⟩
      rw [hu] at hresq
      have hdef := Int.emod_def ((N-1) * a) q
      rw [hresq] at hdef
      left
      have : (q - (N-1)) * a = 1 + q * (a - ((N-1) * a / q) - 1) := by
        linear_combination hdef
      rw [this, Int.add_mul_emod_self_left]
      have : (1 : ℤ) % q = 1 := Int.emod_eq_of_lt (by omega) (by omega)
      omega
  · -- the Dirichlet element is a genuine base element
    exact ⟨u₀, Or.inl ⟨hu1, huN, hu⟩, hu3⟩

/-! ## Non-vacuity: the concrete (N, D) = (31, 4) member — THM-1285's 4/127 -/

example : ((2*4-1) * 55) % ((31+1)*4 - 1) = 4 := by decide

/-- The concrete rung floor at the first D=4 tower member: at `t = 55/127`,
every element of `{1..31}\{30} ∪ {120}` is at circle-distance `≥ 4/127` from
every integer — the kernel-checked floor half of `M({1..29,31,120}) = 4/127`. -/
theorem member_4_127 (v : ℤ) (hv : (1 ≤ v ∧ v ≤ 31 ∧ v ≠ 30) ∨ v = 120)
    (m : ℤ) : (4 : ℝ) / 127 ≤ |(v : ℝ) * ((55 : ℝ) / 127) - m| := by
  have hv' : (1 ≤ v ∧ v ≤ 31 ∧ v ≠ 31 - 1) ∨ v = 4 * (31 - 1) := by
    rcases hv with ⟨h1, h2, h3⟩ | h
    · exact Or.inl ⟨h1, h2, by omega⟩
    · exact Or.inr (by omega)
  have h := rung_floor_witness 4 31 55 (by norm_num) (by norm_num)
    (by decide) v hv' m
  norm_num at h
  exact h

end RungFloor
end LonelyRunner

#print axioms LonelyRunner.RungFloor.cert_dvd
#print axioms LonelyRunner.RungFloor.rung_core
#print axioms LonelyRunner.RungFloor.rung_band
#print axioms LonelyRunner.RungFloor.rung_floor_single
#print axioms LonelyRunner.RungFloor.rung_floor_witness
#print axioms LonelyRunner.RungFloor.exists_binder_multiplier
#print axioms LonelyRunner.RungFloor.Q_dvd_iff_binder_dvd
#print axioms LonelyRunner.RungFloor.rung_dead_of_common_factor
#print axioms LonelyRunner.RungFloor.small_moduli_seal_low
#print axioms LonelyRunner.RungFloor.dirichlet_step
#print axioms LonelyRunner.RungFloor.small_moduli_seal
#print axioms LonelyRunner.RungFloor.packing_seal
#print axioms LonelyRunner.RungFloor.member_4_127
