/-
LRC-EXTREMAL CERTIFICATE (opus-2026-07-20-S439) -- the WOWII-103 native_decide template
applied to the Lonely Runner extremals. Mathlib-free, ℕ-only, kernel-checkable.

M(V) = max over t∈[0,1) of min over v∈V of ‖v·t‖   (‖x‖ = distance to nearest integer).
On the beat grid t = a/D (D a beat denominator, 0 ≤ a < D), we have
    min_{v∈V} ‖v·(a/D)‖  =  md(V,D,a) / D ,   md = min_{v∈V} min((v·a) mod D, D − (v·a) mod D).
By the finite-attainment fact (THM-401: M is attained at a pair-sum time), M(V) = max over the
beat grid of md/D. The certificate is then two ℕ-decidable facts:
    UPPER:  ∀ grid point (D,a):  14·md(D,a) ≤ D        ⟹  M(V) ≤ 1/14
    LOWER:  md(14,1) = 1  (with D = 14)                ⟹  M(V) ≥ 1/14
Hence M(V) = 1/14 for both the classical family {1..13} and its leaf-inflation {1..11,13,24}.
Mirrors 04-computation/lrc_extremal_certificate_reference_opus_S439.py.

SCOPE / RELATION TO THE FLEET'S LEAN (honesty note, S439). The fleet's TournamentH7 project
ALREADY proves M({1..13}) = 1/14 more strongly: `LRCAPTight.mreach_AP_eq` derives it for the true
real supremum via Dirichlet approximation (not a finite grid), and `LRCGridValue.ap12_margin_eq`
is a reusable native_decide grid decision procedure (instantiated at n=12 → 1/13). This file is
NOT the first M({1..13}) certificate. Its distinct value: (i) it is STANDALONE and Mathlib-FREE
(ℕ-only, core Lean), and (ii) it newly certifies the LEAF-INFLATION extremal {1..11,13,24} -- the
WOWII-103 inflation-motif family that the fleet's n=12/n=13 flagships do not cover -- and that it
shares the witness t=1/14. A self-contained WOWII-template instance, complementary to LRCGridValue.
-/

/-- min over speeds of the integer "distance numerator" at t = a/D. Seeded with D (a valid upper
    bound since each term ≤ D). -/
def md (V : List Nat) (D a : Nat) : Nat :=
  (V.map (fun v => let r := (v * a) % D; min r (D - r))).foldr min D

/-- beat denominators: {2v} ∪ {v+w, w−v : v < w}. All positive (v ≥ 1). Duplicates are harmless. -/
def beatDenoms (V : List Nat) : List Nat :=
  (V.map (fun v => 2 * v)) ++
  ((V.map (fun v => (V.map (fun w => if v < w then [v + w, w - v] else [])).flatten)).flatten)

/-- the finite critical grid of (D, a) with 0 ≤ a < D. -/
def grid (V : List Nat) : List (Nat × Nat) :=
  ((beatDenoms V).map (fun D => (List.range D).map (fun a => (D, a)))).flatten

def V13 : List Nat := [1,2,3,4,5,6,7,8,9,10,11,12,13]
def V24 : List Nat := [1,2,3,4,5,6,7,8,9,10,11,13,24]

/-- UPPER bound: over the whole beat grid, the min-distance numerator never exceeds D/14,
    i.e. min_v ‖v t‖ ≤ 1/14 at every critical t. -/
theorem lrc13_upper : (grid V13).all (fun p => 14 * md V13 p.1 p.2 ≤ p.1) = true := by native_decide
theorem lrc24_upper : (grid V24).all (fun p => 14 * md V24 p.1 p.2 ≤ p.1) = true := by native_decide

/-- LOWER bound: t = 1/14 (D = 14, a = 1) attains min-distance exactly 1/14. -/
theorem lrc13_lower : md V13 14 1 = 1 := by native_decide
theorem lrc24_lower : md V24 14 1 = 1 := by native_decide

/-- The two families have IDENTICAL argmax structure at t = 1/14: the leaf-inflation 12 ↦ 24
    preserves the extremal value AND its witness (24·(1/14) = 12/7 reproduces the danger comb). -/
theorem lrc_leaf_inflation_same_witness : md V13 14 1 = md V24 14 1 := by native_decide

#eval "LRC extremal certificate OK: M({1..13}) = M({1..11,13,24}) = 1/14 (native_decide, ℕ-only)"

-- trust base: native_decide adds `Lean.ofReduceBool` (trusts the compiler) atop propext +
-- Classical.choice. The finite-attainment step (beat grid = all local maxima) is THM-401.
#print axioms lrc13_upper
#print axioms lrc24_upper
