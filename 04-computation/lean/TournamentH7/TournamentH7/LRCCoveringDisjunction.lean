/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S129)
-/
import TournamentH7.LRCCoveringReach

/-!
# The covering-disjunction endpoint of branch 2 (HYP-4652)

Branch 2 of the `(C)`-crux (`00-navigation/LRC14-PROOF-MAP.md`) is the finite covering system:
*every non-AP blocker clears at some modulus `q ≤ Q₀`.*  mac-mini's atom
`LRCCoveringReach.reach_ge_of_covering` supplies the **per-`q`** step — an avoid-band rotation `c`
mod `q` forces `μ/q ≤ M(v)`.  What the residue check actually produces is a **disjunction over a
finite covering set** `q ∈ Q`, together with the arithmetic that the chosen band radius clears
`2/25` (`μ/q ≥ 2/25`).  This file packages exactly that:

> **`loose_of_covering_witness`** — a single clearing witness `(q, c, μ)` with `μ/q ≥ 2/25`
> (encoded `2q ≤ 25μ`) and the avoid-band condition ⟹ `M(v) ≥ 2/25` (loose).
> **`loose_of_covering_set`** — the same with the witness packaged existentially, the shape the
> residue enumeration emits: `(∃ q c μ, …) ⟹ M(v) ≥ 2/25`.

So branch 2 reduces cleanly to the **purely arithmetic** obligation *"every compressed non-AP
12-blocker admits such a witness with `q ≤ 39`"* (klein-S144: `q ≤ 31` on compressed families, 0
gaps to height 650 000) — a finite mod-`q` residue check, now isolated behind this Lean endpoint
with no analysis left in the reduction.  `M(v) = sSup (margin v '' [0,1])` throughout, matching the
covering-reach atom and `LRCMreachConcrete`.
-/

namespace LonelyRunner
namespace CoveringDisjunction

open TournamentH7.LRCWitness

/-- **Per-witness clearance.**  A rotation `c` mod `q` (`0 ≤ c ≤ q`, `0 < q`) whose avoid-band radius
`μ` clears `2/25` (`2q ≤ 25μ`, i.e. `μ/q ≥ 2/25`) and which keeps every speed `v i` in the safe band
`[μ, q−μ]` (mod `q`) forces `M(v) ≥ 2/25`.  Threads mac-mini's `reach_ge_of_covering` through the
`μ/q ≥ 2/25` arithmetic. -/
theorem loose_of_covering_witness (v : Fin 12 → ℤ) (q c μ : ℤ)
    (hq : 0 < q) (hc0 : 0 ≤ c) (hcq : c ≤ q) (hμ : 2 * q ≤ 25 * μ)
    (hclear : ∀ i, μ ≤ (v i * c) % q ∧ (v i * c) % q ≤ q - μ) :
    (2 : ℝ) / 25 ≤ sSup (margin v '' Set.Icc (0 : ℝ) 1) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hμR : (2 : ℝ) * (q : ℝ) ≤ 25 * (μ : ℝ) := by exact_mod_cast hμ
  have h2 : (2 : ℝ) / 25 ≤ (μ : ℝ) / (q : ℝ) := by
    rw [le_div_iff₀ hqR]; linarith
  exact le_trans h2
    (TournamentH7.LRCCoveringReach.reach_ge_of_covering v q c μ hq hc0 hcq hclear)

/-- **The covering disjunction (branch-2 endpoint).**  If *some* modulus `q` in the covering system
admits a clearing witness `(c, μ)` — band radius `μ` clearing `2/25` and every speed off the hole
mod `q` — then `M(v) ≥ 2/25`, so `v` is loose (out of the gap `(1/13, 2/25)`).  This is the exact
shape the finite residue check emits per family; discharging its premise for every compressed non-AP
12-blocker (`q ≤ 39`) is all that remains of branch 2. -/
theorem loose_of_covering_set (v : Fin 12 → ℤ)
    (hcov : ∃ q c μ : ℤ, 0 < q ∧ 0 ≤ c ∧ c ≤ q ∧ 2 * q ≤ 25 * μ ∧
        (∀ i, μ ≤ (v i * c) % q ∧ (v i * c) % q ≤ q - μ)) :
    (2 : ℝ) / 25 ≤ sSup (margin v '' Set.Icc (0 : ℝ) 1) := by
  obtain ⟨q, c, μ, hq, hc0, hcq, hμ, hclear⟩ := hcov
  exact loose_of_covering_witness v q c μ hq hc0 hcq hμ hclear

/-- **Mod-25 is the `q = 25, μ = 2` instance** (sanity check that the disjunction subsumes case 1).
An avoid-band rotation `c` mod 25 keeping every speed in `[2, 23]` clears, since `2·25 ≤ 25·2`. -/
theorem loose_of_mod25_instance (v : Fin 12 → ℤ) (c : ℤ) (hc0 : 0 ≤ c) (hc25 : c ≤ 25)
    (hclear : ∀ i, 2 ≤ (v i * c) % 25 ∧ (v i * c) % 25 ≤ 23) :
    (2 : ℝ) / 25 ≤ sSup (margin v '' Set.Icc (0 : ℝ) 1) :=
  loose_of_covering_witness v 25 c 2 (by norm_num) hc0 hc25 (by norm_num)
    (fun i => ⟨(hclear i).1, by have := (hclear i).2; omega⟩)

/-! ## `(C)` ⟸ covering completeness

The crux `(C)` — *the dilated AP is the unique 12-family with `M < 2/25`* — reduces to a single
arithmetic obligation once the covering disjunction is the endpoint: **every non-AP family admits a
covering witness.**  The three-case split (non-blocker mod 25, blocker `q ≤ 39`, mult-of-25 small
`q`) is invisible here — all three are instances of the one existential `HasCoveringWitness`. -/

/-- `v` is a **dilated arithmetic progression**: an AP `a + d·(0,1,…,11)` read in some order `σ`.
The unique tight locus of the 12-speed gap (the AP `{1,…,12}` is `a = 1, d = 1, σ = id`). -/
def DilatedAP (v : Fin 12 → ℤ) : Prop :=
  ∃ a d : ℤ, ∃ σ : Equiv.Perm (Fin 12), ∀ i, v i = a + d * ((σ i).val : ℤ)

/-- `v` admits a **covering witness**: a modulus `q` with a clearing rotation `(c, μ)` (`μ/q ≥ 2/25`,
every speed off the hole mod `q`).  This is the premise of `loose_of_covering_set`, so it forces
`M(v) ≥ 2/25`.  Cases 1/2/3 of the crux are all instances of this one existential. -/
def HasCoveringWitness (v : Fin 12 → ℤ) : Prop :=
  ∃ q c μ : ℤ, 0 < q ∧ 0 ≤ c ∧ c ≤ q ∧ 2 * q ≤ 25 * μ ∧
    (∀ i, μ ≤ (v i * c) % q ∧ (v i * c) % q ≤ q - μ)

/-- **The crux `(C)` as a statement.**  Every 12-family with `M(v) < 2/25` is a dilated AP. -/
def CruxStatement : Prop :=
  ∀ v : Fin 12 → ℤ, sSup (margin v '' Set.Icc (0 : ℝ) 1) < 2 / 25 → DilatedAP v

/-- **The covering-completeness obligation — ⚠ NOT a FINITE check; `(G)`-equivalent (MISTAKE-116,
opus-S130).**  Every non-dilated-AP 12-family admits a covering witness (`∃ q`, no bound).  This was
believed to be a *finite* mod-`q` residue check (`q ≤ 39`), but mac-mini-S36 (verified opus-S130)
showed the clearing modulus is UNBOUNDED: compressed varying-k families `≡ AP mod lcm(2..Q₀)` escape
every `q ≤ Q₀` and clear only at `nextprime(Q₀)`.  So `q` is genuinely unbounded, and this statement
is EQUIVALENT to the analytic `(G)` ("every non-AP 12-family is loose"), NOT a finite reduction of it
(the klein `q ≤ 31` figure was a height-sampling artifact `≤ 650k`; the escape families sit at height
`~10¹⁴`).  Still an honest `Prop` obligation — but discharging it needs a **scale-uniform**
decorrelation/Fourier argument, not a finite pile of certs. -/
def CoveringComplete : Prop :=
  ∀ v : Fin 12 → ℤ, ¬ DilatedAP v → HasCoveringWitness v

/-- **`(C)` ⟸ covering completeness.**  If every non-AP family admits a covering witness (at SOME
modulus `q`, unbounded — see the `CoveringComplete` warning), then the crux holds:
`M(v) < 2/25 ⟹ v` is a dilated AP.  The AP is exempt precisely because it is the unique `M`-minimizer
(`1/13`), with no clearing witness at any `q`.  This is a VALID implication, but `CoveringComplete` is
`(G)`-equivalent (analytic), so this does NOT reduce `(C)` to something finite. -/
theorem crux_of_covering_complete (h : CoveringComplete) : CruxStatement := by
  intro v hM
  by_contra hnAP
  exact absurd (loose_of_covering_set v (h v hnAP)) (by linarith)

#print axioms loose_of_covering_witness
#print axioms loose_of_covering_set
#print axioms loose_of_mod25_instance
#print axioms crux_of_covering_complete

end CoveringDisjunction
end LonelyRunner
