import TournamentH7.TNCDetectionDepth

/-
# The moment-nullcone interface (THM-1775)

A single abstract structure for the template that THM-1775 shows tournament transitivity, TNC,
and GMC(2) all instantiate: a moment sequence `phi : ℕ → α` whose "all moments vanish" nullcone
is detected at a finite depth = the order of its governing recurrence.  The detection engine is
`TNCDepth.zeros_propagate`, already proved sorry-free; this file packages it and gives the
constructor every instance actually uses — building the propagation `step` from a monic linear
recurrence (Cayley–Hamilton for tournaments, the holonomic recurrence for TNC/GMC).

`sorry`-free; the tournament/TNC/GMC specializations feed `ofMonicRec` with their own
recurrences (Cayley–Hamilton, THM-1670, THM-1740 respectively).
-/

namespace MomentNullcone

variable {α : Type*} [Zero α]

/-- The data of a moment sequence together with the forward zero-propagation its recurrence
supplies: order `r ≥ 1`, and the step "`r` consecutive zeros force the next term to vanish". -/
structure Data (α : Type*) [Zero α] where
  phi : ℕ → α
  order : ℕ
  order_pos : 1 ≤ order
  step : ∀ m, 1 ≤ m → (∀ i, i < order → phi (m + i) = 0) → phi (m + order) = 0

/-- **Finite detection depth.**  If the first `order` moments vanish, every moment vanishes:
the nullcone `{phi(m)=0 ∀m≥1}` is cut out by `phi(1),…,phi(order)`.  This is the template's
common conclusion, discharged by `zeros_propagate`. -/
theorem detect (M : Data α)
    (init : ∀ j, 1 ≤ j → j ≤ M.order → M.phi j = 0) :
    ∀ m, 1 ≤ m → M.phi m = 0 :=
  TNCDepth.zeros_propagate M.order_pos M.step init

/-- Contrapositive: escaping the nullcone happens within the first `order` moments. -/
theorem escape_within (M : Data α)
    (hne : ∃ m, 1 ≤ m ∧ M.phi m ≠ 0) :
    ∃ j, 1 ≤ j ∧ j ≤ M.order ∧ M.phi j ≠ 0 :=
  TNCDepth.nonzero_within_depth M.order_pos M.step hne

end MomentNullcone

/-- Build a `MomentNullcone.Data` from a **monic order-`r` linear recurrence**
`phi (m + r) = combo m` where `combo m` is any term that is `0` whenever the `r` preceding
moments all vanish.  This is the shape produced by Cayley–Hamilton (tournaments), the toral
recurrence (TNC, THM-1670), and the moment recurrence (GMC, THM-1740): with leading
coefficient a unit, the recurrence rearranges to `phi (m+r) = − (lower-order combination)`,
which is `0` when the lower terms are `0`. -/
def MomentNullcone.ofMonicRec {α : Type*} [Zero α]
    (phi : ℕ → α) (r : ℕ) (hr : 1 ≤ r)
    (combo : ℕ → α)
    (rec : ∀ m, phi (m + r) = combo m)
    (combo_zero : ∀ m, 1 ≤ m → (∀ i, i < r → phi (m + i) = 0) → combo m = 0) :
    MomentNullcone.Data α where
  phi := phi
  order := r
  order_pos := hr
  step := fun m hm hz => by rw [rec m]; exact combo_zero m hm hz

/-!
## The three instances (recorded; each supplies its recurrence to `ofMonicRec`)

* **Tournament (rational floor).** `phi m = tr(A^m)`; Cayley–Hamilton gives
  `tr(A^{m}) = −∑_{i=1}^{n} c_i · tr(A^{m−i})` (monic, `c_i` = char-poly coefficients), so
  `order = n` and the nullcone `{tr(A^m)=0 ∀m}` = nilpotent = transitive is detected at depth
  `n`.  (THM-1775, THM-895.)
* **TNC (algebraic).** `phi m = CT(Λ^m)`; the toral recurrence (THM-1670) is order `D = M+N`,
  and THM-1720/1725 show its leading coefficient never vanishes at `m ≥ 1`, so `ofMonicRec`
  applies; nullcone = one-sided, depth `D` (THM-1710).
* **GMC(2) (holonomic).** `phi m = E[P^m]`; the moment recurrence (THM-1740) is finite order
  `K ~ charge span`; nullcone = charge-one-sided, depth `K`.

`H` (Hamiltonian-path count) does NOT instantiate this: it is not a moment sequence at all —
it splits within a co-spectral class at `n = 6` (THM-1780), so it has no governing recurrence
in the trace moments.  The interface correctly excludes it.
-/
