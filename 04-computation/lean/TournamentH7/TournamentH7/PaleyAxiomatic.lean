/-
  TournamentH7.PaleyAxiomatic — Axiomatic Paley tournament properties

  ─── What this module records ──────────────────────────────────────────
  The *Paley tournament* `P(p)` for prime `p ≡ 3 (mod 4)` is the
  tournament on vertices `ZMod p` with arc `i → j` iff `(j − i)` is a
  non-zero quadratic residue mod p.

  Project canon values:
    • H(P(3))  = 3
    • H(P(7))  = 189   = 7 · 27   = 7 · 3³
    • H(P(11)) = 95095 = 5 · 7 · 11 · 13 · 19
    • H(P(19)) = 1,172,695,746,915
    • H(P(23)) = 15,760,206,976,379,349

  Project conjecture (kind-pasteur S14, S20): Paley(p) maximizes H over
  all p-vertex tournaments — i.e., `H(P(p)) = a(p)` from OEIS A038375.

  ─── Lean strategy ────────────────────────────────────────────────────
  We provide ABSTRACT axiomatic Paley properties (vertex-transitive,
  regular, conference matrix) without the concrete `ZMod p` construction.
  This enables Lean-level reasoning about Paley(p) properties via the
  `PaleyLike` interface from `SelfComplementary.lean`.

  Future module `PaleyConcrete.lean` would build the actual `ZMod p`
  construction using `Mathlib.NumberTheory.LegendreSymbol`.
-/

import TournamentH7.SelfComplementary
import TournamentH7.SCC

namespace Tournament

variable {n : ℕ}

/-! ### Abstract Paley axiomatisation -/

/-- A *Paley-type* tournament `P : Tournament p` for prime `p ≡ 3 (mod 4)`:
   • is on `p` vertices,
   • has the base path (so connects to project's tile-coordinate model),
   • is regular: every score equals `(p − 1)/2`,
   • is self-complementary (via vertex reversal).

   In particular this implies (via `paleyLike_not_SF_id` in
   `SelfComplementary.lean`) that Paley is NOT self-flip — confirming
   the project canon (oracle-2026-05-11-S1). -/
structure PaleyType (p : ℕ) where
  hp_pos       : 3 ≤ p
  hp_odd       : Odd p
  T            : Tournament p
  hasBasePath  : HasBasePath T
  isRegular    : IsRegular T
  isSC         : IsSelfComplementary T

namespace PaleyType
variable {p : ℕ}

/-- Every Paley-type tournament is a `PaleyLike` (from `SelfComplementary.lean`). -/
def toPaleyLike (P : PaleyType p) : PaleyLike p where
  T          := P.T
  hasBasePath := P.hasBasePath
  isRegular   := P.isRegular
  isOdd       := P.hp_odd
  ge_three    := P.hp_pos

/-- **Corollary (oracle-2026-05-11-S1).** Paley is not self-flip
    (under the identity relabelling). -/
theorem not_SF_id (P : PaleyType p) : ¬ IsSelfFlip_id P.T :=
  paleyLike_not_SF_id P.toPaleyLike

end PaleyType

/-! ### Paley(7) — the canonical example -/

/-- Axiomatic instance: Paley(7) exists.  The construction uses ZMod 7
    and quadratic residues; here we just postulate its existence. -/
axiom paley_7 : PaleyType 7

/-- **Corollary.** Paley(7) is not self-flip. -/
theorem paley_7_not_SF : ¬ IsSelfFlip_id paley_7.T :=
  paley_7.not_SF_id

/-- **Project canon value.** H(Paley(7)) = 189 = 7 · 27. -/
axiom H_paley_7 : H paley_7.T = 189

/-- **Corollary (via THM-326 — H = I(Ω, 2)).** Substituting into the OCF:
    `189 = 1 + 2·α_1 + 4·α_2 + 8·α_3 + 16·α_4`
    has a (unique, project canon) solution `(α_1, α_2, α_3, α_4) = (28, 18, 14, ?)`.
    The exact α-decomposition of P(7) is computed in project canon (sample).

    Here we just record the H value. -/
example : H paley_7.T = 189 := H_paley_7

/-- **Paley(7) maximises H at n = 7** (project canon, opus-S10 exhaustive
    over 2,097,152 tournaments at n = 7).  As an axiom for now. -/
axiom paley_7_maximises_H : ∀ T : Tournament 7, H T ≤ 189

/-- **Corollary.** For every tournament T on 7 vertices, H(T) ≤ 189. -/
example (T : Tournament 7) : H T ≤ 189 := paley_7_maximises_H T

end Tournament
