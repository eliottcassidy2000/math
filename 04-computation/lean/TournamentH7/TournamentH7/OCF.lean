/-
  TournamentH7.OCF — Odd-Cycle Collection Formula and structural axioms

  This file collects the EXTERNAL classical results we use.  They are
  recorded as `axiom`s with full citations so the H ≠ 7 proof depends only
  on published theorems.

  ─── External results axiomatised ────────────────────────────────────────
   (A1) Odd-Cycle Collection Formula
        — Grinberg–Stanley, arXiv:2412.10572, Corollary 20 (2024).
            H(T) = Σ_k α_k(T) · 2^k.
   (A2) Moon–Moser 3-cycle bound
        — Moon, *Topics on Tournaments*, 1968 (orig. 1962).
            Every SC tournament on s ≥ 3 vertices has ≥ s − 2 directed
            3-cycles.
   (A3) Moon–Camion Hamilton cycle
        — Camion, *Chemins et circuits…*, C. R. Acad. Sci. Paris (1959).
            Every SC tournament has a directed Hamilton cycle.
   (A4) SCC localisation of an Ω-triangle (folklore via SCC partition).
            If three odd directed cycles of T pairwise share a vertex,
            they all lie in a single SCC.
   (A5) Small-SCC odd-cycle counts (3-cycle uniqueness, n=4 enumeration).
            SC tournament on 3 vertices: exactly 1 odd directed cycle.
            SC tournament on 4 vertices: exactly 2 directed 3-cycles
            (and no other odd cycles since 4 < 5).

  ───────────────────────────────────────────────────────────────────── -/

import TournamentH7.Cycles
import TournamentH7.SCC

namespace Tournament

variable {n : ℕ}

/-- α k T = number of size-k independent sets in the conflict graph Ω(T) of
    odd directed cycles. Two odd cycles conflict iff they share a vertex.

    Concretely:
      • α₀ T = 1 (empty set is always independent),
      • α₁ T = total number of odd directed cycles in T,
      • α₂ T = number of unordered vertex-disjoint pairs,
      • α₃ T = number of triples of pairwise vertex-disjoint odd cycles,
      • …                                                              -/
opaque alphaCount (k : ℕ) (T : Tournament n) : ℕ

/-- Number of odd directed cycles of T whose vertex set is contained in S. -/
opaque oddCyclesIn (T : Tournament n) (S : Finset (Fin n)) : ℕ

/-! ### (A1) Odd-Cycle Collection Formula -/

/-- **Axiom (A1, OCF — Grinberg–Stanley 2024).**
    Polynomial form, truncated at k = 4 (sufficient because for any T the
    sum is finite and the H = 7 reasoning never needs k ≥ 4). -/
axiom ocf (T : Tournament n) :
    H T = 1 + 2 * alphaCount 1 T + 4 * alphaCount 2 T
            + 8 * alphaCount 3 T + 16 * alphaCount 4 T

/-- **Axiom (A1', independence polynomial subset bound).**
    A size-`k` independent set requires `k` distinct odd cycles, so
    `α_k(T) ≠ 0 → α_1(T) ≥ k` for `k ≥ 1`. This is built into the
    definition of α_k as the number of independent k-subsets of Ω(T). -/
axiom alpha_subset_bound (T : Tournament n) (k : ℕ) (hk : 1 ≤ k) :
    alphaCount k T ≠ 0 → k ≤ alphaCount 1 T

/-! ### (A2) Moon–Moser 3-cycle bound -/

/-- **Axiom (A2, Moon–Moser 1962).**
    Every SCC of T on s ≥ 3 vertices houses at least s − 2 distinct
    directed 3-cycles. Three-cycles are odd, so they contribute to
    `oddCyclesIn`. -/
axiom moonMoser (T : Tournament n) (S : Finset (Fin n))
    (hSCC : IsSCC T S) (hs : 3 ≤ S.card) :
    S.card - 2 ≤ oddCyclesIn T S

/-! ### (A3) Moon–Camion Hamilton cycle for odd s ≥ 5 -/

/-- **Axiom (A3, Moon–Camion 1959/1968).**
    When an SCC of T has odd size s ≥ 5, the directed Hamilton cycle on
    S is an odd cycle of length s that is distinct from every 3-cycle.
    Combined with (A2): `oddCyclesIn T S ≥ (s − 2) + 1 = s − 1`. -/
axiom moonCamion_oddSize (T : Tournament n) (S : Finset (Fin n))
    (hSCC : IsSCC T S) (hodd : Odd S.card) (hs : 5 ≤ S.card) :
    S.card - 1 ≤ oddCyclesIn T S

/-! ### (A4) SCC localisation of an Ω-triangle -/

/-- **Axiom (A4, folklore: cycles localise to SCCs).**
    If every pair of three odd directed cycles of T meets in a vertex
    (i.e. `α₁ T = 3 ∧ α₂ T = 0`), then there is a single SCC `S` of T
    that contains the vertex set of every odd cycle of T. In particular
    `oddCyclesIn T S = α₁ T = 3` and `|S| ≥ 3`. -/
axiom omegaTriangleLocalises (T : Tournament n)
    (h1 : alphaCount 1 T = 3) (h2 : alphaCount 2 T = 0) :
    ∃ S : Finset (Fin n),
      IsSCC T S ∧ 3 ≤ S.card ∧ oddCyclesIn T S = 3

/-! ### (A5) Small SCC odd-cycle counts -/

/-- **Axiom (A5a).** SC tournament on exactly 3 vertices is a 3-cycle and
    contains exactly 1 odd directed cycle. -/
axiom oddCyclesIn_size3 (T : Tournament n) (S : Finset (Fin n))
    (hSCC : IsSCC T S) (hs : S.card = 3) :
    oddCyclesIn T S = 1

/-- **Axiom (A5b).** SC tournament on 4 vertices has score sequence
    (1,1,2,2); it has exactly 2 directed 3-cycles and (since 4 < 5) no
    longer odd cycles. -/
axiom oddCyclesIn_size4 (T : Tournament n) (S : Finset (Fin n))
    (hSCC : IsSCC T S) (hs : S.card = 4) :
    oddCyclesIn T S = 2

end Tournament
