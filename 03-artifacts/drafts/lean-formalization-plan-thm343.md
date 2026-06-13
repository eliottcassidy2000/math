# Lean Formalization Plan: THM-343 (H(T) ≠ 7)

**Session:** opus-2026-05-28-S6
**Target proof:** THM-343 — "For every tournament T, the number of Hamiltonian paths H(T) ≠ 7."
**Cleaner proof version:** see `01-canon/theorems/THM-343-H7-impossible.md` (S6 revision).

This document maps each step of the proof to Lean / Mathlib infrastructure, identifies what needs to be defined from scratch, and provides a proof skeleton.

---

## 1. Mathlib Status (as of May 2026)

What exists:
- `Mathlib.Combinatorics.Digraph.Basic` — `Digraph V` type (adjacency relation).
- `Mathlib.Combinatorics.Quiver.*` — quiver paths and basic morphisms.
- `Mathlib.Combinatorics.SimpleGraph.*` — extensive, but undirected.
- `Mathlib.Combinatorics.SimpleGraph.Connectivity` — walks, paths, cycles for SimpleGraph.

What does NOT exist:
- A `Tournament` type (complete oriented graph).
- Strongly-connected components for `Digraph`.
- Directed cycle / Hamilton cycle for `Digraph`.
- Moon's theorem.
- The Odd-Cycle Collection Formula (OCF) of Grinberg-Stanley.

Thus the formalization needs to build a tower of definitions and lemmas. Estimated effort: substantial (several thousand lines).

---

## 2. Foundational Definitions Required

```lean
import Mathlib.Combinatorics.Digraph.Basic
import Mathlib.Combinatorics.Quiver.Path
import Mathlib.Data.Fintype.Basic

/-- A tournament on a finite vertex type V: a digraph whose adjacency relation
    is irreflexive and "tournament-like": for every pair u ≠ v, exactly one of
    u ⟶ v or v ⟶ u holds. -/
structure Tournament (V : Type*) [Fintype V] where
  Adj : V → V → Prop
  irrefl : ∀ v, ¬ Adj v v
  tournament : ∀ u v : V, u ≠ v → (Adj u v ∨ Adj v u) ∧ ¬ (Adj u v ∧ Adj v u)

variable {V : Type*} [Fintype V] [DecidableEq V]
namespace Tournament
variable (T : Tournament V)

/-- A directed walk is a list of vertices with consecutive adjacencies. -/
def Walk : V → V → Type*  -- via Quiver/Path

/-- A directed cycle of length k is a closed walk of distinct vertices,
    represented as a list (v₀, v₁, …, v_{k-1}) with v_i ⟶ v_{i+1 mod k}. -/
structure DirCycle (k : ℕ) where
  vertices : Fin k → V
  distinct : Function.Injective vertices
  edges : ∀ i : Fin k, T.Adj (vertices i) (vertices ((i+1) % k))

/-- An odd directed cycle has odd length ≥ 3. -/
def OddCycle := Σ k, { c : DirCycle T k // k ≥ 3 ∧ Odd k }

/-- Hamiltonian path = directed walk visiting all vertices exactly once. -/
structure HamPath where
  perm : Equiv (Fin (Fintype.card V)) V  -- a permutation of vertices
  edges : ∀ i : Fin (Fintype.card V - 1), T.Adj (perm i) (perm ⟨i+1, …⟩)

/-- The Rédei H-invariant: number of Hamiltonian paths. -/
def H : ℕ := Fintype.card (HamPath T)

end Tournament
```

## 3. Required Lemmas (in dependency order)

### L1 — Rédei (1934)

```lean
theorem redei_odd (T : Tournament V) : Odd (T.H) := sorry
```
**Status:** classical. Not in Mathlib. Should be formalizable in ~50 lines.

### L2 — Sub-tournament

```lean
def induced (T : Tournament V) (S : Finset V) : Tournament S where
  Adj := λ u v => T.Adj u.val v.val
  …
```

### L3 — Strong connectivity

```lean
/-- Vertices u and v are mutually reachable in T. -/
def stronglyReachable (T : Tournament V) (u v : V) : Prop :=
  ∃ p : T.Walk u v, True ∧ ∃ q : T.Walk v u, True

/-- T is strongly connected if every pair is mutually reachable. -/
def IsSC (T : Tournament V) : Prop := ∀ u v : V, stronglyReachable T u v
```

### L4 — Directed cycle is SC

```lean
lemma DirCycle_SC (T : Tournament V) (k : ℕ) (c : T.DirCycle k) :
  IsSC (T.induced (Finset.image c.vertices Finset.univ)) := sorry
```
**Difficulty:** easy — follow the cycle's arcs to reach any vertex from any other.

### L5 — Union of SC sub-digraphs sharing a vertex is SC

```lean
lemma SC_union (T : Tournament V) (S₁ S₂ : Finset V)
  (h₁ : IsSC (T.induced S₁)) (h₂ : IsSC (T.induced S₂))
  (h_shared : (S₁ ∩ S₂).Nonempty) :
  IsSC (T.induced (S₁ ∪ S₂)) := sorry
```
**Difficulty:** easy — reachability transits through the shared vertex.

### L6 — Moon (1966) Corollary 2.1

```lean
/-- Number of directed k-cycles in T. -/
def cycleCount (T : Tournament V) (k : ℕ) : ℕ :=
  Fintype.card (T.DirCycle k)

/-- Moon's bound: a strong tournament on n vertices has ≥ n-k+1 directed k-cycles,
    for 3 ≤ k ≤ n. -/
theorem moon_subtournaments
  (T : Tournament V) (hSC : IsSC T) (k : ℕ) (hk : 3 ≤ k) (hkn : k ≤ Fintype.card V) :
  Fintype.card V - k + 1 ≤ T.cycleCount k := sorry
```
**Difficulty:** moderate — induction on n via Moon's argument (every SC tournament on n verts contains an SC subtournament on n-1 verts; the removed vertex is in some k-cycle by vertex-pancyclicity). Estimated ~200-400 lines.

**Alternative**: prove only the k=3 special case (needed for the main argument when s ≥ 6) plus the k=5 case at n=5 (needed for s = 5). That's a much narrower commitment.

### L7 — Camion (1959) — special case k=n of L6

```lean
theorem camion (T : Tournament V) (hSC : IsSC T) (h : 3 ≤ Fintype.card V) :
  ∃ c : T.DirCycle (Fintype.card V), True := sorry
```

### L8 — n=4 SC uniqueness

```lean
theorem n4_SC_unique_score (T : Tournament V) (h : Fintype.card V = 4) (hSC : IsSC T) :
  T.scoreSequence = [1, 1, 2, 2] := sorry
```

### L9 — c₃ counting formula

```lean
theorem three_cycle_count (T : Tournament V) :
  T.cycleCount 3 = (Nat.choose (Fintype.card V) 3) -
    (Finset.univ.sum (λ v => Nat.choose (T.outDeg v) 2)) := sorry
```

### L10 — OCF (Grinberg-Stanley)

```lean
/-- Conflict graph of odd cycles. -/
def Omega (T : Tournament V) : SimpleGraph (T.OddCycle) := …

/-- Independence polynomial evaluated at 2. -/
def alphaSum (T : Tournament V) : ℕ :=
  ∑ k, (independenceCount T.Omega k) * 2^k

theorem OCF (T : Tournament V) : T.H = T.alphaSum := sorry
```
**Difficulty:** HIGH. The Grinberg-Stanley proof goes through ham(D̄) and a sign-reversing involution. Substantial effort.

**Alternative**: take OCF as an axiom (acceptable for a finite-case theorem).

### L11 — Chain constraint

```lean
theorem chain_constraint (G : SimpleGraph V) (k : ℕ) (hk : 1 ≤ k)
  (h : ∃ S : Finset V, S.card = k ∧ G.IsIndependent S) :
  (Finset.univ.filter (λ S => S.card = k-1 ∧ G.IsIndependent S)).card ≥ k := sorry
```

## 4. Main Theorem Skeleton

```lean
theorem THM_343 (T : Tournament V) : T.H ≠ 7 := by
  intro h_eq
  -- Step 1: decomposition uniqueness via OCF + chain constraint
  rw [OCF] at h_eq
  obtain ⟨α, h_alpha⟩ := unique_decomposition_3 h_eq chain_constraint
  -- α = (3, 0, 0, …)
  obtain ⟨C₁, C₂, C₃, h_distinct, h_pairwise⟩ := exists_three_cycles_K3 T α
  -- Step 2: T[V'] is SC, where V' = V(C₁) ∪ V(C₂) ∪ V(C₃)
  let V' := (C₁.vertices ∪ C₂.vertices ∪ C₃.vertices : Finset V)
  have h_SC : IsSC (T.induced V') := SC_union_three T C₁ C₂ C₃ h_pairwise
  -- Step 3: case analysis on s = |V'|
  let s := V'.card
  interval_cases s
  · -- s ≤ 2: impossible since cycle has ≥ 3 vertices
    omega_or_simp
  case 3 =>
    -- only 1 directed 3-cycle, contradicting 3 distinct cycles
    exact absurd h_distinct (n3_unique_cycle T V' h_s)
  case 4 =>
    -- max 2 directed 3-cycles in 4-vertex SC tournament
    have := n4_max_cycles T V' h_SC h_s
    omega
  case 5 =>
    -- Moon: ≥ 3 three-cycles + ≥ 1 five-cycle ≥ 4 odd cycles in T[V']
    have h_three : (T.induced V').cycleCount 3 ≥ 3 := moon_subtournaments _ h_SC 3 (by omega) (by omega)
    have h_five : (T.induced V').cycleCount 5 ≥ 1 := moon_subtournaments _ h_SC 5 (by omega) (by omega)
    -- but α₁ = 3 says T has exactly 3 odd cycles ≥ ones in T[V'], contradiction
    have : 3 < T.alpha 1 := …
    omega
  case _ s_ge_6 =>
    -- Moon: ≥ s-2 ≥ 4 three-cycles in T[V']
    have h_three : (T.induced V').cycleCount 3 ≥ s - 2 := moon_subtournaments _ h_SC 3 (by omega) (by omega)
    omega_or_simp
```

## 5. Effort estimate

| Component                      | Difficulty | Est. LoC |
|--------------------------------|------------|----------|
| Tournament type + basic ops    | Easy       | 100-300  |
| Sub-tournament                 | Easy       | 50       |
| Directed cycle / Hamilton path | Easy       | 100-200  |
| H invariant                    | Easy       | 50       |
| Strong connectivity            | Easy       | 100      |
| SC union lemma (L5)            | Easy       | 50       |
| L4 (cycle is SC)               | Trivial    | 20       |
| Moon's theorem (L6)            | Medium     | 300-500  |
| Camion (L7) [or via L6]        | Trivial    | 5        |
| n=4 score lemma (L8)           | Easy       | 50       |
| c₃ formula (L9)                | Medium     | 100      |
| OCF (L10) [or AXIOMATIZE]      | Hard       | 1000+    |
| Chain constraint (L11)         | Easy       | 50       |
| Main theorem THM_343           | Easy       | 100      |
| **Total**                      |            | **~2500-3500** |

If we **axiomatize OCF** (taking Grinberg-Stanley as a black box), the workload drops to ~1500 LoC, all of it standard tournament theory.

## 6. Critical Path

The minimal-formalization-effort version:

1. Define `Tournament`, `DirCycle`, `H`, `OmegaConflictGraph`, `IndPoly`.
2. Take OCF as **axiom**: `axiom OCF : T.H = T.alphaSum`.
3. Prove L11 (chain constraint) — easy combinatorics on independent sets.
4. Prove Step 1 (unique decomposition) by computation: enumerate `(α₁, α₂, α₃)` with sum-3 weighted ≤ 3 and chain constraint, get `(3, 0, 0)` uniquely.
5. Prove L4 (cycle is SC), L5 (SC union).
6. Prove L6 for **k=3 only** (the simpler "Moon-Moser" form): SC tournament has ≥ n-2 three-cycles. Use Bondy's elegant proof: induct on n via "remove an in-vertex" or score-shifting.
7. Prove L7 (Camion via L6 at k=n, OR via direct argument). Note: at n=5 we need a five-cycle; this is Camion + n=5.
8. Prove L8 (n=4 unique SC score) — finite check.
9. Combine in the case-analysis above.

This gives a "good enough" formalization in maybe 800-1500 LoC of Lean.

## 7. Validation of Audit Findings

The S6 audit (`04-computation/audit_thm343_s6.py`) verified each assumption exhaustively at small n:

| Assumption                        | Verification                  |
|-----------------------------------|-------------------------------|
| OCF (H = I(Ω, 2))                  | n ≤ 6 exhaustive, 0 failures  |
| Chain constraint (Kruskal-Katona) | n ≤ 6 exhaustive, 0 violations|
| SCC partition + cycle ⊆ SCC       | n ≤ 5 exhaustive, 0 failures  |
| n=4 SC unique score (1,1,2,2)     | confirmed                     |
| c₃(T) formula                     | n ≤ 5 exhaustive, 0 mismatches|
| Moon ≥ n−2 three-cycles in SC     | n ≤ 6 exhaustive, 0 violations|
| Camion Hamilton in SC             | n ≤ 6 exhaustive, 0 violations|

These give STRONG numerical evidence. Lean formalization would replace this evidence with airtight proofs.

## 8. Next Steps

1. **Decide axiomatization scope**: prove OCF in Lean, or accept it as axiom.
2. **Build a Lean repo skeleton**: `tournament-redei/` with the type definitions above.
3. **Pull-request candidate**: `Mathlib.Combinatorics.Tournament.Basic` — could grow naturally from this work.
4. **Sub-target**: First formalize Camion (1959) — useful standalone result and a stepping stone to Moon.

---

*This plan is a roadmap, not a one-session task. Estimated 3-6 months of focused work for the full formalization; 1-2 months for the axiomatized version.*
