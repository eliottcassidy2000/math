# THM-1445 — The switching-class H-sum is n!

**Status:** PROVED (verified n=3..6 exactly; general proof below).
**Author:** death-star-2026-07-20-S61b (HYP-8300).

## Statement

Let n ≥ 2 and let H(T) denote the number of directed Hamiltonian paths of a tournament
T on the labeled vertex set {1,…,n} (Rédei's number; H(T) is odd). Let the **switching
classes** be the orbits of tournaments under vertex switching (switch at v = reverse all
arcs incident to v); equivalently the cosets of the cut space B = ⟨δ(v) : v⟩ ⊆ F₂^{E(K_n)}
in the space of arc-orientations. There are 2^{C(n−1,2)} switching classes, each of size
2^{n−1}. Then for **every** switching class 𝒞,

  **Σ_{T ∈ 𝒞} H(T) = n! .**

## Proof

Switching at v adds δ(v) (the star of v) to the orientation vector in F₂^{E(K_n)}, so the
switching classes are exactly the cosets of the cut (bond) space B = ⟨δ(v)⟩, dim B = n−1,
|B| = 2^{n−1}.

Fix a directed Hamiltonian path π = (v_{σ(1)} → … → v_{σ(n)}), σ ∈ S_n. Its n−1 consecutive
edges form a spanning tree T_π of K_n (a path is a tree), and π prescribes an orientation
on each of those edges. A tournament D realizes π as a directed Hamiltonian path **iff** D
agrees with π's prescribed orientations on the n−1 edges of T_π (the other C(n−1,2) arcs
are unconstrained).

Restriction to the tree edges, ρ : F₂^{E} → F₂^{E(T_π)}, maps the cut space B
**isomorphically** onto F₂^{E(T_π)}: the fundamental cut-sets of the spanning tree T_π form
a basis of B, one per tree edge, and each fundamental cut contains exactly its own tree edge,
so ρ carries that basis to the standard basis. Hence within any coset D₀ + B (= a switching
class 𝒞) the tree-edge values ρ(D), D ∈ 𝒞, range over all 2^{n−1} combinations, **each
exactly once**. In particular **exactly one** D ∈ 𝒞 agrees with π on T_π — i.e. realizes π.

Therefore
  Σ_{T ∈ 𝒞} H(T) = #{(T, π) : T ∈ 𝒞, π a directed Ham path in T}
                 = Σ_{π} #{T ∈ 𝒞 realizing π} = Σ_{π} 1 = |S_n| = n!,
the sum over the n! directed Hamiltonian paths (vertex orderings). ∎

## Remarks

- **Global corollary:** summing over the 2^{C(n−1,2)} classes recovers
  Σ_T H(T) = n!·2^{C(n−1,2)} (each directed Ham path lies in 2^{C(n−1,2)} tournaments).
- **Odd/even reading.** H is odd-valued (Rédei) and an even function (complement-invariant,
  HYP-534: only even-degree Walsh coefficients). This theorem says the odd values, summed
  over a switching class, give the constant n!. Since switching classes = tilings
  (mac-mini THM-474 / kp THM-1430) and, at **odd n**, correspond to even graphs / two-graphs
  (opus THM-1430), **the even-graph/switching projection of H is constant (= n!)**: all of
  H's tournament-specific (odd-cycle-collection) content lives in the switching-class *fiber*
  (the tiling / 2-adic direction), invisible to the even-graph quotient.
- The mean of H over a class is n!/2^{n−1}, an integer iff n is a power of 2 (v₂(n!) ≥ n−1
  ⟺ s₂(n) ≤ 1); so the identity is genuinely a statement about the *sum*, not a per-element
  average.

## Cross-links
Rédei (H odd) · HYP-534 (H is an even function; even-degree Walsh only) · THM-474 / kp
THM-1430 (tilings = switching classes) · opus THM-1430 (switching graphs = two-graphs = E_n,
even-graph bridge odd-n only) · HYP-3808 (fiber-parity checksum) · S60 two-arithmetics
(the 2-adic/tiling direction carries H's content) · THM-466 (2-adic digits of H = odd-cycle
census) · HYP-8300.
