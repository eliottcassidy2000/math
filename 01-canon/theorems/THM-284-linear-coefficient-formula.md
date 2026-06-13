# THM-284: Linear Coefficient of H in Tiling Bits = 2^(skip-1)

**Status:** PROVED  
**Source:** opus-2026-04-03-S28  
**Proved for:** all n ≥ 3 (verified n=3..8, with combinatorial proof for all n)

## Statement

Let H(x₁, ..., xₘ) be the Hamiltonian path count expressed as a multilinear polynomial in the tiling bits, where the base tournament is transitive (all bits = 0) and tile j = (xⱼ, yⱼ) has skip sⱼ = xⱼ - yⱼ ≥ 2.

Then the linear coefficient of xⱼ is exactly 2^(sⱼ - 1).

Equivalently: H(transitive tournament with arc (x,y) reversed) = 1 + 2^(x-y-1).

## Proof

The transitive tournament on vertices {1, ..., n} has the unique Hamiltonian path n → n-1 → ··· → 1 (giving H = 1).

Reversing the arc from x→y to y→x (where x - y = s ≥ 2) preserves the original HP (which uses only base-path arcs n→n-1, ..., 2→1). The new HPs must all use the reversed arc y→x.

Any new HP has the form:
  [descend through some subset of "between-vertices" to y] → y → x → [descend from x through remaining vertices]

The "between-vertices" are {y+1, y+2, ..., x-1}, of which there are s-1. Each can be visited either before the jump (on the way down to y) or after the jump (on the way down from x). These choices are independent, giving 2^(s-1) new HPs.

Since no HP uses the original arc x→y (it used x→x-1 instead), the total is H = 1 + 2^(s-1). ∎

## Corollary: Cycle Count from Transitive

By the OCF formula H = I(Ω, 2), this implies the single-flip tournament has exactly 2^(s-2) directed odd cycles, all forming a clique in Ω(T) (sharing the lower vertex y).

## Consequences

- The apex tile (n,1) has the largest linear coefficient: 2^(n-2)
- The main effect of each tile on E[H] is proportional to 2^(skip-1)
- The tile effect hierarchy follows the perpendicular distance from the base path
- Grid-reflected tile pairs (x,y) and (n+1-y, n+1-x) have the same skip, hence the same linear coefficient, explaining the effect symmetry
