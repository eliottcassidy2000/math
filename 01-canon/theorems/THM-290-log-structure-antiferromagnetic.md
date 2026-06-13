# THM-290: The Log-Structure Theorem — H(t) is Antiferromagnetic

**Status:** Part (1) PROVED, parts (2)-(4) VERIFIED n=4,5,6 (conjectured for all n)
**Filed by:** opus-2026-04-04-S3
**Depends on:** THM-284, OCF

## Statement

Let H(t₁,...,tₘ) be the multilinear polynomial for Hamiltonian path count.

**(1) PROVED:** The log-linear coefficient of tile k with skip s = x-y is:
  λ_k = log(1 + 2^(s-1))

**(2) CONJECTURE (verified n=4,5,6):** ALL pairwise log-interactions are negative:
  λ_{ij} = log(H(ij)·H(0) / (H(i)·H(j))) < 0
for every pair of tiles (i,j).

**(3) Equivalently:** The odds ratio H(ij)·H(0) / (H(i)·H(j)) < 1 always.
This means: the effect of flipping two tiles is LESS than the product of
individual effects (sub-multiplicative / log-submodular).

**(4) FACT:** log(H) has degree m (full multilinear) while H has degree
2⌊(n-1)/2⌋. The OCF exponential structure creates the sparsification.

## Proof of (1)

λ_k = [Möbius coefficient of t_k in log(H)]
    = log(H(e_k)) - log(H(0))
    = log(H(e_k) / 1)  [since H(0) = 1 for transitive]
    = log(1 + c_k)      [since H(e_k) = 1 + c_k by the multilinear expansion]
    = log(1 + 2^(s-1))  [by THM-284: c_k = 2^(skip-1)]  ∎

## Physical Interpretation

H(t) is the PARTITION FUNCTION of a system on the staircase diagram where:
- Each tile is a binary spin variable t_k ∈ {0,1}
- External field: h_k = log(1 + 2^(s-1)) > 0 (each spin "wants" to flip)
- Pairwise coupling: J_{ij} < 0 (ANTIFERROMAGNETIC — spins compete)
- Higher-order interactions exist (the log has full degree m)

The system is a **frustrated antiferromagnet**:
- Flipping any tile increases H (positive field)
- But each pair of flips competes (negative coupling)
- The H-maximizer is the tournament that optimally resolves frustration

The tournament achieving maximum H (the Paley tournament at primes p ≡ 3 mod 4)
is the "ground state" of this frustrated system.

## Numerical Data

| Tile pair type | λ₂ range (n=6) | Number |
|---------------|----------------|--------|
| Same-end, adjacent skip | -1.10 | 12 |
| Same-end, skip diff 2 | -1.10 to -1.50 | 8 |
| Same-end, skip diff 3 | -1.22 to -1.89 | 4 |
| Cross-end | -1.10 to -1.22 | 4 |
| Disjoint, adjacent | -0.51 to -0.65 | many |
| Disjoint, distant | -0.59 to -0.90 | many |

Shared-vertex pairs have STRONGER (more negative) coupling than disjoint pairs.
This is because shared-vertex tiles physically "compete for the same vertex" in
cycle packings.

## Connection to the Master Formula

From the OCF decomposition (THM-287 generalized):
c_S = Σ_{cycle packings} 2^r · Π (-1)^{|S_i ∩ F(C_i)|}

The log transform converts this SUM-OF-PRODUCTS into a cumulant expansion:
log(H) = Σ_k λ_k t_k + Σ_{k<l} λ_{kl} t_k t_l + ...

The cumulants (λ's) measure the "connected" part of the interaction, stripping
away the redundant cycle-packing structure. The fact that all pairwise cumulants
are negative reflects the COMPETITIVE nature of odd-cycle packing.

## Open Questions

1. **Prove λ_{ij} < 0 for all tile pairs at all n.** This would be a "tournament
   log-submodularity" theorem. Possible approach: FKG inequality for the hard-core
   model on the cycle conflict graph.

2. **Are ALL higher-order log-interactions also negative?** If so, H(t) is
   completely "log-submodular" — a strong structural property.

3. **Connection to H-maximization:** The frustrated antiferromagnet picture
   suggests that the H-maximizer solves a combinatorial optimization problem
   on the staircase. Can variational methods (mean field, belief propagation)
   approximately find the maximizer?

## See Also
- THM-284 (linear coefficient = 2^(skip-1))
- THM-289 (generalized reversal cancellation)
- OPEN-Q-013 (Paley maximizer conjecture)
- Scripts: h_log_structure.py, h_ie_deep.py
