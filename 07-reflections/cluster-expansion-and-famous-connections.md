# Cluster Expansion and Famous Problem Connections

**Session:** opus-2026-05-28-S3b

## The SC Formula IS the Mayer Cluster Expansion

The inclusion-exclusion formula for SC tiling counts:
  SC(n) = Σ_{S ⊆ cuts} (−1)^|S| · 2^{f(S)}

is exactly the **Mayer cluster expansion** (polymer gas partition function) from statistical mechanics.

In the Mayer expansion for a lattice gas:
  Z_conn = Z_total · exp(−Σ_polymers φ(polymer))

Our situation:
- Z_total = 2^m (all tilings)
- Z_conn = SC(n) (strongly connected tilings = "connected" part)
- "Polymer" = a set S of cuts with weight (−1)^|S| · 2^{f(S)} / 2^m

The ratio SC(n)/2^m → 1 exponentially, which means:
  Σ_S φ(S) → 0 exponentially

The leading term: single-cut polymers at boundary positions k=1 and k=n−1 each contribute −1/2^{n−2}. Together: −2/2^{n−2} = −1/2^{n−3}.

So SC(n)/2^m ≈ exp(−1/2^{n−3}) ≈ 1 − 1/2^{n−3}, confirming the non-SC ratio ≈ 8/2^n.

**This is the Mayer first virial coefficient!**

## Bivariate GF as Polymer Partition Function

The bivariate GF F(x,y) = xB(xy)/(1 − x − xB(xy)) is the grand canonical partition function of a 1D polymer gas:

- "Monomers" (vertices not in any SC block): fugacity x/(1−x)
- "Dimers/polymers" (SC blocks of width ≥ 2): fugacity B(xy) per block

This is a **hard-core 1D polymer model**: SC blocks cannot overlap (they cover disjoint cut intervals). The partition function factorizes as:
  F = (monomer fugacity) × (polymer fugacity) / (1 − polymer fugacity × monomer)

This is formally identical to the **Grand canonical ensemble** of a 1D gas with hard-core interactions.

## Connection to Bloom-Sawin Sum-Product (arXiv:2605.28781)

The paper disproves the Erdős-Szemerédi sum-product conjecture for ℝ using:
- **Additive structure**: totally real number field units embed as a lattice in ℝ^d
- **Multiplicative structure**: units form a multiplicative group

For our tournament H-values:
- **H(T) is always ODD** (Rédei's theorem)
- **H(T)+H(T')** is always EVEN (odd+odd)
- **H(T)·H(T')** is always ODD (odd×odd)

This **automatic parity separation** means H-values live in a stronger structure than what Bloom-Sawin constructs. The H-spectrum is "automatically" in different additive vs multiplicative residue classes.

The natural analog of Bloom-Sawin in our setting: can the independence polynomial map T → I(Ω(T), x) exhibit sum-product compression? The "sum" is the additive combination of polynomials, the "product" is the product polynomial (= disjoint union of conflict graphs). The evaluation at x=2 gives H(T).

**New question** (HYP-1750): Is there a sequence of tournaments T_1,...,T_N such that {I(Ω(T_i), x)} is simultaneously small in additive size (few distinct polynomials mod some ideal) and in multiplicative size (few distinct products)? This would be a "sum-product" phenomenon in the polynomial ring Z[x].

## Connection to Permanent Theory (van der Waerden / Bregman-Minc)

H(T) is the count of Hamiltonian paths = evaluation of a permanent-like formula. The IE formula for SC(n) has the same structure as Ryser's formula for permanents:

  perm(A) = (−1)^n Σ_{S⊆[n]} (−1)^|S| Π_i (Σ_{j∈S} A[ij])

Our IE = Ryser applied to the tiling model's "block matrix" structure, where the "row product" becomes 2^{f(S)}.

The **van der Waerden conjecture** (now Egorychev-Falikman theorem) gives a lower bound for permanents of doubly stochastic matrices. The analog for our SC tiling sequence: a "lower bound" for SC(n) in terms of 2^m. We proved SC(n)/2^m → 1, which is actually an UPPER bound (SC ≤ 2^m, equality impossible). The gap is exactly measured by the polymer expansion.

## Zero-Sum Theory (Erdős-Ginzburg-Ziv)

**New discovery**: At n=5, H(T) mod 5 never equals 2!

From the OCF with i₂=0 at n=5: H = 1 + 2i₁ (sum of all odd cycles).
- H ≡ 2 (mod 5) requires i₁ ≡ 3 (mod 5)
- But no 5-vertex tournament has exactly 3 (or 8, 13, ...) odd directed cycles

This is a MODULAR CONSTRAINT on the independence polynomial: I(Ω(T), 2) ≢ 2 (mod 5) for n=5.

The EGZ connection: any sequence of 9 = 2·5−1 five-vertex tournaments contains 5 whose H-values sum ≡ 0 (mod 5). Since H ∈ {1,3,5,9,11,13,15}, the mod-5 values are {0,1,3,4}. By EGZ (in Z/5Z), any sequence of 9 values from {0,1,3,4} contains 5 with sum ≡ 0 (mod 5). This is guaranteed but the missing residue 2 makes the constraint tighter.

**HYP-1751 (Zero-sum)**: For any prime p, the H-spectrum H_p = {H(T) : T is p-vertex tournament} misses some residue class mod p.

At n=5, p=5: missing residue is 2.
Prediction: at n=7, p=7: H mod 7 might miss some residue(s).

## Summary of Famous Problem Connections

| Famous Problem | Connection to Our Work | Status |
|---|---|---|
| Mayer/cluster expansion | SC formula = polymer gas | PROVED (THM-341, C3) |
| van der Waerden / Ryser | SC IE = permanent formula | IDENTIFIED |
| Bloom-Sawin sum-product | H-values have parity separation | OBSERVED |
| EGZ zero-sum | H mod p has gaps (e.g., H≢2 mod 5 at n=5) | COMPUTED |
| Lee-Yang zeros | I(Ω,x) always positive at x=2 | TRIVIAL (coeff positive) |
| #P complexity | H(T) via tiling = poly for circulants | FOLKLORE |
