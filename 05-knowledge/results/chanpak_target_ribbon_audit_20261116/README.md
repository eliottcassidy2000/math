# Target inversion, literal span-3 growth, and conflict graphs

Session 2026-11-16. This packet distinguishes exact construction from residue coverage and a literal lower-growth certificate from the verified thirteen-unit certificate. Chan–Pak remains open.

## Choose the factors from the actual target

[The target-first decoder](PRODUCT_INVERSION.md) intersects the exact digital remainder interval with the target to find candidate products, then tests individual factor pairs by **integer equality**, not just congruence. It recovers the three targets in the earlier seam bank. It also gives a width-free ray normal form: a larger spacer can be traded for padding in a width, and one spacer vertex already suffices for the full width-free spectrum of that architecture. A sharp type of limitation remains: **no target 2^m+3 is reached with t≥max(k,r)+2**, for any cut. Spacers above (1+ε)log₂ n give only O(X/(log X)^ε) distinct values up to X. A flexible exact inverse is not a universal selector.

## Bounded growth and rank potential

- **Transitive prefixes:** a final vertex can affect only the last six old labels, as all earlier labels are excluded from every odd directed cycle. Exact finite replay proves the uniform **sharp factor 39** on this entire subfamily. The all-prefix factor-39 bound is a separate question on arbitrary nontransitive ribbons.
- **All literal prefixes, different observable:** an explicit strictly positive row on the 119 raw fragment states has one-step upper growth **39/10**, and dominates H, giving H≤5·3.9^n. The sum of all H readings over extensions of length at most three is another positive row. Uniform one-step growth for H is equivalent to this smoothed readout being uniformly comparable to the actual H on reachable prefixes. The smoothing removes the raw zero but does not establish comparison. [Verified envelopes](POSITIVE_ENVELOPE.md).
- **Rank state:** exact rational rank witnesses have the same suffix (7,0,2) and rank 3 but take different ranks (3 and 2) when appended with the identical letter 0. A literal cone certificate cannot use only those tags. The smaller thirteen-unit cone does replay exactly. [Status and certificate](RANK_CERTIFICATE_STATUS.md). An equal-length strong alias and a counterexample to transitive cone generation are in [Rank budget and exact tests](RANK_BUDGET.md).
- **Density:** the cited wave used a spectral-radius lower bound as an upper bound for a word-count inference. Even the repaired, rigorous word-count lower bound does not imply a lower bound on distinct integers without control of collisions. [Exact power-saving criterion](DENSITY_GATE.md). Literal span-3 power-saving density is not established here.

## One apex versus two

By [Grinberg–Stanley](https://arxiv.org/abs/2307.05569), the conflict-graph independence polynomial at 2 is exactly H, when the vertices are all individual nontrivial odd directed cycles and conflicts mean vertex intersection. Thus universal coverage for physically realizable conflict graphs would be the full conjecture, not a simpler permission to choose an abstract graph. A one-apex clique already reaches every fixed admissible residue (prior Session 25). [A narrower explicit mirrored tape construction](CLIQUE_RESIDUES.md) does this with one common tape per modulus.

[Two-apex surgery](TWO_APEX.md) reverses one arc between two one-apex blocks. It realizes the exact conflict graph (K_a ⊔ K_b) ∨ K_c, with a physical formula for c. This repairs **29, 35, 49, 59, 71** in the *absolute* one-apex spectrum, not previously missed fixed congruence classes. The padding formula controls precisely the linear coefficient while freezing ab. Complete finite replay through 10,000 is kept separate from any density claim. [Independent audit](TWO_APEX_AUDIT.md).

The next proof burden is a genuine absolute selection or collision bound on a growing family. Matching a residue of the target plays the role of finding an admissible carry, not of paying it.
