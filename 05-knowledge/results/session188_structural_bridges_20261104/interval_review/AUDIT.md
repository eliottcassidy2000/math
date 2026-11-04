# Independent audit of arbitrary interval reversals

**Verdict: PASS. No discrepancy found.**

I reviewed `interval_modules/PROOF.md` and `verify.py`, and inspected the inherited S187 named transitive-order proof used by the arbitrary replacement theorem. Producer normal and optimized runs agree and match the supplied normal receipt.

## All-order proof

The induced-subset injection is valid. A subset of the source ordered sum has one induced piece per occupied block and exactly the product of their Hamilton counts. Choose a fixed Hamilton order in the target induced quotient, concatenate one internal Hamilton path from each piece, and reverse internal paths when that block was conversed. Distinct choices remain distinct after the concatenation. Applying the resulting nonnegative inequalities to every part of every ordered partition proves coefficientwise domination of all ordered path-cover counts.

If the target quotient is transitive, it is another ordered join of those same pieces, so every subset product and every cover count agrees. If it is nontransitive, a directed triangle's three Hamilton paths survive deterministic vertex insertion, giving at least three quotient Hamilton orders. The h=1 profile gain is therefore at least twice the source block product. Nonnegative ambient coefficients and c_1=H(Q)>0 give the strict ambient scalar bound. This handles arbitrary fixed ambient quotients, not only ordered contexts.

The two-dimensional finite-difference formula correctly inverts the sum over intervals containing each pair. Boundary terms have increasing indices unless explicitly outside the range. There are exactly binom(k,2) nonsingleton intervals, so the map is a bijection between interval parity patterns and named quotient tournaments. Exactly k! labeled tournaments are transitive. This count includes the empty pattern among the scalar-preserving cases, and correctly excludes singleton-block intervals that would create additional internal choices.

For arbitrary replacement of a transitive source module, a nontransitive replacement has a strict scalar increase. A different transitive order is handled by the inspected S187 universal named rank-query theorem, whose full cover profile increases coefficientwise and has a strict h=1 contribution in every quotient. Thus the two cases exhaust every different same-name replacement. This is not extended to an arbitrary nontransitive source module.

## Independent finite checks

`independent.py` does not import the producer verifier. It checks interval inversion for every pattern through k=6 and independently recovers the k! scalar-preserving quotient counts.

It computes induced-subset Hamilton counts by literal permutations, then ordered partition profiles. Across 320 quotient/block/converse controls, including nontransitive internal blocks, it verifies profile domination, equality precisely for transitive quotients, and the quotient-Hamilton multiplier.

For both strong-ambient named examples, a separate literal Hamilton-word engine reconstructs the entire Boolean polynomial and all 1,024 actual query values. The crossing contrast has 144 nonzero coefficients and degree four; the nested contrast has 136 nonzero coefficients and degree four. Every saved coefficient matches, with all unlisted coefficients zero. The whole-converse control retains all 124 nonzero coefficients and has no quartic term. No comparison uses a sorted coefficient multiset or drops high-degree coordinates.

Normal and optimized independent outputs match byte for byte. Artifacts are `independent.py`, `independent.json`, `independent.optimized.json`, and the paired producer replay receipts.
