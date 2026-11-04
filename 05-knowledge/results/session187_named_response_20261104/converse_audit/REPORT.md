# Independent audit: SCC-local response reduction

**PASS, within the stated scope.** Reviewed `../converse/PROOF.md`, both author scripts, and both author receipts. The theorem is an exact specialization of the full named response, not a projection onto endpoint moments or an assertion of intrinsic descent through a nontrivial strong quotient.

## Proof review

The quotient condensation of a tournament is a transitive tournament. With every earlier SCC pointing into both ports and both ports pointing into every later SCC, the literal extended tournament is an ordered sum with the active extension as its middle block. Hamilton paths concatenate uniquely across ordered-sum blocks. The factor is therefore exactly the product of Hamilton counts of the outside SCCs. Positivity suffices for cancellation; the additional oddness statement is the standard tournament Hamilton-path parity theorem and is not needed here.

The same lifted assignments work for both same-named substituted cores because the outside quotient and active physical names are unchanged. Hence full ambient equality implies full active-block equality. When the root's quotient SCC is a singleton, the active block is precisely the core; otherwise it is the whole strong quotient substitution, not the core alone. The proof explicitly leaves this distinction intact. Neither intrinsic strongness nor strongness of the extended block is needed. Arbitrary exterior ordered padding cannot cancel an active witness, but arbitrary nonordered exterior contexts are not canceled by this argument.

The inherited C2 homogeneous-pair obstruction is credited rather than re-proved or expanded to the prime residual. The packet does not claim that a prime twelve-core is excluded. The full cubic theorem is used only for whole converse; arbitrary internal reversal retains the full quartic basis.

## Independent computation

`verify.py` in this audit directory uses an independently written predecessor subset DP and explicit adjacency-matrix extension construction. It recomputes quotient reachability with graph search, selects the root SCC independently, and reconstructs physical lifted masks. It checks five queries in each of the 150 author cases, and every one of the sixteen literal eleven-vertex hostile terms: **766 independent checks**.

For the active hostile, a separate literal Hamilton-word polynomial expansion from the sibling collision checker verifies every coefficient in both the arbitrary internal-reversal quartic and the whole-converse cubic. This is independent of the author's finite-difference interpolation. Directed reachability independently verifies that the active five-vertex tournament is strong. Normal and optimized audit receipts agree.

## Corrected reporting issue

During audit, the original `sccs` reporting function sorted a list with a key that iterated that same list. Python makes a list appear empty during in-place sorting, so the original component order was not topological. Component membership, chosen active block, lifted masks, and factor computations did not depend on that order. The author independently caught and fixed this by freezing a copy before sorting, and reran the receipts. This was a reporting-order issue, not a counterexample to the theorem or arithmetic.

The all-order conclusion rests on the ordered-sum proof. These tests are corroboration, not an unbounded reconstruction census or a hidden-converse exclusion.

The final added S185 gate-padding section was also reviewed: both scalar and complete localized-response factors follow from the same ordered-sum argument, and the active cubic is correctly distinguished from unrestricted ambient reversal.
