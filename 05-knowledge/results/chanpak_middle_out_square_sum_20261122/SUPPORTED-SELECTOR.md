# Joint cuts from the middle: exact supported slices

**PROVED.** In the notation of the [two cut gauges](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md), write \(q_-,q_+\) for the two normalized cut increments, so that \(H\) increases by \(2q\). Each \(q_\epsilon\) is a directed cut with exactly one positive polarity per unordered pair. A common one-vertex extension at the next exponent needs the single contact set \(S\) to satisfy
\[
 (q_-(S),q_+(S))=(D,D),\qquad D=(E_{p+1}-E_p)/2.
\]
The proofs also allow unequal target coordinates, in particular \((D+1,D-1)\) if the new lower/upper orientation is allowed to exchange the old one. Start by pinning the contacts at the two endpoints of the old reversal in each of their four possible patterns. Then work outward along a chosen common spanning path, when available, fixing contacts nearest the midpoint of the two marked positions first. This is an ordering for search, not an assumed path reflection on tournaments; correctness is independent of the order.

## Compare both ports, not their midpoint alone

The two positive scalarizations
\[
 U=8q_-+9q_+,\qquad V=9q_-+8q_+
\]
are again directed cuts, even though the normalized signs of the host change at the marked arc. They preserve the full two-coordinate selector: \(U+V=17(q_-+q_+)\) and \(V-U=q_--q_+\). Thus the common target becomes \(U=V=17D\). The exact ratio here is a useful pair of probes, not a theorem privileging those weights. It is the reversal of which coordinate carries the higher weight that matters: a single scalar 8:9 loses information, just as a single square-sum reflection loses the alternating sheet. Other positive probe weights can strengthen rejection.

Here is a polynomial exact rule on a supporting boundary. Fix any pins, take any positive rational \(\alpha,\beta\), and let \(F=\alpha q_-+\beta q_+\). Its minimum over the pinned chamber is one source–sink min-cut (attach forced vertices to the source or sink with capacities exceeding all finite weights). If the minimum is above the projected target, reject. **If the minimum equals the projected target, solve the simultaneous problem exactly in polynomial time:** all exact solutions must be in the minimizing family, and that family is a chain.

Indeed, for incomparable contact sets \(A,B\), the cut defect
\[
 F(A)+F(B)-F(A\cap B)-F(A\cup B)
\]
is the sum of the symmetric capacities between \(A\setminus B\) and \(B\setminus A\). It is strictly positive because each host has a positive capacity on every unordered pair. The two chambers still contain the union and intersection, so two incomparable cuts cannot both minimize. The inclusion-least minimum cut is given by residual reachability. To get the next cut in the chain, retain that entire inside set and force one additional free vertex inside in turn. Among the new minima of the same old value, take the inclusion-least (equivalently, the one with smallest cardinality); if none exist the chain is exhausted. This takes at most \(O(n^2)\) min-cuts overall. The same proof works for scalarization with just one positive coefficient. Merely finding an optimum of the projected target is not enough: check both coordinates on every member of the chain.

A simple complete selector combines these short-chain exits with outward branching. Prune with the four positive tests \(q_-,q_+,8q_-+9q_+,9q_-+8q_+\), or any suitable enlarged bank. Also project the target on signed directions, such as \(q_+-q_-\). For each unordered pair get the minimum and maximum possible contributions under its two pins and sum those independent bounds; the signed interval is valid, though not necessarily tight. When no test resolves a chamber, branch on the next outward contact. This terminates, returning exactly the shared contact sets (take both children if enumerating). It need not be polynomial on interior targets. The code is written with a single witness as its return value.

## The best positive probe can also be selected exactly

Let \(\phi(\lambda)=\min_S[\lambda q_-(S)+(1-\lambda)q_+(S)]-(\lambda d_-+(1-\lambda)d_+)\) on \([0,1]\) for target \(d\). Its maximum gives the sharp lower convex bound available from positive weighting. It can be found by min-cut and rational bisection. If \(B\) bounds the integer range of the coordinates and the target magnitudes, all breakpoints have denominator at most \(4B\): each is an intersection of two integer-affine cut lines. At an interior rational query the least and greatest active slopes are found by lexicographic perturbations: scale the integer capacities for the given probe by more than \(4B\) and add or subtract one unit of \(q_--q_+\) in the secondary objective. Both perturbed capacities remain nonnegative. At endpoints use the inward perturbation. A bracket of width less than \(1/(32B^2)\) contains at most one such reduced-rational breakpoint. Concavity supplies the correct bisection direction; a bracket of active slopes around zero already certifies a maximizer. If neither endpoint nor a queried point is optimal, the maximizing set inside the final bracket contains a breakpoint, which can be recovered by bounded-denominator rational reconstruction. This uses \(O(\log B)\) min-cut calls. If the maximum is zero, its minimizing chain is an exact shared-cut selector; positive rejects, negative leaves genuine discrete work.

## Literal controls

All following weights are actual tournament increments (normalized to half counts), not arbitrary test matrices.

* Rows \((16,1,3,7,14)\), marked reversal \((0,3)\), host counts \(9,11\). Pin 0 in and 4 out. The two-cut minimum sum is 8 and the unique minimizing mask is 15, with \((q_-,q_+)=(4,4)\). The shared contact set raises the pair to \(17,19\).
* Rows \((4,1,2,7)\), marked reversal \((0,3)\), counts \(3,5\). Pin 0 in and 1 out. The sum minimum is 4, but the unique optimizer is mask 1 with cut pair \((1,3)\). There is no common contact for \((2,2)\) in that chamber.
* Rows \((4,1,2,7,15)\), same marked reversal and counts \(3,5\). Each marginal contains \(6\), but **there is no shared contact** for \((6,6)\), even without pins. The witnesses hitting either coordinate are masks \(8,23,27,29,30\), with cut pairs \((6,4),(3,6),(6,9),(6,7),(6,5)\). Here both weighted projections have their maximum bounds above and their minimum below the target; the integer problem is not an independent pair of min-cuts.

See [the replay](check_supported.py). For complementary structural exact reductions see [the split–quotient selector](SPLIT-QUOTIENT.md) and [the affine-lattice selector](LATTICE.md).
