# All-length exclusion for mixed outward edits across a source gap one

20 October 2026. Status: **PROVED**, with independently audited smooth-content input, cited analytic dependencies, and finite-exact certificates. The complete proof application and all 904 physical residual paths passed a [separate root review](root-review/ROOT_AUDIT.md). No assertion about arbitrary disjoint edits is made.

## 1. Precise statement

Start from a primitive lower zero-intercept mechanical word with p>=3, q=p+r, 1<=r<p, and gcd(p,r)=1. Consecutive source ones have residue labels u+r and u, where

    max(0,p-2r)<=u<min(r,p-r).

Move the earlier one left and the later one right. Their physical adjacent edges are disjoint, both moves are source-legal, and either order realizes the local change 0110 -> 1001.

**Theorem.** Every resulting branch word is nonintegral, for either sign of 2^q-3^p and at all lengths.

The p=2,q=3 residue boundary is not a path in this family: the two outward edges target the same zero and cannot both be executed. It is excluded by physical legality, not by discarding a genuine negative cycle.

## 2. Large-weight reduction

For the actual content g=gcd(|2^q-3^p|,B), the [smooth common-root content theorem](../collatz_two_swap_outward_gap_one_content_20261020/OUTWARD_GAP_ONE_CONTENT.md) gives

    g<=22^(4p^(2/3)).

Its proof is specific to this legal family: it uses both 2-adic and 3-adic valuations, algebraic norms, and the inequalities u<p-r and r>p/3 to rule out a common algebraic component of two independent lattice relations. This is not the external-prime argument, which does not apply to the smooth edit constant -4.

For p>=2744=14^3,

    log(g/3^p)
      <=4p^(2/3)log22-p log3
       <(19/21-13/12)p=-5p/28<-p/9.                (2.1)

Here log22<19/6 is certified by 21 rational positive Taylor terms for exp(19/6), and log3>13/12 by the first two positive artanh terms.

If the output were integral, then g=|D| for D=2^q-3^p, so the nonzero multiplicative form would satisfy

    |2^q3^(-p)-1|<exp(-p/9).

The [audited absolute-gap certificate](../collatz_two_swap_endpoint_exclusion_20261020/ENDPOINT_TWO_EDIT_CLASSIFICATION.md), with its cited explicit real Matveev theorem and exact both-sign continued-fraction checks, excludes this inequality whenever p>=100 and p<q<2p. That arithmetic lemma depends on the absolute gap, not on which edit family produced it. Therefore every possible integral output has p<2744.

## 3. Complete small-weight candidate reduction

The [independently verified signed-descent certificate](../collatz_two_swap_overlap_exclusion_20261020/OVERLAPPING_TWO_SWAP_CLASSIFICATION.md), Section 3, proves that every positive magnitude 1<=m<=N=2^24 enters the cycle 1,2, and every negative magnitude in the same inclusive range enters one of the magnitude cycles 1; 5,7,10; or 17,25,37,55,82,41,61,91,136,68,34. Both a checked unsigned-128-bit accelerated implementation and a materially different checked unsigned-64-bit full-half-step implementation executed that complete range, with fatal overflow and termination checks. The present proof imports that frozen finite theorem; it does not claim a new global verification range or rerun an unreported larger range.

Any other integer cycle has integer minimum absolute magnitude M>=N+1. Multiplication of the odd-step ratios and log2>2/3 give, for either sign,

    |q/p-log3/log2|<1/((3M-1)log2)
                  <3/[2(3N+2)]<1/(2N).

Since p<2744<4096 and N=4096^2, this is strictly below 1/(2p^2). Coprimality forces q/p to be an ordinary continued-fraction convergent of log3/log2.

Exact rational intervals give the complete candidate list below 2744:

    (p,q)=(2,3),(5,8),(12,19),(41,65),
          (53,84),(306,485),(665,1054).

The next denominator is 15601, from 24727/15601. Eighty-term positive artanh series certify every needed digit. The (2,3) shape has no legal paths here and is kept as a boundary control.

The known cycles are treated separately, not subjected to the false minimum premise. Their primitive shapes are (1,2),(1,1),(2,3),(7,11). Coprime p,q excludes multiple traversals. Only (7,11) adds a shape to the above list; it is the known negative -17 cycle and is not a continued-fraction candidate.

## 4. Complete actual-path residual

`verify_small_outward_candidates.py` constructs both orders of every legal pair and independently enumerates physical gap-one outward edge pairs. It includes original-cut crossings, computes the carry in two ways, and tests the full gcd rather than a resultant surrogate.

| p | q | Ordered actual paths | Maximum actual content |
|---:|---:|---:|---:|
| 2 | 3 | 0 | not applicable |
| 5 | 8 | 4 | 1 |
| 7 | 11 | 6 | 1 |
| 12 | 19 | 10 | 1 |
| 41 | 65 | 34 | 551 |
| 53 | 84 | 44 | 11 |
| 306 | 485 | 254 | 1 |
| 665 | 1054 | 552 | 37 |
| total | | 904 | |

All 904 ordered outputs are nonintegral. In particular, the six paths at the known -17 shape do not produce that cycle. The nontrivial partial factor 551 at p=41 is retained as actual content, not confused with full divisibility. Normal and optimized Python give identical receipts with 3,653 explicit checks.

Together with (2.1), the finite signed minimum proposition, the known-cycle handling, and complete rational candidate exhaustion, these checks close the exact family stated in Section 1.

## 5. Dependencies, reproduction, and limitations

The separate content theorem remains independently meaningful. The present exclusion also requires the cited Matveev estimate, the frozen absolute-gap interval certificate, the frozen signed-descent theorem through 2^24, and the new complete 904-path residual. A failure of any stated dependency would affect the corresponding implication; none is replaced by finite pattern evidence.

Reproduce the new exact checks from this directory:

    python3 verify_small_outward_candidates.py > receipt.replay.json
    python3 -O verify_small_outward_candidates.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

The result concerns primitive k=1 mixed outward swaps across a source gap one. It does not assert all disjoint edits, nonneighboring source labels, or general Collatz convergence. Any broader structural classification must state and prove its own partition and retain repeated-source exceptions.
