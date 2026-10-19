# Independent proof audit: fixed-count contact atlases

**Verdict: PASSED.** The full draft
`square_sum_fixed_count_contacts_20261019.md` was checked against the
both-port bounded-section notes, the rational recurrence-bound note, and
the inherited signed-PPT bridge theorem. No theorem-blocking gap was found.
This audit concerns the proof; the companion program's execution and check
counts are separate verification tasks.

## Proof-critical checks

1. **The cusp signs and linear coefficients agree on both ports.**
   Both small translating sections satisfy `w(S)=-h+O(h²)`, where
   `w=-X/Y`. At the plus anchor, `du/omega=-2` at the cusp. At the minus
   negative-sheet anchor, `du/omega=+2`, and the positive output parameter
   is `v=-u`. Thus translation by `[j]S` changes either positive parameter
   by `2jh+O(h²)`, giving `v_(2j+1)=1+(2j+1)h+O(h²)`.
   The limiting anchors `(3,4)` and `(-1,0)` and the identity are smooth
   points of the nodal fiber. The formal law and local translations used
   in the argument therefore remain valid, and both inverse coordinates
   are regular at the cusp.

2. **The common cutoff is effective.** After cancellation, the finite
   rational functions for consecutive differences, the upper endpoint,
   and the output denominators have the stated nonzero values at zero.
   The coefficient-tail estimate in (10) yields a positive rational
   cutoff, including its endpoint. It certifies distinct admissible
   parameters simultaneously for both ports. The square roots are
   nonzero because their radicands are in `(1/2,3/2)`.

3. **The recurrence amplification gives r genuinely different petals.**
   The draft correctly chooses the direction of greater spread width,
   obtaining an arc length `ell>1/20`. A signed displacement less than
   `ell/r` has its first r multiples strictly inside that same arc, with
   no wraparound. Strict spread monotonicity excludes duplicates and the
   anchor. The distinction between the small step and its multiples is
   reflected correctly in the atlas budget
   `M_r=2r floor(20rW_r)`. The positive-sheet switch for the minus-port
   fallback is retained.

4. **Height and forbidden-set claims have the proper quantifiers.**
   The atlas and operation count are finite for each fixed r. Rational
   height inequalities give effective r-dependent polynomial bounds in
   the anchor hypotenuse, with no claimed polynomial dependence on r.
   Primitive normalization and unique spread identify the unordered
   petals correctly. Avoiding fewer than r forbidden petals does not
   introduce their heights into the output bound; the draft separately
   acknowledges their encoded-input comparison cost.

5. **Finite gluing preserves exactly the stated seed data.** Repeated
   fixed-r extension gives the displayed iterated-power bound. Taking
   `H=g lcm(C_1,...,C_k)` makes the original scale-g seed an
   integer-square-scaled copy, by the factor `(H/(gC_1))²`. It does not
   preserve the original integer center or labels. The inherited port
   order yields the stated simple `(2k+1)`-cycle and AP-paired fan as a
   subgraph. The bound in the largest original seed label follows using
   `g>=1` and `S_k>=1`.

## Presentation clarification implemented

The rooted-tree observation is valid by choosing one new petal per
parent constraint and forbidding all existing petals. The draft now
explicitly says it realizes all prescribed tree edges and that additional
contacts may occur, implementing the audit's recommendation. No
induced-contact-tree assertion is proved or needed.

No spread-direction, graceful-rank, fixed-center, interval-coverage, or
multiple-parent contact conclusion was used in the proof or inferred
from it.
