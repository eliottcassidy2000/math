# Independent audit of complementary offsets and exact observer additions

Status: VERIFIED. The principal proofs and computations in COMPLEMENTARY_OFFSETS.md, the uniform-substitution and full-observer sections added to NOTE.md, and CLOSED_OBSERVER.md are valid. Scope limits concerning the spectrum conjecture and ternary affine source closure are correctly stated.

## Complementary offsets

The one-vertex C/D transfers follow from the two exact starting-path counts and the deletion defects. Their temporal compositions CDDC and CDDD give the stated profiles and two consecutive equal scalar offsets. The guard assumptions and strongness conditions suffice, including the zero-defect order-three boundary.

Protected left/right pairs remain protected under both operations. Distinctness of q,b,s,c is used correctly to retain the opposite witness. The old unique-in/unique-out guard incidences also survive. L and R commute up to the two new labels. Separate affinity in 2^k−1 and 2^l−1 proves the bilinear formula; its four corner counts determine the coefficients.

The named source rows, final rows and observers are correct. The family 2(2^(k+1)+1)(2^(l+1)+1)−1 and the induced strong chain 5→17→53→161 follow. The +4-family warning is correct: its double-sided (2,2) count is 307, not 295.

Source attribution was checked against S163 FIRST_STEP.md at d9d9f3e483c54d83ebd569b91a85b7d015923cfe. The prior obstruction concerns repetition of its same two-vertex attachment; the new second attachment genuinely differs and does not contradict it.

Independent finite checks in extended.py cover:
- 814 unique-in-neighbor guard cases among all labeled strong tournaments of orders three through five.
- 40 protected-pair sources with five additional lattice points per source.
- Both named two-parameter families at all 0<=k,l<=3.
- Both named complementary offset and fixed-observer examples.

## Uniform substitution and full-observer valuations

For an ordered cover by k nonempty paths, freely permuting the k slots proves k! divisibility. Fixing a quotient word of module runs makes the internal choices factor as the product of these ordered cover counts. Thus the gate becomes a singleton of weight 2t modulo three, exactly as stated. The same argument applies to both observer-query graphs when the port is outside all replaced modules, because each query vertex sees each gate uniformly.

Each inactive gate contributes a factor divisible by three to every run-word summand, giving 3^r divisibility of all query counts. Differences and division by two retain it; the exact first-coordinate valuation prevents additional common factors of three. Independent checks cover both wrapper ports for all s=1,2 settings and all sixteen strong two-gate substituted observers.

## Exact wrapper compression and closed operations

All six listed collapsed-run word lists were independently regenerated and match. Their counts yield query values P+R, P+2R, 3P+3R, hence observer (P+R,R/2,P+R/2). The exact formula was independently tested on every one of the 75 nonempty labeled internal tournaments through order four.

Appending B followed by T1 scales (P,R) by (9,168) or (11,200). Both displayed observer matrices are algebraically correct. The sources remain strongly connected induced extensions with the original physical port retained, and every output remains in the same wrapper class. These facts justify repeated operations, not merely one-step coincidences.

The operations commute on their induced exact profile/observer maps. They need not yield isomorphic fixed-port tournaments when the different strong gate modules are placed in the opposite ordinal order. Publication wording should distinguish these statements.

The closed-family count of s+1 observer states, assignment multiplicities and gcd valuation law follow as stated. Independent literal two-step tests give:

- inactive/inactive: (112977,56448,56529);
- inactive/active or active/inactive: (134499,67200,67299);
- active/active: (160121,80000,80121).

Their strongness, induced original four-vertex source, and valuations all pass.

## Packaging and wording notes

- COMPLEMENTARY_OFFSETS refers to GUARD_ALGEBRA.md, which was absent from its submitted temporary directory. Ensure the published target exists, or remove that reference; the formulas used in this document are self-contained.
- Replace the ambiguous phrase “Because 2 is invertible modulo 3^r, it divides each observer coordinate” with “Since 2 is invertible modulo 3^r, 3^r divides each observer coordinate.”
- Describe CLOSED_OBSERVER's commutation at the exact-profile/observer level, rather than claiming unproved commutation of the literal fixed-port tournaments.

These are clarification/packaging issues, not failed mathematical claims.

## Reproduction

Place independent.py, extended.py, and closed.py together. Run `python3 extended.py` and `python3 closed.py`, and repeat with `python3 -O`. No submitted core/check modules are imported. The updated original check.py and offset_check.py also pass both modes.
