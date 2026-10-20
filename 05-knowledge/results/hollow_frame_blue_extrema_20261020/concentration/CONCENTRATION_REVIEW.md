# Concentration proof review

Status: **PROVED; INDEPENDENT PROOF REVIEW PASSED**, 2026-10-20.

## Reviewed mathematics

The author and a separate mathematical reviewer independently obtained
the same coefficientwise majorization

    (1+t)/(1-t) <= (1-2t)^-1

in the order on nonnegative formal-series coefficients. The reviewer then
read and approved the complete proof, including:

- preservation of coefficientwise order under positive integer powers;
- the bound for every allowed j, not just fixed j;
- the exact j=1 contribution and the summed O(k^-2) remainder;
- the L^2 and probability-concentration conclusions;
- strictness of E[h^2]/E[h] < max h for k>=2, using positivity and
  nonconstancy of h;
- the least-odd-integer-strictly-above rounding, including Q_3=7.

Approval applies to the proof in `HOLLOW_FRAME_CONCENTRATION.md` and its
bounded consequences.

## Exact validation

The author's verifier computes coefficients by both repeated series
multiplication and the independent composition formula. It tests the
termwise bounds, finite and rational tail bounds, inherited exact
variances through k=7, and odd rounding. Runs through k=128 and k=256
passed; the latter used `python3 -O` and checked all 16,640 allowed
summands. Only the larger run is retained as
`concentration_checks.json`.

The verification uses exact integer/rational arithmetic. Its displayed
decimals are descriptive only. The finite computations support the
implementation; the reviewed proof establishes the all-k theorem.
No independent computational rerun is claimed by this review record.

## Scope and preservation

This supplement is add-only. It modifies neither the main extrema
packet nor its second-moment supplement. It makes no new exhaustive-frame
claim, central limit theorem claim, exact-maximum conjecture, or
literature-priority claim.
