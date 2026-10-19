# Independent review: conditional paired halving and the stopping point

2026-10-19. **PASS after two narrow boundary/scope corrections.** No all-length annihilator theorem, integer-period theorem, paired-halving theorem, or sharp-horizon theorem is proved by this review.

Reviewed artifacts:

- `GILBREATH_PAIRED_HALVING_STOPPING_NOTE.md`: SHA-256 `8532757d648ee0a97d5efa512be608c6a5fa5b77a25d82679153f546322f2431`.
- `verify_even_front_annihilator.py`: SHA-256 `2018d0688b22fa2f050effa138826c6730a8b136357573625f05b4724e8ab9ae`.

## Corrections applied with root authorization

1. The age-zero parent conclusion is now explicitly conditioned on the even-front annihilator hypothesis. It is false for arbitrary ternary vectors without that property.
2. The bound `finite tau <= Q-2` is now restricted to L >= 1. At L = 0, the initial constant-3 vector has Q = 1 and tau = 0. The uniform conditional bound is `2 floor((Q-1)/2)`.

No other mathematical text or code was changed.

## Phase criterion and counting implication

The annihilator lemma is correct for finite, empty, and one-sided infinite words. Before the first nonbinary output, the running value is `b+h(j)`. A 2 is safe exactly when that value is one; safe passage through a 2 leaves the running bit unchanged. Thus the required phase is exactly `b=1+h(j)` at every 2. If there are no 2s, both bits work. The argument is a local induction and needs no periodicity assumption.

For reachable fronts containing a one, the vector is ternary: before the first one the alphabet is {0,3}, and its first one folds running values into {1,2}; later binary boundary steps preserve the bound two. At an even front length, the open annihilator hypothesis therefore applies to every such parent, including age-zero parents.

For a parent whose first 2 is at r, the binary outputs before r are complementary for the two boundary choices. Exactly one child retains a 2 at r. The annihilator hypothesis makes the other child entirely binary. When r=0 this gives ages zero and infinity, so no new positive-age child enters the count. Binary parents have two binary children. The all-zero front has first entry three, and both children remain dangerous at age zero. These cases exhaust all parents. Hence the conditional prefix-extension bijection `A_r(L+1)=A_r(L)` for r>0 and even L is valid, including L=0.

The normalization is also correct. Exactly F=L+1 fronts have tau zero, because their length-L front word contains at most one one. Combining this with the zero-run probability formula gives `c6(F)=2 sum_(r>0) A_r(F-1) 2^(-r)`. Thus equality of these weighted counts at odd F is precisely `E(F+1)=E(F)/2`. No finite horizon or interchange of an infinite probability sum is needed: each front length has only finitely many words.

## Generating-function indexing and conditional horizon

Write p for the input parity series and q for the output under T_b. The age-zero convention gives

    q(z) = (p(z)+b)/(1+z),       p_0(z)=1/(1+z).

Unrolling in chronological bit order yields exactly

    p_L(z) = [1 + sum_(s=0)^(L-1) b_s (1+z)^(s+1)]/(1+z)^(L+1).

In particular the earliest bit has exponent one, not exponent L. The numerator has degree at most L and value one at z=1, so there is no cancellation and the denominator pole has exact order L+1. With Q the least power of two at least L+1, multiplying numerator and denominator to denominator `1+z^Q` gives a numerator of degree below Q, proving pure parity period Q. Conversely, a period d would require `(1+z)^(L+1)` to divide `1+z^d`; the latter has root multiplicity `2^v_2(d)` at z=1. Hence Q divides d. The least parity period is exactly Q, and the binary linear complexity is L+1.

The full integer period remains a separate hypothesis. If it holds, the first defect lies in the first Q entries. Combining this with the inherited even-delay theorem gives the corrected uniform bound above, or Q-2 for L>=1. It gives an upper bound only; it does not prove attainment.

## Hostiles and verifier audit

The displayed phase hostiles and same-parity examples are valid. The period-eight word has 2s at ages four and six with different cumulative parity phases, giving child first defects six and four. Its parity is the genuine binary trace of front 101100. The isolated age-eight spike preserves that parity and the pair invariant while violating the proposed horizon, so those sidecars cannot establish reachability.

For precision about the two-step counterexample, appending zeros to the initial finite word gives these exact infinite transforms:

    21012 followed by zeros
    -> 21102 followed by constant twos
    -> 21002 followed by repeating 02.

Their prefixes suffice for the claimed failure of annihilation. The later two tails should not themselves be read as zero tails. All three sequences obey the stated pair-language restrictions.

The producer's eventual-period algorithm is sound: its repeated state includes both the input-cycle position and the running value, so cycle detection is exact. Reducing a repeating block and pulling a matching suffix of the preperiod backward preserves the represented infinite sequence. Two cycles suffice for the phase test because a cycle changes cumulative parity by at most one bit. The bounded minimality loop establishes only its stated dyadic-period universe. Normal and optimized producer reruns match the saved receipt.

## Independent controls

A separate verifier, `session6-review/verify_paired_halving_independent.py`, imports no producer code. It certifies a proposed output period by checking two whole blocks while the input period divides the block length. Equal blocks force equal running states at their endpoints, so determinism proves the infinite continuation; this is not a finite-preview guess.

Its 785,592 explicit checks cover:

- all 9,841 finite ternary words of lengths zero through eight, including the exact annihilating bit;
- all 32,767 front words through L=14, exact integer cycle lengths, the displayed generating-function coefficients, parity complexity L+1, and finite first-defect bounds;
- even-front annihilation with age-zero, positive-age, and binary parents distinguished;
- equality of positive-age histograms for every tested adjacent even/odd front-length pair;
- the exact tails of the two-step hostile.

Normal and optimized independent receipts are byte-identical:
`session6-review/paired_halving_independent_receipt.json` and `session6-review/paired_halving_independent_receipt_optimized.json`.
The receipt SHA-256 is `096adb872f54e177f283494965805663b5f37b0ed131343f1887275bf47c5aed`.

The note now states a valid conditional route and a genuine unresolved reachable-state property. The finite controls do not close either open hypothesis.

## Publication clarification

After this review, the main note was clarified to print the exact infinite tails of the two-step hostile, exactly as verified above. Its current SHA-256 is 93d5b3ddd629df9d530bdb3dd0a968e96fee25bbdbcf7be3c2f5645323b9bc29. The historical reviewed hash is retained above. Independent scripts and receipts originally in session6-review are copied into the local review/ directory. No algorithm or theorem changed.
