# Independent audit of the ribbon second checkpoint

Status: PASS for the explicit infinite family and the stated finite-box obstruction. No all-target coverage proof, general no-affine theorem, or independent replay of the full census is claimed.

## Source grammar and boundary

The literal tournament has one vertex per letter. For i<j, the edge points j→i exactly when j−i is within the span and the corresponding bit occurs in letter j; otherwise it points i→j. Backward bits aimed before the first vertex are ignored. The matrix order is consequently A_last ... A_first acting on the stored initial column, followed by the stored output row. The audit checks this order rather than reversing the word or treating its cyclic rotations as equal physical seeds.

The source files ribbon.py, intrep.py, span2_intrep.txt and span3_intrep.txt are byte-identical to commit acbad6c8636cf8cbb675b4c947fa19e33ac0df23 under 04-computation/experiments/ribbon_spectrum_20261104/. Their SHA256 values are recorded in SOURCES.json. The exact representation itself is inherited, not newly proved by a few sample counts.

The representation's empty-word count is zero. The displayed family formula would give one at k=0, so excluding k=0 is essential and correct. The proof uses the genuine initial vector for a fresh tournament, not an arbitrary eigenvector or a steady-state boundary.

## Positive infinite family

For each block 02032333 and 02133232, independently loaded exact integer matrices give rank(Q)=2, characteristic polynomial x^12(x−265)(x−2), and Q^3=267Q^2−530Q. A separate predecessor-based subset dynamic program, with literal adjacency built directly from the grammar and no imported graph-count routine, gives H(B)=65 and H(B^2)=16825, for both blocks (8 and 16 vertices).

The matrix identity implies h Q^(k+2) v=267 h Q^(k+1) v−530 h Q^k v for all k>=1. Its distinct roots are 265 and 2. Solving the two initial conditions yields exactly

    H(B^k)=(63*265^k+200*2^k)/263, k>=1.

The audit additionally verifies divisibility and exact matrix evaluation for k=1,...,12. Both coefficients are nonzero; the physical seed therefore does not isolate the doubling eigenspace. This family is a positive exact construction, not a doubling induction or a coverage statement. No strong-connectivity claim is inferred.

## Finite census and arithmetic safety

The producer's verifier enumerates every nonempty block in the stated boxes. Their cardinalities are sum(4^j,j=1..9)=349524 and sum(8^j,j=1..5)=37448. Prefix products use left multiplication and reduction modulo the screening prime. A nonsingular Q−2I modulo a prime is necessarily nonsingular over the rationals. Candidates are screened at a second prime and surviving candidates receive exact determinant/characteristic-polynomial checks. Comparison against saved first-prime candidate lists is asserted after fresh enumeration, so the verifier does not simply trust those lists.

The audit reviewed that loop and independently rechecked all 50 and 30 stored first-prime candidates with separately implemented modular Gaussian elimination at primes 1009 and 1013. This leaves respectively 16 and 0 second-prime survivors. Each of the sixteen listed span-2 cyclic rotations was independently checked to have eigenvalue 2 and not eigenvalue 1. The complete noncandidate enumeration was not independently rerun; its finite-exact status relies on the retained full producer run and inspected verifier. This limit is explicit in RESULT.json.

Signed-int64 safety was checked from the actual source matrices. The span-2 dimension/max absolute letter entry are 14/5, and the span-3 values are 62/4. With residues <=1012, raw matrix-product sums are bounded by 70840 and 250976. Modular elimination intermediate subtraction/product magnitudes are bounded by 1025156. Thus neither prime's screen can silently overflow int64. The exact family calculations use arbitrary-precision SymPy integers.

## Scope of the no-affine consequence

If H(u B^k)=A 2^k−c for every k>=0, with A=H(u)+c nonzero and c nonzero, its scalar minimal recurrence polynomial is exactly (x−2)(x−1). Any matrix Q representing this sequence has its characteristic polynomial as an annihilator by Cayley–Hamilton. Hence both 1 and 2 must be eigenvalues of Q. This necessity is independent of the prefix's reachability details and remains true even if one only starts the recurrence after a fixed finite transient.

The finite screen finds no such spectral pair in the enumerated block boxes: all sixteen eigenvalue-2 span-2 survivors lack eigenvalue 1, and there are no span-3 eigenvalue-2 survivors. Therefore the claimed nonconstant fixed-block affine doubling sequence is excluded there, even on a constrained set of physical prefix states. The assumptions c!=0 and A!=0 cannot be dropped: pure doubling needs only eigenvalue 2, and constant sequences need only eigenvalue 1. Larger blocks, state-dependent choices between blocks, nonlocal edits, and collective reachable-set arguments are not excluded. The earlier eigenvalue-1 simplicity obstruction to arithmetic progressions is not a proof of this new finite-box conclusion.

## Reproduction and controls

Run `python3 audit/ribbon/verify_audit.py` in a published packet containing the sibling ribbon directory. The checker also detects the original source-workspace layout. RESULT.json records independent tests; SOURCES.json records inherited-file identity. The full producer census is deliberately not part of this audit replay. The empty-word discrepancy, missing 1-mode in all sixteen rotations, and nonzero 265-mode in the physical seed are the pertinent hostile controls.

The first-checkpoint assembled extrema proof/verifier hashes match the audited final freeze. Its packet checksum manifest and portable cylinder checker were checked again; no portability issue was found.

Reviewed source freeze:
- RESULT.md: 7c961e810ee46c25bd03faf75a41f9a1447f04392875ad4c1b5dce906b6b2894
- verify.py: e77a8675f557351539b8cd102af6334769571838c9e237072f1d3f84b45762ca
- screen.py: c1d30fc6bbcb1c7bcd72ca2b58a73f8b6046f84d2787e18b7bb32e853c5b7d69
