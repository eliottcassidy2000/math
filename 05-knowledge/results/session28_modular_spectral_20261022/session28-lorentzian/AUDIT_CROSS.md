# Independent cross-audit: prime-density moment obstruction and Hessian bridge

Status: PASS for the parent prime-density theorem and the peer Hessian theorem, with two minor peer-report scope qualifications requested below. Audited2026-10-22. No parent/peer file was edited by this lane.

## Parent artifacts

Audited in full:
- session28-root/PRIME_DENSITY_SPECTRAL_OBSTRUCTION.md
- session28-root/check_prime_density_square.py
- session28-root/check_prime_density_square.json

`check_cross_audit.py/json` records exact input SHA256 values. The audit independently multiplies the five factors of q, squares it, recomputes all decimal prime counts through10^7 by a fresh bytearray sieve, and incorporates the four larger counts from the primary source table. Granville, *A good new millenium for the primes*, Table1 on printedpage2, was independently opened via the web tool and its10^8,10^9,10^10,10^11 rows checked:

    5761455,50847534,455052511,4118054813.

Source: https://dms.umontreal.ca/~andrew/PDF/PrimesICM1.pdf . This is a source verification beyond10^7, not a claim of a local sieve through10^11.

The independently obtained polynomial coefficients are(0,−1,3,−1,−3,2); evaluation gives exactly−2447633/22500000000; the uniform-measure control gives1/396. The density index is indeed m_n=delta_(n+1), so square degree10 consumes delta_11, not delta_10. The report correctly enumerates which nonzero indices are used.

### Quantifier audit

- A positive measure on any subset of the real line would make every polynomial square nonnegative. Therefore this hostile excludes the stated full scalar Hamburger-moment realization, not merely a nonnegative-spectrum realization.
- Self-adjointness plus the stated vector power domains makes the spectral-moment implication valid. The negative square, not an assumed bounded spectrum, drives the contradiction.
- The report does not extend the result to arbitrary transformed observables, off-diagonal matrix elements, signed measures, or every late shifted tail. Those exclusions are correct and important.
- The resolvent claim is explicitly restricted to a genuine moment-generating germ with these moments; it does not silently identify a branch singularity with an eigenvalue.
- The PNT block asymptotic is sound: subtracting the two fixed-ratio scales leaves a nonzero leading difference. Abelian comparison gives the logarithmic divergence. This is a separate classical asymptotic input, not needed for the exact negative-square conclusion.
- The rational no-go follows from the local meromorphic behavior at1. The positive-weight level-one holomorphic modular no-go follows because the zero constant term forces cusp behavior at infinity, while the S-transformation forces exponentially small behavior (up to a power) at tau=iy→0. The actual divergence contradicts it. This is properly scoped to the native series and modular class named, not all modular transforms/levels/vector-valued models.

### Logarithmic-integral model

The displayed positive model density is correct:

    delta^li_(b,n+2)=∫_0^1 x^n (b x−1)/[(b−1)log(bx)] dx.

Indeed inserting the elementary reciprocal integral and then integrating u from1 to b gives exactly that weight. Its numerator and denominator have the same sign; its removable value at x=1/b is1/(b−1). At x=0 the limiting density is0; at x=1 it is1/log b. It is thus a finite positive measure on[0,1]. The model uses k≥2 so no logarithmic-integral singularity at1 is crossed. The actual sequence's nonpositivity and this smooth model's positivity are compatible.

No correction to the parent theorem or script is required by this audit.

## Peer Weil REPORT.md

The signed-incidence Hessian proof, its coefficient c_n=(n−2)!/2^(n−2), all eigenvalue multiplicities, the score/circulation interpretation, and its quadratic Taylor term agree with this lane's independent derivation and permutation-based checks. It is consistent with GEOMETRY.md: the latter emphasizes that the score-norm extension is not the unique multilinear independent-edge expectation. These are complementary statements, not a contradiction.

The report's n≥2 Hessian statement is valid at n2 (the Hessian is zero). Its later nondegenerate-saddle statement correctly assumes n≥3.

Two minor scope clarifications were sent directly to the peer worker:
1. The unsigned two-subset module formula 1⊕V⊕S^(n−2,2) assumes n≥4. At n3 it is1⊕V, and at n2 it is1. This does not affect the signed-module/Hessian theorem.
2. The eigenvalue tail comparison lemma should explicitly say finite Hermitian matrices, or state discrete-spectrum/min-max hypotheses. The exact Schur-complement criterion with C>0 is automatic in finite dimension; a general operator version needs bounded invertibility/domain hypotheses. The surrounding unbounded-form warning already points in the right direction.

The eventual-polynomial moment obstruction, linear63 example, distinction between vector-resolvent and determinant poles, Cayley sign correction, positive comparator, and positive-diagonal/not-positive-block distinction are algebraically consistent on the finite-dimensional scopes stated. This cross-audit does not certify every external analytic source or any proposed full Weil positivity theorem. Those are outside its assigned evidence boundary.

## Final packaging and peer revision

The two minor peer scope requests above were applied and directly rechecked before publication: n≥4 is now explicit for the three-term unsigned Specht decomposition; the tail lemma says finite Hermitian matrices; the Schur condition requires positive boundedly invertible C. Original full-audit peer hash was8df1ca404080386cdcf38ccaec39f2b3ab6726acf48e1e0b2414046b25322881; the corrected peer report hash is e565c74e672f37ad02a962193018da81fc62ca315d51b902912aeaf953aee4e0. The original four source hashes are preserved in ORIGINAL_AUDIT_INPUTS.json; current replay hashes are in check_cross_audit.json. This explicitly distinguishes the original audited snapshot from the minor corrected revision.

The portable default reproduction command is `python3 session28-lorentzian/check_cross_audit.py`, run anywhere: sources are located relative to the script's packet parent. For separate source directories use `--root-dir ROOT --weil-report REPORT`. There are no absolute default paths. A temporary copy laid out exactly as the publication packet was tested and its default invocation produced byte-identical JSON to the explicit-source invocation.
