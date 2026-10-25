# Independent audit: finite-field six-factor threshold

Status: VERIFIED algebraic proof, with CITED Weil input; FINITE-EXACT hostile controls. Reviewer did not import producer arithmetic. Target: `session55-finite-field-geometry/THEOREM.md`, read 2026-10-25.

## Verdict

The uniform odd-field theorem with q >= 14400 is correct. The earlier q >= 14641 claim is also safe but weaker. The complete classification is now verified by the independent census replay documented in the closure below. The characteristic-two obstruction is independent of this threshold.

## Independent derivation

For v outside {0,1,-1,2,-2}, the six affine functions 1+u, 1-u, 1+u+v, 1+u-v, 1-u+v, 1-u-v have distinct roots. Conversely, every member of that excluded set causes a collision among their roots. This statement is valid also in characteristics three and five, where some excluded values coincide.

The product weight W(u)=product(1+chi(L(u)))/64 has empty-subset contribution q, six singleton contributions zero, and fifteen quadratic contributions each at least -1. For degrees 3,4,5,6 the respective total absolute bounds are 40,45,24,5 times sqrt(q). Hence sum W >= (q-15-114 sqrt(q))/64.

At a root exactly one factor vanishes, so W <= 1/2; deleting the six roots costs at most 3. The nine offsets are {su+tv:s,t in {-1,0,1}}. For nonzero v, they are pairwise distinct exactly when u avoids {0, v, -v, 2v, -2v, v/2, -v/2}. This removes at most 7 further units. Repeated exclusions only make the estimate safer. Thus the good-u count is at least (q-655-114 sqrt(q))/64. At q=14400 the numerator is 65. It is increasing for sqrt(q)>=120. This proves the claimed threshold, once v exists.

For selecting v, the two-factor weighted sum is (q-chi(-1))/4. The root loss at v=±1 is at most 1, and excluding the five values above loses at most 5 further units. The producer's deliberately loose lower bound (q-25)/4 is sound and positive for q>25. A sharper bound is unnecessary.

No division by three occurs. In characteristic three the bad projective directions collapse to four and in characteristic five to six. There is no condition degree < characteristic in the squarefree quadratic-character bound.

## Character-sum dependency, including an explicit proof reduction

The exact imported assertion is: for odd q and a squarefree nonconstant f in F_q[X] of degree d, with chi(0)=0, |sum_x chi(f(x))| <= (d-1)sqrt(q). Nonzero leading coefficients are allowed. Products of distinct affine factors meet all its hypotheses.

Here is the standard proof reduction, making clear which deep theorem is imported. For the smooth projective model C of y²=f(x), squarefreeness and odd characteristic give genus g=floor((d-1)/2). The affine point count is q+S, where S=sum chi(f(x)). If d is odd there is one rational point at infinity, so Hasse–Weil |#C(F_q)-(q+1)|<=2g sqrt(q) gives |S|<=(d-1)sqrt(q). If d is even and a is the leading coefficient, there are 1+chi(a) rational points at infinity; hence |S+chi(a)|<=(d-2)sqrt(q). The triangle inequality and sqrt(q)>=1 finish the displayed bound. For d=1,2 the same calculation has genus zero, giving the familiar exact formulas. This is a proof from Hasse–Weil, not an independent proof of Hasse–Weil itself.

A researcher-authored statement is Tony Forbes, *Applications of Weil's Theorem on Character Sums*, Theorem 1, page 2 (talk notes, 2007), https://www.theoremoftheday.org/MathsStudyGroup/ADF32A.pdf . It states the distinct-root version for arbitrary multiplicative character order, with constant multiples of powers excluded. It discusses the geometric source but does not supply a full proof of the deep Weil theorem. The producer's Babai URL could not be retrieved in this audit. Thus a primary original proof of Weil itself was NOT independently read; the dependency remains explicitly CITED, not newly proved.

## Optional center-zero corollary

For q≡1 mod4, set c=0,u=1. It is sufficient that v-1,v,v+1 are nonzero squares, since -1 is a square. The three-factor weight has sum at least(q-3-2sqrt(q))/8. Deleting its three roots costs at most 3/2; deleting the other collision ratios ±2,±1/2 costs at most 4. The good-v count is therefore at least(q-47-2sqrt(q))/8, positive for q>=65. This verifies the producer's optional corollary independently.

## Exact hostile tests

Run `python3 session55-audit-finite-field/audit.py` from the workspace root. Output: `audit-results.json`. The audit implements its own tuple-polynomial field arithmetic, independently chosen irreducible presentations, and checks x^(q-1)=1 for every nonzero field element.

Fields tested: 3,5,7,11,29,9,25,27,49,81,125,243. For every v it checks the six-root criterion, including excluded v=0,±1,±2. For every nonzero v and every u it checks the entry-collision iff. It also checks admissible-v existence for q>25. On up to two admissible slices per field, it evaluates all 63 nonempty subset sums, checks the individual Weil inequality by integer squaring, the exact singleton/quadratic statements, the expanded weight identity, and the root-weight loss. All tests PASS. These controls validate the implementation and boundary mechanisms, not the infinite theorem by themselves.

## Literature and novelty boundary

Cain (2019), https://arxiv.org/html/1908.03236v2 , Section 7, states conjectures rather than a full finite-field classification. Its strict-prime-power conjecture must be limited to odd characteristic, since its own Corollary 5.1 gives infinitely many characteristic-two obstructions.

David Lai's author-maintained source was located and read, not merely its README. Pinned source: https://github.com/jiji7879/magic_squares_fields/blob/0677087714b7dd446095db3e7020c813e6eef7fb/paper/main.tex . Commit date: 2026-09-29. The paper imports Lidl–Niederreiter Theorem 5.41, derives an eight-factor bound (q-28-741sqrt(q))/256, removes at most nine parameter values, and obtains 553736. Its finite-field classification section explicitly lists the sixteen odd Parker orders. Its complete computational certificate has not been independently rerun here. Therefore our defensible comparison is a smaller sufficient threshold and independent finite verification, not discovery of a previously unclaimed classification. No inference about the integer magic-square-of-squares problem follows.


## Closure: complete classification independently verified

Status: VERIFIED proof plus FINITE-EXACT complete remaining universe. This closure supersedes the earlier pending-census qualification. Producer census SHA-256: 37bfadf8dc7b223a82699e6c29e8207d58fab220acf8c8cce8b03e3130fa5150.

Run `python3 session55-audit-finite-field/verify_census.py` from the workspace root. This uses `field_arithmetic.py`, copied from the audit's independently written tuple-field implementation, and never imports producer functions. The script completed successfully; machine-readable output is `census-audit-results.json`.

- An independent Eratosthenes sieve and prime-power expansion confirms every one of the 1,730 odd prime powers below 14,400 appears exactly once, with correct characteristic and extension degree.
- All 1,714 positive witnesses have nine distinct entries, all squares (checked by Euler exponentiation), and all eight line sums equal, using their supplied polynomial field presentation.
- All extension moduli are independently proved irreducible: successive Frobenius powers are tested against the modulus by polynomial gcd for every degree through half the modulus degree, and the final Frobenius identity is checked. This is different from the producer's trial-division implementation. A reducible polynomial has an irreducible factor of degree at most half its degree, detected by the corresponding x^(p^k)-x gcd.
- All sixteen negative fields are searched exhaustively in both normalized planes, with no witness found: orders 3,5,7,9,11,13,17,19,23,25,27,31,43,47,67,243. Center-one enumeration applies only the necessary opposite-pair square condition before enumerating every ordered pair. The normalization proof in the producer theorem establishes completeness.
- All 71 records in `small_counts.json` agree with independent exhaustive counts. The smallest field allowing a nonzero center is F59, which has eight labeled center-one witnesses. The earlier non-Parker fields F29,F37,F41,F49,F53 have normalized center-zero counts 4,4,2,2,4 respectively, and zero center-one counts.

Combining this independently checked finite census with the audited q>=14400 theorem proves the complete odd-field classification. Adding the characteristic-two argument gives all Parker finite fields: every field of characteristic two, together with the sixteen listed odd orders. This is independent verification of a classification already explicitly claimed by Lai; it does not assert priority or settle the integer problem.
