# Independent audit of the unequal cyclic cross-profile criterion

**PASS. No mathematical correction found.** One bibliographic correction was requested and made by the author: the linked arXiv v1 PDF numbers the odd-cycle expansion Theorem 1.39, whereas the HTML conversion uses 1.38. The proof now specifies the PDF version. The audit did not edit the producer. The theorem is conditional on the six ordered cross-profile inequalities at every relevant index; it does not assert universal profile log-concavity or full response injectivity.

## Exact kernel and finite certificate

The old-run matrix is correct: a walk in the directed three-cycle is forced once its start and length are known, and each additional full circuit multiplies its monomial by abc. Reversing the quotient transposes the matrix because weights commute. In the two-port formula, the two summands represent the two relative orders of x,y. Only x before y allows an empty middle segment; the fixed edge x->y forbids the other empty middle segment. Each port incidence and each empty end segment is accounted for.

Multiplying old-run series creates ordered occurrences of each physical block. Applying L fills these occurrences with a single ordered spanning cover of that block, so no old label is repeated or omitted. Setting q_0=0 is essential because all three blocks are nonempty. The finite profile supports justify coefficientwise evaluation of the rational series.

The independent checker does not reuse the producer's polynomial algebra. It counts literal quotient-run words by multiplicity, for all 64 exponent vectors in {0,1,2,3}^3, and compares the cyclic contrast with the displayed S/(1-abc)^2 coefficients. This is a complete rational identity certificate, not a numerical fit: after multiplying by (1-abc)^3, both numerators have degree at most three in each variable. Since the denominator has constant term one, agreement on that coefficient box determines every possible numerator coefficient.

## Sign argument and zero supports

For k>=1, each term x_i x_j z_l is at most x_i y_j y_l by applying the assumed inequality to the ordered pair (j,l). Each of the six x_i y_j z_l terms is at most y_i y_j y_l by the ordered pair (i,l). Hence the first three positive terms are absorbed by the three first negative terms, and the other six by -6y_1y_2y_3. No division occurs, so vanishing support entries cause no gap.

For k=0, all z_i vanish, leaving the stated strictly negative term because every child has positive Hamilton count. Thus the weighted sum is strictly negative. Orders two and three have precisely the three listed q profiles, and their only possibly nonzero cross constraint is indeed at k=1; all nine profile pairings satisfy it. Equal log-concave profiles give the same condition by substitution, but separate log-concavity of unequal profiles does not. The TT2/TT3/TT5 sum +4540 is a valid sign-boundary control, not a collision.

## Physical names and full response

The mixed-seven source uses fixed blocks {0,1}, {2,3}, {4,5,6}. The audit reconstructs the simultaneous XOR of the entire vertex set and all three child sets, and checks that this equals the saved target with unchanged interiors. Both scalar counts are 159. The named query (3,12) and all 64 coarse queries are independently recounted.

Every full two-port Hamilton word contains at most four port-to-old incidences, so its Boolean polynomial has degree at most four. The independent checker evaluates all 1,471 assignments of total Hamming weight at most four using backward subset recursion, then performs its own subset Möbius inversion. It compares every nonzero coefficient and every implied zero with the saved complete contrast. Cubic and quartic terms are retained, not inferred from marginal agreement. This check uses neither the producer's forward DP nor its literal Hamilton-word polynomial routine.

## Reproduction

`audit.py` rejects producer changes during its run and binds both source hashes in `receipt.json`. It also reruns the unchanged producer normally and with optimized Python, comparing both outputs and the saved record. Independent normal/optimized receipts agree. All checks use explicit exceptions. The exact kernel coefficient bound is a proof step; the finite named seven-vertex data remain a concrete control rather than a general collision classification.
