# Independent audit: zero slack and bounded-deletion sparsity

Session 195, 5 November 2026. **PASS.** This audit binds the original zero-slack packet, including its universal structural theorem and its **finite-bank** sparsity theorem. No stronger candidate refinement is incorporated into this review.

Reviewed proof SHA-256: `03681c3d613ca0284a3d00b5086dd5abda934c535c1072fad96fb1fae7c00cd1`. Exact companion source hashes appear in `reviewed_sources.json`. Suggested publication destination: `session195/zero-slack-audit/`.

## 1. Exact slack identity and structural equivalence

The identity

    f=R+e-h-d=sum_(distinct x,y in N+(r)) B_yx(T-r)

is correct for every nonempty rooted tournament, including singleton and weak cases. R uses an optional empty first path and a nonempty second path starting at r; B uses two nonempty paths. These conventions must not be interchanged.

The R-covers that concatenate are in bijection with Hamilton words, split immediately before r, and contribute h. Nonconcatenating covers with second path exactly(r) contribute d-e. Deleting r from every remaining cover gives precisely the boundary sum, with two distinct root out-neighbors serving as its marked endpoints. The correspondence is invertible. No signed cancellation or presumed nonnegativity is used.

The two-endpoint extension injection is valid: scan a Hamilton word of the other vertices until one named endpoint beats its next vertex, put the other endpoint at the end of the prefix, and put the selected endpoint at the start of the suffix. Earlier vertices beat both endpoints. The no-hit and empty-word cases are correct. Deleting the named endpoints and concatenating recovers the original word, proving injectivity.

It follows that f=0 exactly when the root has at most one out-neighbor. The theorem needs no strongness, grammar, or graph-order hypothesis. In particular, it does not rely on the false assertion that every specified arc lies in a Hamilton path.

## 2. Universal unit-deletion classification

If d=1, then T-r is transitive; the proof of this fact via an induced triangle and deterministic Hamilton-word insertion is valid. Zero slack leaves at most one root out-neighbor. If that neighbor is the i-th vertex of a chain of order n, the exact count is

    1+2^(n-i)-[i=n].

Thus the complete scalar set is exactly {1} together with {2^m+1:m>=1}. The case i=n contributes1, not2. A root with no out-neighbor and the singleton boundary also contribute1.

The strong representatives are valid: let the root beat only the first vertex of a chain of order m+1. Root plus the full chain form a directed Hamilton cycle. The stated H15 graph was independently counted and is strong; because15 is outside the displayed set, no representation of that scalar can have both d=1 and f=0, even outside the grammar.

## 3. Finite-bank sparsity proof

The O((1+log X)^2) result is valid with the stated fixed D and fixed finite gadget bank. Its proof does not infer universal bounded-deletion sparsity for arbitrary unbounded gadget banks.

For every non-D operation, h'>=max(h,R) and R'<=C0*h', with C0 depending only on the fixed bank. Therefore a reachable state with h=t<=D is a bounded state followed by a tail of neutral D operations. Its only unbounded coordinate is R=r0*2^k, with r0<=C0*D and e=0 when k>=1. This is an upper envelope; the proof correctly does not assume every envelope tuple is physical or reachable.

The final K/D suffix is also classified correctly. The coordinates b=d-e, c=h-d and f are nonnegative on physical states. K preserves f and D sends it to2b+c+2f. Thus positive slack cannot disappear after the final root switch. A D in a zero-slack suffix is possible only from(d,d,d,d), which K fixes; at most one such D can occur. The ensuing counts have the stated geometric form.

A D-free suffix K^j gives h_j=d+2^j(h0-e0) for j>=1, while j=0 contributes h0 separately. This boundary distinction is preserved in the proof.

The last root-switch deletion law d=d_G*t bounds its input count by D. Combining the bounded prefix with the final suffix gives the two displayed coefficient cases: g-d_G for j=0 and g-E for j>=1. Both are nonnegative by induced monotonicity and endpoint counting. The positive dyadic coefficient bounds j+k by O(log X), yielding the claimed O(log^2 X) envelope. The finite exceptional and single-geometric families cause no problem.

The concluding prose was explicitly scoped to a fixed finite bank during review. No other correction was required.

## 4. Infinitely many already-attained normalization failures

The actual D^k A graph macro, k>=1, sends a physical state with e=d and R=h to

    ((1+2^k)h,h,(1+2^k)h,h).

It returns to the same equalities. The singleton satisfies the premise, so k=1,2,4 yield every product3^a*5^b*17^c. Distinct exponent triples give distinct counts by unique prime factorization. Restricting each exponent to a constant multiple of log X gives Omega(log^3 X) values at most X.

Comparing that lower bound with the proved O(log^2 X) upper bound establishes infinitely many **already attained** counts with no bounded-d zero-slack representation in the specified full grammar. It does not claim those counts are absent from the grammar or from tournaments.

The graph strongness assertion is correct. After adding the sink chain, the new A root reaches the old root; the old strong graph reaches the sink chain, which returns to the new root. The first macro from a singleton is a directed cyclic extension and is strong. Thus all nontrivial product witnesses can be chosen strong.

The quantified conclusion is: for each fixed finite bank containing A and each fixed D, infinitely many such values lack that normalization **within that bank's grammar**. The proof does not exclude alternative unrestricted representations when D>1.

## 5. Independent exact controls

The standalone C++ verifier imports no producer code. It enumerates all33,867 labelled tournaments on orders1 through6 and checks:

- all202,013 rooted slack identities and root-outdegree equivalences;
- all38,615 zero-slack roots;
- all5,039 cases with d=1 and f=0, whose counts are exactly1,3,5,9,17 in this order-limited control;
- all502,170 unordered named endpoint pairs for the positive boundary bound;
- all10,650 endpoint-pair instances through order5 by literal injection, involving15,770 distinct source Hamilton words;
- complete rooted profiles through order5 by literal permutation/cut enumeration.

The explicit strong dyadic representatives are checked for m=1 through8. The H15 witness is independently counted and checked for strongness.

A separate seven-gadget physical bank has C0=134. Its complete bounded-state exploration through H<=127 retains601 states. The verifier checks134 zero-slack normal-form traces for deletion bounds1,3,5, reconstructing the last root switch, the bounded prefix, the complete D tail and every final suffix formula. This is an independent finite-instance control of the envelope proof, not a substitute for its all-order argument.

Thirty physical product words using factors3,5,17, with total graph order at most11, pass exact profile and strongness checks at every macro boundary. The125 exponent triples with exponents0 through4 also give distinct arithmetic products. The infinite density statement rests on the proof, not this sample.

All tests use explicit exception gates. Optimized and unoptimized bounds-checked runs agree byte for byte. Compiled binaries are not release artifacts.

## Reproduce

From this directory:

    g++ -std=c++17 -O2 -Wall -Wextra -pedantic verify.cpp -o /tmp/zero-slack-audit
    /tmp/zero-slack-audit > /tmp/zero-slack-receipt.txt
    cmp receipt.txt /tmp/zero-slack-receipt.txt
    g++ -std=c++17 -O0 -D_GLIBCXX_ASSERTIONS -Wall -Wextra -pedantic verify.cpp -o /tmp/zero-slack-audit-debug
    /tmp/zero-slack-audit-debug > /tmp/zero-slack-debug.txt
    cmp receipt.txt /tmp/zero-slack-debug.txt

No shared repository checkout was edited. Any stronger proposed bound requires its own separately reviewed proof and is not promoted by this audit.
