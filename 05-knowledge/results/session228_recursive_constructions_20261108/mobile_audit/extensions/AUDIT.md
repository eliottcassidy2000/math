# Independent audit: strong and exponentially large physical fibers

Verdict: **PASS**, no corrections required, on frozen `inputs/STRONG.md` and `inputs/MULTI.md`. The initial audit packet is unchanged; these are separately frozen addenda. Exact source and dependency hashes are in this directory's SHA256SUMS.

## Strong single-triangle theorem

For r≥3, 1≤p≤2r−4, deleting the original leading0 and trailing1 removes an actual global source and sink, respectively, since neither is replaced. Forced Hamilton endpoints preserve D, and the same argument preserves d in the deletion graph. A global source cannot occur in the root-starting member of an ordered two-path cover, so its removal preserves U. A global sink can terminate either path, including being the sole vertex of the other path; these two disjoint possibilities each biject to covers after removal. Thus its removal halves U. The root itself is neither removed vertex, avoiding degenerate exceptions.

The remaining last deletion vertex is the unreplaced0 at original position2r−2; all deletion Hamilton paths end there and can continue to the root, giving E=d=3. The first deletion block is either singleton1 or, when p=1, the triangle of neighborhood size2. In both cases the root reaches that first strong block, and hence all later blocks. Every vertex reaches the final singleton0 and then the root. This proves strong connectivity, including the first-triangle boundary.

Consequently the claimed state `(B_r+(-2)^p−2,3,2^(2r−1)+1,3)` and order2r are correct. Powers give pairwise distinct D. The r3 witnesses23 and29 both have sidecar(3,33,3). This is a genuine strong-fiber strengthening, not merely reusing reducible witnesses.

## Separated multiple triangles

Let A=S1*S0 and F=T1−3A, using independently derived matrices from the earlier audit. Exact multiplication independently confirms

    A(A−4I)(A−I)^2=0,
    F A^h F=0 for h=1,2,3,4,
    F^2 != 0.

Multiplying the degree4 identity on the left and right gives the same forward recurrence for F A^h F. Four consecutive zeros at h1,...,4 prove vanishing for all h≥1. No vanishing at h0 is claimed or needed.

In the expansion of the r pair factors, each replaced pair contributes3A+F. A term with at least two F factors has two successive F factors separated by at least one unreplaced A, since the allowed replacement indices differ by at least2. Any intervening replacements selected as their3A term add more A factors and scalar factors, not a new interaction. Thus all terms with two or more F vanish as matrices. The remaining zero- and one-defect terms give exactly

    D = 3^t H_r + 3^(t−1) sum_(j in P)(4^j−2),
    H_r=(4^r+6r−1)/9.

The base sequence starts1,3,9 at r1,2,3. The cubic recurrence for r≥1 follows from the matrix annihilator and proves the displayed H_r at all orders. Its integrality is therefore guaranteed; independently, 4^r modulo9 cycles4,7,1 and the numerator is divisible by9 in all three residue classes. The theorem restricts t≥1, so every displayed factor3^(t−1) is integral. The verifier includes t0 only as a separate base control using H_r, not by integer truncation of3^−1. The r0 empty-word boundary is not part of the H_r claim.

The multiple-replacement U identity supplies the stated common sidecar. The first and last pairs are unmodified, so source/sink trimming works exactly as above. The trimmed first and last deletion blocks are singleton1 and singleton0; the intervening blocks are strong. This proves strongness without any uniform root orientation inside a triangle. The orders are2r+t+1 before trimming and2r+t−1 after trimming, and E becomes d=3^t.

For fixed t, subtracting2t from the digit sum is constant. Uniqueness of base4 digits0/1 proves distinct D for distinct subsets. There are N=floor((r−1)/2) eligible indices. The central binomial coefficient is at least2^N/(N+1), while trimmed order is O(N); hence exponential growth with graph order is literal and asymptotic, not an exponential in a smaller unrelated parameter. For example, n≤5N+3 along the stated central-subset sequence is already enough to establish an explicit exponential lower bound up to a polynomial factor.

## Hostile interaction and scope

On effective word0101, replacing both pairs yields D41. The two single-replacement values are9 and13, while the unmodified scalar is3. The naive additive expression is39. This reproduces the claimed hostile and explains precisely why the separation hypothesis is necessary. F² is the missing interaction, rather than a failure of the U sidecar identity.

The construction produces exponentially many different scalar Hamilton counts at one full sidecar for each fixed(r,t), even among strong tournaments. It does not produce infinitely many counts at one fixed numerical sidecar, and does not provide an interval, a prescribed target, a recursive inverse, or Chan–Pak completeness. The growing deletion level3^t is explicit; these exponential lower bounds should not be confused with the earlier logarithmic classification at fixed deletion level3.

## Independent receipts

The standard-library `check.py` imports only `../verify.py`, the earlier independent audit module; no producer program is imported. It reconstructs multi-triangle adjacency matrices, checks literal profiles by subset dynamic programming, and verifies strongness by reachability from every vertex. All conditions use explicit runtime checks active under optimization.

Commands from this directory:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

Both executions PASS with byte-identical output:

- 20 literal single-triangle strong profiles, r3–6 and every allowed p;
- 8 literal multi-triangle and trimmed-strong profile pairs, r3–5, all eligible subsets including the separately guarded t0 controls;
- 3,068 scalar transfer cases, r3–21 and every subset of eligible indices;
- uniqueness for every fixed(r,t) group in that range, largest tested fiber252;
- exact annihilator, initial base recurrence values, all required defect-product identities, and the adjacent41-versus39 hostile.

These finite checks corroborate the all-order matrix/combinatorial proof; they are not used as extrapolated evidence for the theorem.
