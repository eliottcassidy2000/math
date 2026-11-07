# Prove the ordered-path half-bound

**PROVED, independently reviewed.** A forced two-colour shuffle closes the ordered-path conjecture and the general positive deletion-increment bounds in HYP-9469(1),(2),(6). Chan–Pak completeness remains open.

For a tournament X and a nonempty partition A,B with any arc A->B, write h=H(X), p=H(X[A])H(X[B]), and q for the Hamilton words listing B before A. The new injection gives

    h >= p+q,       q <= (h-1)/2.

Adding a vertex v->A, B->v therefore gives

    H(X+v) >= max(h+p,2h-p) >= (3h+1)/2.

A cut with no A->B arc preserves H. Thus **every induced extension, of any size, either preserves H or raises it to at least (3h+1)/2**. The factor 3/2 is asymptotically optimal.

A positive deletion increment 2k forces the parent count to be at most 4k-1. Combining this with SCC products, Busch's strong-order bound, and a complete strong-core census through six proves that the parent counts at increments 2,4,6,8 are exactly

- 2: {1,3} (already known in all orders);
- 4: {1,5};
- 6: {3,5,9};
- 8: {1,5,9,11}.

The last three lists were previously finite observations. The new finite cutoff proves that no later-order value is missing from the lists.

## Files and verification

- `PROOF.md`: all-order injection, growth theorem, sharpness family and exhaustion argument.
- `verify.py`: standalone standard-library verifier with explicit exception gates.
- `RESULT.json`, `RESULT.optimized.json`: equal deterministic receipts.
- `SOURCES.json`: source bytes checked against GitHub at main `4c63fd58b566f3e006eb87eaedcf48824741cc84`.

Run `python3 verify.py` and `python3 -O verify.py`, directing output outside the packet and comparing it with the receipt. Checks cover all 74,938 nonzero cross matrices with both dimensions at most four, all 33,867 labelled tournaments through order six and every cut (2,044,814 active cuts), plus literal and fresh extension controls. The independent audit uses a separate permutation-based implementation and verifies the full pair injection through order six, rather than reusing the author's dynamic program.

This is a structural advance and an induced-growth obstruction, not an all-odd construction or a new ambient forbidden path count. Generic path merging is inherited and credited to Busch via THM-4745; the new step forces a nonordered image while retaining both complete input subsequences.
