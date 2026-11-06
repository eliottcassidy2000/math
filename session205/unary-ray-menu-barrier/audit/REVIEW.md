# Independent audit of the unary-ray rebuilding obstruction

Status: PASS in the stated restricted rewrite scope.

Reviewed `session205-rebuilding/UNARY_RAY_BARRIER.md`, `verify_unary_barrier.py`, and `RESULT.json`. The normal and optimized verifier runs both completed and their output exactly matched the author receipt. The verifier SHA256 is `52910d2089fec49ae19a0e435b93b17128ba82631cb1b83a168935708280f3f1`; the receipt SHA256 is `844f78264d91727962bc99114677af8262927d00af3d52942a1d0cdb0554eb98`.

## Mathematical checks

- The 28-state presentation counts unordered spanning path forests. Components beginning at bulk can never receive a predecessor from a later vertex, so the one-bulk-start truncation is valid for Hamilton-path readout. Forgetting a newly introduced vertex has a unique inverse; multiplicities from endpoint-equivalent choices are explicitly summed.
- The characteristic polynomials and ranks are computed with exact SymPy integer/rational arithmetic. Rank(A)=rank(A^2)=18 implies that the zero-eigenvalue generalized space is killed after one multiplication. Therefore context sequences have the claimed exponential-polynomial form for every n>=1, not merely eventually. This precise starting index is essential to the argument.
- The vector recurrence implies the scalar recurrence from n=1 onward. The initial counts 3 and 5 give H(1^n)=2^n+1. The empty word is correctly excluded: its count is 1, not 2.
- For fixed end contexts and a nonerased unary letter, the only nonzero eigenvalues are positive 1,2,3. The eventual dominant-term argument is valid and rules out infinitely many zeros of a nonzero defect. No Skolem-Mahler-Lech theorem or unproved recurrence-zero assertion is needed.
- If a defect vanishes identically for n>=1, n=1 forces a literal tournament with H=7. I checked the corrected all-order statement and proof route in `math/01-canon/theorems/THM-029-alpha1-gap.md` and its independent strong-component alternative in `THM-1370-h-spectrum-omits-7-21-all-n.md`. The proof does not use the retracted common-intersection argument. It is legitimate to evaluate at n=1 even though the proposed universal witness language only contains this ray for n>=5: the exponential-polynomial identity has already been established on the entire n>=1 domain.
- Erasing the unary letter gives a constant output and hence at most one successful exponent. Taking a finite union gives the claimed eventual failure for every branch simultaneously.
- Finite composition preserves the permitted rewrite form. Thus the finite normalization-branch corollary is valid, including count-preserving normalizations.

## Scope and limitations

This is an all-order impossibility theorem for a fixed finite menu of letter-or-erasure maps with fixed end contexts on the full eligible two-apex language. It is not an impossibility theorem for word-valued morphisms, stateful transducers, freely selected insertion positions, an unbounded normalization routine, or a smaller closed witness language. It does not rule out an eventual odd-count threshold and supplies no new scalar coverage interval. The finite 341-word literal checks are regression controls, not the all-order proof.

No mathematical changes are required. Status may be promoted to independently audited in precisely this scope.
