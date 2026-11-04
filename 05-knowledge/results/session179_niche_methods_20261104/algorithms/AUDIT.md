# Independent audit of the algorithmic lane

Verdict: **SOUND, with one software-interface clarification recommended.** No mathematical correction is required.

## Proof review

1. **Semilinear uniform gap:** A period with positive h eventually leaves the strict contracting halfspace; a nonzero integral period with zero h would imply equality of a power of 3 with a power of 2. Thus all nonzero periods and bases have strictly negative h. A finite minimum of their length-normalized gaps gives the claimed uniform delta, including prefactor 1. Empty language and zero periods are correctly handled. The near-critical family has unbounded length and h between -log 2 and 0, so the non-semilinearity and non-context-free conclusions follow.
2. **Conservative grammar:** A signed unary-stack counter can implement the bounded increments +2 and -3 using finitely many internal states. The independently implemented counter agrees with `2a<=3z`. Its uniform contraction bound is exact. If Collatz holds, appending repetitions of the cycle at 1 reduces `2a-3z` by one per repetition and can end at 1. Conversely the stated strict smaller-source coverage gives strong induction. The equivalence is valid and does not solve the coverage premise.
3. **Affine compiler:** Composition order and carry formula are correct. Final divisibility inductively forces prefix integrality and the corresponding parity. Exact slope/threshold distinction, positive intermediate states, the least positive source formula, and the multiplier-5 threshold hostile are all correct.
4. **Quantifiers and resources:** The note does not infer a complexity lower bound from non-context-freeness; it distinguishes all sharp contracting words from a sufficient sublanguage. It acknowledges that a compact straight-line expression can require exponentially large expanded integers. No universal Collatz conclusion is claimed.

## Independent numerical controls

The original `check_certificates.py` was run normally and under `python -O`; both outputs byte-match the supplied `check_certificates.json`. Copies are `audit.normal.json` and `audit.optimized.json`.

An independent checker, `audit_extra.py`, imports none of the author's functions. It uses homogeneous matrix multiplication to derive branch triples; checks every compatible and incompatible residue for every word through length 9; independently walks each contracting arithmetic progression to its first descent witness; and implements the explicit signed pushdown-counter updates through word length 12.

Result: **365,890 checks passed**, including **348,502 incompatible-residue controls**. Normal and optimized outputs byte-match: `audit_extra.json`, `audit_extra.optimized.json`.

## Sources

The primary paper abstract at https://arxiv.org/abs/1006.3825 explicitly states the required Parikh semilinearity theorem and confirms the listed authors/title. The original ACM DOI returned HTTP 403 in this independent check, so no claim is made to have opened its full text. The cited older repository files exist at the inspected snapshot; their headings/status summaries support the scoped inspiration map. Relevant affine-section, finite-lift, valuation-potential, and Mahler-contraction passages were inspected in more detail. This audit does not independently reprove all results in those source files, and no new theorem depends on transferring their speculative claims.

## Recommended implementation clarification

The helper `certified(t,n)` is sound when t is the output of the trusted branch-word compiler. It does not authenticate arbitrary externally supplied triples. For example, a forged triple `(1,2,1)` is accepted at n=3 even though it is not the triple of either one-letter Collatz branch. This does not affect the mathematical theorem or existing checks, which all use compiled triples. State the trusted-compiler precondition, or provide a wrapper that validates the word and recomputes its triple if presenting the code as an untrusted-input certificate checker. The author was notified.

No installations were made. The forbidden MISTAKES file was neither read nor scanned.

**Clarification resolved:** The author added the trusted-compiler precondition to RESEARCH.md and the `certified` docstring. Untrusted submissions must include the word or straight-line expression for recomputation. No executable behavior changed. Final verdict: SOUND with the interface scope now explicit.
