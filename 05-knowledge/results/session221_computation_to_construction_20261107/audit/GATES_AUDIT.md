# Session221 reduced O^m J^2 language audit

**Verdict: PASS** on the mathematical claims and normal-mode literal verification in `session221-gates/REPORT.md`. This is a restricted program-language analysis of an inherited physical macro, not a new primitive gate or absolute spectrum constructor.

## Proof review

The universal formulas follow by splitting paths according to the new vertex with sole outneighbor and by assigning appended sinks independently to the two ordered cover paths. For `m=0`, the `-e` term in the second coordinate is essential. For `m>=1`, `q=2^(m+1)` gives exactly

    (h,R) -> (h+qR, 2h+qR).

The physical order has the stated disjoint reverse-arc endpoints, consecutive crossing intervals, cutwidth at most two, and strongness. These claims were independently checked against the literal construction, not inferred from the arithmetic alone.

The relation `R_j=h_j+h_(j-1)` gives the variable continuant. At each inverse step, `(h-p)/p` belongs to the half-open interval `[q,2q)`, so its highest power-of-two scale identifies q uniquely. Integer divisibility, `0<=earlier<p`, and termination at `(0,1)` are sufficient as well as necessary. The two-adic decoder is also sound: after the first gate the prior cover count is 2 modulo 4. The first gate requires its separate case.

The mod-eight lane is fixed by the first q: 5 when q=4, otherwise 1. Later increments are divisible by eight. Also `h_k>product q_i`; compositions of the exponent sum into parts at least two give the exact Fibonacci upper count `F_(floor(log_2 X)+1)`, hence `O(X^(log_2 phi))` distinct outputs. This applies only to the positive-m reduced language, not all O/Z/J programs or all cutwidth-two tournaments.

The physical count-89 collision preserves order but not R: programs (2,2) and (3,1) give R=98 and 106, so the same next q=4 produces 481 and 513. This proves that scalar H is not a continuation state. Pair decoding does not imply scalar injectivity or all-target occupancy.

## Independent executed checks

`audit_gates.py` imports no producer code. It builds all pair orientations afresh and counts subset paths backwards, using opposite adjacency to obtain root-start counts and exact convolution for R.

Completed independently:

- 1,140 literal gate checks: all labeled tournaments on 1–4 vertices, every root, m=0,1,2,3;
- 30 literal programs over m=1,2,3 through length 3 and physical order 13;
- complete all-program monotone enumeration through H=1,000,000: 4,047 programs, 2,737 scalar values, no pair collisions, least scalar collision 89;
- both count-89 states and their 481/513 continuations counted on actual graphs.

The separate mathematical proofs, not the finite census, establish the all-input gate and infinite-language decoder/density conclusions. No external originality claim was assessed. THM-4904 and earlier O/Z/J work are inherited context, and the general macro is expressly not promoted as new.

## Verification guard and release

The initially reviewed producer verifier used Python `assert`, so optimized mode would skip its checks. The author replaced all 17 assertions with explicit `require` guards, reran normal and optimized verification, and obtained byte-identical receipts. The auditor checked the resulting guard definition, absence of assertion statements, and receipt equality. Final verifier SHA256: `d2b83637f679cee5e7bf404ea6683f8bed70de0ddc552da8fcad9fd37bbaf4b8`. The independent audit script also uses explicit exceptions and was separately rerun successfully under optimization. The robustness issue is closed.

Release `GATES_AUDIT.md`, `audit_gates.py`, `gates_receipt.json`, and `GATES_SHA256SUMS`. Run:

    python3 -B audit_gates.py

The audit script uses only the standard library. No mathematical correctness blocker remains.
