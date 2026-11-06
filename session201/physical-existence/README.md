# Physical adjacent-cut inverse and diagonal Mersenne repair

Session 201. The main proof and verifier passed independent audit without correction. Universal {7,21} coverage remains open.

## Results

* **PROVED, independently audited in the stated scope:** the exact physical pair inverse for full C1, balanced, and synchronized joint membership needs cover heights only through floor(sqrt(3(n-1))), improving the previous floor(2 sqrt(n-1)). Every minimum J-budget statement and physical graph reconstruction is retained.
* **PROVED, independently audited in the stated scope:** coefficient 3 is sharp for the adjacent-letter inequality and for cutting a specified word. Global optimality among all presentations is not claimed.
* **PROVED, independently audited in the stated scope:** for every even a>=2, (2^(a+1)-1)^2+2 has a zero-J balanced precursor and at least four distinct numerical source profiles. The balanced graph uses one more vertex than the original diagonal recipe. No graph-isomorphism multiplicity conclusion is made.

The mechanism is an exact nonnegative 4-by-4 residual matrix for each physical letter. At the crossing of the two physical source heights, choose the better of the two adjacent cuts. Transport both boundaries together in the joint problem.

## Reproduction

    python3 verify.py --cap 5001 > VALIDATION_NORMAL.json
    python3 -O verify.py --cap 5001 > VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json

The verifier uses only Python's standard library. Both receipts match. The exact finite universe has 126,910 physical profiles through U=5000, containing all F targets through 5001. The smaller pair inventory has cutoff122, 514 physical source profiles, and 398 synchronized pairs. Full C1, balanced, joint and all minimum J costs match the forward universe. It checks 19,682 words/cut families through word length8 and literal Hamilton-path/cycle counts for the two first even diagonal repairs.

At target cap5001 the previous cutoff141 has 641 physical sources and 499 synchronized pairs. Thus the raw full pair Cartesian product falls from410,881 to264,196, and the raw joint product from319,859 to204,572. These are counts before target-dependent pruning, not a runtime or asymptotic complexity claim.

## Files and scope

* PROOF.md contains all-order statements, seed cases, cost accounting and sharpness limits.
* verify.py is the verification entry point; it retains the inherited synchronization hostiles and literal graph reconstruction checks.
* dependencies/ freezes the three inherited proofs used here.
* exploratory/ records bounded probes that prompted the final theorem. They are not proof dependencies. The affine probe uses NumPy, unlike the verification entry point.

The odd-index diagonal chamber remains open. Neither the improved inverse nor the explicit even-index family proves that every other odd target has a nonempty physical fiber.

The author verifier and its receipts preserve their frozen pre-audit status string; audit/REPORT.md provides the subsequent independent PASS. Companion corollaries are not part of this checkpoint.
