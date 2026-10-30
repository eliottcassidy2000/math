# Independent audit: fixed-core resonant tail

**Verdict: PASS, scoped.** The proof in `../lrc/RESONANT_TAIL.md` establishes the asserted maximum, complete maximizer sets, and terminal mass slopes for every positive integer tail appended to C={1,...,11,13}. It does not prove unrestricted LRC(14).

## Proof audit

1. The AP subcore bounds clearance by1/12. Equality forces equally spaced points0,t,...,11t, so precisely the four primitive twelfth-grid points remain. Speed13 preserves them.
2. At radius1/14 the finite endpoint arrangement contains154 distinct walls including0 and1. Independent rational wall/open-cell evaluation reconstructs exactly four intervals and six singletons, matching the displayed certificate.
3. On each interval, every core speed remains in a single safe branch from both bottom endpoints to the corresponding top center. Each safety inequality is affine in(time,radius). Checking all three vertices therefore proves the entire triangular region. The active bottom endpoint owners are exactly(13,11),(5,7),(7,5),(11,13). Their walls give the reverse containment. Nesting eliminates every other candidate above1/14.
4. For a resonant tail12m, its phase at every center is zero. The inequality ||12m x||<=12m|x| gives the claimed global upper bound throughout the complete top chamber. The owner5 endpoints attain it because the phase displacement equals the proposed radius, less than1/2. Equality isolates precisely those two endpoints.
5. For m1,2 the same upper bound excludes all radii above1/14. The four nondegenerate bottom intervals cannot contain a safe tail point, while all six primitive fourteenth-grid points do. The two exceptions are necessary and correctly stated.
6. For a nonmultiple of12, multiplying any primitive twelfth-grid point leaves a nonzero multiple of1/12. All four survive, and the core excludes higher clearance. Direct active-owner calculations verify the stated derivative formula, including tied owners resolved by their maximum speed.
7. THM-1230's file explicitly leaves the ladder beyond m10 conjectural. This proof closes that specific obligation; it does not inherit the historical file's broad gap claims.

## Independent controls

Run `python independent_core.py` and `python -O independent_core.py` in this directory. Both outputs agree exactly:1077 explicit exception-based gates. The computation checks all154 core walls and every intervening open cell, all vertex/branch conditions of the four affine triangles, and exact observer/owner slopes for tails1..96 plus997,1009,1012,1200. No producer functions are imported. `receipt.json` records the result.

The parent separately verified maxima from opposite-slope candidate times for a larger tail range. Those checks are separate from this audit. The all-tail conclusion rests on the proof and finite core certificate, not extrapolation from tested tails.
