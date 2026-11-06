# A complete all-depth target-67 hostile for triangle-root amalgamation

Status: **FINITE-EXACT, independently audited, with all-depth completeness proved below.** This is a failure of one physical grammar, not a failure of tournament existence.

The alphabet consists of the O/Z singleton operations and the four rooted triangle cores j=0,1,2,3 from `CONTENT.md`. The root is retained throughout. Every possible incidence pattern of the root on a triangle is included up to an independent cyclic relabelling of that triangle. Start at the rooted singleton. The total number of letters is unrestricted.

The common-content theorem with c=3 shows that a target not divisible by three can use at most two triangle cores, irrespective of the number of binary singleton letters. Retain the exact state (h,x,y,d,m), where m is the number of cores so far. Use the explicit six physical transitions in `verify_triangle_67.py`, with core letters permitted only when m<2.

Here is a complete target inverse for odd outputs coprime to three through N:

1. Start at (1,0,0,1,0).
2. Record every child value h<=N.
3. Expand and merge a child exact state only when h+x<=N.

Every transition is nonnegative and never lowers h. O preserves h and increases x. The core j=3 multiplies h by three and increases x. Every other letter that can change h has h'>=h+x: this is immediate for Z, and for j=0,1,2 the h row has coefficient at least three on both h and x. Thus once h+x>N, no future use of Z or a j<=2 core can lead to a value <=N. A suffix consisting solely of O and j=3 either leaves the recorded h unchanged or makes it divisible by three. It cannot yield a new value in the requested coprime-to-three target class.

The cutoff therefore loses no target in that class. Every expanded state has 1<=d<=h<=N, 0<=y<=d, 0<=x<=N-h and 0<=m<=2, so only finitely many states occur. Exact equal-state merging preserves every future scalar continuation by the physical transfer. This proves completeness at arbitrary total word length, not merely a sampled depth.

The exhaustive inverse with N=67 finds precisely the coprime-to-three odd omissions

    7,55,67.

The first additional prime omission is 67. Therefore no number of triangle core operations, and no placement or number of O/Z singleton spacers, realizes H=67 in this grammar. The general target 67 is known to be realizable by other tournaments; nothing here contradicts that.

The verifier constructs actual adjacency matrices for four positive controls, H=13,41,53,65, and counts their paths by subset DP. The normal and optimized Python outputs must agree:

    python3 verify_triangle_67.py > TRIANGLE67-normal.json
    python3 -O verify_triangle_67.py > TRIANGLE67-optimized.json
    cmp TRIANGLE67-normal.json TRIANGLE67-optimized.json

This provides a concrete stopping reason for the triangle-deletion route: more triangle cores cannot repair the obstruction. A genuinely different deletion type, a different physical operation, or rebuilding outside this retained-root grammar is required. A strong four-vertex deletion has (H,P2)=(5,32) and gcd one, so the common-content obstruction does not apply to it; that observation alone does not prove that it supplies all missing targets.
