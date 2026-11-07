# One old-arc reversal escapes each one-exit obstruction

**Status: FINITE-EXACT physical constructions, pending independent audit.** This is a positive surgery control, not an all-target recurrence.

Use exactly Session213's comb B^s(m,q,tau): a transitive chain x_0,...,x_(m+1), the apex a beating only x_1 among chain vertices, b beating x_0, x_1->b iff tau=1, x_(i+2)->b iff the i-th bit of q is one, and apex arc a->b for s=+, b->a for s=-.

Two single-arc changes suffice:

1. Start with **B^-(8,14,0)**, which has H=3,187. Reverse only x_0->a to a->x_0. The resulting twelve-vertex tournament has **H=7,477**. Its outgoing bit rows, in chain order followed by a,b, are

       [1022,1020,2040,4080,4064,4032,1920,1792,1536,1024,3,1991].

2. Start with **B^+(9,179,1)**, which has H=24,693. Reverse only x_5->a to a->x_5. The resulting thirteen-vertex tournament has **H=27,371**. Its outgoing rows are

       [4094,6140,8184,8176,4064,1984,8064,7936,3584,7168,2048,4130,1329].

The standalone checker constructs each original comb from (m,q,tau,s), checks that exactly the two corresponding row bits change, and counts parent and child both by subset DP and by the three-run transfer. It also compares the original parent's count with Session213's closed quadratic comb formula, implemented literally and credited in the source manifest.

The distinction from exit relocation is exact: each move **retains** a->x_1 and adds a second chain exit. The one-exit class exclusions in `PROOF.md` and `ADDENDUM.md` therefore do not apply to the children. Since every original comb belongs to that class, zero old-arc reversals cannot realize either target from any original comb. Hence one is the exact minimum if the permitted operation is to choose any original comb and reverse old arcs on the same vertex set. This statement does not cover vertex additions, deletions, or more general rebuilding.

The smaller eleven-vertex 7,477 witness in the initial packet is a different positive control. The current twelve-vertex witness is useful because its relation to an original comb requires only one flip; no claim is made that it has minimum order.

    python3 one_flip_verify.py > ONE_FLIP_RESULT.replay.json
    python3 -O one_flip_verify.py > ONE_FLIP_RESULT.optimized.json
    cmp ONE_FLIP_RESULT.json ONE_FLIP_RESULT.replay.json
    cmp ONE_FLIP_RESULT.json ONE_FLIP_RESULT.optimized.json
