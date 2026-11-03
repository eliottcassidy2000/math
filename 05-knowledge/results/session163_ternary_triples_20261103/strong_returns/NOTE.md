# Two strong source returns after ternary offsets

**Scope.** These are literal strongly connected tournaments with marked ports and all-order common suffixes. They realize the branch patterns `1,0` and `2,0` from an infinite base family, with exact source-state equality. They do **not** provide the other branches at the returned sources, a closed state bank, or a proof of the tournament spectrum. The state-transfer lemma and q2 atom are inherited from the [S156 physical return analysis](https://github.com/eliottcassidy2000/math/blob/3f58f67ceef35168766783b0fa6a837a7153aeb7/05-knowledge/results/session156_spectrum_20261103/returns/NOTE.md); no literature-priority claim is made.

For a tournament X with port p, let X^a add a new vertex beating p and losing to all other old vertices, and mark the new vertex. Define its observer

    v(X,p) = (h,x,y),
    h=H(X),  H(X^a)=h+2x,  H(X^aa)=h+2(x+y).

Let t=(2,-1,0). The following graphs are specified by adjacency row masks, bit j of row i meaning i→j, with port the final vertex:

    U  (7 vertices): (112,97,115,103,74,16,32)
    U2 (8 vertices): (252,249,194,228,204,148,32,64)
    U4 (9 vertices): (496,497,499,487,456,400,288,64,128)
    R2 (10 vertices): (1008,1017,1011,997,776,912,816,592,128,256)
    R4 (10 vertices): (996,993,994,999,975,912,800,576,128,256).

Literal Hamilton-path subset enumeration, including the a and aa queries, gives

    U  = (31,14,25),  U2 = (95,41,75), U4 = (97,40,75),
    R2 = (285,123,225), R4 = (291,120,225),
    U2 = 3U+t, U4 = 3U+2t, R2=3U2, R4=3U4.

Each graph is strong. The first three graphs and their offset identities are inherited; the two R graphs give exact marked-source returns for the 0 branch after each offset. Their inherited atom words in processing order are R2=`q1 q1 a a` and R4=`r32 a a a a a`.

## Infinite mixed-suffix strong families

Append either block `q2=(8,1,3,6)` or `b=(0,1,3)` (adjacency row masks) to a marked source X as follows: all old vertices beat all new vertices, except the first new vertex beats the old port; the last new vertex is the new port. In each block every vertex reaches the first vertex, hence either operation preserves strong connectivity: every new vertex reaches the first new vertex and thereby the old port, while old vertices reach every new vertex. The inherited source-transfer lemma gives

    v(X q2)=v(X) P,             v(X b)=v(X) B,
    P=[[5,3,4],[8,7,8],[6,6,7]],  B=[[1,2,2],[4,3,4],[4,2,3]],
    tP=t,                       tB=−t.

Let W be any finite word in b and q2 with an even number of b's, and put N_W=H(U W). Then all five graphs U W, U2 W, U4 W, R2 W, R4 W are strong and retain their physical newest ports. The exact three-coordinate identities above transport to

    H(U2 W) = 3N_W+2,      H(R2 W) = 3(3N_W+2),
    H(U4 W) = 3N_W+4,      H(R4 W) = 3(3N_W+4).

In fact `v(R2 W)=3 v(U2 W)` and `v(R4 W)=3 v(U4 W)`, which is stronger than equality of first coordinates. Infinitely many distinct counts occur already for W=q2^k: x,y≥0 by path-insertion monotonicity, so the transfer has h′=5h+8x+6y≥5h. For the empty W the actual counts are 31→95→285 and 31→97→291. At W=q2 they are 417→1253→3759 and 417→1255→3765. The equalities compare explicitly constructed literal source graphs, not an abstract unmarked count rule. These are fixed finite seed replacements followed by replay of the retained common suffix, not operations on an arbitrary source graph. The seed replacements increase the respective vertex counts by two and one. For odd-b words the same physical returns hold, but the two offset identities become 3N_W−2 and 3N_W−4 instead.

## Universal nonstrong return

For any nonempty marked source X, prepend a directed triangle C3 uniformly, producing C3⊕X and retaining the old marked port. All successive a-query vertices are beaten by the triangle, so their graphs are C3⊕X^a and C3⊕X^aa. Consequently the entire state scales by 3, and this 0-branch operation repeats indefinitely. But C3⊕X is not strong. This elementary observation does not replace the strong examples above or yield the missing offset branches.

## Checks and limitations

Run `python3 strong_return_check.py` and `python3 -O strong_return_check.py`. It independently checks strong connectivity and calculates all three observer coordinates by endpoint/subset dynamic programming on each of the five explicit graphs, and on each b and q2 continuation. It checks the seed relations, both transfers on those examples, and the stated gains on t. The all-word identities rely on the inherited general transfer lemma, not on the finite controls. No source closure is asserted for additional 1 or 2 branches from R2, R4, U2, or U4. In particular, an observed count is not silently treated as a reusable source state.
