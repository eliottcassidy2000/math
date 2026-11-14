# Absolute local combs and the geometry of one-return stars

Session 327. **Proved; literal replay passed.** The local completeness statements concern the specified 99 chart or its one-return-star implementation, not the full tournament spectrum. The retiming theorem concerns the complete low-clock atlas, not new profiles at larger clocks.

## 1. A universal slope/intercept cone, with cardinality costs

In a one-return star I→U→W plus apex v, let S be the nonempty out-cut of v in I and let Y be its in-cut in W. Use the exact S325 state notation

    h=H(I), k=H(W), d=H(W+v),
    A=αβ, c=hd−kλ,
    H(output)=A Q_U(2)+cH(U).

The right block may be empty, with k=d=β=1 and Y empty. Write s=|S|≥1 and y=|Y|.

**Theorem 1.** Every such positive-slope star satisfies

    −A/(y+1) < c ≤ A/s.                                (1)

In particular −A<c≤A. For positive c, a slope-99 phase template A=99·2^r must pay

    99·2^r ≥ s c.                                      (2)

Thus every c>99 requires at least one phase, regardless of which literal templates are chosen.

**Proof.** For every x∈S, split each Hamilton path of I immediately before x. This injects its h paths into two-covers whose second slot starts at x. Different x give disjoint covers, so α≥s h. Dropping the outgoing guard at v injects paths of W+v into the β incoming-relaxed covers, hence β≥d. Also λ≤α. Finally β contains the k covers with empty first slot. For each z∈Y, splitting every path of W immediately after z gives k distinct additional covers ending at z in the first slot. Thus β≥(y+1)k. All h,k,d are positive. The two exact slack expressions are

    A−s c = β(α−s h)+s h(β−d)+s kλ ≥0,
    A+(y+1)c = α[β−(y+1)k]+(y+1)k(α−λ)+(y+1)hd >0.

These prove (1), without a graph-size cutoff or a sampling assumption. In the special upper equality c=A, necessarily S is a singleton, α=h, β=d and λ=0; these conditions also suffice. Equation (2) is a lower bound, not a claim that every allowed phase can be realized. For the published 99 templates, 67 of99 attain the cardinality lower bound on their phase; the other32 carry additional cost in that selected bank. This count is finite and does not prove their phases minimal.

## 2. The complete low-clock atlas

The following is the full list of nonempty physical pairs (H,Q2) with Q2≤304:

| Q2 | possible H |
|---:|---|
| 2,4,8,16 | 1 |
| 12,24,48,96,192 | 3 |
| 32,64,128 | 1,5 |
| 84,144,168 | 9 |
| 100,200 | 11 |
| 104,208 | 13 |
| 108,120,216,240 | 15 |
| 224 | 17 |
| 256 | 1,5,19 |
| 272,296 | 23 |
| 288 | 9,27 |
| 304 | 25 |

There are 35 pairs. This slightly extends the S326 cutoff303 by including (25,304).

Completeness is structural. The odd-cycle formula gives Q2≥2^n plus 2^(n−1) for each directed triangle. Every vertex in a strong tournament lies in a directed triangle: otherwise its out-neighbors have no edge to its in-neighbors, preventing return to it. Hence a strong order-seven letter has at least three triangles and Q2≥128+3·64=320. Order eight with one triangle already gives384; orders at least nine have baseline at least512. Thus every nontrivial SCC at budget304 has order at most six. Exhausting all one-vertex extensions through order six, with canonicalization only within degree classes, gives class counts1,1,2,4,12,56. Retain the strong letters, include the singleton, and take their multiplicative (H,Q2) closure through budget304. This generates exactly the table, independent of the number and placement of singleton SCCs. The replay performs this complete finite calculation from scratch; no larger-tournament census is assumed.

## 3. Exact shape around 30,000

Let J be the odd integers in [29327,30119], a set of397 targets. First consider the **relaxed** chart

    N=99Q+cH,  c∈{1,3,…,197},

with any nonempty physical source (H,Q), giving all digits phase zero even when a star cannot do so.

**Theorem 2 (exact local comb).** Its targets in J are exactly

    N≡9 (mod18), or N≡25,41 (mod46).                   (3)

There are75 of them. The periodic comb underlying (3) has39 of207 odd residue classes in each period414. No three consecutive odd integers lie in this comb, whereas pairs do occur; this is an absolute shape statement, not a modular sufficient condition inferred without a clock.

**Proof.** Positivity forces Q≤floor(30119/99)=304, so the complete atlas applies. Only four pairs can enter J:

    (23,272), (9,288), (27,288), (23,296).

The lower-clock pair (19,256) stops at29087; the higher (25,304) begins at30121; the table excludes every other possibility. The first and last pairs give the two residues41 and25 modulo46. The H=9 middle pair gives9 modulo18. The H=27 middle pair is redundant on J: its digit c ranges31,…,59, so the same target is supplied by H=9 and digit3c≤177. Conversely all residues in (3) inside J give valid odd digits in the required range for those pairs. This proves equality, not just an upper bound.

Modulo414, the18-progression contributes23 classes and the two46-progressions contribute9 each; each intersects the first once and they do not intersect each other. Thus23+9+9−2=39. In a run of three consecutive odds, the18-progression can supply at most one point and the two46-progressions together at most one, since their least positive separation is16. So no such run is covered.

### A whole relaxed lane is physically phase-obstructed

Now allow **any** one-return-star templates with slope99·2^r, r≥0, and offsets c=1,3,…,197, not just the published bank. Absorb phase r into Q in the relaxed formula. The (23,272) lane needs digits105,107,…,137 inside J, all greater than99. Theorem1 forces r≥1, requiring an actual source pair (23,272/2^r). The complete atlas contains none. This lane cannot be implemented by any such star.

Consequently all targets from this entire one-return-star class in J are confined to

    N≡9 (mod18), or N≡25 (mod46),                       (4)

which contains60 of the397 targets (31 of207 odd classes per abstract period414). This is an upper bound, not a claim that all60 have such a star. Exactly15 targets of the relaxed comb are excluded by this universal phase obstruction. The two points29619 and30033 on the discarded41-mod46 lane remain eligible through the9-mod18 lane.

### Exact effect of the published phases

For the particular S325 bank, a candidate relaxed triple (H,Q,c) is executable iff (H,Q/2^r_c) is physical and nonempty. The atlas makes this an exact check. On J:

- (23,272) contributes nothing, for the universal reason above.
- (23,296) requires phase zero and contributes16 targets.
- (9,288) permits phases zero or one and contributes26 targets.
- (27,288) requires phase zero and contributes7 targets.

After six overlaps, the bank realizes exactly43 targets in J. Its largest successive gap there is44, and it has the adjacent odd pair29833,29835. The replay supplies the exact list and witness triples. These are complete chart statements; excluded points can have other strong realizations, as S326 already demonstrates at30013,30029,30031 and30047.

## 4. Common clock doubling separates rather than widens these bands

Use **every** profile in the Q≤304 atlas. Retiming all clocks by the same t≥0 gives the relaxed sets

    B_t = {99q·2^t+c h : (h,q) in the atlas, c=1,3,…,197}.  (5)

These are inherited physical profiles: append t singleton SCCs. For a fixed q, all its targets lie in the closed interval

    [99q·2^t+min H(q), 99q·2^t+197 max H(q)].             (6)

This is only a convex hull; it does not assert all intervening odds occur.

**Theorem 3 (sharp unbraiding).** For every t≥3, the intervals (6) belonging to distinct clocks q are pairwise disjoint and ordered by q. The threshold t=3 is sharp even for the actual published phases.

For the finite atlas, the exact maximum over adjacent clock pairs is

    max (197 max H(q)−min H(q'))/[99(q'−q)] = 662/99 <8,

attained at(q,q')=(288,296). Thus2^t≥8 separates every adjacent interval, and hence all of them. This finite maximum certifies the all-t statement; a bounded scan of retiming exponents is not the proof.

At t=2 there is a genuine published-bank collision:

    99·104·4+123·13 = 99·108·4+1·15 =42783.              (7)

The source (13,104) has literal rows(2,12,9,17,7), and(15,108) has rows(16,5,9,3,14). Digit123 has phase two, using the first source without pads; digit1 has phase zero, using the second source with two pads. The resulting strong cores have orders15 and14, so are nonisomorphic. The replay directly counts both cores at42783. This shows the threshold is sharp for actual realizers, not just overlapping convex hulls.

### An explicit growing omitted interval with realized endpoints

At the adjacent clocks104 and108, for every t≥3 the entire odd interval strictly between

    L_t=10296·2^t+2561   and   R_t=10692·2^t+15          (8)

is absent from B_t. These are respectively the last possible point at clock104 (H=13,c=197) and the first at clock108(H=15,c=1); Theorem3 excludes every other clock. The number of omitted odd targets is exactly

    (R_t−L_t)/2−1 =198·2^t−1274.

Both endpoints are strongly realized by the published bank: digit197 has phase one and digit1 has phase zero, so t≥3 pays both costs. At t=3 the endpoints are84929 and85551, with310 omitted odd targets between them. These are gaps in the retimed low-clock family, not gaps in the full spectrum.

This is the mechanism behind the distinction between an accidental near-30,000 cluster and growing coverage. Fixed-H digit widths stay constant while differences between their clocks double. To fill the growing gaps one must introduce new profile clocks, change the cut, or use a different construction; retiming this complete finite atlas does not do it.

## Reproduction and scope

The replay independently regenerates the small complete profile atlas, verifies all inequalities and exact comb memberships, checks the published-phase subset, regenerates the universal-star slack certificate on literal templates, and directly counts the two strong42783 witnesses. The all-order parts are the split-path inequalities, the algebraic phase obstruction, and the monotone interval separation proved above. No claim of global tournament nonexistence is made. Chan–Pak remains open.
