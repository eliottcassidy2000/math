# Exact equal-chain window capacity and the remaining matching gap

Date: 2026-10-19. Status: proved capacity formula and proof-method obstructions; the all-size support-`m+2` perfect-matching statement remains open here.

This is a continuation of [Chain support thresholds and rank windows](../gold_chain_support_windows_20261019/CHAIN_SUPPORT_THRESHOLDS_AND_RANK_WINDOWS.md). It does not change that frozen checkpoint or its finite certificates.

## 1. Setup and the precise open statement

Let `R_m` be the disjoint union of the chains `A1<...<Am` and `B1<...<Bm`, with `m>=2`. A binary shuffle determines its labeled extension uniquely. Write `C=binom(2m,m)`. The extensions beginning with A and B each number `C/2`.

For extensions `u,v`, their support distance is the number of positions carrying different **labeled** elements. This is not binary Hamming distance. Define `H_m(k)` to join opposite first-letter orientations when their support distance is at most `k`.

The previous checkpoint proves, for every `m`, that `H_m(m+1)` has no perfect matching, even though every vertex has a neighbor. A two-to-one Hall bottleneck supplies the proof. A global exchange of the two chains supplies a perfect matching at support `2m`. Thus the exact threshold `s_m` satisfies

    m+2 <= s_m <= 2m.

Explicit finite certificates establish `s_m=m+2` for `2<=m<=7`. The remaining proposed theorem is

    For every m>=2, H_m(m+2) has a perfect matching.          (OPEN)

Neither the formula below nor the obstacles in Section 4 prove this proposed theorem or refute it.

## 2. Closed formula for the canonical rank-window construction

Both minima have possible positions `1,...,m+1`. The canonical fixed-outside window is therefore the first `m+1` positions. Freeze the entire suffix of length `m-1`. In one fiber, any two extensions have support distance at most `m+1`, so the opposite orientations form a complete bipartite graph.

**Theorem.** The maximum matching restricted to these fixed-suffix fibers has size

    nu_window = C/2 - Delta_m,
    Delta_m = binom(m-1,floor(m/2))
              binom(m-1,floor(m/2)-1).                     (1)

The same value is the optimum fractional matching inside these fibers. It is a lower bound for the capacity of the complete graph `H_m(m+1)`, not a formula for that complete graph.

**Proof.** Suppose the prefix contains `a` A's and `m+1-a` B's, where `1<=a<=m`. There are

    binom(m-1,m-a) = binom(m-1,a-1)

possible suffix words. The orientation counts in each such fiber are

    A_a = binom(m,a-1),    B_a = binom(m,a).

Thus the exact fixed-suffix capacity is

    sum_(a=1)^m binom(m-1,a-1) min(A_a,B_a).                (2)

Its deficiency relative to the total A-first side is the sum of the positive fiber imbalances. Put `j=a-1`, `q=floor(m/2)`, and `c_j=binom(m-1,j)`, taking `c_j=0` outside `0,...,m-1`. Pascal's identity gives

    binom(m,j)-binom(m,j+1) = c_(j-1)-c_(j+1).

The difference is positive precisely when `j>(m-1)/2`; when `m` is odd, the additional term at `j=q` is zero. Consequently the deficiency telescopes:

    sum_(j=q)^(m-1) c_j [c_(j-1)-c_(j+1)]
      = c_q c_(q-1).

This is exactly (1). Complete bipartite fiber matching also proves fractional optimality. ∎

The first values provide useful scale:

| m | C/2 | fixed-window capacity | Delta_m |
|---:|---:|---:|---:|
| 2 | 3 | 2 | 1 |
| 3 | 10 | 8 | 2 |
| 4 | 35 | 26 | 9 |
| 5 | 126 | 102 | 24 |
| 6 | 462 | 362 | 100 |
| 7 | 1716 | 1416 | 300 |

For example, the complete support-four graph at `m=3` has capacity nine, whereas the fixed-window construction has capacity eight. Reusing the local rank bound does not account for exchanges between its fibers.

## 3. An explicit Gold threshold for this positive construction

Suppose an independent ordinal prefix block has `g` extensions and a prescribed first-comparison child has `h` extensions. Append `R_m` and compare its minima second. The initial and child counts are `E=gC` and `b=hC`. Keeping the prefix-block word fixed, the construction above supplies

    nu = h(C/2-Delta_m).

The Gold matching requirement `nu>=2b-E` is therefore equivalent to

    h/g <= 2C/(3C+2Delta_m).                               (3)

For `m=2,...,7`, these sufficient first-branch thresholds are respectively

    3/5, 5/8, 35/57, 42/67, 462/743, 143/227.

This is an exact all-size positive theorem for the specified construction. It does not claim necessity among other matchings or comparison choices. By the standard central-binomial asymptotic, `Delta_m/C` is asymptotic to `1/(2 sqrt(pi m))`, so these thresholds approach `2/3` from below.

For the previous checkpoint's degree-three prefix `T_t`, where `h=3^(t-2)` and `g=(3h+1)/2`, condition (3) becomes `h Delta_m<=C/2`. Increasing `t` defeats this fixed-window certificate while leaving the actual strict-Gold margin equal to `C/2`. This identifies the quantitative burden on any cross-fiber repair.

## 4. Obstacles to short proofs of the open statement

### 4.1 No proper fixed-prefix bijection can work

Fix `k<2m` and require every exchange to preserve positions `k+1,...,2m`. There is no perfect matching between the first-letter orientations under this restriction, regardless of how the first `k` labels are rearranged.

If `k<=m`, choose the suffix of the word `A^m B^m`. The corresponding prefix contains only A's, so this fiber has no B-first extension. If `m<k<2m`, write `k=m+r`, with `0<r<m`, and fix the all-B suffix `B^(m-r)`. The fiber's two orientation counts are

    binom(m+r-1,m-1),    binom(m+r-1,m).

Their difference is

    (m-r)/(m+r) binom(m+r,m) > 0.

Thus that fiber cannot be perfectly paired. In particular, when `m>=3`, a proof at support `m+2` cannot simply rearrange the first `m+2` positions of every extension. Some exchanges must move labels farther to the right while keeping the total number of changed positions small.

### 4.2 A standard adjacent compression can increase a neighborhood

At `m=3`, support five, consider

    u = ABBBAA,    u' = ABBABA.

The second word is obtained by an adjacent `BA -> AB` shift, moving an A earlier in the usual path-dominance order. Nevertheless, its B-first neighborhood has size eight, whereas the first word's neighborhood has size seven. The complete distance table is small:

| B-first word v | distance(u,v) | distance(u',v) |
|---|---:|---:|
| BAAABB | 6 | 6 |
| BAABAB | 6 | 6 |
| BAABBA | 5 | 4 |
| BABAAB | 5 | 4 |
| BABABA | 4 | 2 |
| BABBAA | 2 | 4 |
| BBAAAB | 6 | 5 |
| BBAABA | 5 | 3 |
| BBABAA | 3 | 5 |
| BBBAAA | 4 | 5 |

This rules out the assertion that this individual standard shift always decreases opposite-side neighborhoods. It does not rule out a more sophisticated compression, a coupled transformation of both sides, or a theorem about globally compressed extremal families.

Nor are lexicographic initial segments automatically neighborhood minimizers. At `m=3`, the first four A-first words in lexicographic order have all ten B-first neighbors. The four-word family

    {AAABBB, AABABB, AABBAB, ABAABB}

has only five. An isoperimetric approach needs a genuine structural theorem rather than either of these simple monotonicity assumptions.

### 4.3 Hall testing reduces to independent path families

Let `c` exchange the two chain names. Identify the B-first side with the A-first side through `c`. Because

    distance(u,c(v)) = distance(v,c(u)),

the resulting adjacency is a symmetric graph `G_m(k)` on A-first words. For `k<2m` it has no loops, since a word and its chain-name complement disagree everywhere. Its bipartite double cover is exactly `H_m(k)`. At the proposed bound `k=m+2`, the loopless formulation therefore applies for `m>=3`; the already-settled boundary `m=2,k=4` has loops.

**Lemma.** That double cover has a perfect matching if and only if every independent set `I` of `G_m(k)` satisfies `|N(I)|>=|I|`.

**Proof.** One direction is Hall's condition. For the other, if an arbitrary set `S` violates Hall, put `I=S\N(S)`. No vertex of `I` has a neighbor in `S`, so `I` is independent and `N(I)` is contained in `N(S)\S`. Hence

    |I|-|N(I)| >= |S|-|N(S)| > 0.

Thus an independent Hall obstruction exists whenever any obstruction exists. ∎

This is an exact reduction of the open problem, not its solution. For example, all nonnegative Dyck paths form an independent family for every `k<2m`: in such a word `A_i` has position at most `2i-1`, while `B_i` has position at least `2i`. Hence no labeled position in one Dyck word agrees with the chain-name complement of another, and their distance is `2m`. This family has `C/(m+1)` members. At `m=2`, its size two already exceeds half of the three-vertex graph, proving that support below four fails. At `m=3`, support five gives an independent family of size five whose neighborhood also has size five, as the small replay verifies. These examples suggest that extremal independent path families, rather than arbitrary Hall sets, may be the right objects to classify.

### 4.4 Other candidates remain candidates

Initial-run rotation is locally short but can map two source words to one target. First-return reflection is a genuine bijection, but a first excursion of length `2m` moves all `2m` labels. A global chain exchange also moves `2m` labels. None supplies the missing all-size `m+2` bound.

An unrestricted monotone coupling is also available. Remove the first position: an A-first word corresponds to an `(m-1)`-subset `S` of the remaining `2m-1` positions, while a B-first word corresponds to an `m`-subset `T`. The inclusion graph `S subset T` is `m`-regular on both sides and has a perfect matching. Under the labeled-word interpretation, all A labels move weakly right and all B labels weakly left. But if the one added position is `j` in the remaining word, every full-word label through position `j+1` changes. Restricting `j+1` to at most `m+2` fixes a proper suffix and fails by Section 4.1. Thus monotonicity alone does not provide the support bound; a successful short monotone coupling would have to redistribute paths beyond these one-element inclusion moves.

A possible next approach is to characterize Hall-extremal shuffle-path families under a suitable coupled or weighted compression, or to build an explicit bijection that redistributes the differently biased suffix fibers. Those are research directions only. The current closed conclusions remain the all-size bounds, the exact window capacity (1), the sufficient Gold threshold (3), and the finite sharp certificates through `m=7`.

## 5. Replay

Run `python verify_chain_window_capacity.py`, and independently run it with `python -O`. The verifier enumerates only equal-chain orders through `m=7` to check the fiber formula and the proper-prefix obstruction; it checks the binomial telescoping identity for a larger bounded range using exact integers. It also recomputes the entire ten-word compression table, the four-word lexicographic comparison, and the independent Dyck-family control at `m=3`.

These bounded checks validate the implementation and displayed examples. The all-size claims rest on the proofs above. No new large matching census or unproved extrapolation is used.
