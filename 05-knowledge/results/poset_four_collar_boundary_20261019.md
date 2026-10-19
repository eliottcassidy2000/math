# Four-collar transport, an explicit 7-thin family, and the finite-boundary frontier

**2026-10-19. Status: PROVED scoped; finite exact controls; CITED current preprints; OPEN general conjectures.**

The principal new result is an exact local-to-global certificate for a core poset attached to arbitrary contexts through two single-incomparability interfaces. It gives an explicitly proved infinite family of ordinal-indecomposable, width-three, 7-thin posets with a designated pair of balance at least `101/236`. It also identifies precisely why an arbitrary balanced core cannot replace the certificate. Neither the classical poset 1/3–2/3 conjecture nor the separate tournament conjecture is solved.

## 1. Inheritance and current literature

### The inherited mechanism

- `05-knowledge/hypotheses/HYP-9277-tournament-one-third-balance-narayana-equality.md` concerns uniform Hamiltonian paths of tournaments. Its strong-majority reduction applies to any random linear order, but its full conjecture is distinct from the poset conjecture.
- `01-canon/theorems/THM-4655-strong-tournament-balance-lemmas.md` supplies an exact cut-vertex transfer chain. Its local-certificate viewpoint is the useful inheritance here, not an identification of the two probability laws.
- `05-knowledge/results/golden_order_bridge_20261014j.md` gives the exact thin-order/Fibonacci comparison model, the marked cyclic-seam obstruction, and correctly weighted ideal flow. Its warning is essential: restrictions of uniform extensions are usually not uniform extensions of the induced core.
- `05-knowledge/results/strong_balance_pisano_golden_eleven_20261014.md` distinguishes proved local injections from failed average and stochastic-domination shortcuts. We retain that discipline.

The closest positive mechanism is source-bound prefix/suffix counting. The canonical hostile is the cyclic closing constraint on Fibonacci swap words, which preserves all common comparisons but is not a full poset-extension family. The corrected near miss is replacing a core's induced law by uniform measure. The least-used coordinate exploited below is the pair of endpoint flags.

### Verified literature boundary, as of 2026-10-19

For a finite nonchain poset, put `δ(P)=max_(x≠y) min(Pr(x<y),Pr(y<x))` and `π(P)=max_x |{y:y∥x}|`.

Aires–Chan–Pak–Panova's September 2026 preprint proves `δ(P) ≥ (5−√5)/10+ε`, with admissible `ε=3^(−3^(16·10^18))`. Its Theorem 1.4 gives an absolute K such that width greater than K implies `δ(P)>e^(−1)−10^(−100)>1/3`. Its small-range argument selects BFT triples using maximal displacement from mean-rank order. This is a preprint result, not a general 1/3 proof. The cited Kahn–Saks resolution is an *in-preparation* companion, so we do not use that announcement as an independently checked theorem. [Primary source](https://arxiv.org/html/2609.30888v1).

Aires–Kahn prove: for every K and ε>0, some L(K,ε) has the property that width below K and range above L imply `δ(P)>1/2−ε`. This follows through variance/range estimates and balancing in bounded width. [Primary source](https://arxiv.org/abs/2510.26134).

**CITED COROLLARY / deduction.** Combining the preceding two statements with ε=1/6 shows that a counterexample has universally bounded width and range. This is not a cardinality bound: arbitrarily long finite posets can have bounded width and bounded range.

Gupta's v2 preprint reports exact verification of the Gold Partition and 1/3–2/3 conjectures through fourteen elements, with 128 equality classes, all ordinal sums of singletons and `2+1`. The previous repository census is treated as inheritance, not rerun here. [Primary source](https://arxiv.org/abs/2607.23926v2).

## 2. The two-interface construction

Let A,B,C be disjoint finite posets, allowing A or C to be empty. Choose distinct ports `a,b∈B`, with a minimal and b maximal. If A is nonempty, choose a maximal `u∈A`; if C is nonempty, choose a minimal `v∈C`.

Start with the ordinal sum `A ⊕ B ⊕ C` and delete exactly the cross comparisons

- `u<a`, when A is present;
- `b<v`, when C is present.

Keep every other cross comparison. Call the resulting poset P. It is indeed transitive: a removed comparison cannot be re-created through an intermediate element because its left endpoint is maximal in the preceding part and its right endpoint is minimal in the following part.

Write

`r=e(A−u)/e(A)` if A is present, otherwise r=0;
`s=e(C−v)/e(C)` if C is present, otherwise s=0.

Thus `0≤r,s≤1`. Let `L(σ)=1{σ starts with a}` and `R(σ)=1{σ ends with b}` for `σ∈Ext(B)`.

### Theorem 1. Exact two-interface marginal law

The restriction of a uniform extension of P to B has probability proportional to

`w_(r,s)(σ) = (1+r L(σ))(1+s R(σ)).`                    (1)

In particular, for any event E determined by the B-order, define

`N=e(B), N_a=e(B−a), N_b=e(B−b), N_ab=e(B−{a,b})`,

and let `C,C_a,C_b,C_ab` count E among B-extensions with no flag, the left flag, the right flag, and both flags, respectively. Then

`Pr_P(E) = (C+r C_a+s C_b+rs C_ab)/(N+r N_a+s N_b+rs N_ab).`    (2)

The flag counts involving E must be interpreted as events in B, including when E itself uses a or b.

**Proof.** If a precedes u in an extension of P, every element of `A−u` precedes a and u precedes every element of `B−a`. Hence a,u are consecutive. Swapping them restores the A-before-B boundary. Conversely, this swap is possible precisely when the chosen A-extension ends in u and the chosen B-extension starts in a. The same statement holds for v,b at the other boundary.

The two swaps have disjoint endpoints because a≠b. Thus they commute, including when B has only two elements. Every P-extension has a unique normal form consisting of an A-extension, a B-extension, a C-extension, and a choice of the enabled boundary swaps. For fixed σ, the number of left choices is `e(A)+e(A−u)L(σ)` and the number of right choices is `e(C)+e(C−v)R(σ)`. Multiplication and cancellation of `e(A)e(C)` prove (1) and (2). Empty contexts contribute one normal-form choice. ∎

This proof is an exact extension bijection, not a heuristic replacement of a context by an independent random environment.

## 3. The four-collar theorem

The four *collar posets* are B itself, B with a singleton left context, B with a singleton right context, and B with both singleton contexts, using the construction above. They realize `(r,s)=(0,0),(1,0),(0,1),(1,1)`.

### Theorem 2. Necessary and sufficient finite certificate for a fixed comparison

For any `0≤β≤1/2` and a fixed pair x,y in B, the following are equivalent:

1. That pair is β-balanced in every two-interface context P described above.
2. That same pair is β-balanced in all four collar posets.

More generally, the exact range of its probability over all allowed contexts is the interval between the least and greatest of its four corner probabilities, in the sense that these extrema are attained. The actual set of attained intermediate values need not fill the whole interval.

**Proof.** For fixed E, write numerator `F(r,s)` and denominator `Z(r,s)>0` from (2). Both are bilinear. Each is the bilinear interpolation of its four corner values:

`F(r,s)=Σ_(i,j∈{0,1}) λ_ij(r,s) F(i,j)`,

with nonnegative coefficients `λ_00=(1−r)(1−s)`, `λ_10=r(1−s)`, `λ_01=(1−r)s`, `λ_11=rs`, summing to one. Consequently

`F(r,s)/Z(r,s)=Σ_ij [λ_ij Z(i,j)/Z(r,s)] · F(i,j)/Z(i,j)`.

This is a convex combination of the four corner probabilities. Every corner is realized by an empty or singleton context, proving both directions and sharp extrema. ∎

A fixed finite decision-tree certificate whose correctness is given by linear inequalities in extension-event counts can be transported by the same argument. In particular, a fixed adaptive two-comparison Gold strategy transfers when its four path deficits are nonnegative at every corner. This is a statement about one fixed strategy, not separate existential Gold certificates chosen independently at the corners.

## 4. A completely explicit asymmetric core

Take B on labels 0,…,8 with cover relations

`0<2<7<8`, `3<4<5<8`, `0<6`, `1<6`.

Its full predecessor masks, least significant bit = label 0, are

`(0,0,1,0,8,24,3,5,189)`.

Choose ports `a=0`, `b=6` and event `E={2 precedes 4}`. The exact flag counts are

| flag | total | event E |
|---|---:|---:|
| none | 664 | 344 |
| starts at 0 | 280 | 196 |
| ends at 6 | 160 | 80 |
| both | 70 | 49 |

### Hand derivation of every count

Remove 1 and 6. The remaining poset consists of two chains `(0,2,7)` and `(3,4,5)`, followed by 8. Its extensions are the twenty shuffles of two three-letter chains, followed by 8. Let t be the position of 0 in such a seven-element extension.

For `t=1,2,3,4`, the numbers of shuffles are `(10,6,3,1)`. Among these, the numbers with 2 before 4 are `(7,3,0,0)`.

Now insert 1 in one of eight gaps, then insert 6 after both 0 and 1. For a fixed underlying shuffle with first-chain start t, the number of choices is

`Σ_(j=0)^(t−1)(8−t)+Σ_(j=t)^7(8−j)=36−t(t+1)/2`.

These weights are `(35,33,30,26)`. Thus

`N=10·35+6·33+3·30+1·26=664`,
`C=7·35+3·33=344`.

After deleting 0, shuffle a two-element chain `(2,7)` with `(3,4,5)` and append 8. There are ten shuffles, seven having 2 before 4. The independent chain `(1,6)` has 28 possible placements, so `(N_a,C_a)=(280,196)`.

After deleting 6, the original twenty shuffles each have eight positions for the isolated 1, and ten shuffles satisfy E. Hence `(N_b,C_b)=(160,80)`. Deleting both 0 and 6 leaves ten shuffles, seven satisfying E, and seven positions for 1, giving `(N_ab,C_ab)=(70,49)`.

### Corollary 3. A sharp context-independent designated-pair bound

For every two-interface context around this B,

`53/103 ≤ Pr_P(2 precedes 4) ≤ 135/236`.

The four corner values are

`43/83`, `135/236`, `53/103`, `669/1174`.

Therefore the designated pair has balance at least

`101/236 > 1/3`.                                      (3)

Both probability extrema are attained by singleton/empty contexts. This bound concerns the designated pair; the maximum balance δ(P) can be larger.

The core has width three, height four, maximum incomparability degree seven, and no nonidentity automorphism. The selected pair has incomparable strict predecessor sets and incomparable strict successor sets, so its proof is not a twin exchange. The hand calculation proves the count certificate; independent exact enumeration checks it.

The core's cover graph is a tree, so its bare balancedness lies in the known forest-cover class. The content here is the exact arbitrary-context certificate and its repeated-block consequence, not a claim that this nine-element core was previously outside known balance results.

## 5. An infinite ordinal-indecomposable 7-thin family

Take m≥1 ordered copies of B. Put every element in an earlier copy below every element in a later copy, except delete exactly the comparison

`6 in copy i < 0 in copy i+1`, for i=1,…,m−1.

Call this P_m. Every copy has the two-interface form against its entire preceding and following contexts, not merely against adjacent copies. Hence **every copy contributes the designated pair with the uniform bound (3)**, for all m.

The incomparability graph of B is connected, and the missing cross comparisons connect these graphs by bridges. Consequently P_m is ordinal-indecomposable. Every antichain spanning different blocks has at most two elements, so its width is exactly three. The incomparability degrees inside B are

`(4,7,5,5,5,5,6,5,2)`.

The left port gains at most one extra incomparable element and the right port gains at most one. Thus `π(P_m)=7`, exactly. Also `|P_m|=9m`, and concatenating the chains `3<4<5<8` gives height at least 4m. This supplies arbitrary size and height with a proved balance certificate at the first range beyond the general 6-thin result. It does **not** settle all 7-thin posets, and no claim of an exhaustive classification of previously known subclasses is made.

For m≥2 the cover graph is no longer a forest. Labeling the second block 9,…,17, the first two blocks already have the undirected Hasse cycle `0–6–10–8–7–2–0`. Thus the core's forest-cover observation alone does not prove this repeated-family statement.

### Exact transfer law

Partition a block's extensions by the two endpoint flags. The count table `(left flag, right flag)` is

`W=[[294,90],[210,70]]`.

For E it is

`W_E=[[117,31],[147,49]]`.

A boundary contributes `1+(right flag of preceding block)(left flag of next block)`. Therefore

`M=[[504,160],[714,230]]`, `M_E=[[264,80],[411,129]]`,

and, with `e=(1,0)` and `1=(1,1)^T`,

`e(P_m)=e M^m 1`,
`#Ext(P_m; E in block j)=e M^(j−1) M_E M^(m−j) 1`.

Since `tr(M)=734` and `det(M)=1680`,

`e(P_m)=734 e(P_(m−1))−1680 e(P_(m−2))`,
`e(P_0)=1`, `e(P_1)=664`.

This transfer matrix is a context-preserving analogue of the Fibonacci matrix in the inherited thin-order model. Its eigenvalues are `367±√133009`; no golden-ratio identity is being asserted for this different family.

## 6. Why the certificate cannot be weakened to “the core is balanced”

Let B be `a<c` plus an isolated y, with ports a,c. Its three extensions have internal minority probabilities 1/3. Give it singleton contexts on both sides. The endpoint-flag totals are `(3,2,2,1)`, so the new total is eight.

Both `Pr(a<y)` and `Pr(y<c)` become 3/4. Thus **no pair entirely inside the core remains 1/3-balanced**. The full five-element poset is the Fibonacci thin order Q_5; its two new boundary pairs have minority probability 3/8.

This is genuine boundary migration. An internal-core witness need not survive even when the enlarged family is the full extension set of a poset. The four-collar hypothesis is substantive, and a general proof must permit witnesses to leave the initially chosen core.

## 7. A failed extension of the maximal-displacement selector

Consider the seven-element poset with predecessor masks

`(0,0,0,3,1,21,55)`.

It has 55 extensions. In increasing mean-rank order its labels are

`(0,1,2,4,5,3,6)`.

The mean-displacement numerators, with denominator 55, are

`(35,34,−18,−28,19,−27,−15)`.

Vertex 0 uniquely maximizes absolute displacement. Yet it has no balanced partner: it is comparable to labels 3,4,5,6, and its minority probabilities against 1 and 2 are respectively `18/55` and `17/55`, both below 1/3. The pair 1,2 has probability `28/55`.

So a route that promotes “choose a BFT triple touching the maximal-displacement vertex” to “that vertex itself has a 1/3-balanced partner” fails. This does not contradict the cited BFT-triple theorem: its witnessing pair may be the other two vertices.

## 8. Exact finite-window compiler for the remaining bounded-range problem

The following elementary statements are proved here as a useful compiler. No priority claim is made for the individual bandwidth or dynamic-programming observations.

Fix any reference extension `x_1,…,x_n`, and suppose `π(P)≤D`.

### Proposition 4. Bandwidth and bounded displacement

1. If `x_i∥x_j`, i<j, then `j−i≤d(x_i)+d(x_j)−1≤2D−1`.
2. In every extension f, `|f(x_i)−i|≤D`.
3. Every k-element ideal I obeys

   `{x_1,…,x_(max(0,k−D))} ⊆ I ⊆ {x_1,…,x_(min(n,k+D))}`.

4. At each rank there are at most `binom(2D,D)` such ideals.

**Proof.** For (1), every element strictly between x_i and x_j in the reference order must be incomparable with at least one endpoint; being comparable to both would give `x_i<z<x_j`. Excluding the endpoints leaves at most `d(x_i)+d(x_j)−2` positions.

For (2), moving an element past any other element reverses that pair relative to the reference extension, so the crossed element must be incomparable. The absolute displacement is at most the total number of crossed elements, at most D. Every ideal is a prefix of some extension, so (3) follows. The middle interval has at most 2D elements, of which D are selected away from the ends. Boundary padding by already-selected dummy labels to the left and unselected dummy labels to the right embeds every boundary case into the D-subsets of a 2D-window, proving (4). ∎

Consequently exact uniform extension counts and comparison probabilities have a bounded-width integer transfer realization. A state records which D positions of the active 2D-window have been used. A transition inserts an available element, shifts the window, and checks that the outgoing oldest label has already been used. It has weight one. Padding handles the first and last D ranks.

If `F_k(I)` counts prefixes and `G_k(I)` counts completions, then the contribution of an adding-x transition `I→I∪{x}` to `Pr(x<y)` is

`F_k(I) G_(k+1)(I∪{x})/e(P)` when y∉I.

This is precisely the inherited source-bound flow, now with a rankwise state bound independent of n. All pair-event numerators can be accumulated while traversing the transitions, keeping the at most D incomparable partners of the added element. A straightforward implementation uses `O(n D^2 binom(2D,D))` exact-integer arithmetic operations and stores `O(n binom(2D,D))` prefix/suffix integers. This is an arithmetic-operation statement, not a bit-complexity claim.

Combined with the cited universal range bound for counterexamples, this makes the unsolved general problem an unbounded-length problem over a fixed finite interface alphabet. **Finite state does not itself imply a bounded counterexample size or a finite proof.** The integer prefix/suffix weights are unbounded, and bulk stationary limits can discard the balancing boundary, as the Fibonacci example demonstrates.

## 9. What remains and how the themes actually connect

A concrete unresolved question is the **one-sided common-corner witness property**: for every finite nonchain A and maximal u, does some single pair remain balanced under all weights `1+s·1{u last}`, `0≤s≤1`? Equivalently, does one pair work both at s=0 and s=1? The unrestricted property includes the original conjecture at s=0 and is not being sold as an easier solved reduction.

The next genuinely new proof obligation for bounded range is to combine robust-core certificates with **witness migration** through non-certifying interfaces. A possible strategy is to prove a decreasing boundary potential or classify allowed transfer products that evade every local comparison certificate. An arbitrary stationary transfer cycle cannot simply be declared a counterexample: it must be realizable by a full finite poset extension family with its true boundary vectors.

The cross-domain map is explicit:

- Tournament cut-vertex decomposition and this poset decomposition both produce nonnegative transfer kernels with endpoint states. They preserve path/extension counts only in their own carriers. Their states and weights are different.
- The Fibonacci/golden bridge supplies the exact warning about seam closure and the infinite-volume limit. The present theorem supplies a finite set of genuine boundary tests, rather than a match of constants.
- Square-sum/graceful and signed-Collatz work can reuse the methodological idea “retain the admissibility/endpoint state before composing transitions.” No map carrying their global existence or termination predicates to the poset balance predicate is proved here.

## 10. Exact reproducibility and limits

From the repository root run `python3 05-knowledge/experiments/poset_four_collar_20261019.py`. The JSON records:

- the analytic core count table, checked independently by explicit extension enumeration and subset recurrence;
- actual enumeration of all four collar posets;
- 81 nontrivial context/endpoint combinations, checking the exact marginal formula;
- bridge chains of 1 through 12 cores (9 through 108 elements), comparing the 2×2 transfer count with a sparse ideal recurrence and the bounded-window recurrence;
- independent designated-pair counts for every block through three cores;
- the boundary-migration control and seven-element displacement-selector hostile.

All decisions and bounds use integers or rational numbers. The all-m theorem rests on the bijection and four-corner argument, not on the finite chain checks. The control script performs no repository edits.
