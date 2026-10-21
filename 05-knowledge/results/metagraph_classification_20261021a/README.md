# Three presentations in every order, with literal metagraph ports

Research date: 2026-10-21. **PROVED from the explicitly cited inherited inputs** below; independent proof/API audit passed. **FINITE-EXACT:** the accompanying source-bound controls. No claim of complete metagraph adjacency or black-core connectivity.

The only strongly connected tournament with three Hamiltonian-path automorphism orbits is the regular cyclic tournament on five vertices, denoted $R_5$ here. Consequently, every tournament with three presentations is an ordered SCC stack over $K_1,C_3,R_5$ containing exactly one $R_5$.

This settles the first spectrum-floor clause of [HYP-9352, presentation-count laws](../../hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md). Its stronger monotonicity and $2n-1$ assertions, the forbidden presentation counts $7,21$, and the adjacency/detour questions do not follow. The hypothesis's pendant-center question was already settled by the separate audited deficit-two packet; that is an input here, not a new result.

## 1. Inheritance, conventions, and the missing case

For a labelled tournament $T$, let $H(T)$ be its directed Hamilton-word count, $a(T)=|\operatorname{Aut}(T)|$, $h(T)=H(T)/a(T)$, and $b(T)$ its blue Hamilton-word orbit count. A word is blue when its positional reversal is an anti-automorphism. Automorphisms act freely on words; $h$ is odd, and on a self-converse tournament $b$ is odd. A non-self-converse tournament has $b=0$.

We use these current sources:

* [THM-4741, presentation counts and the corrected rotation alternatives](../../../01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md). The rotation inequality is inherited. Its corrected index-two alternative is essential.
* [The audited uniform prefix-rerouting classification](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md): the nontrivial strong pure-blue atoms are exactly $C_3,R_5$. We do **not** use the incomplete all-order argument in THM-4730(B).
* [The audited deficit-two classification](../deficit_two_centers_20261020/README.md): the unique strong self-converse tournament with $h-b=2$ is $S_4$, with $(h,b)=(5,3)$. Its all-order cutoff and finite certificates are inherited, not rerun here.
* [THM-2477, rigid strong classification, with its all-order supersession](../../../01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md): $h=1$ iff the SCC word uses only $K_1,C_3$.
* [THM-2454, center law](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md), [the exact deficit-divisor extension](../deficiency_three_20261020/centers.md), and [THM-4731, fundamental merged degree convention](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md).

The external input is **CITED** Benhocine–Wojda, *On the existence of specified cycles in a tournament*, Journal of Graph Theory 7 (1983), 469–473, [publisher abstract and DOI](https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190070413). Every tournament of order $n\ge7$ contains a directed Hamilton path whose first vertex beats its last: the $D(n,2)$ instance of their oriented-cycle theorem. The primary abstract was read for this exact statement; no new proof of that theorem is claimed. The same input was independently checked in the [rigid-closure audit](../merged_metagraph_finite_pure_atoms_20261020/AUDIT.md).

The closest corrected near miss is [MISTAKE-648](../../../01-canon/MISTAKES.md): odd automorphism order does not make the cycle-cut index odd. Index two occurs when $n=2\pmod4$. The new ingredient below handles that case rather than dropping it. MISTAKE-647 records why re-expanding the already known Narayana and pendant series is not new work.

Our concept board is: cycle-cut orbits; intrinsic score classes; closed versus open Hamilton words; ordered SCC reconstruction; converse-merged degree; literal source/word/target ports.

## 2. A new four-word obstruction for two equal rotation orbits

**Lemma.** Suppose a tournament has a marked Hamilton cycle
$0\to1\to\cdots\to2m-1\to0$, where $m\ge3$ is odd, and $i\mapsto i+2$ is an automorphism. Then it has at least four Hamilton-word automorphism orbits.

Put $A=\{0,2,\ldots,2m-2\}$ and $B=\{1,3,\ldots,2m-1\}$. The step-two rotation acts transitively on each set, so both induced tournaments are regular of degree $(m-1)/2$. Let each vertex of $A$ beat $k$ vertices of $B$. Transitivity and the total cross-edge count show that each vertex of $B$ beats $m-k$ vertices of $A$. Their full scores are therefore

$$s_A=(m-1)/2+k,\qquad s_B=(m-1)/2+m-k.$$

These differ because $m$ is odd. Thus **every** automorphism preserves the two sets, even if it does not preserve the originally marked cycle.

A regular tournament is strong: otherwise a source SCC $X$ of size $x<n$ would have average full outdegree $(x-1)/2+n-x>(n-1)/2$. By Hamiltonicity, $T[A]$ and $T[B]$ have directed Hamilton cycles. The marked alternating cycle supplies at least one cross arc each way. For an arc $a\to b$ with $a\in A,b\in B$, cut the internal $A$-cycle immediately after $a$ and the internal $B$-cycle immediately before $b$. Concatenation gives a Hamilton word of pattern $A^mB^m$. The opposite-direction arc similarly gives $B^mA^m$.

Together with the two alternating cuts, the four patterns are

$$ABAB\cdots AB,\quad BABA\cdots BA,\quad A^mB^m,\quad B^mA^m.$$

They are distinct for $m>1$ and invariant under every automorphism, proving the lemma. The inherited parity theorem upgrades the bound to at least five, but the proof below needs only four. The case $m=1$ is excluded explicitly: there is no strong tournament on two vertices.

This is a source-preserving map from Hamilton words to binary **score-class words**. It preserves a sufficient obstruction to orbit equivalence. It does not claim that equal binary words are equivalent, or that the binary word reconstructs actual endpoints or metagraph neighbours.

## 3. Complete strong classification at presentation count three

First suppose $T$ is strong, self-converse, and $h=3$. Since $b$ is odd and $0<b\le3$, either $b=3$ or $b=1$. In the first case $T$ is pure-blue, so the inherited strong-atom classification gives $R_5$ (the other atoms have $h=1$). In the second case $h-b=2$, so the inherited classification gives $S_4$, whose $h=5$ is a contradiction.

Now suppose $T$ is strong with $h=3$, without assuming self-converse. Choose a Hamilton cycle, let $R$ be its successor permutation, and let

$$q=[\langle R\rangle:\operatorname{Aut}(T)\cap\langle R\rangle].$$

The $n$ cuts of this cycle occupy exactly $q$ path orbits, so $q\le3$. The intersection has odd order: a nonidentity involution would reverse the edge of a transposed pair.

* If $q=1$, $R$ is an automorphism. The tournament is an odd circulant. Negating its cyclic coordinates is an anti-automorphism, so the preceding self-converse argument gives $R_5$.
* If $q=2$, write $n=2m$. The intersection is generated by $R^2$ and has odd order $m$. Strongness excludes $m=1$. The four-word lemma contradicts $h=3$.
* If $q=3$ and $n\ge7$, the cycle cuts already give three **closed** path orbits. The cited bypass gives an **open** path. Closedness is invariant under automorphisms, so this is a fourth orbit, again a contradiction. If $q=3$ and $n<7$, oddness of $n/3$ forces $n=3$; the only strong tournament is $C_3$, with $h=1$.

The regular cyclic $R_5$ has $(H,a,h,b)=(15,5,3,3)$, verified directly by all Hamilton words and automorphisms. This completes the iff theorem. Together with the inherited $h=1$ classification and oddness, every strong tournament other than $K_1,C_3,R_5$ has $h\ge5$.

**Small hostile.** The nonregular five-vertex bypass exception also has $H=15$, but $a=3,h=5,b=1$. All its Hamilton words close. It is not $R_5$ and is not covered by the $n\ge7$ bypass input. The checker retains its literal ten arcs and distinguishes both objects.

## 4. Exact enumeration and reconstruction

Ordered SCC decomposition is unique. Every Hamilton word traverses its SCCs in that order, and both $H$ and $a$ multiply. Since $h=3$ is prime, exactly one SCC has $h=3$ and all the others have $h=1$. The strong theorem therefore proves

$$h(T)=3\quad\Longleftrightarrow\quad
T\text{ is a stack over }\{K_1,C_3,R_5\}\text{ with exactly one }R_5.$$

Let $A(z)=1/(1-z-z^3)$. Counting **ordinary isomorphism classes**, not labelled graphs or converse pairs, gives

$$F_3(z)=z^5A(z)^2,\qquad S_3(z)=z^5A(z^2).$$

The second series counts self-converse classes: their SCC word is a palindrome, with $R_5$ in the center. All these classes are pure-blue. Every non-self-converse class has $b=0$. Thus the converse-merged count is $(F_3+S_3)/2$, whereas its non-self-converse part is $(F_3-S_3)/2$.

Writing $f_n=[z^n]F_3$, with negative indices zero, the exact recurrence is

$$f_n-2f_{n-1}+f_{n-2}-2f_{n-3}+2f_{n-4}+f_{n-6}=\mathbf1_{n=5}.$$

The first ordinary counts at $n=5,\ldots,15$ are
$1,2,3,6,11,18,30,50,81,130,208$.

The public `decode_three(source)` returns the literal ordered SCC vertex sets, with directed-cycle coordinates for each atom. `audit_three(source,receipt)` reconstructs all internal and cross edges and rejects a substituted source. It does not merely accept a multiset of component sizes. In particular, $K_1\to R_5$ and $R_5\to K_1$ are different ordinary sources with the same unordered atoms.

`port_receipts` emits three actual Hamilton words representing the complete path-orbit fibre. It retains source adjacency, the word, its ordered endpoints, the exact nonpath-complement target adjacency, and blue status. Each $C_3$ contributes one representative; the $R_5$ representatives are its three Hamilton words starting at cyclic coordinate zero. Its rotation group has order five and is its full automorphism group. Product action of the SCC automorphisms proves completeness. `audit_port` replays the literal word and target without searching.

The SCC type word is a complete unmarked reconstruction invariant **on this classified stratum**. It forgets source vertex labels, local cycle gauges, actual endpoint pairs, and the selected target. Those are the reason for retaining the marked frames and port receipt. No claim is made that presentation count alone, an unordered atom list, or a degree sequence reconstructs a tournament.

## 5. The next trivalent layer has one remaining finite center polynomial

Use $Q$ for THM-4731's **fundamental converse-merged paired-black graph**, whose degrees are $h$ on a non-self-converse pair and $(h-b)/2$ on a self-converse class. This is not the literal network with doubled paired-black edges.

Let $R_6(z)$ count strong self-converse types with $h-b=6$. The previously proved center-divisor law says that a self-converse deficit-six type is either rigid palindromic padding of one such center, or has a half of presentation count three and a center in $\{\varnothing,K_1,C_3\}$. Substituting the now proved $F_3$ gives the exact series for **all degree-three nodes** of $Q$:

$$B_3(z)=\frac{z^5A(z)^2-z^5A(z^2)}2
       +A(z^2)R_6(z)
       +z^{10}A(z^2)^2(1+z+z^3).$$

The three terms are disjoint: non-self-converse $h=3$ pairs; rigid-half strong deficit-six centers; nonrigid $h=3$ halves. A self-converse $h=3$ type has degree zero, not three.

The [effective black-defect order bound](../effective_black_defect_order_bound_20261020/README.txt) makes $R_6$ a finite polynomial of degree at most $4096$: take $d=6,s=8,r=256$ in $n\le2sr$. This is an inherited coarse bound, not a practical proposed census or a completed center classification. Thus $B_3$ is rational, with a denominator dividing

$$(1-z-z^3)^2(1-z^2-z^6)^2.$$

The displayed denominator yields an eventual recurrence of order at most eighteen; it is not asserted minimal. Unknown initial data still include the coefficients of $R_6$. Determining that finite center alphabet is the precise remaining classification obligation. Neighbour identities, edge multiplicities, bridges, and connectivity require their port data and do not follow from this node-count formula.

**An inherited enumeration hostile.** HYP-9352(4)'s proposal that every fixed-degree count is eventually quasi-polynomial is incompatible with the already proved degree-zero alphabet. Choose any word of $k$ blocks, each block either $(1,3)$ or $(3,1)$, and mirror it about an empty center. This gives $2^k$ distinct pure-blue classes of order $8k$, hence at least $2^k$ degree-zero nodes. A quasi-polynomial has polynomial growth. The repaired expectation here is a rational series with possible exponential growth, not a quasi-polynomial. No shared hypothesis file is edited by this packet.

## 6. Finite controls, reproducibility, and status limits

Run from the repository root:

```text
python -B 05-knowledge/results/metagraph_classification_20261021a/check.py --verify
python -B -O 05-knowledge/results/metagraph_classification_20261021a/check.py --verify
```

The checker uses no third-party packages. Its finite universes are:

* all $1099$ labelled tournaments on one through five vertices, with independent Hamilton-word and score-compatible automorphism enumeration;
* all $2184$ marked step-two-invariant Hamilton-cycle frames with two orbits of sizes $3,5,7$: respectively $8,128,2048$ cases. All four actual words and their distinct intrinsic score patterns are checked. This is a complete finite symmetry-normalized universe, not a random sample;
* the inherited complete ordinary-class TSV files at orders three through eight, checking their class counts and independently replaying every three-presentation survivor;
* all one-$R_5$ atom words through order fifteen, the three source-bound port receipts for each, and coefficient recurrences through order forty;
* the genuine $q=2$ frame, $S_4$, both order-five names, reversed SCC order, source replacement, swapped ports, a fabricated target, and float/bool alias controls.

The receipt stores literal example source rows, internal cycles, four-word certificates, complete marked port receipts, and count rows. The finite census only supports the analytic proof; it does not replace the cited bypass or the already audited all-order atom/center inputs.

Pre-freeze implementation controls caught two ordinary representation mistakes: the inherited census uses big-endian triangular bits, while the local exhaustive generator uses little-endian bits; and an SCC sort key must inspect a separate source list, because Python's in-place `list.sort` temporarily empties its input. Both are repaired and the reversed two-vertex SCC test is retained. No theorem was inferred from those preliminary outputs.

The independent review additionally caught an internal-cycle sidecar that accepted a twice-traversed cycle because it compared vertex sets alone. The consumer now requires each literal cycle to contain its owned vertices exactly once; the doubled-cycle hostile is retained. The four actual Hamilton words already certified the mathematical lower bound, so this was a receipt-authentication correction.

Final receipt: **21,964 always-active checks**, normal and optimized execution equal; LF SHA256 `33839183bbb4285e524818847c6f10cc581e27827f2b8be6dcb3eb28acc73d02`.

Independent family-agent audit: **PASS** on the complete proof, source decoder, literal target/port API, and optimized-versus-saved replay, including the cycle-length correction. Root's independent proof read also passed. The external bypass theorem and inherited audited all-order classifiers remain named inputs rather than newly re-proved results.
