# The smallest stable refinement of the full slot observer

**Status: PROVED IN NOTES after independent proof/code audit; see [AUDIT.md](../AUDIT.md).** The finite-group and linear-algebra facts below have elementary proofs. The complete computations are FINITE-EXACT for n=3..7. No all-order class separation, metagraph generation recurrence, Hamiltonian-spectrum closure, or Chan–Pak proof is asserted.

## Inheritance, target, and live concepts

[THM-796, three-sorted recursive incidence](../../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md) already supplies actual dart fibres and the weighted node kernel. [THM-5153, exact slot pairings](../../../../01-canon/theorems/THM-5153-exact-slot-pairings-of-merged-tournament-classes-and-self-converse-structure-to-order-nine.md) already gives the counts of full-Q classes through order9; our Q counts are independently reproduced controls, not newly discovered counts. [The prior slot/wiring note](../../joint_metagraph_dyadic_20261113/SLOT_WIRING.md) proves the first full-Q closure failure at5. Its equal-Q H9 pair, codes8 and10, is the canonical hostile: equal Q, scores, H, Aut, converse status and degree still leave different actual partners.

The corrected near miss is the prior confusion between degree-sorted canonical codes and unrestricted canonical minima; all codes in this packet are explicitly unrestricted minima. The least-used relevant sidecars here are the source-dart degree and blue-dart degree. Anchor: the smallest extra information needed to close full-Q transitions. Niche: distinguish literal incidence from normalized probability. Wildcard: optimize a small set of fixed physical anchor columns.

Concept board: the full slot polynomial; automorphism/converse source mass; colored partner incidence; refinement depth; literal anchor measurements; linear versus algebraic observation closure. The order-seven delayed pair connects all six: normalization hides its different masses, Q hides different targets, color repairs the loss, and one extra refinement round restores the target distinction.

## 1. Objects and source weights

Fix n≥3. The ordinary class space consists of all tournaments up to relabeling; R is converse. The merged nodes are ordinary classes modulo R. A fixed-spine presentation has arcs i→i+1. Native complement flips every **nonpath** arc and retains these path arcs. Reflection of the spine realizes converse. A dart is one actual fixed-spine presentation; native complement pairs darts into lines. Parallel lines are retained, and a loop contributes two diagonal darts.

Write A=A_B+A_K for the merged symmetric integer dart matrix, split into blue and black. A presentation is blue exactly when reversing its path order is an anti-automorphism. For a merged node u represented by T,

    d_u = sum_v A(u,v) = |R(T)| H(T)/|Aut(T)|,

where |R(T)| is1 or2 according as the ordinary class is self-converse or chiral. The ordinary native matrix has degree H(T)/|Aut(T)|. Both matrices, and both blue parts, are exported in [data.json](data.json); this preserves the unmerged gauge for the pairing lane. Black is always A minus blue.

The normalized node walk is P=D^(-1)A, where D=diag(d). Colored subkernels P_B=D^(-1)A_B and P_K=D^(-1)A_K use the same total source mass, not separate color denominators. Uniform physical path sampling in one ordinary representative induces the same merged-node probability law; the free path-orbit and converse factors cancel in that normalized ratio. Those factors remain necessary for raw incidence counts and literal darts.

For T let F_j(T) count ordered covers into j nonempty directed paths. The complete slot polynomial is

    Q_T(t)=sum_j F_j(T) binom(t,j).

Q is invariant under converse. The binomial and monomial coefficient vectors determine each other, so the initial partition Π_Q is equality of the full polynomial, not of a few evaluations. H=Q(1) supplies no extra coordinate once Q is retained.

## 2. The coarsest stable partition is canonical

Let V be any finite state set, let M_1,...,M_s be specified matrices indexed by V, and let Π_0 be an initial observation partition. Define Π_(r+1) by refining each Π_r cell according to the signatures

    u ↦ ( current cell of u,
          (sum_(v∈C) M_c(u,v))_(c=1..s, C∈Π_r) ).       (1)

Stop at the first unchanged partition Π_*. Then Π_* is the **unique coarsest refinement of Π_0** for which every displayed cell sum is constant within each source cell. There are at most |V|−|Π_0| strict rounds.

**Proof.** Each strict round increases the finite number of cells. At termination the defining row sums are constant, so the partition is stable. If Θ is any stable refinement of Π_0, it refines every Π_r by induction: every Π_r cell is a union of Θ cells, whose row sums are separately constant on a Θ source cell; summing them makes (1) constant there. Thus Θ refines the limit. This proves the minimality and uniqueness claimed. This is the standard partition-refinement argument, not a new general algorithmic theorem.

For M=A this is an equitable partition preserving exact dart-weighted walk counts after projection. If L is the cell-indicator lift, the stable quotient matrix B satisfies AL=LB, hence A^kL=LB^k for every k. For M=P it is strong lumpability: every source state in a cell has the same probability of entering each target cell. Taking M_1=P_B,M_2=P_K preserves the joint color and target distribution; the same intertwining works for every colored word. Thus the closure is for every starting state, not merely a stationary average.

There is a useful all-size normalization rule. **If d is constant within every initial cell, then the A and P refinement sequences coincide at every round.** Within a source cell, dividing every row-sum signature by the same positive number does not change equality. Every subsequent cell retains that constant d, proving the assertion inductively. The same statement holds for the two colored matrices. Without retaining d, normalization can merge signatures which literal row sums separate; Section4 gives the first example in this census.

This theorem classifies the smallest extra partition information once the actual matrices are known. It is not a recurrence constructing those matrices from Q alone. Counting “one scalar” is not a meaningful alternative minimality notion, since an arbitrary scalar could encode every class. Sections5–6 therefore restrict the allowed additional measurements explicitly.

## 3. Complete closure census through order seven

The exact universe comprises every labelled tournament at n=3..7, partitioned into ordinary relabeling orbits; every fixed-spine presentation; and every physical Hamiltonian path in every ordinary representative. There are no tournament filters. Refinement counts include the initial partition and omit the final duplicate stable count here.

| n | merged nodes | Q cells | augmented Q cells | A and colored-A closure | P closure | colored-P closure |
|---:|---:|---:|---:|---|---|---|
|3|2|2|2|2|2|2|
|4|3|3|3|3|3|3|
|5|10|8|9|8→10|8→10|8→10|
|6|34|24|29|24→34|24→34|24→34|
|7|272|156|204|156→272|156→271→272|156→272|

“Augmented Q” means (Q,folded scores,Aut,SC,d). Its closure is discrete after one strict round at5,6,7 under each of the four matrix choices. Full Q plus d begins with9,28,190 cells at5,6,7 and also reaches the discrete partition in one round.

Consequently every exact memoryless partition model which refines Q needs all merged nodes at these orders. This is a finite result about partitions of the existing class space. It does not rule out history-dependent models, arbitrary encodings, changing-size recurrences or a different algebraic representation. It says nothing about n≥8.

## 4. Why uncolored normalization costs a second round at seven

The unique nonsingleton cell after one P-refinement of Q at n=7 is {5,11}. Its exact data are:

| source code | H | Aut | SC | d | blue degree | complete native row |
|---:|---:|---:|:---:|---:|---:|---|
|5|3|3|no|2|0|A(5,301)=2|
|11|3|3|yes|1|1|A(11,291)=1|

Both have Q=t^7+2t^5. Their target nodes291 and301 have the same polynomial

    Q=t^7+12t^5+24t^3+4t,

and H41, Aut1. Node291 is self-converse with degree41 and blue degree5; node301 is chiral with degree82 and blue degree0. Therefore both source P rows send probability1 to the same initial Q cell. Their raw A signatures instead record masses2 and1, and their colored P signatures record a black move versus a blue move.

The two target states291 and301 split during the first refinement. Their deterministic predecessors consequently split during the second. The mechanism is a lost intermediate class distinction, not a failure of the refinement algorithm. Retaining d or color repairs the delay at the first round.

Detection of non-lumpability and completion of a stable refinement are distinct depths. The prior reversible Gram identity detects a nonzero omitted-coordinate error at walk length2 as soon as a partition is not lumpable. It does not say that one partition-refinement pass restores all missing states. Nor does a count of distinct row signatures equal a linear Krylov rank.

## 5. One intrinsic bit closes the order-five observer

The only Q collisions at5 are {2,5} and {8,10}. Their blue degrees are respectively (0,1) and (1,3). Hence the single Boolean

    b(T) = 1 exactly when the merged node has one blue dart

splits both collision pairs. Thus (Q,b) distinguishes all10 nodes, and is closed for every actual colored or uncolored walk. Zero additional bits fail, while one succeeds: **one binary coordinate is minimum at this order in the ordinary information sense**, and this witness is an existing physical observable rather than an arbitrary class encoding.

This does not continue as a separation claim. The same observer has24 cells at6 and162 at7. Even (Q,blue degree,d) has only29 cells at6 and199 at7, so it is not stable there, since the coarsest stable refinement of Q is discrete.

A second explicit bit uses a fixed actual anchor class: b_8(u)=1 if A(u,8)>0, where code8 is the specified order-five tournament, not merely its ambiguous Q label. The four relevant entries are A(2,8)=2, A(5,8)=0, A(8,8)=0, A(10,8)=2. This one bit also closes Q at5. The transitive-anchor bit alone does not: it separates only the H9 pair. Q together with SC status and the transitive-anchor bit does separate all10 classes.

## 6. Minimum fixed anchor measurements: 1,2,9

Restrict now to measurements from a **fixed list of actual target classes**. For an anchor v, measure A(u,v), or P(u,v) in the normalized regime. In a colored regime the one anchor yields the two-component vector (A_B(u,v),A_K(u,v)), or its total-degree-normalized version. The cost here is the number of named anchors; a colored anchor has two numeric outputs.

At the audited orders, the exact minimum number of anchors needed in addition to Q is:

| n | A anchors | P anchors | colored A anchors | colored P anchors |
|---:|---:|---:|---:|---:|
|3|0|0|0|0|
|4|0|0|0|0|
|5|1|1|1|1|
|6|2|2|2|2|
|7|9|9|9|9|

Because the coarsest stable Q refinement is discrete, “distinguishes all nodes within each Q cell” and “is a stable refinement” are equivalent for these orders. This equivalence is a consequence of the closure census, not an assumption in the optimization.

Witnesses are {8} at5; {81,145} for A at6 and {81,24} for P at6. At7 the witnesses are

    A: {301,291,161,4722,677,288,4708,4642,5749},
    P: {301,291,161,4722,549,288,4728,805,5368}.

The same witnesses work in the corresponding colored regime. These codes identify unrestricted canonical tournament representatives and hence literal physical anchors; none is identified only by Q.

**Exact optimization proof.** The universe consists of unordered pairs of distinct nodes in the same Q cell. A candidate anchor covers exactly the pairs on which its measured outputs differ. An anchor set distinguishes the initial cells iff it covers this universe. There are184 such pairs at7. Duplicate coverage masks may be identified, and a mask contained in another may be discarded: replacing it by its superset never increases cost or loses a covered pair.

The checker branches on an uncovered pair, trying every remaining candidate covering it. Every cover must choose one of these candidates, so the recursion is exhaustive. It memoizes uncovered masks and remaining budgets. A family of pairs with mutually disjoint candidate sets gives a valid lower bound; another bound divides the number of uncovered pairs by the largest one-column coverage. These pruning conditions only discard impossible budgets. A packing proves lower bound5 at7; exhaustive branches reject budgets5,6,7,8, and the displayed nine-anchor sets certify the upper bound. The separate raw/normalized and colored/uncolored searches visit roughly49,000–54,000 distinct states each. Full witnesses, failed budgets and lower-bound pair lists are in [analysis.json](analysis.json).

## 7. Partition closure versus linear closure

The block-constant function space of a partition is a unital algebra under pointwise multiplication. The coarsest stable refinement can equivalently be described as the smallest such algebra which contains the initial cell indicators and is invariant under the selected transition operators. Indeed finite pointwise algebras are precisely spaces of functions constant on a partition, and (1) enforces operator invariance.

The linear Krylov space generated by initial indicators under operator words need not be closed under pointwise products. Thus its dimension is a different minimality question. The prior packet computed exact minimal linear moment dimensions for **score** observations at5,6. Neither that result nor the number271 in Section3 is being relabelled here as a full-Q linear rank calculation. The partition argument above is sufficient for the claimed finite memoryless classification.

## 8. Verification, pins, and connection contract

[orbits.cpp](orbits.cpp) constructs the complete ordinary orbit map by adjacent-transposition generators, starting each orbit at its least unvisited code. It scans all2^21 labelled orientations at7. The Python checker independently verifies unrestricted minimality of every ordinary representative by all n! permutations. Stabilizer orders follow from n!/orbit size.

[check.py](check.py) builds both ordinary and merged colored matrices in two ways: all fixed-spine words with native complement, and all physical Hamiltonian paths of each ordinary representative followed by native complement and division by the free automorphism action. They agree entrywise. Ordered covers are counted by all permutations with break counts, and independently by subset-path dynamic programming followed by ordered subset convolution. The break method uses the exact binomial count of optional breaks after every backward adjacency is forced to break. These are complete polynomials, not sampled evaluations.

Positive controls include class totals, full-Q class counts, Q parity/integrality, H=F_1, source degrees, doubled loops, black twin multiplicities, and the discrete orders3,4. Hostiles include the equal-Q H9 pair, the delayed order-seven deterministic pair, and the eight-anchor impossibility. The independent pairing lane's complete four matrices agree with ours at3..6; its later audit owns any additional independent order-seven comparison.

The snapshot pin is a8b6d8dbf0c6b601813785efb0b5904b46c23d43. These instructions compile only the local standard C++17 orbit helper; all Python calculations use the standard library. Reproduce:

```bash
python3 05-knowledge/results/incidence_closure_20261113/closure/check.py
python3 -O 05-knowledge/results/incidence_closure_20261113/closure/check.py --output data.optimized.json
python3 05-knowledge/results/incidence_closure_20261113/closure/analyze.py
python3 -O 05-knowledge/results/incidence_closure_20261113/closure/analyze.py
```

Source→target map: actual path-marked tournaments project to merged nodes and then Q cells. Q preserves all uniform-exterior scalar fillings but destroys path-to-partner incidence. The canonical restoration is the coarsest stable partition for the specified operator palette; named anchor columns give a separately optimized physical measurement scheme. The cheapest decisive probe is a pair of equal-Q classes with unequal row sums to one current cell. No named finite quotient here is promoted to an unrestricted class-wiring recurrence.
