# Replacing a core through one outlet: exact two-cover transfer and a norm obstruction

**Scope.** The formulas below hold for arbitrary nonempty tournament cores, not just transitive spines. The missing values are missing from the stated constructions only. No claim is made about the full tournament spectrum.

Write h(U) for the number of Hamiltonian paths, with h(empty)=1, and
```
Q(U) = sum over S ⊆ V(U) of h(U[S]) h(U[V(U)−S]).
```
Thus Q counts two ordered path slots (either possibly empty) covering U. Put P(U)=Q(U)−h(U). For nonempty U, h is odd by Rédei, Q is even by complementary pairing, and P is odd.

## 1. Arbitrary replacement at a single outlet

Suppose a vertex v of a tournament T has exactly one outneighbor q. Replace v by any nonempty tournament U whose vertices inherit all of v's exterior arcs. Write R=T−v, A=h(R), and let E and O be the following quantities for T: sum h(T−C) over directed cycles C containing v, of even and odd length, respectively. Cycles are counted up to rotation (not reflection). Then the new count is exactly
```
h(T[v←U]) = (A−2E)h(U) + (E+O)Q(U).
```
The converse statement holds for a unique in-neighbor by reversing all arcs.

By the [odd-cycle packing formula](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20), h is the sum of 2 to the number of directed odd cycles over vertex-disjoint packings, including the empty packing. A crossing directed cycle exits U through q. Since a simple cycle visits q at most once, it has exactly one contiguous nonempty U-run; there is at most one crossing cycle in a packing. Contract that run to v to recover an old cycle C. The outside remainder contributes h(T−C). For an old odd C, the U-run must have odd order, and for old even C it must have even positive order (length changes by the run order minus one).

The signed sum of the two-cover contributions by the parity of |S| is zero. Expand both h factors by odd-cycle packings. For each uncolored packing, its cycles and unused singleton vertices give nonempty collections of odd blocks; assigning each block freely to one of two slots cancels the signed sum. Thus each parity has total Q(U)/2. Subtract h(U) from the even total to remove its empty slot. Counting packings not crossing the seam gives A h(U); adding the single crossing cycle with weight 2 proves the formula. This generalizes the two-cover mechanism of the [earlier padded cyclic wrap](../chanpak_cross_problem_theorems_20261114/PRIME_WRAP.md).

## 2. A genuinely two-source identity and what it still cannot reach

Take independent nonempty cores U,V and a singleton q in cyclic orientation U→q→V→U, complete between blocks. Then
```
h(C3[U,q,V]) = P(U)P(V) + 2 h(U)h(V).
```
For the old source T with U replaced by v, R=q→V. The nonempty V-run has even order for an old even crossing cycle and odd order for an old odd one, so E=Q(V)/2−h(V), O=Q(V)/2. Substitute into §1, or directly count the opposite-parity pairs of U- and V-runs in odd cycles using q. U and V can be changed independently.

Two useful lower bounds for a core W of order n are
```
P(W) ≥ n h(W);
P(W) ≥ 2^n−1, improved to P(W) ≥ 3·2^(n−1)−3 if W is cyclic.
```
For the first, all n+1 splits of each Hamiltonian path give distinct ordered covers. For the second, in the two-slot packing expansion a packing with k odd cycles covering s vertices contributes 2^k·2^(n−s+k) to Q: every cycle and every remaining singleton chooses its slot. Subtraction of h removes 2^k. Each term is nonnegative and the empty packing gives 2^n−1. Any cyclic tournament has a directed triangle; the singleton packing of this triangle contributes another 2^(n−1)−2. Equality in the cyclic bound holds for a single triangle plus singleton strongly connected components.

**Exact small all-order spectrum of the two-core singleton wrap through 98.**
```
3, 5, 9, 11, 15, 17, 23, 27, 33, 37, 47, 51, 65, 69, 91, 93, 95.
```
The bounds force each core to have order at most 6; if a core has order 6 it must be transitive, since 93+2·3=99. An independent subset-endpoint path DP exhausts every labelled core through order 5 plus the order-six transitive profile; [replay_one_outlet.py](replay_one_outlet.py) lists exactly the stated set and directly checks literal glued tournaments.

A short symbolic exclusion explains the earlier 63/67/69/71 contrast. For two transitive cores, use (2^a−1)(2^b−1)+2; through 71 the possibilities are 3,5,9,11,17,23,33,47,51,65. Two cyclic cores already contribute at least 9·9+2·3·3=99. If only U is cyclic and V is transitive of order t, then t≥4 overshoots. At t=3 the only possible count in the specified four-target set is 69, attained with U=C3; for U of order at least 4 the lower bound gives 7·21+6>71. At t=2 only cyclic orders 3 and 4 can fit; the respective profiles are (h,Q)=(3,12), and (3,24) or (5,32), giving 33,69,91. At t=1, cyclic orders 3 and 4 yield 15,27,37; order ≥6 overshoots. With n=5, at most one odd cycle can occur in a packing, so N=Q+h=33+18 times the number of triangles +6 times the number of directed 5-cycles. This excludes 67 and 71 modulo 6. For N=63, a cyclic U must have exactly one triangle and two 5-cycles. But if there is only one triangle, every outside vertex is either a uniform winner or uniform loser against it (a mixed vertex creates a second triangle). No pair of such outside vertices can make a cycle with the triangle: a loser beating a winner would create another triangle. The triangle is therefore its entire SCC and there is no 5-cycle. Hence 63 also fails.

The count 69 has two different ways through this exact formula: (h,Q)=(3,12) paired with (1,8), or (3,24) paired with (1,4). At 23 one sees a simpler use of asymmetry: transitive orders 2 and 3 give 23, impossible for the synchronized subset in §3.

## 3. Equal normalized profiles force a quadratic norm

Synchronize the two ratios P(U)/h(U)=P(V)/h(V)=x/y with x,y coprime. Both are odd because the original coordinates are odd. Write h(U)=a y and h(V)=b y. The formula becomes
```
N = a b (x² + 2y²).
```
If N=p is prime, a=b=1. If N=7p with p odd prime, the coprime norm is never divisible by 7 because −2 is a nonsquare modulo 7. It must be p and the two multipliers are 1 and 7 in either order. In particular the cofactor p for either lane must be 3 mod 8.

For p≡3 mod 8 the positive representation p=x²+2y² is unique. Indeed the ring Z[√−2] is norm-Euclidean: rounding the real and √−2 coordinates of a quotient incurs norm at most 1/4+2/4<1. The prime splits by the quadratic-residue test and its two prime factors are associates of conjugates. The only units are ±1, so absolute positive coordinates are unique. They are odd by reducing modulo 8. This gives an **exact target-first criterion**:
- prime p: seek profile (h,Q)=(y,x+y);
- target 7p: seek profiles (y,x+y) and (7y,7(x+y)).

In both cases, any sought core has order n≤floor(x/y) by P≥nh. If the profiles exist, substituting them suffices, and repeating the same witness gives the prime target. Thus this particular inverse search is genuinely finite with an explicit order cutoff read from the target's norm, and fails immediately when that cutoff is too small. Example: p=59 gives (x,y)=(3,5) and fails. Neither criterion is claimed for unsynchronized ratios.

Even the equal-profile subclass misses infinitely many globally attainable norms. Consider N=3·5^(2j), j≥1. Since −2 is a nonsquare modulo 5, any representation x²+2y²=N forces 5 to divide both x,y. Iterating gives (x,y)=(5^j,5^j) for positive coordinates, not necessarily reduced. The bound x≥ny then forces n=1, contradicting h=5^j. Yet N has an ordinal witness: one directed triangle, then 2j copies of a four-vertex h=5 core, since h multiplies over ordinal sums. This tests exact coverage inside the norm congruence, not merely the absence of three other odd residues modulo 8.

## 4. Quantitative loss when synchronizing the two sources

For the equal-profile subclass, the number of distinct outputs at most X is O(X log log X/log X), by an elementary lattice bound. With y≤X^(1/4) there are O(X^(3/4)) choices of x,y bounded by the norm. Otherwise y=h(W)≤n!≤n^n forces n≫log X/log log X, while x≥ny and x≤√X. For each x there are O(1+x log log X/log X) choices of y. This is a capacity bound, not a count of distinct norms after collisions.

There is also an upper bound O(X log log X/log² X) for prime cofactors p≤X in **both** synchronized lanes p and 7p (equal ratios, not necessarily equal profiles). In either lane there is a core with h=y, hence the same angular restriction. In a dyadic shell Z≤p≤2Z, discard the small-y case y≤Z^(1/4) at cost O(Z^(3/4)). The other positive norm generators have angle arctan(√2 y/x)≪log log Z/log Z. Every generator of prime norm gives a distinct prime ideal for a different p. Apply [Coleman, Theorem 3](https://matwbn.icm.edu.pl/ksiazki/aa/aa65/aa6515.pdf#page=21) to the angle sector about the real axis in the fixed field Q(√−2); its width parameter l is a constant times log log Z/log Z, within the stated range Z^(−3(1−ε)/10)<l<1/2. Angle is taken modulo π because the units are ±1. The shell bound is O(Z log log Z/log² Z) plus the discarded O(Z^(3/4)). Sum dyadically down to √X and bound all smaller cofactors by √X. The synchronized family therefore misses a density-one proportion even within the congruence-compatible prime and 7p lanes. This use of angle equidistribution is independent of the exact low-count exclusions.

**Design implication.** The right inverse of a two-core singleton wrap asks for a physical factorization of N−2h(U)h(V). Using independent profiles helps at 23 and 69 but provably cannot reach 63,67,71. An exact induction toward Chan–Pak would need an additional contact mode, a longer boundary, or a different wrap. Synchronizing normalized source profiles imposes a stronger norm obstruction even with unrestricted growth.

## 5. Removing the singleton restriction still does not close the small gaps

For three arbitrary nonempty cores, write L_i(j) for the number of ordered j-tuples of **nonempty** vertex-disjoint paths covering core i (zero for j>its order). A Hamiltonian path in a cyclic three-block wrap is exactly a cyclic run word, together with these ordered path covers on the three sides. Thus the full count is the sum over r≥1 of
```
3 ∏_i L_i(r)
 + sum over i of L_i(r+1)∏_(j≠i) L_j(r)
 + sum over i of L_i(r)∏_(j≠i) L_j(r+1).
```
The three terms count cyclic words of lengths 3r, 3r+1 and 3r+2. In the equal case any of the three starting sides is allowed. In the other cases the unique side, or unique pair of sides, with excess visits determines the starting side. The coefficient is nonnegative for every product. For fixed side order n, L_i(j) is minimized by the transitive tournament: each ordered partition of its vertices into nonempty sets admits a path per set by Rédei, and for a transitive tournament each choice is unique. Those minimum counts are nondecreasing in n for every fixed j by adding an extra vertex to the first slot, except for the harmless impossible j=0 case (never used).

If at least one side is a singleton, §2 already covers the construction. If all three sides have exactly two vertices, the count is 45. If all have at least two and some have at least three, the run formula is bounded below by the transitive model of side sizes (3,2,2), whose count is 137. Thus through 98 the exact spectrum of **arbitrary** cyclic triple wraps is precisely the list in §2 with 45 added. It still misses 63, 67, and 71 irrespective of how freely all three cores are selected. This is a restriction on the cyclic quotient, not on core size, local congruences, or the unrestricted problem.

## 6. A positive repair: change the host to an actual non-cyclic-triple boundary

A single-outlet host is more general than the cyclic triple. Let U be transitive of order m≥1. Then h(U)=1 and Q(U)=2^m, so §1 gives a straight dyadic family
```
N(m) = A−2E + (E+O)·2^m.
```
The following literal hosts exhibit all four small comparison values. R is the exterior tournament, listed as its outgoing bit masks (bit j means an arc to vertex j, indexed from 0). All of U beats the outlet q, and every other vertex of R beats all of U.

| Exterior rows R | q | (A,E,O) | Exact family N(m) | Small comparison |
|---|---:|---:|---:|---:|
| (38,0,2,7,15,30) | 0 | (5,9,10) | 19·2^m−13 | 63 at m=2 |
| (20,1,2,7,14) | 3 | (9,7,11) | 18·2^m−5 | 67 at m=2 |
| (4,1,2,7) | 3 | (3,3,6) | 9·2^m−3 | 69 at m=3 |
| (18,0,3,7,14) | 2 | (5,3,6) | 9·2^m−1 | 71 at m=3 |

In particular, by using width m−1 in the 18-slope host and width m in the two 9-slope hosts, **every** m≥2 gives all three consecutive odd counts
```
9·2^m−5,     9·2^m−3,     9·2^m−1.
```
At m=3 these are 67,69,71. The arithmetic target can choose an actual host (the additive offset) and source width by checking exact divisibility and whether the quotient is a positive power of two. The 63 construction gives the all-order arm 19·2^m−13; the literal witness at m=2 has rows (4,5,152,3,11,31,63,123). Changing the exterior host rather than retuning the three cyclic blocks is the resource which reaches the excluded targets. A finite set of such dyadic families is still sparse: this is an exact changing-boundary construction on infinitely many counts, not an arbitrary-target compiler.
