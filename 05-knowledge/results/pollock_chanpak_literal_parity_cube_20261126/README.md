# Pollock and Chan–Pak share a literal parity cube

2026-11-26. This continues the [figurate four-contact carry](../figurates_local_global_carry_and_norm_20261126/), [independent tournament responses](../independent_responses_and_two_holes_20261126/), and the [fixed-address versus absolute-source distinction](../five_core_web_20261111/fivefold/DYADIC_SHELLS.md). We construct two shared **exact endpoint pairs** with three vertex-disjoint tournament switches, give an infinite strong family of the same switch type, and explain the minimum order. This is a duality of the marked response, not an equivalence of the conjectures or a proof of either.

## 1. Pollock's cube: two fibers, same lower moments

Let $F(x)=(x^3-x)/6$, $x\ge1$. Write the four-index Walsh coordinates as $S,U,V,W$, with $S$ the ordinary sum and the other three the pair cuts, as in the prior note. Three independent sign flips of $U,V,W$ leave $S$ and $R=U^2+V^2+W^2$ unchanged. The four-term total is
\[
 (S^3+3SR-16S)/96+UVW/16 .
\]
When $S$ is even and $S\ge4+|U|+|V|+|W|$, all eight sign choices have positive integer indices provided the initial one is integral. Thus the total is exactly a constant plus an all-three Walsh character. The action of even numbers of sign flips on each four-state parity side is a regular Klein-four action, and merely permutes the four indices.

Equivalently take positive integers $t,a,b,d$ and a fifth, fixed index $z$. The low and high tuples are respectively
\[
 (t,t+a+b,t+a+d,t+b+d,z),\qquad
 (t+a,t+b,t+d,t+a+b+d,z).
\]
Their gap is $abd$ and their index sum and index square sum agree. This is the earlier product-of-three-directions carry, now viewed as a pure parity response on a cube. Walsh contrasts of degrees one and two vanish **over the full cube**; first differences from a fixed corner do not vanish.

## 2. A literal seven-vertex tournament cube

Use vertices $0,\ldots,6$. The base out-neighbor rows are
\[
(123;\ 45;\ 1346;\ 156;\ 0356;\ 026;\ 01).
\]
For clarity, the third row (for vertex 2) is $\{1,3,4,6\}$. Keep all arcs fixed except the disjoint arcs $01,23,45$, forward at the indicated base. The following three permutations, together with the identity, preserve the rest of the marked digraph while permuting the marked pairs:
\[
 (2\,3)(4\,5),\qquad (0\,1)(2\,4)(3\,5),\qquad (0\,1)(2\,5)(3\,4).
\]
They form a Klein four-group acting transitively and freely on each of the two parity sides of the switch cube. In particular they are actual isomorphisms between tournaments on the same side, not just equality of path counts. They fix vertex 6. The correspondence with the Pollock cube is an isomorphism of the **two regular parity orbits** after choosing one representative on each side. It need not identify a named individual cut sign with a named tournament arc permutation.

The two sides have Hamiltonian path totals **141 (odd switch parity) and 143 (even)**, and all states are strong. Counts and strong connectivity are verified directly by the independent subset-path recurrence in the replay. The path-count function has just its constant and its degree-three Walsh terms.

**Minimality for disjoint arcs.** Fewer than six vertices cannot support three disjoint arcs. On six vertices the third mixed response is always zero for such a triple. Indeed every path not using all three pairs as consecutive vertices cancels from the signed sum. Contract each pair to a block with its entrance and exit chosen by the arc direction. For a fixed order of the three blocks a contributing path has two connectors with indicators $e,f$. Pair this order with the reverse order and reverse all three internal directions. This negates the Walsh sign and complements each connector, giving the signed sum of $ef-(1-e)(1-f)=e+f-1$. Each of these three terms omits at least one internal direction, so each sums to zero. Thus order seven is minimal for a nonconstant pure parity cube on disjoint arc flips. Our exhaustive fixed-matching control at orders six and seven confirms zero nonzero third responses at six and finds 24 parity cubes at seven (all strong); matching labels can be fixed without loss by relabeling. This minimality does not concern triples whose arcs share vertices.

## 3. A Stirling family of strong cubes, and two common absolute pairs

Replace vertex 6 by any nonempty tournament module $M$, with identical external directions for its vertices. Let $A_r(M)$ count ordered covers of $M$ by $r$ nonempty directed paths: equivalently sum the product of internal Hamiltonian counts over ordered partitions into $r$ nonempty subsets. A Hamiltonian path of the enlarged tournament consists of an ordering of the six outside vertices, permissible slots for the module runs, and an ordered path cover. So on either parity side the count is $\sum c_r A_r(M)$.

The exact coefficients on the two sides, computed by the slot definition and replayed independently, are:

| Side | $c_1$ | $c_2$ | $c_3$ | $c_r,\ r\ge4$ |
|---|---:|---:|---:|---:|
| Odd | 141 | 136 | 38 | 0 |
| Even | 143 | 137 | 38 | 0 |

To reproduce the slot coefficient, enumerate permutations of the six outside vertices. Every backward adjacent outside step requires a module slot, and that slot and any discretionary slot must have the requisite arcs to and from vertex 6. If $z$ slots are required and there are $f$ permissible optional slots, add $\binom{f}{r-z}$ to $c_r$. At most three runs can occur, since all but possibly the final run must exit to one of the two out-neighbors of vertex 6. The permutations above extend by fixing the whole module. The strong connectivity of the base state extends too: outside vertices can reach any module vertex and vice versa using the quotient paths, and two module vertices can use a walk through the outside and then return if necessary. The original three switches remain individual arcs.

For a transitive $k$-vertex module, $A_r=r!{k\brace r}$, the number of onto assignments to ordered run labels: each used subset has one Hamiltonian path. In particular the two totals are
\[
 L_k=38\cdot 3^k+22\cdot2^k-17,\qquad L_k+2^k .
\]
This is an all-height, three-switch pure parity response, not a set of all unit flips. The gap grows dyadically while its starting count grows on a ternary clock; the quotient is still a constraint. For a general module the gap is $2A_1(M)+A_2(M)$, always positive and even by swapping the two ordered parts in $A_2(M)$.

Now normalize a tournament total by $H=2N+1$, so its increments get the same units as tetrahedral increments. The following exact endpoint pairs agree **in both problems**, rather than just sharing an amplitude:

| Module | Strong tournament counts | Pollock targets | Two lists of positive centered-convention indices |
|---|---|---|---|
| Four vertices, reverse bit 2 of the transitive order | 5065, 5097 | 2532, 2548 | $(11,15,17,17,9)$ and $(13,13,15,19,9)$ |
| Five vertices, reverse mask 40 | 20247, 20355 | 10123, 10177 | $(1,8,11,16,38)$ and $(2,7,10,17,38)$ |

Masks mark reversed pairs in lexicographic order, starting at bit zero. Their ordered cover profiles $(A_1,A_2,A_3)$ are $(5,22,36)$ and $(15,78,198)$. The Pollock parameters $(t,a,b,d,z)$ are $(11,2,2,4,9)$ and $(1,1,6,9,38)$. Each row can be checked by substituting all five $F(x)$ and directly counting all eight literal enlarged tournaments. For comparison, ordinary ordered join of the seven-vertex cube with two 3-count and two 5-count factors gives another common pair $(15862,16087)$ after normalization (not a strong enlarged cube). The corresponding Pollock parameters are $(5,3,3,25,28)$.

No inference of unrestricted coverage follows. In particular using *only* transitive modules and equating their normalized endpoints to the mixed-difference Pollock family need not succeed: an exhaustive scan of that exact subfamily has no match through $k=26$. The nontransitive module is a useful genuine knob, changing the absolute target and its companion together.

## 4. Triangular numbers in other families: a compulsory layer in every tournament

The module transform supplies a triangular-number family independent of the ordinary polygonal/polyhedral definitions. For every $k$-vertex tournament $M$,
\[
 \frac{A_{k-1}(M)}{(k-1)!}=\binom{k}{2}
\]
for $k\ge2$. The partition into $k-1$ parts merges exactly one pair, and the pair has one directed path. The next layer, for $k\ge3$, retains orientation:
\[
 \frac{A_{k-2}(M)}{(k-2)!}
 =\binom{k}{3}+3\binom{k}{4}+2c_3(M).
\]
Here there is either one triple (one path if transitive, three if cyclic) or two pairs (three ways to match four selected vertices). Thus triangular numbers sit on the first subdiagonal of a **weighted Stirling triangle**; the tetrahedral binomial count and the actual directed-triangle correction appear one step later. This is a map from micro choices to macro Hamiltonian path counts, not a dimensional analogy. The triangular number also appears as the dimension of the independent-arc cube $\binom{k}{2}$ and as the maximal inversion number of a permutation (degree of the Mahonian polynomial $\prod_{j=1}^k(1+q+\cdots+q^{j-1})$). Those are other natural filtrations containing the triangular sequence, but do not inherit the Stirling run law.

## 5. The Collatz comparison: legal source versus formal response

For the shortcut map, fix a length-$n$ parity word $(\epsilon_j)_{j=0}^{n-1}$. The composition on a legal source is $(A x+B)/2^n$, where
\[
 A=\prod_j(1+2\epsilon_j),\qquad
 B=\sum_j 2^j\epsilon_j\prod_{i>j}(1+2\epsilon_i).
\]
Each formal word does have one legal source residue mod $2^n$: induct on length, since changing the source by $2^n$ after fixing length $n$ changes the next state by the odd number $A$, so it toggles the next bit. This is the address guard behind the [weighted Collatz passages](../collatz_weighted_passages_first_crossing_20261126/).

As formal polynomials, holding unmarked bits zero, the mixed response of the numerator $Ax+B$ on a marked set $I$ of size $r>0$ is $2^{r-1}(2x+2^{\min I})$: in $B$ the only monomial contributing to that coefficient has its additive bit at $\min I$. This is an exact small-to-large multiplicative/additive ledger, but the eight formal outputs do **not** use the same legal initial integer. By contrast the tournament eight states are literal, and the Pollock paired endpoints agree in lower moments. The comparable problem in all three settings is to retain the source guard along with the response; the shared cube does not transfer a convergence or coverage theorem.

Replay: `python replay.py`; the optional fixed-matching census is `g++ -O3 cube.cpp -o cube && ./cube 6 && ./cube 7`. The optional transitive-module no-match scan is `python common.py`.
