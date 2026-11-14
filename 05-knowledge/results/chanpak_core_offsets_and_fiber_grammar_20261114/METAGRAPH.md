# The SCC rewrite grammar and the bounded-jump metagraph

**Proved:** the ordered rewrite, pointed and converse inventory, and
fixed-jump cutoff below. **Finite-exact:** the gap-four alphabet, by the
attached replay. The vertex-fiber formula is inherited, not a novelty
claim. This is the **vertex-deletion** metagraph, not the fixed-order
blue or black arc-flip graph.

Write H(empty)=1. Every tournament is uniquely an ordinal word of
strongly connected components (SCCs); all arcs point forward between
components, and H is their product. These are actual isomorphism-class
letters, not just counts. This is the order-join reduction in
[THM-1862](../../../01-canon/theorems/THM-1862-order-join-reduction-principle.md).
The turning point is to mark the child vertex before passing to classes.

## 1. Exact local rules, at all orders

For each strong class C and each Aut(C)-orbit of vertices, pick a
representative v and record the rule

    SCC-word(C-v)  -->  C             [marked by the orbit of v].

Allow C to be singleton: its left side is the empty word. All growth
edges of the vertex-deletion metagraph are precisely the applications
of these rules to contiguous subwords. Indeed, in any pointed child
T,v, there is a unique SCC C containing v. All other SCCs are untouched
by deletion: the fixed forward arcs between the old SCCs will not
merge anything, and refining C replaces it with a contiguous subword.
Conversely, any rule application gives a physical vertex extension by
a forward join with its exterior. The pair consisting of the resulting
child word and its chosen marked letter and vertex orbit is a complete
edge key. Forgetting the key can identify parallel edges, for example
distinct deletions in one singleton run. Do not count those as
distinct edges in a simple class graph. This is an exact structural
recursion; the strong pointed alphabet itself is not a finite alphabet
globally.

On weights, the two local endpoint counts multiply by the same
exterior product m. Weight fibers alone cannot drive this recursion.
The two self-converse order-five H=9 births with outgoing rows
(14,28,24,16,1) and (10,28,25,16,1) have labeled deletion decks
(1,5,5,5,1) and (1,3,5,5,3), respectively. Hence one has 3 as an
immediate predecessor weight and the other does not, even with
order, SCC size and self-converse status fixed. Labels here are only
for replay, not additional metadata in the counterexample.

## 2. Uniformly finite alphabet at any bounded positive jump

The [odd-cycle-packing identity](https://arxiv.org/html/2412.10572v3#S3.SS4)
gives, on pointing v in C,
    
    H(C)-H(C-v) = 2 sum over directed odd cycles Q through v of H(C - V(Q)).

In a packing containing v choose its unique cycle through v and then
any packing on the remainder. Every summand is at least 1. In a
nontrivial strong tournament v lies on a directed cycle of every length
3 through |C| (Moon). For completeness: some out-neighbor of v beats
some in-neighbor, or the out-neighbor set cannot reach v. To extend a
v-cycle of length at least 3, insert any outside vertex with mixed
incidences to it at a suitable transition. If no such vertex exists,
every outside vertex is in A (loses to the cycle) or B (beats it).
Strong connectivity gives an arc A->B, since otherwise A is terminal
(if either set is empty, strong connectivity also fails). In the cycle
replace one non-v vertex by those two vertices in order. Thus

    H(C)-H(C-v) >= 2 floor((|C|-1)/2) >= |C|-2.

In the whole graph a strict jump J equals m times this local jump.
Consequently all jumps are even and nonnegative. Neutrality is
**equivalent** to inserting a singleton at an SCC boundary. Different
orders with the same count may have quite different neutral
backbones. For each fixed ordered backbone of r nontrivial SCCs, its
equal-H inter-order graph is exactly the directed Pascal lattice
of gap lengths (b_0,...,b_r) in N^(r+1); an edge increments precisely
one coordinate. Singleton insertions at several positions in the
same run descend to that same simple edge.

For fixed J>0, the only candidate exterior products are the odd
divisors m of J/2 *that are themselves physically available*. For each
candidate m, enumerate pointed strong C through size J/m+2 and retain
only local jump J/m. The only possible endpoints are then
    
    (m H(C-v), m H(C)).

For each retained C and attainable m these endpoints occur: put an
actual m-witness outside C. Moreover the complete pointed edges are
obtained by using every possible split of every exterior SCC word at
the marked C, with arbitrary singleton padding. For a finite alphabet
of exterior nontrivial blocks it suffices to use the
[Busch bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3),
which gives |D|<=1+floor(log_5(m^3)) per nontrivial exterior block D.
No n-dependent search is needed once J is fixed. In particular
**no strictly positive bounded one-vertex jump can connect unbounded
absolute counts**. This is compatible with exact fixed-gap examples
built by arc switches or changing a core; those are different moves.

* J=2: the only endpoints are 1->3 and 3->5 (enumerate just the
  triangle and the unique strong order-four class). At child order n
  their pointed edge counts are respectively max(n-2,0) and
  2 max(n-3,0). The first has one local root orbit and the second
  has two. Quotienting pointed edges by converse gives
  ceil(max(n-2,0)/2) and max(n-3,0).
* J=4: only 1->5 and 5->9. The order bound is 6; the replay records
  2 local pointed rules at size 4 and 5 at size 5, respectively,
  with none of that jump at size 6.

## 3. Coefficient recurrences and converse

Use formal multiplication [a][b]=[ab] for H labels, with z for order.
Let S be the weighted series of strong classes including the singleton.
The inherited vertex recurrence is V=1/(1-S): the coefficient of
z^n[h] is the number of unlabelled order-n classes of count h and
forgetting H gives A000568. This is the formula in the repo's
[actual Hamiltonian multiplicities packet](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md).

The edge version uses pair labels [a,b] multiplied componentwise.
Let D count the pointed strong rules above, with z-degree |C| and label
[H(C-v),H(C)]. Let Dc be the subset fixed by rooted converse
(there is an anti-isomorphism of C fixing the chosen vertex).
Define Delta V by [h] -> [h,h], and Sigma V by z^n[h] -> z^(2n)[h²,h²].
Then the total and the converse-fixed pointed-edge series are

    E = (Delta V)^2 D,              Ec = (Sigma V) Dc.

The two arbitrary exterior words in E explain its factors. A fixed
pointed child has its marked SCC in the center; right exterior must
be the reversed sequence of converses of the left exterior; the
center must be rooted self-converse. That gives Ec without introducing
false fixed points in an unpointed self-converse class. The series of
converse-orbits of pointed edges is (E+Ec)/2, coefficientwise. This
quotient keeps the marked-edge inventory; it need not be the simple
graph obtained by subsequently collapsing parallel edges. For the
unpointed fixed-class series, define S_c using self-converse strong
classes; the identical palindrome decomposition, now with optional
center, gives sum(V_n[h] z^(2n)[h²])*(1+S_c). These formulas retain exact
multiplicity, rather than averaging over possible extensions from a
count-only state.

The existing finite-inventory formula for vertex fibers can be seen
geometrically here: for a nontrivial backbone of r letters and total
size b, its order-n contribution is binom(n-b+r,r) for n>=b. For
pointed strict growth with r exterior nontrivial blocks, distinguish
the one active block and there are r+2 gaps, giving
binom(n-b+r+1,r+1) with b including that active block. At fixed J the
series over all endpoints is eventually polynomial of degree at most
1+Omega(J/2), because each exterior nontrivial component contributes
a factor at least one prime to m. The alphabet is finite for every
bounded H or positive J, not uniformly finite as the cap grows.
This is not an unrestricted A000568 speedup or a proof of Chan–Pak.

Replay: python3 replay.py 6. The code builds the complete layers by
vertex addition and degree-refined canonicalization, recomputes H
with subset DP, and compares literal pointed orbit enumeration to
the transfer formulas. The all-order statements use the SCC argument
and packing identity, not the six-layer test.
