# Factorial frozen-tail fibers and their signed-permutation geometry

Date: 2026-10-19. Status: **PROVED**, with [independent review](review-algebra/REVIEW.md) and exact finite controls below.
No external-priority claim. The Graceful Tree Conjecture and the repository's global one-color exchange conjecture remain open.

## 1. Inheritance and the neglected operation

A graceful tree on the fixed labels `0,...,m` chooses exactly one edge of each absolute-difference color `1,...,m`. The prior note [the frozen-fiber obstruction note](../graceful_fiber_obstructions_20261019.md) freezes the upper colors and constructs a fiber with at least k components under one-color exchanges: a distinguished isolated fixed-forest component must retain the color of its unique attachment. The subsequent [two-color buffer](../graceful_two_color_buffer_20261019/GRACEFUL_TWO_COLOR_BUFFER.md) repairs that family using at most two colors. The general [upper-tail connectivity](../graceful_upper_tail_two_color_connectivity_20261019/UPPER_TAIL_TWO_COLOR_CONNECTIVITY.md) proves width-two connectivity of every upper-tail fiber, with an exponential sufficient diameter bound.

Here the underused operation is to **duplicate the isolated attachment at separated numerical positions**. The full attachment assignment becomes a permutation. This amplifies the earlier k-component obstruction to exactly k! components, while making the entire width-two metric explicit.

The canonical hostile is the old single-leaf attachment-color invariant. The corrected near miss is to extrapolate its small counts 2,4,8: those are finite observations, not an exponential theorem. The least-used sidecar is the assignment of colors to individual fixed-forest components, not merely the set of available color contacts.

Concept board: fixed-forest contraction; one-color tree exchanges; perfect matchings; signed permutations; word metrics. All are connected by the explicit maps below, with no imposed tournament orientation.

## 2. An explicit upper-tail forest

Fix k>=1 and m>=k²+2k. Put L=k+1 and

    b_i=iL,       i=1,...,k,
    S={b_1,...,b_k}.

Freeze the edge of every color d=k+1,...,m as follows:

    color d not in S:       {0,d};
    color b_i:             {i,b_i+i}.

Call the resulting graph F. Each replacement really has the indicated color. Its upper endpoint lies in the label interval because `b_i+i <= k(k+2) <= m`, and is not in S because its residue modulo L is i in {1,...,k}.

**Lemma 1.** F is a forest whose components are the k singleton vertices b_i and one tree H containing every other label.

**Proof.** Initially the retained star edges join 0 to all labels greater than k outside S. Each label i in {1,...,k} is attached by exactly the edge `{i,b_i+i}` to that star. Thus the nontrivial part is a tree and includes precisely the complement of S. No selected edge touches S. There are m-k selected edges and k+1 components. QED.

For every i and every free color d in {1,...,k}, both labels b_i-d and b_i+d lie in the interval and in H: b_1-k=1, b_k+k=k²+2k<=m, and distinct b_i are separated by at least k+1. Thus the exact nonloop color-d contacts after contracting F are the k spokes H--b_i, each with exactly two actual edge realizations

    {b_i-d,b_i},       {b_i,b_i+d}.

Every other color-d edge has both endpoints in H and is a loop after contraction.

## 3. Exact state and component theorem

**Theorem 2.** The graceful completions containing F are in bijection with

    S_k x {-,+}^k.

They number `k! 2^k`. Under single-color exchanges the state graph is the disjoint union of exactly k! copies of the k-dimensional cube Q_k.

**Proof.** A completion contracts to a tree on H,b_1,...,b_k. As every allowed nonloop contact is a spoke, each b_i must receive exactly one edge. There are k edges and k colors, so their assignment to the leaves is a permutation: let pi(d)=i mean color d attaches b_i. Independently choose a sign epsilon_d selecting its left or right realization. Conversely every such assignment and sign vector makes a tree and uses each color once, hence gives a graceful completion.

A single-color exchange cannot change pi(d): that would leave its former leaf unattached and give its new leaf two parallel contacts to H. It can change the sign epsilon_d, and this always remains a tree. Therefore each fixed-permutation fiber is exactly Q_k, and different permutations lie in different components. QED.

The smallest interval supplied by this construction has `(k+1)^2` vertices. Hence there are examples with exactly

    floor(sqrt(n))-1 factorial

components on every interval with n>=4 vertices, by taking k=floor(sqrt(n))-1 and including all extra labels in H as above. More explicitly the component count is `(floor(sqrt(n))-1)!`; this is a construction, not an extremal upper bound. Along n=(k+1)^2 its logarithm is asymptotic to sqrt(n) log(n)/2 by Stirling's formula.

## 4. Exact width-two distance

Allow a move to change one or two color coordinates, requiring the endpoint to remain a graceful tree. In the signed-permutation description, a move either:

1. keeps pi fixed and changes at most two signs; or
2. transposes the leaves assigned to two colors and independently chooses their two new signs.

These exhaust the moves: two bijections differing in at most two positions either agree or differ by a transposition.

For states x=(pi,epsilon) and y=(tau,eta), let r=k-c, where c is the number of cycles, including fixed points, of the relative permutation `tau^{-1} pi`. Let f be the number of colors d for which `pi(d)=tau(d)` but `epsilon_d != eta_d`.

**Theorem 3.** Their width-two distance is exactly

    dist_2(x,y) = r + ceil(f/2).

Consequently the full width-two graph is connected and has diameter k-1 for k>=2. For k=1 the diameter is one.

**Upper bound.** Sort every nontrivial relative cycle of length ell by ell-1 transpositions. This uses r transpositions in total and never touches initially fixed colors. Each moved color participates in some transposition; assign its desired sign on its final participation. Then correct the f sign errors on fixed colors, two at a time, using ceil(f/2) additional moves.

**Lower bound.** Consider any path with q transposition moves and h other moves. Relative to the target, each transposition changes the permutation cycle count by one. If b of these transpositions decrease that count, then

    q = r + 2b.

Let t count the initially fixed, sign-incorrect colors that are touched by at least one transposition along this path. Immediately before its first such touch, each is still a singleton cycle of the current relative permutation. Its first-touch transposition therefore merges that singleton with another cycle and decreases the cycle count. Such a transposition first-touches at most two of these colors. Thus `2b>=t`, and `q>=r+t`.

The other f-t sign-incorrect fixed colors never participate in transpositions. They require non-transposition moves, each of which changes at most two signs, so `h>=ceil((f-t)/2)`. Hence

    q+h >= r+t+ceil((f-t)/2) >= r+ceil(f/2).

This proves the formula.

For the diameter, a k-cycle gives distance k-1. If s>=2 colors are nonfixed and their permutation has c'>=1 nontrivial cycles, then r=s-c'<=s-1 and f<=k-s, yielding `r+ceil(f/2)<=k-1`. If s=0 the distance is at most ceil(k/2)<=k-1 for k>=2. QED.

This linear metric is specific to the constructed star-contact family; it does not improve the general rank-staircase diameter theorem.

## 5. Connection contract and scope

- Source: gracefully labeled trees containing the explicitly frozen upper-tail forest F.
- Target: signed perfect matchings between the k colors and k isolated fixed components; equivalently `S_k x {-,+}^k`.
- Map: contract H, retain each spoke's color and its left/right realization.
- Preserved predicates: completion validity; exact one-color and width-two adjacency; shortest width-two distance.
- Lost data: the internal shape of H and actual vertex labels. The formulas for F and b_i restore them. Retaining only the uncolored contact star would lose the permutation and signs and erase the obstruction entirely.
- Transfer: permutation transposition length computes reconfiguration distance; the hypercube sign factor computes every one-color component exactly.

The state set has the cardinality of the hyperoctahedral group, but the theorem only requires the explicit signed-permutation coordinates and moves. No assertion about a standard Coxeter generating set is made.

The construction changes the abstract tree shape under exchange. It neither gracefully labels every abstract tree nor proves global one-color connectedness when the high colors are released. It strengthens the obstruction for frozen upper tails and fully resolves the geometry of this particular enlarged family.

## 6. Verification and literature

Run `python3 verify_factorial_fibers.py` and `python3 -O verify_factorial_fibers.py` from this directory. The checker contracts the actual arithmetic forest, independently enumerates its nonloop edge choices and all acyclic rainbow completions, computes one-color components, checks the hypercube description, and compares every width-two distance at k<=4 with the formula. Larger construction checks do not substitute for the proof.

The relevant established setting is rainbow spanning-tree and independent-transversal reconfiguration. D. Fon-Der-Flaass, *Local equivalence of transversals in matroids*, Electronic Journal of Combinatorics 3(1), R24 (1996), https://doi.org/10.37236/1248, studies one-coordinate equivalence and supplies sufficient conditions and diameter questions. Its full-rank strict-Hall theorem is already imported in [three-color release](../graceful_three_color_release_20261019/THREE_COLOR_RELEASE_DICHOTOMY.md). The present construction is an explicit arithmetic realization, and no priority claim is made for the abstract star/perfect-matching mechanism.
