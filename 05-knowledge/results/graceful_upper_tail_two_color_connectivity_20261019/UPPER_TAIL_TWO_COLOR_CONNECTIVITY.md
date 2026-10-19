# Rank-staircase reconfiguration and arbitrary graceful upper-tail fibers

## 1. General matroid rank-staircase theorem

Let M be a finite matroid with rank function rho and closure cl. Let A_1,...,A_s be pairwise disjoint subsets of its ground set satisfying

    rho(A_i) >= s-i+1,  i=1,...,s.

Any two independent selections choosing exactly one element from each A_i are connected by at most 2^s-1 operations of these kinds:

1. Replace one selected element within its color class, preserving independence.
2. Replace two selected elements within their respective classes, preserving independence and the closure of the entire selected set.

Loops may belong to a class, but never to an independent selection. Parallel elements and distinct realizations are retained as distinct ground elements. The statement concerns actual elements, not simplification classes. For a graphic matroid, closure preservation is equivalent to preservation of the connected-component partition on the fixed vertex set.

### Proof by induction

For s=0 there is one empty selection and no moves. For s=1, any two permitted elements are nonloops and a single replacement suffices.

For s>=2, apply induction in M to the last s-1 classes, whose required rank inequalities are unchanged, and retain a current element a in A_1. A recursive closure-preserving double replacement I -> I' preserves compatibility with a, since cl(I)=cl(I'). It also preserves the closure of I union {a}.

For a recursive single replacement e -> f, write I=J union {e} and I'=J union {f}, where J is independent of size s-2. Work in N=M/J. Both e and f are nonloops. Since A_1 is disjoint from J,

    rho_N(A_1) = rho_M(A_1 union J)-rho_M(J) >= s-(s-2)=2.

The current pair {a,e} is independent in N. If {a,f} is independent, change e to f directly. Otherwise a and f are parallel nonloops. Choose a' in A_1 such that {a',f} is independent, possible because rho_N(A_1)>=2. If {a',e} is independent too, first change a to a', then e to f. Both intermediate selections are independent. Otherwise a' and e are parallel nonloops. The old pair {a,e} and new pair {a',f} occupy the same two distinct parallel classes, and so have equal closure in N and are independent. Replace both simultaneously.

Equality of closures lifts from N back to M after adjoining J: for X disjoint from J, cl_(M/J)(X)=cl_M(J union X) minus J, while both closures in M contain J. Thus this last operation has the required full closure invariant.

Each recursive operation costs at most two operations. Once the last s-1 colors are at their exact targets, replace the A_1 element by its target; independence follows from the prescribed full target selection. With L_0=0 and L_s<=2 L_(s-1)+1 we obtain L_s<=2^s-1. Omit replacements that make no change.

### Deterministic construction and existence

Fix a total order on each finite class. Recursively produce the path on classes A_2,...,A_s. For each single step, use the preceding tests in their stated order, and choose the first a' with {a',f} independent in M/J. Independence in M/J is tested by independence of J union the proposed pair in M. This is a terminating independence-oracle algorithm; no enumeration of the reconfiguration graph is needed. The same rank inequalities also guarantee that some independent selection exists: choose from A_s,A_(s-1),...,A_1 in that order. When choosing from A_i, the already chosen set has size s-i, strictly smaller than rho(A_i), so it cannot span A_i.

## 2. Lifting partial operations using one connected color

Suppose r>=1, and F has r+1 components and r free colors. Suppose Q_1, the color-one contact graph on those components, is connected. Suppose colors 2,…,r meet the staircase

    rank(Q_d) ≥ r+1−d,  d=2,…,r.

Then every two rainbow spanning-tree completions of F are connected by at-most-two-color moves, with diameter at most 2^(r−1).

Remove color one from the source and target and apply the partial theorem with s=r−1. Every partial state is a forest with two components after adjoining F. Lift each partial operation to a full completion.

- A closure-preserving two-color operation keeps those two components unchanged, so the current color-one edge still reconnects them. Perform the operation directly.
- For a single-coordinate replacement e -> f, the proposed partial state is an independent forest with exactly two components. Retain the current color-one edge if it crosses those components. Otherwise choose any color-one edge crossing them; one exists because Q_1 is connected. The resulting full tree differs in at most the changed partial color and color one.

After at most 2^(r−1)−1 such moves, every color except one agrees with the target. A final single-color replacement of the color-one edge reaches the target exactly. For r=0 the fixed forest is already a spanning tree and the fiber is a singleton of diameter zero. For r=1 this is simply a nonempty one-color fiber, a clique, and the same bound 1 applies.

## 3. Arbitrary integer-difference upper-tail fibers

Let m>=0 and 0<=r<=m. On vertices 0,…,m, fix any acyclic selection F containing one edge of each difference color r+1,…,m. The completion fiber consists of all gracefully labeled trees containing F, with colors 1,…,r free.

**Theorem.** This entire fiber is connected under moves replacing at most two free-color edges simultaneously, with diameter at most 2^(r-1) for r>=1, and zero for r=0.

For r=0 the assertion is immediate. Otherwise the forest F has r+1 components. The full distance-d graph on0,…,m has exactly d components, its residue chains modulo d. Adding F cannot increase the component count. Thus the color-d contact graph Q_d on the r+1 contracted F components has at most d connected components, giving

    rank(Q_d) ≥ r+1−d.

For d=1 this also says Q_1 is connected. All hypotheses of the connected-buffer-color theorem are therefore satisfied, proving the conclusion. Nonemptiness is independently guaranteed by the existing descending upper-tail completion construction; alternatively the same staircase supports greedy reverse-order selection of an acyclic rainbow completion.

The case r=m has no frozen upper tail and proves width-two connectivity of the full space of gracefully labeled trees on the fixed label interval. The fixed-forest conclusion is stronger: every allowed upper-tail choice remains unchanged throughout the reconfiguration.

## 4. What is and is not settled

The original one-color upper-tail fibers can be disconnected; the previous hostile family proves that width one is not universally sufficient. The theorem shows that width two suffices for every upper-tail fiber. The exponential bound is sufficient, not claimed optimal. Previously proved small-r or special-family bounds remain sharper.

One-color global GT1 remains open here. The operations may change the abstract tree type: this theorem connects existing gracefully labeled trees and does not prove that every unlabeled tree has a graceful labeling. It does not settle the Graceful Tree Conjecture.

The proof uses elementary forest rank, parallel contacts after contraction, and connected-cut repair. No historical-priority claim is made without comparison with the established rainbow-tree and matroid reconfiguration literature.

## 5. Deterministic graphic algorithm

Order each distance-color class by its smaller endpoint. Apply the matroid algorithm to colors 2 through r in the graphic matroid contracted by F, implementing independence by connected-component tests on F plus the proposed edges. Lift each returned single step by retaining color one whenever possible and otherwise taking the first color-one edge crossing the proposed two-component cut. Lift closure-preserving double steps without changing color one. Finally install the target color-one edge. A direct verifier checks the fixed tail, exact difference colors, and connectivity after every full move. With m edges on m+1 vertices, connectivity proves the tree invariant.

The accompanying controls exercise this deterministic procedure on global fibers and nonempty frozen tails. These finite controls supplement the proof; the universal claim rests on the rank-staircase induction.
