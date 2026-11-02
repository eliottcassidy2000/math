# Connected cross-parent ear atlases, exact two-ear seams, and unit-move obstructions

**Session137 status:** PROVED elementary consequences of inherited results; the parent-graph theorem also uses the CITED strong-subtournament theorem already imported in THM-4097. The small atlases and substitution controls are FINITE-EXACT. Full odd H-spectrum completeness except 7 and 21 remains OPEN.

## 0. Inheritance, scope, and concept board

Source revision: `bb7cb2d1f61a14632389814469317131d897b640`, read-only at `/tmp/session136-source` (the later `/tmp/session137-source` mirror has the same revision). No source files were changed. All research deliverables are in this packet.

Load-bearing inherited results, cited by ID **and** slug:

* `01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md`: nonconstant ears of strong parents, the endpoint/exposed-gap carrier, and Moon's nested strong-subtournament input. Its finite order-nine interval is not extrapolated.
* `01-canon/theorems/THM-4111-uniform-ear-average-and-recursive-selected-bank-growth.md`: unbounded selected-image maxima do not imply interval growth. This is a scope guard, not a new theorem here.
* `01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md`: response lattices, unit flip/exchange components, and scalar interval gluing. Strong order-five code 8 is the canonical gcd-saturated hostile.
* `01-canon/theorems/THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion.md`: the orientation-supported positive cut representation `F_T(S)=H(T)+sum_(i in S,j outside S,i->j)c_ij(T)`.
* `01-canon/theorems/THM-4177-root-split-johnson-current-and-source-sink-boundary-census.md`: the affine skew root current and odd-path capacity decomposition. This least-used relevant sidecar is essential: neither the root current nor its path atoms is claimed new here.
* `01-canon/theorems/THM-082-deletion-contraction-ham-paths.md`: directed exposed-block contraction, used for a second seam proof.
* `05-knowledge/results/session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md`: positive integer path kernel and `c_ij>=2H(T-{i,j})`. Its one-parent holes remain one-parent holes.

Live concepts: (1) all-parent coverage, (2) literal common-child overlap, (3) root-current seam transport, (4) exact arithmetic versus feasible moves, (5) automorphism divisibility. The anchor is cross-parent coverage; the niche is a constructive overlap atlas; the wildcard is symmetry forcing large local seams. No LRC connection is asserted.

**What is new in this packet:** a constructive all-order connected overlap graph with a linear radius bound; an exact common-core inter-ear seam identity and its exact lattice; an unbounded symmetry obstruction to small seams; and the sharp order-six failure of even the full reparented unit-state atlas to connect one globally solid scalar interval.

## 1. Conventions and the common-child square

For a tournament `T`, let `H(T)` count directed Hamiltonian paths, with `H(empty)=1`. An ear `T+x_A` has `x->i` exactly for `i in A`. Write `F_T(A)=H(T+x_A)` and

    E(T)={F_T(A): empty != A != V(T)}.

Take a base tournament `B` on `V`, and two ears

    U_A=B+x_A,       U_C=B+y_C.

Let `D_1` be the completion with `x->y`, and `D_0` the completion with `y->x`. Then, with the new vertex explicitly identified in each cut,

    H(D_1)=F_(U_A)(C)=F_(U_C)(A union {y}),
    H(D_0)=F_(U_A)(C union {x})=F_(U_C)(A).             (1)

These are literal common children, not only equal numerical counts. If `B` is strong and `A,C` are nonempty proper, both parents and both children are strong. To see the last assertion, each new vertex has an incoming and outgoing connection to the strong base, so every pair of vertices is mutually reachable regardless of the inter-ear arc.

Thus all nonconstant ears of one strong base form a **clique in the common-child overlap graph**. Repetition/isomorphism can merge vertices or the two common-child values; no distinctness claim is made without further hypotheses.

## 2. A connected all-order parent-overlap graph

Let `P_n` have as vertices the isomorphism classes of strong tournaments of order `n`. Join two vertices when representatives have a common strong `(n+1)`-child obtained by a nonconstant one-vertex ear from each.

### Theorem A — constructive connectedness

For every `n>=4`, `P_n` is connected. More precisely, if

    J_n=C3[TT_(n-2),K1,K1],

then every vertex is at graph distance at most `n-4` from `[J_n]`. Consequently

    radius(P_n)<=n-4,       diameter(P_n)<=2(n-4).       (2)

The distances are on **isomorphism classes**, not fixed labelled tournaments.

### Proof

Use the strong-subtournament theorem in THM-4097 to choose a nested chain

    T_4 induced in T_5 induced in ... induced in T_n=T,

with every `T_k` strong. No preselected triangle is required. Every strong tournament of order four is `C3[TT_2,K1,K1]`: it contains a directed triangle, and each of the six nonconstant ears over that triangle is this class. Let `{r,t}` be its homogeneous pair.

Starting with `T_n`, successively replace the vertex of `T_k-T_(k-1)` for `k=n,n-1,...,5` by a clone of `r`. All previous replacements are retained. Orient the clone block transitively, for instance put each new clone immediately after `r` in its block order. The current tournament is the strong `T_k` with a transitive block substituted at `r`. Such a substitution is strong: paths in the quotient lift between different blocks, and a path leaving and returning to the block supplies reachability in the reverse direction between its vertices.

Before replacing the next vertex, deleting it leaves a transitive blowup of the strong `T_(k-1)`. Hence the old and new parents share a **strong** `(n-1)`-core. Their attachment cuts to that core are nonconstant, so (1) gives common strong children and an edge of `P_n`. There are `n-4` replacements. At the end the old pair `{r,t}` together with all clones is a transitive module of size `n-2`; the two other vertices form the cyclic quotient with this module. This is `J_n`, proving (2). QED.

The hub count is

    H(J_n)=1+2^(n-2).                                  (3)

This familiar nearly-transitive count is included only as an audit control. By the inherited odd-cycle formula, every cycle uses the two singleton blocks, so at most one can occur in a disjoint collection. Directed odd cycles correspond to odd-cardinality nonempty subsets of the transitive block. Their number is `2^(n-3)`, giving (3).

### Corollary A1 — where a global gap must live

The union of the real intervals `[min E(T),max E(T)]`, over all strong order-`n` parents, is one interval. If `a<b` are consecutive members of the global strong order-`n+1` spectrum, then some parent has consecutive ear values `u<=a<b<=v`. In particular, with `gap(X)` the largest difference of consecutive elements of a finite set (zero for a singleton),

    gap(S_(n+1)^str) <= max_(T strong, |T|=n) gap(E(T)). (4)

Proof: common-child overlap joins the corresponding numerical intervals. Connectedness therefore makes their union an interval. Any point in `(a,b)` belongs to one such parent interval. Its consecutive parent values bracketing that point must bracket all of `(a,b)`, because a parent value is a global value. This proves (4).

This is **not** a solidity theorem. It says a global gap cannot be blamed on disconnected parent ranges: at least one individual parent already crosses it while missing it.

## 3. Exact inter-ear seam from the inherited root current

For the base `B`, let `s_i,e_i` count its Hamiltonian paths starting and ending at `i`. Let `Q_ij` count old-vertex orders in which `i` is immediately followed by `j` and every other adjacency is valid; the marked `ij` gap itself may be invalid. Define

    b_i=s_i+e_i+sum_j Q_ji,
    beta_ij=Q_ij-Q_ji,                beta^T=-beta.      (5)

For a cut indicator `a=1_A`, put

    alpha(A)=b+beta a.                                  (6)

By THM-4177 Section 2,

    alpha_i(A)=c_ix(B+x_A).                             (7)

These are the integer-scale symmetric exposure capacities. In particular they are positive even integers; from the odd-path decomposition,

    alpha_i(A)>=2H(B-i).                                (8)

The right side is independent of `A` but the capacity generally is not.

### Theorem B — exact seam reciprocity

For arbitrary base `B` and arbitrary cuts `A,C`, including constant cuts,

    Delta_B(A,C):=H(D_1)-H(D_0)
       =(1_C-1_A)^T alpha(A)
       =(1_C-1_A)^T alpha(C)
       =(1_C-1_A)^T b + 1_C^T beta 1_A.                (9)

### First proof: the ear marginal

In the parent `U_A`, move the observer's cut from `C` to `C union {x}`. Its exact marginal is

    F_(U_A)(C union {x})-F_(U_A)(C)
      =sum_(i in A)alpha_i(A)-sum_(i in C)alpha_i(A).

This is the orientation-supported cut formula: all old-old cut terms cancel, and only arcs incident with `x` remain. Combining with (1) gives the first equality in (9). Since

    (c-a)^T beta(c-a)=0,

one can replace `alpha(A)` by `alpha(C)`. Expansion gives the last equality. QED.

### Second proof: exposed-block contraction

The paths changed by reversing `xy` are exactly those using that adjacency. Contract `x->y`: the block inherits incoming edges from `x` and outgoing edges from `y`. Its count is

    sum_j c_j s_j + sum_i(1-a_i)e_i
      +sum_(i,j)(1-a_i)c_j Q_ij.

The opposite contraction interchanges `a,c`. Their difference is precisely (9). This proof is independent of the cut formula and uses only directed Hamiltonian-path contraction.

### Consequences and boundary

1. `Delta(A,C)=-Delta(C,A)` and `Delta(A,A)=0`. Equal attachment cuts are a homogeneous clone pair, so swapping their inter-ear orientation cannot change H.
2. If `A` is properly contained in `C`, then

       Delta(A,C)=sum_(i in C\A)alpha_i(A)
          >=2 sum_(i in C\A)H(B-i)>0.                  (10)

   The two common children are then distinct numerical overlaps between the same two parents. Reversing the inclusion reverses the sign.
3. For incomparable cuts, positive and negative capacity contributions can cancel. The exact criterion for a unit seam is

       |(1_C-1_A)^T(b+beta 1_A)|=2.                     (11)

   Neither the individual parent H values nor an unoriented capacity multiset determines this condition.
4. Along a single-coordinate attachment change `C=A union {i}`, the seam gap is exactly `alpha_i(A)`. Thus positive capacity gives an obstruction, not a guarantee, to a unit seam: if `H(B-i)>1`, such a seam cannot have magnitude 2.

The transported object is a common-core pair of parent networks with retained cut labels. The information discarded by replacing it by scalar H or by an unlabelled overlap edge is exactly the feasible seam value (9).

## 4. Exact seam lattice and an unbounded all-order obstruction

For `m=|B|>=4`, let

    g_B=gcd({b_i} union {beta_ij:i<j}).                 (12)

### Theorem C — exact common-core seam lattice

The additive subgroup generated by all `Delta_B(A,C)` with `A,C` nonempty proper is exactly `g_B Z`.

Every seam is a combination of the coefficients in (12), so one inclusion follows from (9). Conversely choose distinct `i,j,k`. The adjacent-cut seams

    Delta({j},{i,j})=b_i+beta_ij,
    Delta({k},{i,k})=b_i+beta_ik,
    Delta({j,k},{i,j,k})=b_i+beta_ij+beta_ik

are all permitted when `m>=4`. The first two minus the third recover `b_i`; subtracting `b_i` recovers every `beta_ij`. Thus both generated groups agree. QED.

All generators are even: `b_i=alpha_i(empty)` is an even exposure capacity, and `beta_ij=alpha_i({j})-alpha_i(empty)` is the difference of two such capacities. Moreover some nested seam is positive by (10), so `g_B` is a positive even integer.

### Theorem D — symmetry forces arbitrarily large seam divisors

Let `Q` be any strong tournament on `k>=3` vertices, and form

    B=Q[C3,C3,...,C3],             |B|=3k.

Then

    2*3^(k-2) divides g_B.                              (13)

In particular, no bound independent of the base order can guarantee a nonzero seam of bounded magnitude among the co-amalgamations of all nonconstant ears of one strong base.

Proof: independent rotations of the `k` triangle blocks are automorphisms. For a path with one prescribed endpoint, rotations of the `k-1` other blocks act freely, because no nonidentity vertex permutation fixes a linear word. Hence `3^(k-1)` divides every `s_i,e_i`. For an exposed word with specified adjacent endpoints `i,j`, at least `k-2` blocks are unmarked; their rotations act freely while preserving the marked word condition, so `3^(k-2)` divides every `Q_ij`. Formula (5) gives that divisor for every `b_i,beta_ij`. Their evenness and coprimality of 2 and 3 prove (13). The base and all its nonconstant ears are strong. Nested cuts provide a nonzero seam, so the assertion is not vacuous. QED.

This does **not** show an obstruction in the global H-spectrum. Changing the base, its old arcs, or other child coordinates can bypass this particular common-core seam lattice.

### Exact lattice is still not a feasible-step theorem

For `B=C3[C3,C3,C3]`, using base masks

    (450,452,449,23,39,15,184,312,120),

one has

    H(B)=3159,       g_B=6,
    min_(A,C nonconstant, Delta!=0)|Delta|=30.           (14)

The scan covers all `510^2=260100` ordered nonconstant attachment pairs and produces 507 distinct signed gaps. A minimal witness is `A=9`, `C=68`; both parent counts are 11463, and the two literal common children have counts

    H(D_1)=50265,       H(D_0)=50295.                    (15)

Thus even the **exact** seam lattice may suggest step 6 while no actual seam realizes 6,12,18,24. This is the cross-parent counterpart of THM-4118's arithmetic-saturation warning, not a global missing-value claim.

For the strong quotient with masks `(8,1,3,6)`, substitution of four C3 blocks gives an order-twelve base with `H=241137` and `g_B=18`, attaining the lower divisor in (13).

## 5. Even the full reparented unit-state atlas misses a scalar seam

Consider all labelled strong child tournaments of a fixed order `N`. Permit:

* changing the distinguished deleted vertex without changing the child (reparenting);
* a one-bit cut flip over any strong deletion parent; or
* a one-in/one-out exchange over any strong deletion parent.

Retain a physical move only when its H change is `0` or `+/-2`, as in THM-4118. After quotienting only the bookkeeping of the distinguished root, the vertices are labelled strong children and the allowed edges are the two corresponding incident-arc changes. The verifier builds **all** such edges, not a selected-parent bank.

FINITE-EXACT results:

| child order | labelled strong children | unit components | component H images |
|---:|---:|---:|---|
| 5 | 544 | 1 | `{9,11,13,15}` |
| 6 | 22320 | 5 | `{15,17,19}`; `{23,25}`; `{27,29,31,33}`; `{37}`; `{41,43,45}` |

At order six the globally attained H values contain the solid odd interval `23:2:33`, but no allowed physical unit path crosses from 25 to 27. Reparenting cannot repair this: it was already included over every strong deletion parent. Relabelling cannot repair it either, since the five components have disjoint H images, which relabelling preserves.

This is sharper than failure of one-parent solidity. It refutes the prospective claim that a globally solid H interval must be witnessed by one connected component of the entire all-parent unit atlas. The survivor is exactly THM-4118's **scalar interval graph**: the intervals `[23,25]` and `[27,33]` touch by 2 and therefore certify their solid numerical union without pretending a physical edge exists.

The first failed implication is

    literal common-child connectedness + saturated small-step arithmetic
        -/-> feasible unit-state connectedness
        -/-> an all-order interval induction.

Conversely, scalar adjacency remains a lawful completeness certificate. The remaining task is to prove enough such scalar adjacent interval coverage uniformly, not merely to prove that parent presentations communicate.

### A useful repair: nonstrong parents retain missing short seams

Object coverage by strong parents does not imply coverage of short moves by strong parents. A complete additional finite search of Hamming-distance-two pairs between the labelled order-six H=25 and H=27 fibres finds exactly **1,440** such bridges (each counted from its H=25 endpoint). Every bridge flips two arcs incident with one common vertex, that vertex's deletion parent is **nonstrong**, and the two cut bits change in the **same** direction. No disjoint-two-arc bridge and no one-in/one-out exchange bridge exists.

A minimal explicit nonstrong parent is

    B=(4,1,2,7,15),       H(B)=3,       B=TT_2 => C3.

The following entire Boolean square consists of **strong children**:

| cut mask | 16 | 17 | 24 | 25 |
|---:|---:|---:|---:|---:|
| H of child | 27 | 23 | 33 | 25 |

Thus the diagonal `25 -> 16`, removing two cut vertices, supplies the missing unit change `25 -> 27`. None of the square's ordinary one-bit edges is a unit change. Both restrictions matter: merely adding nonstrong parents while retaining only flips/exchanges still misses this diagonal; merely adding same-direction two-bit moves while retaining strong parents also misses it.

Adding this one diagonal edge to the full order-six strong-parent unit atlas merges the `{23,25}` and `{27,29,31,33}` components, giving exactly the four maximal odd intervals in the global spectrum. This is a finite repair, **not** an all-order induction.

There is a simple all-order carrier for admitting nonstrong parents safely. If a tournament B has strong components

    C_1 => C_2 => ... => C_k,

then `B+x_S` is strong **if and only if**

    S intersects C_1,       C_k \ S is nonempty.        (16)

Necessity: otherwise the source component has no incoming arc in the child, or the sink component has no outgoing arc. Sufficiency: x can enter C1, traverse all component levels, and return from Ck; strongness inside each component supplies every required endpoint. For k=1 this is the familiar nonconstant-cut criterion. For k>1 the feasible cuts form the product of the nonempty subsets of C1, unrestricted subsets of the middle components, and proper subsets of Ck.

The useful new question is therefore to organize scalar interval certificates using **all strongly completing parent cuts**, with the SCC-endpoint conditions (16), rather than assuming that Moon's object-complete strong-parent cover is also move-complete.

## 6. Reproduction, validation, and honest frontier

Run in this directory:

    python3 -B audit.py > audit.out
    python3 -B -O audit.py > audit_optimized.out
    cmp audit.out audit_optimized.out

The verifier uses only Python's standard library and local `explore.py`, never inherited code. It compares literal exposed-word enumeration with an independently organized subset endpoint convolution; checks the two-ear identity by direct child H-DP for every base and every attachment pair through base order four; uses literal child permutation controls; verifies all clone paths for every labelled strong parent through order six; constructs the entire labelled unit atlas at orders five and six; checks every H=25 to H=27 Hamming-distance-two pair at order six, and checks the order-nine/twelve symmetry examples. Runtime and storage are modest. Explicit checks, not Python `assert`, survive `-O`.

The separate parent-appointed auditor independently reconstructed the order-six graph and the seam identity; see its report in the Session137 audit packet when integrated by the parent. That audit is not simulated by importing this verifier.

**Direct progress:** all-order cross-parent overlap is now constructive and connected, with radius `n-4`; exact seam values and their lattice are available from the inherited root current.

**Obstructions:** local co-amalgamation seams can have arbitrarily large obligatory divisibility, the lattice can understate the smallest feasible step, and all-parent unit-state connectivity can fail inside a globally solid scalar interval.

**Open:** uniform scalar interval coverage/overlap, unbounded right endpoints of solid selected-bank intervals, all-order strong-prime witnesses, and the full positive-odd H-spectrum except 7 and 21. None is claimed here.
