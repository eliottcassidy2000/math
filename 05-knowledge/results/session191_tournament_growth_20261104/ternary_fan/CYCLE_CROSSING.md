# A one-birth, one-arc ternary fan crossing from one cycle to three

**PROVED, all positive clocks.** For every k>=1 let N_k=9+11·2^k. There is an explicit strong tournament with H=N_k and exactly one Hamiltonian cycle. Adding one vertex, in either of two specified ways, produces strong tournaments with H=3N_k+2 and H=3N_k+4, each with exactly three Hamiltonian cycles. The two children differ in one arc. All existing vertex names, all old arcs, the physical attachment root, and every exterior incidence are retained.

Unlike the first two-port fan, the children deliberately leave the uniquely Hamiltonian class. The same named old vertex has start counts 1,3,4 at the parent and two children. This is an explicit bounded increase in the endpoint budget, not a claim of repeatable ternary closure.

## 1. The three physical cores

On names 0,...,4 take

    P=(14,20,8,2,13),     r=1, v=3.

Add y=5. For the lower child let y beat precisely old vertices {0,3}; for the upper child let y beat precisely {0}. These give

    A=(14,52,40,2,45,9),
    B=(14,52,40,34,45,1).

In both cases every other old vertex beats y. Thus P is the identical induced parent of both children. Reversing the single arc y->3 in A gives B. Each core is strong: P is strong, y reaches old vertex0, and old vertex1 reaches y.

With F_r and J(u,r) as the ordered two-path-cover counts in `PROOF.md`, the exact profiles are:

| Core | H | F_r | e_r | J(u,r), in physical vertex order |
|---|---:|---:|---:|---|
| P | 9 | 11 | 1 | (3,1,1,1,5) |
| A | 29 | 33 | 3 | (4,3,4,3,15,4) |
| B | 31 | 33 | 3 | (3,5,4,4,14,3) |

The full start/end vectors and actual core cycles are generated in the receipt. Literal permutation/cut enumeration and a separate subset recurrence independently verify every displayed coordinate.

## 2. The shared exterior and exact path/cycle laws

For a nontrivial strong core D, put W_k(D,r) on D, a transitive block Z=(z_1,...,z_k), and w, with

    w -> Z -> D,       r -> w -> D-r,       k>=1.

The already proved cover bijection gives

    H(W_k(D,r)) = H(D)+2^k F_r(D),
    s_u(W_k(D,r)) = J_D(u,r) for u in D,
    s_w(W_k(D,r)) = H(D),
    s_z_i(W_k(D,r)) = 2^(k-i) F_r(D).

For completeness, the corresponding cycle law is

    hc(W_k(D,r)) = e_r(D).                           (1)

Here cycles are oriented and counted modulo cyclic rotation. Every cycle has a unique representation starting at w. Since D has no arcs to Z and only r beats w, that representation is w, all of Z in its forced transitive order, a Hamiltonian path of D ending at r, then w. Conversely every such core path gives exactly one cycle. This proves (1) without fitting finite data.

The wrapper is strong because every core vertex reaches r and then w, while w reaches Z and the core. Its outer vertex w always has sole in-neighbor r and s_r(W_k-w)=0. The exact zero-response mark therefore persists even though the Hamiltonian-cycle count changes.

## 3. All-clock ternary fan and its structural change

Using the core table, the three counts are

    N_k=9+11·2^k,
    H(W_k(A,r))=29+33·2^k=3N_k+2,
    H(W_k(B,r))=31+33·2^k=3N_k+4.

The parent has k+6 vertices and each child has k+7. The first four parent/children triples are

    (31;95,97), (53;161,163), (97;293,295), (185;557,559).

Construct either child by adding y to the physical parent wrapper. The new exterior incidences are Z->y and w->y; all previously present incidences remain unchanged. Switching between the children reverses only y->3. There are no intermediate arc-reversal states on this one-edge route.

Equation(1) gives cycle counts exactly 1,3,3, independently of k. At the fixed old vertex v=3, the start counts are exactly 1,3,4. Every start count of either child is at least three: this follows directly from the displayed J-vectors and the exterior start formulas. Thus neither child has a hidden unit-start port; the departure from the unit-start interface is genuine. The three-cycle count agrees with the universal lower bound s_u>=hc from cutting cycles.

This explicit construction illustrates why allowing a larger endpoint budget can make the physical operation simpler: the first fan used two births and two flips while preserving a unit port; this fan needs one birth and one flip while crossing from the one-cycle to the three-cycle class. The statements concern these two explicit examples, not an optimality theorem across all constructions.

## 4. Reproduction and limits

    python3 verify_cycle_fan.py
    python3 -O verify_cycle_fan.py

The verifier is standalone Python with no input data or third-party dependencies. It independently checks all three core profiles by literal permutations/cuts versus subset recurrences, the physical common-parent relation, and the exact single changed arc. It verifies 18 physical wrappers for k=1,...,6 and their complete start vectors. It directly enumerates every Hamiltonian cycle at k=1,2 for all three cores, six cycle checks in total. Both modes give identical stdout and receipt hash. The all-order statements follow from the path-cover and cycle bijections, not from finite extrapolation.

The exterior and two-bin calculus are inherited from the reviewed Session189/190 work, as cited by `PROOF.md`. No novelty against the literature or unexamined incoming work is asserted. This packet supplies a concrete all-clock fan that crosses the cycle strata. A continuation that returns both children to admissible parent states, a universal bounded-cycle grammar, and arbitrary odd-target realization remain open.
