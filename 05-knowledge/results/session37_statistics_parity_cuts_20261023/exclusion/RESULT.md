# Tournament exclusion, mixed statistics, and the corrected spectral boundary

Status: PROVED elementary algebra for Sections 1–4; FINITE-EXACT with a complete finite universe for Section 5; explicit counterexample independently certified by permutation and cycle enumeration. Repository-relative findings, not a literature-priority claim. This note does not claim a quantum model, a polynomial-time Hamiltonian algorithm, or a connection between arithmetic transcendence and Bose/Fermi statistics.

## Executive result

The natural unsigned Hamiltonian carrier is a **commutative square-zero algebra**, not an exterior algebra on one odd generator per vertex. It has a faithful even-sector exterior realization using **two** odd generators per vertex. The familiar determinant/permanent distinction can be expressed exactly as a two-species cycle filter: equal formal fermionic and bosonic multiplicities cancel even directed cycles and double odd directed cycles.

But discarding vertex labels before combining these sectors destroys compatibility data. A complete census and an explicit hostile sharpen the old repository boundary:

- The pair `(det(zI−A), per(zI+A))` determines `H(T)` on all tournaments of each fixed order `n <= 8`.
- It fails at `n=9`, and hence at every `n>=9` by adding sources.
- The old THM-506 conclusion “failure from eight” confused failure of **affine reconstruction** with failure of **arbitrary functional determination**. An affine-rank test does not imply a collision in the actual invariant fibers.

The `n<=8` assertion is finite-exact, not an all-order structural proof. The complete universe at eight has `2^21=2,097,152` fixed-Hamiltonian-path presentations and `3,217` distinct determinant/permanent keys, with no key splitting H. Every eight-vertex tournament is isomorphic to one of these presentations by Rédei's theorem. Negative conclusions at smaller orders follow by source extension, rather than by extrapolation.

## 0. Inheritance, source boundaries, and concept board

Closest proved mechanisms:

- Session35 `recursions/CONTEXT-AND-PORTS.md`: the square-zero path resolvent and marked determinant quotient, with the explicit six-vertex adjacency-cospectral `H=29/25` hostile.
- `THM-002-ocf.md`: the established odd-cycle formula, not a new formula here.
- `THM-506-permanental-companion.md`: signed and unsigned cycle-packing coefficients, but its rank-only non-determination inference needs correction.
- `THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`: reversal cancels odd-edge signed paths; this is different from odd-vertex directed cycles in OCF.
- Session36 `context/RESULT.md`: the path-cover and substitution carrier retains more than a single scalar and is interface-relative.

Primary source checked: John Irving and Mohamed Omar, *Revisiting the Rédei–Berge Symmetric Functions via Matrix Algebra*, arXiv:2412.10572v3, https://arxiv.org/html/2412.10572v3 . Equation (2) is the squarefree MacMahon permanent identity; Lemma 1 is the determinant quotient for walks; Proposition 2 is the determinant/permanent Hamiltonian formula. Therefore those identities, including their tournament specialization, are inherited, not new priority claims. The related OCF source is Grinberg and Stanley, https://arxiv.org/abs/2307.05569 .

Live board:

1. Hard exclusion: commuting nilpotents; retain labels and disjointness.
2. Fermionic alternation: exterior signs; retain a choice of top-degree orientation if using one generator per label.
3. Bosonic symmetric powers: arbitrary occupation before squarefree projection; no analytic convergence or positive-state assumption.
4. Cycle parity: parity of each cycle length is distinct from parity of the number of cycles.
5. Scalar det/per: retains occupied-size and cycle-count-parity marginals, but can lose cycle-length and overlap data.
6. Linear forests: run reversal uses edge-length parity, not a particle-statistics assignment.

## 1. Exact carrier and genuine exterior realization

Work over

    R_n = Q[x_1,...,x_n]/(x_1^2,...,x_n^2),  X=diag(x_i).

Its basis is `x_S` indexed by vertex subsets. Variables commute, and `x_S x_U=0` exactly when `S` and `U` intersect. For the 0/1 adjacency matrix A,

    Phi_A = 1 + 1^T X(I−AX)^(-1)1

has `[x_S]Phi_A=H(T[S])`, including the empty-subset value one. The inverse is finite because the positive-degree ideal is nilpotent. Walks with a repeated vertex vanish; surviving paths retain **positive** multiplicity. This is hard exclusion without antisymmetrization.

Let Lambda be the exterior algebra on the ordered odd generators

    theta_1, eta_1, theta_2, eta_2, ..., theta_n, eta_n.

The map `iota(x_i)=theta_i eta_i` is a well-defined injective algebra homomorphism `R_n -> Lambda_even`. Indeed each image squares to zero, even-degree images commute, and distinct subsets map to distinct exterior basis elements. The top coefficient in the displayed interleaved order satisfies

    [theta_1 eta_1 ... theta_n eta_n] iota(Phi_A) = H(T).

Swapping two vertex labels exchanges two degree-two blocks and contributes sign `(-1)^(2*2)=+1`; the paired realization retains unsigned relabeling invariance. If one instead chooses the grouped convention `theta_1...theta_n eta_1...eta_n`, a global factor `(-1)^(n(n−1)/2)` must be included.

This is an algebraic realization of commuting nilpotents **inside the even sector**. It does not turn x_i into odd fermionic variables or assert a physical system of paired fermions.

### Optimal number of odd generators

For n>=1, the paired realization is dimensionally optimal in a precise sense.

**Even-sector minimality theorem.** A faithful unital algebra embedding `R_n -> Lambda(Q^m)_even` exists if and only if `m>=2n`.

Proof: any nilpotent exterior element has zero scalar part. In the even sector every image of x_i therefore has degree at least2. The nonzero element `x_1...x_n` must have nonzero image, of degree at least2n. Since exterior degrees above m vanish, `m>=2n`. Disjoint pairs provide sufficiency. No gradedness assumption on the embedding is needed, only that its image lies in the even sector. QED.

**Ungraded minimality theorem.** If the image may use odd as well as even terms, a faithful unital embedding `R_n -> Lambda(Q^m)` exists if and only if `m>=2n−1`.

Proof: write l_i for the degree-one part of the image of x_i. Because the images commute, the degree-two part of their commutator is `2 l_i wedge l_j=0`. Over Q this means all nonzero l_i span at most one common line. Thus every term of the n-fold product involving two degree-one parts vanishes. Each surviving term has degree at least `1+2(n−1)=2n−1`; injectivity forces this product to be nonzero. For sufficiency choose `x_1 -> theta_1` and `x_i -> theta_(2i−2) theta_(2i−1)` for i>=2. These commute, square to zero, and their subset products are linearly independent. This construction singles out one odd generator and is not an all-even, label-symmetric realization. QED.

The distinction is essential: the lower bound2n is for the **even-sector** bridge actually used in the determinant/positive-count interpretation. These bounds concern an algebra representation, not computational complexity or physical particle numbers.

### Minimal hostile to one-generator “fermionization”

Replacing each path monomial by the ordered product `theta_(v1)...theta_(vk)` in an exterior algebra on just n odd generators gives the signed Hamiltonian sum

    H_alt(T) = sum_(Hamiltonian permutations pi) sign(pi),

not H. It transforms by the sign of a relabeling, whereas H is invariant.

For the four-vertex tournament with outneighbor bitmasks `[8,1,3,6]`, equivalently upper-pair-lex mask4, the five Hamiltonian paths are

    (0,3,2,1), (1,0,3,2), (2,0,3,1), (2,1,0,3), (3,2,1,0).

Their signed sum is −1. Thus even taking the absolute value gives1 rather than5. This is minimal for the absolute-value failure: for n<=3 the transitive tournaments have one path and the cyclic triangle has three paths of the same permutation parity.

## 2. Precise mixed-statistics identity

In R_n the following finite identity is the squarefree MacMahon identity:

    det(I−AX)^(-1) = per(I+AX).

It is special to squarefree projection. It is not an identity of ordinary univariate polynomials after identifying every x_i with a common scalar z.

Introduce formal parameters f,b. Nilpotence defines arbitrary powers by finite binomial series. Put

    Z_A(f,b;X) = prod_i(1+x_i) det(I+AX)^f det(I−AX)^(-b).

The finite logarithm and trace identity give

    log Z_A = sum_i x_i
              + sum_(directed simple cycles C)
                  ( b + (−1)^(|C|−1) f ) x_(V(C)).

Here a simple cycle is counted modulo cyclic rotation. The rooted trace contributes |C| copies, canceled by the logarithm's denominator. There are no loops or two-cycles in a tournament.

Consequently

    [x_V]Z_A(f,b;X)
      = sum_(vertex-disjoint cycle collections L)
          (b+f)^(o(L)) (b−f)^(e(L)),                 (1)

where o(L), e(L) count the odd- and even-length cycles. Remaining vertices are singletons supplied by `prod(1+x_i)`.

At f=b=1 the even cycles cancel and odd cycles have weight2. The marked tournament determinant lemma gives

    Z_A(1,1;X)=Phi_A,
    H(T)=[x_V]Z_A(1,1;X).

This is the established OCF in explicit mixed-statistics coordinates, not a new OCF claim.

The terminology is algebraically precise only as follows. The determinant generates exterior-power traces, whereas its inverse generates symmetric-power traces, in the usual formal trace identities. The matrix AX is not a Hermitian Hamiltonian, and nilpotent vertex activities are not physical probabilities. `per(I+AX)` is the **squarefree projection** of the bosonic series, not the unrestricted bosonic partition function. For integer nonnegative f,b, the factors may be read as species-colored cycles with signs in the exterior sector. On the diagonal f=b, even cycles cancel between the species and odd cycles reinforce; the coefficient2 does not count a two-dimensional truth object or a universal quantum state.

The coefficient form also recovers the inherited determinant/permanent convolution:

    H(T) = sum_(S subset V) det((I+A)[S]) per(A[V\S]).

Transpose invariance identifies `(I+A)[S]` with the tournament complement's principal submatrix in Irving–Omar Proposition2.

## 3. Even-cycle response is an exact new sidecar, not a new OCF

Define independent variables u=b+f and v=b−f, and

    P_T(u,v)=[x_V] Z_A((u−v)/2,(u+v)/2;X).

By (1), P is the full cycle-packing polynomial graded by numbers of odd and even cycles. In particular `P_T(2,0)=H(T)`. Differentiating at v=0 gives the exact even-cycle insertion identity

    (d/dv) P_T(2,v)|_(v=0)
       = sum_(directed simple even cycles C) H(T[V\V(C)]).       (2)

Proof: differentiation selects one even cycle and requires all remaining cycles to be odd; its complement contributes OCF. Higher derivatives divided by j! select j pairwise-disjoint even cycles and the Hamiltonian count of their complement. This is an elementary consequence of the master packing polynomial, recorded as an explicit response identity, not a claim of literature novelty.

For the order9 collision below, the first-order response agrees at 1567, but the second derivative divided by2 differs: it is 96 versus 95, the disjoint-four-cycle pair count. Thus the hierarchy detects a coordinate discarded by the unmarked characteristic/permanental pair; not every response level detects every collision.

## 4. The linear-forest parity is not the same parity

For n>=1 and skew adjacency S=A−A^T, the inherited forest identity is

    2^(n−1) H(T) = sum_(directed linear forests F) (n−|E(F)|)! prod_(ij in F) S_ij.

Reversing one path with an odd number of edges changes its weight by−1 and is a sign-reversing involution. Thus only forests all of whose path components have an even number of edges survive. Such a path has an **odd number of vertices**. The empty forest remains valid.

This shares a reversal-based parity filter with OCF, but the objects, signs, and operations differ:

- OCF: actual directed cycles of A, odd vertex/cycle length, unsigned weight2.
- Forest expansion: paths in the complete digraph weighted by S, even **edge** length, cancellation among opposite orientations.
- Raw exterior algebra: permutation signs under exchanging individual odd generators.

Calling all three “odd” without naming the carrier loses a sign or shifts parity by one.

## 5. Corrected determinant/permanent boundary

### 5.1 Complete finite universe through eight

The primary exhaustive program fixes `0 -> 1 -> ... -> 7` and varies the remaining21 arcs. Every eight-vertex tournament admits a Hamiltonian path and hence has a relabeling in this universe. The programs calculate:

- per(zI+A) by row-assignment subset DP over integer coefficient arrays;
- det(zI+A) by the same assignment order with the exact permutation-inversion sign;
- H by independent end-vertex subset DP.

The dictionary is keyed by both complete polynomials. All `2,097,152` presentations were examined; there are `3,217` keys, each with a unique H. The positive controls include the transitive presentation and nontransitive paths. There is no isomorphism-class count claim: presentations are not classes.

The usual characteristic polynomial is equivalent data, since

    det(zI−A)=(-1)^n det(−zI+A).

If a collision existed at smaller n, adjoining `8−n` sources would multiply both polynomials by `z^(8−n)` and preserve H, contradicting the eight-vertex census. Thus the result covers all fixed orders n<=8, including the empty tournament convention.

### 5.2 Explicit failure at nine

Use upper-pair-lex bits `(0,1),(0,2),...,(7,8)`, with bit1 meaning i->j for i<j. The two masks and equivalent outneighbor rows are

    T: mask29759988540,
       [120,397,337,228,362,198,386,21,169], H=2843;
    U: mask18875161993,
       [274,116,193,135,300,77,281,115,174], H=2841.

Both have, in ascending powers z^0 through z^9,

    per(zI+A): [384,692,629,380,154,68,28,0,0,1],
    det(zI+A): [20,−82,167,−192,154,−68,28,0,0,1].

Directly enumerating all9! permutations gives those two polynomials and the Hamiltonian counts. A second route enumerates rooted simple cycles, then all disjoint cycle packings; it agrees exactly with the permutation outputs and OCF.

The common odd-cycle totals below nine are c3=28,c5=154,c7=398, and both have D33=94. Their differing odd-cycle packing coordinates are

    T: c9=191, D35=209, T333=11;
    U: c9=192, D35=210, T333=10.

The H difference is therefore

    2*(−1)+4*(−1)+8*(+1)=2.

The invisible trades are compensated by D44=96/95 and by the constant `c9+T333=202`. Thus the scalar fermionic/bosonic summaries keep cycle-count parity and total occupied size but not individual cycle lengths and compatibility.

There is an especially concise exact kernel witness. Let `C_T(z;u,v)` sum all vertex-disjoint cycle collections with weight `z^(number of occupied vertices) u^(odd cycles) v^(even cycles)`. Then for these actual tournaments,

    C_T(z;u,v)−C_U(z;u,v)
       = z^9 (u^3−u) + z^8 (v^2−u^2).             (3)

The unmarked permanent specialization is `(u,v)=(1,1)` and the `det(I+zA)` specialization is `(u,v)=(1,−1)`. Equation(3) vanishes identically in z at both points. The Hamiltonian specialization is `(z,u,v)=(1,2,0)`, where it equals2. This supplies a **realized** polynomial kernel element, not merely a formal carrier-vector or rank obstruction. No unrealized cycle-count trade is being substituted for an actual tournament.

Adding a new source makes the adjacency block upper triangular. Both det(zI+A) and per(zI+A) gain a factor z: a non-fixed source cannot belong to a cycle cover. Every Hamiltonian path starts at the source, so H is unchanged. Repeating this gives failure at every n>=9.

### 5.3 What the rank computation did and did not prove

The independent audit supplies a three-tournament certificate of genuine affine failure **without** invariant-fiber failure. The order-eight masks172254133,172256153,172238829 have the same characteristic polynomial and H values177,197,209. Their permanent coefficient vectors differ only in the three lowest entries, which are respectively `(12,31,40)`, `(16,39,44)`, `(18,43,46)`. Thus the last-minus-first permanent difference is3/2 times the middle-minus-first difference, but32 is not3/2 times20 for H. The invariant vectors are distinct, so this proves no affine reconstruction, not non-determination. These three examples were independently checked by all8! permutations in `session37-inheritance-audit/n8_affine_obstruction_direct.json`.

A rank increase after appending H to a chosen invariant-coordinate matrix disproves an affine formula for H in those coordinates. It does not disprove some nonlinear function on the finite set of realized invariant values. The eight-vertex census is an exact example of this distinction. We make no claim of a simple closed reconstruction formula at eight or efficient reconstruction at variable n.

## 6. Reproduction and scope

From the workspace root:

    g++ -O3 -std=c++17 session37-exclusion/exhaustive_n8.cpp -o /tmp/session37_n8
    /tmp/session37_n8 2097152
    python3 session37-exclusion/verify_witness.py
    python3 session37-exclusion/verify_grassmann.py

Receipts: `exhaustive_n8.out`, `witness_verification.json`, `witness_verification.out`, and `grassmann_verification.json`. The exterior test checks every basis product in both explicit embeddings for n=1..7, the signed-count hostile/minimality for n<=4, the exact realized kernel polynomial, and the first three even-cycle responses. `search_n9.cpp` records the deterministic discovery search with seed370021; `export_n8.cpp` and `n8_classes.csv` provide one representative and H per invariant key. The csv columns are mask,H, nine ascending permanent coefficients, nine ascending determinant coefficients. File ordering reflects the C++ hash container and has no mathematical meaning.

Independent audits confirmed the nine-vertex witness by their own 9! permutation enumerations. `session37-inheritance-audit/exhaustive_n8_cycle_audit.cpp` independently checks the full 2^21 eight-vertex universe using rooted-cycle dynamic programming,651 disjoint-cycle-pair placements, and OCF rather than the primary row-assignment/path algorithm. It confirms 3217 keys and no H collisions. The audit receipt also reports H range 1..661 and total presentation-H sum 815136768. No Git mutation was made by this worker.

No claim here identifies parity, rationality, Dedekind cuts, or transcendence. All finite tournament matrix eigenvalues are algebraic because their characteristic polynomials have integer coefficients. A transcendental analytic limit, if considered separately, requires a specified sequence, convergence proof, and number-theoretic argument; determinant/permanent vocabulary supplies none of these.
