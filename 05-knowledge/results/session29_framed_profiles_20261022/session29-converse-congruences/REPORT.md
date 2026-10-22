# Converse symmetry, odd-cycle overlap, and Hamiltonian-count congruences

Status: proof-complete deductions from the cited odd-cycle collection formula, with independent finite-exact checks; prepared for parent review. These are structural consequences, not a claim of literature novelty or a resolution of the global Hamiltonian-count spectrum.

## Inheritance and attribution

The closest established mechanism is `01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md`: if alpha_k counts collections of k pairwise vertex-disjoint directed odd cycles, then H(T)=sum_k 2^k alpha_k. Cycles are cyclically ordered directed cycles, not merely odd vertex subsets. Distinct cycles can have the same support.

The external primary source was checked directly: Irving and Omar, *Revisiting the Redei–Berge Symmetric Functions via Matrix Algebra*, Corollaries 19–20, https://arxiv.org/html/2412.10572v3 . Corollary 20 attributes the underlying formula to Grinberg–Stanley Theorem 1.39 and Lemma 6.5. Its displayed count is for the complement, including loops; loops do not affect Hamiltonian paths, and reversal identifies that count with H(T).

The nearest symmetry sidecar is `01-canon/theorems/THM-552-c3-parity-dichotomy-self-converse.md`, which proves the analogous triangle-parity distinction. The present argument extends the action to all directed odd cycles and to their disjoint packings. `THM-409-self-converse-perspective-flip-involution.md` carefully distinguishes arbitrary anti-automorphisms from involutions. We choose an involutory anti-automorphism rather than assuming every anti-automorphism is one. A targeted search of physically present canon files found no exact statement of the overlap congruence below; this is not an exhaustive novelty search.

The hostile example is a directed triangle followed by a sink: its H=3 at even order four, and it is not self-converse. Thus the symmetry hypothesis cannot be dropped.

## 1. Anti-involutions exist and all choices are conjugate

Every tournament automorphism group has odd order. Indeed, an automorphism of order two would interchange some pair of vertices, reversing their arc, which is impossible; Cauchy's theorem excludes even group order.

Let a be an anti-automorphism of a self-converse tournament. The automorphism a² has odd order m. Then sigma=a^m is an anti-automorphism, sigma²=1, and sigma is an involution when there are at least two vertices. An anti-involution fixes at most one vertex, since an arc between two fixed vertices would have to reverse itself. Thus it has no fixed vertices at even order and exactly one at odd order.

If sigma and tau are two anti-involutions, put g=tau sigma in Aut(T), of odd order m. We have sigma g sigma=g^(-1). Hence h=g^((m+1)/2) satisfies

h sigma h^(-1)=g^(m+1) sigma=g sigma=tau.

Therefore all anti-involutions are conjugate by an automorphism. Any cycle-orbit statistics defined below are independent of the choice, up to canonical isomorphism type, although no particular vertex permutation is distinguished.

## 2. Exact fixed-packing polynomial

Let Omega(T) be the conflict graph whose vertices are directed odd cycles, adjacent when their supports meet. The anti-involution acts on cycles by

rho(C)=sigma(reverse(C)).

This preserves directedness, cycle length, and support intersection. Let F be the number of rho-fixed cycles. Let A be the number of two-element cycle orbits. Partition these A orbits into:

- D admissible pairs: the two cycle supports are disjoint;
- I overlapping pairs: their supports meet.

Thus A=D+I. All rho-fixed odd cycles contain the unique fixed vertex of sigma, when it exists, so no two fixed cycles can lie in a disjoint packing. At even order F=0.

Define Q to have one vertex for each admissible two-cycle orbit, with adjacency when the unions of their supports intersect. For each fixed cycle C, let Q_C be the induced subgraph on admissible orbits whose union of supports avoids C. Write J_T(x)=I(Omega(T),x), the ordinary independence polynomial.

The polynomial counting rho-fixed packings by number of individual cycles is exactly

Phi_T(x)=I(Q,x²)+x sum_{rho(C)=C} I(Q_C,x²).

Proof: every fixed packing consists of whole two-cycle orbits and possibly one fixed cycle. A whole orbit is eligible exactly when internally disjoint; different orbit choices must have disjoint union supports. The converse construction is unique.

All remaining packings occur in two-element orbits. Consequently there exists a polynomial M_T(x) with nonnegative integer coefficients such that the exact identity

J_T(x)=Phi_T(x)+2 M_T(x)

holds coefficientwise. M_T counts the nonfixed packing orbits; this is not a claim that Q alone recovers J_T. At even order it simplifies to

J_T(x)=I(Q,x²)+2 M_T(x).

In particular every odd-indexed packing coefficient alpha_(2j+1) is even at even order, and alpha_(2j) is congruent modulo two to the coefficient of x^j in I(Q,x). This is a quotient-level strengthening of triangle parity, valid at every packing size.

## 3. The reflection-overlap congruence

For every self-converse tournament,

H(T) ≡ 1+2F+4I (mod 8).

Indeed OCF gives H≡1+2 alpha_1+4 alpha_2 modulo eight. Here alpha_1=F+2A. A rho-fixed two-cycle packing cannot contain two fixed cycles, since those intersect at the unique fixed vertex. It is therefore exactly one admissible two-cycle orbit, so alpha_2≡D modulo two. Hence

H≡1+2F+4(A+D)≡1+2F+4I (mod 8).

At even order this specializes to

H(T) ≡ 1 (mod 4),
H(T) ≡ 1+4I (mod 8).

Thus the next binary digit measures the parity of overlapping mirror-pairs of odd cycles, not just the total number of cycles. By Section 1, F and I are genuine invariants of the self-converse isomorphism class, not artifacts of the chosen anti-involution.

At odd order, F records the fixed-center channel. The formula does not say F is determined by triangles alone.

## 4. Consequences for multiplicities and merged metagraphs

Let a_n(h) count ordinary tournament isomorphism classes of order n with Hamiltonian count h. Converse preserves H and acts as an involution on these classes. At even n, every class with h≡3 modulo four is non-self-converse by Section 3. Therefore

a_n(h) is even whenever n is even and h≡3 (mod 4).

The same conclusion holds when restricted to strongly connected classes, or any converse-invariant class family. In the merged metagraph, each such node necessarily represents a pair of distinct ordinary classes; no singleton converse-fixed node can have that count. Its number of fixed-spine tilings is 2h/|Aut(T)|, with the usual action on ordered Hamiltonian paths free.

Equivalently, the ordinary-class count polynomial sum_T z^(H(T)), reduced coefficientwise modulo two, has support only at exponents 1 modulo four when n is even. This is a useful necessary census checksum. It does not classify wiring, identify which permitted counts occur, or prove the global {7,21} exclusions. In particular 21 is not excluded by this congruence.

## 5. Port parity and the Lorentzian seam

This independently explains the mod-four restriction in the even-palindrome construction of the Session29 transfer report.

For any nonempty tournament X with a marked portal a, let H=H(X), S count Hamiltonian paths starting at a, E count those ending at a, and

Q=sum_{nonempty proper U containing a} E(X[U],a) H(X-U),
R=E+Q.

Make a copy Y=X^op with corresponding portal b. Start with all cross arcs from Y to X, except replace b->a with a->b. The elementary two-block path decomposition gives

H(new tournament)=H²-S²+R².

The no-seam paths contribute H²-S²; paths using the seam have block patterns XY, YXY, XYX, YXYX and sum to (E+Q)². This construction is self-converse by swapping the two copies. Since every tournament has odd H, reduction modulo two gives R≡S modulo two for every X and every a. Consequently the square-minus expression is 1 modulo four, because H is odd and squares of equally-paritied integers agree modulo four.

This argument works for arbitrary X, not merely the eighteen five-vertex letters. It explains why a Lorentzian quadratic expression can be an odd nonnegative Hamiltonian count without making its ambient quadratic form positive definite.

## 6. Independent checks and limitations

`check.py` uses no producer or repository modules. It enumerates directed odd cycles rooted at their least vertex, lists disjoint cycle packings, builds the involution and quotient-polynomial coefficients, and compares OCF and all congruences against separate Hamiltonian-path subset DP.

The universe is 663 orientations: every tournament invariant under the fixed reversal anti-involution at orders one through six (1,2,4,16,64,512 cases), and 32 deterministic samples at each of orders seven and eight. This is not a count of isomorphism classes and not a complete census at orders seven or eight. For every tested portal, the independent subset counts also satisfy R≡S modulo two. A non-self-converse order-four H=3 hostile example confirms the hypothesis boundary.

Run `python3 check.py` or `python3 -O check.py`. Outputs match exactly; explicit RuntimeError checks cannot disappear under optimization. The proof, rather than the finite universe, supplies all-order scope.

The quotient Q retains which mirrored cycle-pairs can coexist and forgets all nonfixed packing orbits unless M_T is retained. It is a genuine orbit/overlap representation, not a complete invariant of tournaments or metagraph wiring.
