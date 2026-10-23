# An exterior representation of tournaments and exact component statistics

Status: PROVED elementary character and component identities; FINITE-EXACT independent permutation/graph controls below. This is a new connection in this session, not a claim of literature priority. “Bosonic” and “fermionic” refer only to symmetric versus exterior occupation rules in an explicitly defined algebra. No physical particles, spin, dynamics, or spin-statistics theorem are asserted.

## 1. Recovery, scope, and the missing sign

Inherited:

- Unique ordered SCC decomposition identifies the rational ordered-join algebra of tournament isomorphism types with the free associative algebra on strong tournament types. Session36 `CURL_COCENTER_AND_LIE_BOUNDARY.md` identifies its ordinary cocenter with SCC necklaces, but not with a compositional quotient.
- The old SCC/Lyndon codec, `05-knowledge/results/tournament_even_bijection_20261011c.md`, sections 2–4, identifies the ordinary free-Lie dimensions with the numbers e_n of connected automorphism-even graph types. The older note is present in the preserved `math` checkout; its link is absent from the materialized sparse Session35 publication checkout, so it was read without modifying that checkout.
- THM-4721, `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`, gives the signed Walsh representation and component lemma. Its “even” means edge-reversal sign trivial on every graph automorphism, NOT even vertex degrees, even edge number, or even vertex permutations.
- Royle–Praeger–Glasby–Freedman–Devillers, [Tournaments and Even Graphs are Equinumerous](https://arxiv.org/html/2204.01947v2), Theorem 1.1 and Corollary 3.2, provide the published graph/tournament count equality and the odd-order tournament-automorphism fact. The exterior proof below does not import the former equality: it proves its own different graph enumeration by a direct character calculation.

Canonical hostile: the path P3 is Royle-even but fails the exterior sign condition introduced below; the triangle C3 does the reverse. Thus replacing the Walsh squarefree monomials by exterior monomials is NOT a basis-preserving identification.

Live board: ordered SCC composition / Royle reversal signs / exterior edge determinants / component occupation / primitive curl periods. The exact source of the additional sign, rather than a verbal analogy, is the anchor.

## 2. A genuine exterior representation has the tournament character

Let V_n=Q^n with its usual permutation action of S_n. Put W_n=Λ²V_n, and let P_g be the permutation matrix of g. A basis of W_n is x_ij=e_i∧e_j for i<j. The group permutes these unordered-pair basis lines, with a minus sign exactly when the endpoints reverse order.

Let T_n be the set of labelled tournaments on [n]. Then, for EVERY n≥0 and EVERY g∈S_n,

    trace(g | Λ(W_n)) = det(I + Λ²P_g) = |Fix_Tn(g)|.          (1)

Here Λ(W_n) means the whole finite exterior algebra, with ordinary trace rather than supertrace.

Proof. A signed permutation cycle of length l and sign product s contributes 1−s(−t)^l to det(I+tR).

If g has odd order, every pair orbit has odd length and sign product +1: a negative sign product would give even order to the action on that orbit. Each orbit therefore contributes 2 at t=1. The result is 2^{c_2(g)}, exactly the number of choices of an invariant tournament orientation, one choice per pair orbit.

If g has even order, it has an even vertex cycle. For a vertex cycle of length two, its two vertices give a −1 eigenvector in Λ²V_n. For an even cycle of length l≥4, the l distinct cyclically oriented consecutive pairs

    e_1∧e_2, e_2∧e_3, …, e_l∧e_1

span an ordinary l-cycle permutation representation inside W_n. It has eigenvalue −1. Therefore det(I+Λ²P_g)=0. The same even vertex cycle makes a fixed tournament impossible: its half-turn reverses an antipodal pair. Both sides of (1) vanish. The empty and singleton cases are immediate. ∎

Consequences:

1. Λ(Λ²Q^n) and Q[T_n] are isomorphic S_n representations. This is an isomorphism of representations deduced from their characters in characteristic zero, not a specified map of basis sets or of algebras.
2. Writing a_n for the number of tournament isomorphism types,

       a_n = dim Λ(W_n)^{S_n}.

3. The precise exterior-degree refinement is the nonnegative polynomial

       P_n(q) = (1/n!) Σ_(g∈S_n) det(I+qΛ²P_g).

   It has P_n(1)=a_n. Its coefficients count a graph family specified next. They do NOT yet grade individual tournament types canonically.

The equality in (1) is stronger than an equality of invariant dimensions: it gives every subgroup-fixed invariant dimension after restriction, and each irreducible multiplicity. Choosing an intertwiner in each degree can be extended to an isomorphism of ordinary linear species under bijections. The proof does not specify a canonical species intertwiner or one compatible with ordered join, deletion, converse, or other extra operations; it does not give a natural set bijection.

## 3. Exactly which graphs survive exterior antisymmetrization?

For a graph G on [n], let C_1(G) be the rational vector space on oriented edges, with opposite orientations identified by a minus sign. For g∈Aut(G) define

    δ_G(g) = det(g | C_1(G)).

Call G exterior-admissible if δ_G(g)=+1 for every automorphism. This condition is independent of orientation choices and edge ordering. It is deliberately given a different name from Royle-even.

The line generated by the exterior monomial ∧_(ij∈E(G))x_ij has stabilizer character δ_G. More explicitly,

    δ_G(g) = ε_G(g) · sign(g acting on E(G)),                (2)

where ε_G is the inherited edge-reversal sign. The second factor is the exterior/Koszul sign from reordering the selected edge generators.

For a monomial orbit, its invariant subspace is one-dimensional if the stabilizer character is trivial and zero otherwise. One can see this by averaging the monomial: an automorphism with sign −1 makes the average equal to its own negative; if all stabilizer signs are +1, all contributions to any chosen orbit coefficient agree and the average is nonzero. Therefore

    [q^m]P_n(q)
      = number of exterior-admissible graph types with n vertices and m edges.   (3)

Equations (1)–(3) give an exterior graph partner equinumerous with tournaments, distinct from the inherited Royle-even partner.

The determinant test has a homological form. The equivariant exact chain sequence

    0 → H_1(G;Q) → C_1(G) → C_0(G) → H_0(G;Q) → 0

gives, by multiplicativity of determinants in exact sequences,

    δ_G(g)=det(g|H_1(G;Q)) · sign(g|V(G)) · sign(g|components).       (3a)

The last inverse determinant equals itself because it is a sign. For a connected graph, the last factor is 1. This is a real bridge to the orientation action on cycle space, not a claim that an arbitrary graph complex or physical homology theory has been constructed.

Small decisive tests:

- P3: its nontrivial automorphism has reversal sign +1 and swaps its two edges, so δ=−1. It is Royle-even but not exterior-admissible.
- C3: a reflection reverses an odd number of oriented edges and swaps two edge basis lines; the two minus signs cancel. Every automorphism has δ=+1. It is exterior-admissible but not Royle-even.
- The empty graph is admissible and contributes degree zero for every n.

Thus the polynomial at n=3 is P_3(q)=1+q³, not the Royle edge-count polynomial 1+q².

There is an even sharper duality boundary. Tournament converse fixes both tournament types at n=3, so its trace on the S3-invariant subspace is 2. Exterior-degree parity has trace 1−1=0 on the exterior invariant subspace. Thus there is NO S3×C2 intertwiner taking tournament converse to exterior parity. Character equality for S3 cannot silently be extended to this extra involution.

Three concrete infinite families make the statistics tangible:

- Complete graphs K_n, n≥1: the edge action is Λ² of the vertex action, with determinant sign(g)^(n−1). They are exterior-admissible exactly when n is odd. Those with n≡1 mod4 are even-edge/bosonic components; those with n≡3 mod4 are odd-edge/fermionic components.
- Paths P_n, n≥1: H_1=0 and reversal has vertex sign (−1)^floor(n/2). They are admissible exactly when n≡0 or1 mod4. The first nontrivial examples are fermionic P4 and bosonic P5.
- Simple cycles C_m, m≥3: an even m-cycle rotation already has δ=−1. For odd m, rotations have δ=+1 and reflections have δ=(−1)^((m+1)/2), using their action −1 on H_1. Thus exactly m≡3 mod4 survive, and every such cycle is an odd-edge/fermionic component.

## 4. The component theorem produces exact boson/fermion occupation

For a connected exterior-admissible graph C, assign parity |E(C)| mod 2. Let b_(n,m) count such types with even m and f_(n,m) count those with odd m.

Component theorem. A graph G is exterior-admissible if and only if:

(i) every connected component is exterior-admissible; and

(ii) every connected component type with an odd number of edges occurs at most once.

Proof. An internal bad automorphism extends by the identity, proving necessity of (i). A transposition exchanging two isomorphic m-edge components has determinant (−1)^m on their edge spaces, proving necessity of (ii). Conversely, every automorphism is a product of internal component automorphisms and permutations of identical component copies. Internal determinants are +1. A copy permutation τ of an m-edge component has determinant sign(τ)^m; this is +1 when m is even, and there is no nontrivial copy permutation when an odd-edge type appears at most once. ∎

Set

    P(z,q)=Σ_(n≥0)P_n(q)z^n.

Then the EXACT bigraded Euler product is

    P(z,q)
      = ∏_(C connected exterior-admissible, |E(C)| even)
          (1−z^{|V(C)|}q^{|E(C)|})^{-1}
        · ∏_(C connected exterior-admissible, |E(C)| odd)
          (1+z^{|V(C)|}q^{|E(C)|}).                         (4)

All products are formal: finitely many factors affect a specified vertex degree. This is the requested load-bearing statistics law. Even-edge component types have arbitrary occupation; odd-edge types have occupation zero or one. Two disjoint triangles are the first repeated fermionic-component hostile: each triangle individually survives, but their exchange has determinant −1, and their union does not survive.

There is an actual superalgebra, not only a numerical product. Take

    B = ⊕_(n≥0) Λ(W_n)_(S_n)

with coinvariants under vertex relabelling, and multiplication given by disjoint union of vertex blocks followed by exterior multiplication of edge generators. Internal relabellings make this well-defined. Exchanging two blocks contributes (−1)^(m_1 m_2), so B is supercommutative for edge parity. Graph component decomposition proves

    B ≅ Sym(B_conn,even) ⊗ Λ(B_conn,odd).                  (5)

The connected spaces consist of their one-dimensional surviving graph-orientation lines. Choosing one generator per line makes (5) a standard symmetric/exterior basis statement. No canonical sign choice is needed for the dimensions or occupation criterion.

Important boundary: B is not the ordered-join algebra of tournament types. The latter is noncommutative; B is supercommutative. Equation (1) does not supply a compatible multiplication or identify these connected graph factors with tournament SCCs.

## 5. Exact cancellation, connected balance, and a dyadic recurrence

For every n≥3 and every permutation g, W_n has an eigenvalue +1. If g has at least two vertex cycles, V_n has at least two independent eigenvectors of eigenvalue 1, whose wedge works. If g consists of one n-cycle, use two distinct inverse eigenvalues ζ and ζ^{-1}; these are distinct for n≥3 and their wedge eigenvalue is 1. Consequently

    det(I−Λ²P_g)=0      (n≥3).

At n=0,1,2 the invariant exterior algebra has only its degree-zero vector. Therefore

    P_n(−1)=0 (n≥3),     P(z,−1)=1+z+z².                 (6)

Thus, among all exterior-admissible graph types of every order n≥3, exactly half have an even edge number and half an odd edge number: each half has a_n/2 members. This is cancellation in an exterior-degree grading, not a new theorem that tournament automorphisms have odd/even parity classes.

Write b_n=Σ_(m even)b_(n,m) and f_n=Σ_(m odd)f_(n,m). From (4),

    A(z)=Σa_nz^n=∏(1−z^n)^(−b_n)(1+z^n)^(f_n),          (7)

    ∏(1−z^n)^(−(b_n−f_n))=1+z+z²=(1−z³)/(1−z).

The exponents of an integral formal Euler product are unique, by comparison at the first differing degree. Hence the connected balance law is

    b_n−f_n = 1_(n=1)−1_(n=3).                           (8)

Recover e_n from the inherited Royle-even/Lyndon Euler product A(z)=∏(1−z^n)^(−e_n). Replacing 1+z^n by (1−z^{2n})/(1−z^n) in (7) gives

    e_n = b_n+f_n−1_(2|n) f_(n/2),

or, using (8),

    2f_n = e_n + 1_(2|n)f_(n/2) − 1_(n=1) + 1_(n=3).    (9)

This is a new-to-this-session exact connection to the older even/odd counting mechanism: ordinary connected-even counts determine the exterior boson/fermion component counts by repeated halving. It is not a correspondence between individual connected graph types. The e_n are neither all bosonic nor all fermionic: the doubling correction matters.

Values (n=1,…,8):

    e_n: 1, 0, 1, 2, 8, 43, 398, 6413
    b_n: 1, 0, 0, 1, 4, 22, 199, 3207
    f_n: 0, 0, 1, 1, 4, 22, 199, 3207.

The first six b/f entries are direct complete graph-orbit controls. Entries seven/eight here follow rigorously from (9) and the inherited exact e_n values, not from an exhaustive graph-atlas run in this file.

For odd m≥5 and k≥0, iterating (9) gives the explicit dyadic-ray formula

    f_(2^k m) = Σ_(j=0)^k e_(2^{k−j}m)/2^{j+1}.

The m=1 and m=3 rays carry the exceptional terms in (9). Integral graph counts prove integrality of these dyadic combinations. This is arithmetic of counting sequences; it has no asserted connection to rational versus irrational real numbers.

## 6. A further check: complement symmetry is parity-sensitive

The determinant character of W_n is

    det(Λ²P_g)=sign(g)^(n−1).

For odd n it is trivial. Exterior duality therefore gives

    P_n(q)=q^{binom(n,2)}P_n(q^{-1})  (n odd).

Equivalently, graph complement preserves exterior admissibility at odd vertex order. At even order it introduces the vertex-sign twist and need not preserve admissibility. This is a separate parity from edge parity and from the Royle reversal sign. The full polynomials in the computation receipt verify the odd-order palindromicity through seven vertices.

## 7. Controls, reproduction, and the connection ledger

Run from the session workspace:

    python3 session37-parity/check_exterior_statistics.py
    python3 -O session37-parity/check_exterior_statistics.py

The checker uses only Python's standard library and explicit exceptions, not removable assert statements. It checks:

- all 5,914 permutations across n=0,…,7: exterior/tournament character equality and the supertrace cancellation;
- every graph isomorphism type through n=6, constructed from all labelled edge masks and all relabellings: 1,1,2,4,11,34,156 types;
- the complete exterior-degree polynomial by two independent methods, signed permutation determinants and graph-stabilizer testing;
- the bigraded connected-component product, not only its q=1 specialization;
- Royle-even counts, their distinction from exterior-admissible types, connected balances, and the recurrence (9);
- the explicit P3, C3, and two-triangle hostiles.

Output: `exterior_statistics.json`; optimized output is required to match.

Independent audit: `session37-statistics-audit/graph_signs.json` enumerates all graph types through n=7 and independently confirms the b/f rows; `session37-statistics-audit/characters.json` checks all 272 permutation conjugacy types through n=12. The independent proof audit accepted (4)–(9), the coinvariant-algebra construction, and the complement statement, and sharpened the species/converse boundary above. These independent files are not inputs to this lane's checker.

Connection ledger:

Source → target: labelled tournament permutation representation → full exterior representation Λ(Λ²Q^n). Map: equality of characters gives an abstract linear intertwiner, not an explicit basis map. Preserved: all S_n character data and resulting invariant dimension. Lost/not supplied: individual tournament identity, SCC order, Hamiltonian-path count, multiplication, and exterior-degree interpretation of a particular tournament. Sidecar needed: an explicit, operation-tested intertwiner or a deliberately chosen type codec. Cheapest hostile: P3 versus C3 rules out the naive Walsh-to-wedge monomial map; disjoint triangles verify the occupation rule.

Source → target: exterior graph superalgebra → its Hilbert series. Map: vertex/edge grading. Preserved: multiplicities by two gradings and the occupation rule. Lost: graph identity, component incidence, automorphism group beyond its determinant test. Sidecar: actual connected graph-orientation lines. No relation to the physical spin-statistics theorem follows.
