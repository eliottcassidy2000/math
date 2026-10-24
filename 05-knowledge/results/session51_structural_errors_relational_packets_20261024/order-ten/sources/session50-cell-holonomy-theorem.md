# Equitable-cell holonomy and structural tournament unpointing corrections

24 October 2026. **PROVED:** the equitable-cell restrictions, the signed
triangle-cell automorphism theorem, the strong-order-eight vanishing theorem,
the SCC error-series law, and the all-r matching-family count. **FINITE-EXACT:**
the exhaustive regular-seven input and explicit finite controls.
**OPEN:** an all-order scalar recurrence for A000568, and exact refinement-error
prediction for unrestricted higher-order strong tournaments. No literature
priority or asymptotic speed claim is made.

## 1. Inheritance and the exact question

The frozen Session50 growth packet proves the parent-only unpointing trace
identity and defines stable directed color-refinement upper bounds. It also
certifies the first regular-seven hostile and total discrepancies 14 and 28
at orders seven and eight. This packet does not claim those finite observations
anew. It explains the order-eight value structurally and extracts the next
signed compatibility contribution to the refinement error.

The closest inherited mechanism is Session39's `AFFINE_CONNECTIONS.md`,
which already proves the complete affine-F3 connection decoder, the signed
linear kernel of size3^kappa, the corresponding rotational automorphism
group, and the 216-presentation/16-class three-triangle example. Session40's
`marked-ports/PATH_TRANSPORT_AND_FOREST_RULE.md` already states the odd-cell
small-core classification, forest/joint-holonomy decoder and coprime crossing
lemma. Session40's `PRIME_CELL_BURNSIDE_FORMULAS.md` already treats two
matching-coupled cyclic regular prime cells, exact class counts and their
automorphism denominator. The repository also uses directed refinement in
THM-4131's exact canonicalizer, without claiming cells equal vertex orbits.

Those mechanisms and family counts are inherited, not new here. The new
extraction is their exact contribution to root-refinement error, its additive
SCC transfer, structural vanishing at strong order8, the resulting corrected
upper-trace prediction, and family contributions to later error coefficients.
We reproduce the necessary inherited proofs for a self-contained correction
argument. Source bytes, paths, commits and hashes are archived in `pinned.txt`.

The anchor is the A000568 root-transport correction; the niche is an exact
small-cell symmetry classifier; the wildcard is the signed-graph phase-space
count. The live concepts are pointed root transport, equitable cell sizes,
matching ownership, cycle signs, and ordered strong-component composition.
The missing coordinate is compatible rotations, not another cell-size scalar.

Throughout, T is a finite ordinary tournament and Aut(T) its ordinary
orientation-preserving automorphism group. Stable directed color refinement
starts with one color and repeatedly records the old color and the number of
outneighbors in each current color. In a tournament the corresponding
inneighbor counts are determined by these counts and cell sizes. Cells are
individual invariant color classes, not an unlabelled block system that
an automorphism may permute.

For each Aut(T)-vertex orbit O let C(O) be its stable color cell. Define

    e(T) = sum_(O in V(T)/Aut(T)) (|C(O)| - |O|).

Each term is nonnegative. If D is the growth packet's root-changing matrix,
then e(T) is exactly this child's contribution to the excess of the
refinement upper trace over tr(D). This is not the different symmetry defect
Delta=sum_p(D[p,p]-1).

## 2. Odd cells and coprime crossing rigidity

**Lemma 1 (inherited equitable-cell restrictions).** Every stable cell C induces a regular tournament, so |C| is odd.
For distinct stable cells A,B of sizes a,b, let r be each A-vertex's number
of outgoing neighbors in B, and s each B-vertex's number of incoming neighbors
in A. Then

    a r = b s.

If gcd(a,b)=1, the crossing is uniform: all arcs go A to B or all go B to A.

**Proof.** Equitability applied to C itself makes the internal outdegree a
constant d. Counting its internal arcs gives |C|d=|C|(|C|-1)/2, so
2d=|C|-1. For distinct cells, outgoing equitability and the tournament
complement give constant r and s; count the same A-to-B arcs from either
side. Coprimality makes b divide r, and 0<=r<=b, leaving only r=0 or r=b.
End proof.

In particular every singleton has a uniform relation to every other cell.
This statement does not require an observer, chosen Hamiltonian path, or
a prior automorphism computation.

## 3. Exact signed holonomy for cells of sizes one and three

Assume every stable cell has size one or three. Every size-three cell is a
directed triangle. Choose a cyclic coordinate x in F3 on each triangle i,
with its arcs x to x+1. Its automorphisms are precisely rotations
x to x+k_i, k_i in F3.

For two triangle cells i,j, equitability makes their 3 by 3 crossing matrix
have constant row and column sum d in {0,1,2,3}. For d=0 or 3, the crossing
is uniform and imposes no restriction. For d=1 it is a permutation matrix;
for d=2 its complement is a permutation matrix. In either case let sigma_ij
be the matching permutation from i to j. Every permutation of F3 has the
unique affine form

    sigma_ij(x) = epsilon_ij x + c_ij,
    epsilon_ij in {+1,-1}, c_ij in F3.

Thus a nonuniform crossing is preserved by rotations exactly when

    k_j = epsilon_ij k_i.                              (1)

Make a signed graph H whose vertices are the triangle cells and whose edges
are their nonuniform crossings, with sign epsilon_ij. A connected component
is balanced if every cycle has product of edge signs +1; an isolated vertex
is balanced. Let b be the number of balanced connected components, and let u
be the total number of triangle cells in unbalanced connected components.

**Theorem 2 (inherited signed kernel; root-error corollary).** There is an exact faithful identification

    Aut(T) = {k in F3^(triangle cells): (1) on every edge}
           isomorphic to C3^b.

A triangle in a balanced component is one automorphism orbit. Every vertex
of a triangle in an unbalanced component is fixed by every automorphism.
Every singleton is fixed. Consequently

    e(T) = 6u.                                         (2)

**Proof.** Every automorphism fixes each stable cell setwise. On a triangle
it must be a rotation, and it fixes singletons. Uniform crossings and
singleton crossings impose no condition. On a nonuniform crossing, preserving
the permutation matrix or its complement is exactly
sigma_ij(x+k_i)=sigma_ij(x)+k_j, which is (1). These conditions are also
sufficient: every internal and every crossing arc is then preserved. No
additional global compatibility condition remains.

On a connected signed graph, choose a root and propagate its value along
paths. In a balanced component the propagated sign is path-independent and
the root value is arbitrary in F3. In an unbalanced component, a negative
cycle gives k=-k, forcing k=0, and hence all values vanish. This proves the
group structure. A nonzero balanced-component rotation acts transitively on
each of its triangle cells, while unbalanced components act trivially.
A balanced triangle contributes 3-3=0 to e; an unbalanced triangle has three
singleton orbits and contributes 3(3-1)=6. End proof.

**Coordinate invariance.** Changing the origin of an oriented cyclic
coordinate changes c_ij but not epsilon_ij. The sign is intrinsically whether
the matching preserves or reverses the two directed triangles. Reversing an
edge replaces sigma by sigma inverse and leaves its sign unchanged. More
generally, if one chooses either generator of each C3 independently, write
the new generator as the old one to power tau_i in {+1,-1}. The new edge
sign is tau_j epsilon_ij tau_i; every cycle product is unchanged. Thus
balancedness, b, u and the actual permutation group are coordinate-independent.

The signed kernel and its rotational automorphism interpretation are
already in Session39. Its specialization to individually invariant stable
cells makes that group the full ordinary automorphism group; formula(2)
extracts the exact error relevant to unpointing. Crucially, translation
phases cancel from the automorphism constraints. Full class/wiring decoding
still needs the joint affine phase holonomies, as Session39 proves. Thus
exact isotropy and root-orbit sizes need less information than complete
class identity. Even two differently chiral triangle crossings can have
different class identities and the same exact root-orbit sizes.

The theorem is special to sizes one and three. A larger odd cell need not be
vertex-transitive, and permutations of larger prime sets need not be affine.
The regular-seven hostile is already outside the theorem's domain.

## 4. Structural exactness through six and at strong order eight

**Lemma 3.** A regular tournament on five vertices is the unique cyclic
regular tournament C5 and is vertex-transitive.

**Proof.** Fix v. Name its outneighbors a,b with a to b and its inneighbors
c,d with c to d. Vertex b has lost to v and a, so b beats both c,d. Vertex c
already beats v,d, so a beats c. Vertex a then already beats b,c, so d beats a.
All remaining arcs were specified. In the cyclic order (v,a,b,c,d), each
vertex beats its next two vertices, proving the claim. End proof.

**Theorem 4.** Every tournament of order at most six has e(T)=0. Every
nonregular tournament of order seven has e(T)=0. Every strongly connected
tournament of order eight has e(T)=0.

**Proof.** Stable cell sizes are odd. At order at most eight, if there is no
seven-cell, either every cell has size one or three, with at most two
triangle cells, or there is one five-cell and all other cells have sizes
one or three. With at most two triangle cells the signed interaction graph
has no cycle, so Theorem2 gives exact orbits. A five-cell is internally
vertex-transitive by Lemma3. Its crossing to every other cell is uniform
by Lemma1 because the remaining sizes are one or three. Its automorphisms
therefore extend independently, and any remaining triangle is transitive.
Thus its cells too are the exact vertex orbits.

A seven-cell cannot occur below order seven; at order seven it means one
regular cell. At order eight it leaves a singleton, and their crossing is
uniform by Lemma1. That makes the tournament nonstrong. These cases exhaust
all odd partitions; an eight-cell is impossible. End proof.

This proof is independent of an order-eight census. It also proves that a
size-one/three-cell refinement failure needs at least three triangle cells,
hence order at least nine. It does not assert that all strong tournaments
below nine are exact: the regular-seven hostile is an explicit exception.

## 5. Strong-component transfer of the error

Let C1,...,Ck be the strongly connected components of T in their unique
order, with all arcs Ci to Cj for i<j. Put L_i=sum_(j>i)|Cj|. For v in Ci,

    outdegree_T(v) = L_i + outdegree_Ci(v).

The possible outdegree intervals [L_i,L_i+|Ci|-1] are disjoint, so the first
degree refinement never mixes two components. Inside a component, outside
counts are constant across all its vertices. Induction through the refinement
steps shows that its global color partition is exactly its standalone
partition, up to renaming colors.

Every automorphism fixes each SCC individually because it preserves the
unique total order of components. Conversely component automorphisms extend
independently across uniform crossings. Therefore Aut(T) is the direct
product of the component automorphism groups, with their own vertex orbits,
and

    e(T) = sum_i e(Ci).                                 (3)

Let a_n be the ordinary tournament class count, including a_0=1; let s_n
be the strongly connected class count, counting a singleton as strong. Define
formal ordinary generating series

    A(z)=sum_(n>=0) a_n z^n,
    S(z)=sum_(n>=1) s_n z^n,
    E(z)=sum_[T] e(T) z^|T|,
    E_s(z)=sum_[T strong] e(T) z^|T|.

Unique ordered SCC decomposition gives A=1/(1-S). Marking the component that
contributes to the additive sum (3) gives the exact identity

    E(z) = A(z)^2 E_s(z).                               (4)

Indeed every contribution is one marked strong component, preceded and
followed by arbitrary ordered lists of strong components. This counts
ordinary classes and an additive weight, not labelled objects or a Burnside
fixed-tournament sum. The series are formal; no convergence is asserted.

## 6. The order-seven finite input and order-eight structural correction

Only regular order-seven tournaments can contribute error at that order by
Theorem4. `check.py` enumerates all 2,640 labelled regular tournaments by
assigning every pair orientation and pruning only impossible residual degree
requirements. It partitions them into full permutation orbits, and independently
tests every 7! vertex map for automorphisms of each representative. There
are exactly three classes:

- Aut order7, labelled orbit720, vertex-transitive, e=0.
- Aut order3, labelled orbit1680, vertex orbits3,3,1, e=14.
- Aut order21, labelled orbit240, vertex-transitive, e=0.

The complete labelled sizes sum to2640, and each orbit passes
orbit-size times automorphism-order equals7!. Representatives and literal
orbit data are in `receipt.json`. The inherited hostile
(74,84,97,22,37,11,56) is also checked directly. All three regular classes
are strong. Thus the finite input is E_s,7=14, and Theorem4 supplies
E_s,n=0 for n<=6 and E_s,8=0.

Equation(4), with a_0=a_1=1, therefore yields

    E_8 = E_s,8 + 2 a_1 E_s,7 = 28.                     (5)

The only positive-error order-eight classes are the order-seven hostile
followed by a singleton and a singleton followed by the hostile. They are
distinct because their unique ordered SCC size lists differ.

For the growth packet's pointed extension states P8, let

    U8 = sum_(p in P8) |stable cell of the marked root in U_p|.

Then (5) proves

    A000568(8) = (U8 - 28)/8.                            (6)

U8 requires only the order-seven parent representatives, their actual
permutation automorphism actions to form attachment-subset orbits, and
ordinary color refinement of each extension. It requires no order-eight
canonicalization and no deletion/re-rooting isomorphism tests. With the
separate verified parent-only value U8=55068, (6) gives6880. This is a
structural finite-order improvement over full root transport, not a formula
from n alone, a proof of bounded higher-order error, or closure of an unmarked
blue/black line graph. A standalone predictor is maintained separately from
this proof/control packet.

## 7. A native strong nine-vertex unbalanced cycle

Take three directed triangles A,B,C with cyclic labels x in F3. Put
A_x to B_y exactly when y=x, A_x to C_y exactly when y=x, and B_x to C_y
exactly when y=-x. All other cross arcs are reversed. Adjacency bitsets are

    (74,148,289,86,293,139,182,285,107).

Vertices in A,B,C have degrees3,4,5 respectively, so degree refinement
separates the three cells. Each crossing is equitable and each internal
triangle regular, so these are exactly the stable cells. The matching
signs are +,+,-; their product is negative. Theorem2 gives Aut(T)=1,
nine fixed vertices, and e(T)=18. All three cells are strong and each
pair has arcs both ways, so the whole tournament is strong. The checker
independently tests all9!=362880 vertex maps and finds only the identity.

Changing the last matching to y=x makes the sign triangle balanced and
produces a positive control with Aut=C3 and e=0. This exhibits exactly
what cell sizes and quotient crossing degrees forget: compatible matching
signs around a cycle. No singleton pins or external observer are used.

## 8. All-r matching-triangle class counts and the next error bound

For any r>=3, make r directed triangles indexed i=0,...,r-1. For every i<j,
choose one permutation sigma_ij in S3 and put (i,x) to (j,y) exactly when
y=sigma_ij(x). All other crossing arcs are reversed. Let m=binom(r,2).
This specializes the inherited affine-connection carrier to a complete
support with intrinsic cell names. The degree in cell i is

    1 + 2i + (r-1-i) = r+i.

Thus cells have distinct degrees and are intrinsically ordered. The cells
are exactly the stable classes, all isomorphisms preserve them, and internal
isomorphisms are rotations. Every such tournament is strong, because the
internal triangles are strong and every pair of cells has crossing arcs
in both directions.

Write sigma_ij(x)=epsilon_ij x+c_ij. Each of the2^m sign assignments is
unchanged by rotations. The rotations k_i change the phases by

    c_ij -> c_ij + k_j - epsilon_ij k_i.

The linear map from F3^r to F3^m has kernel dimension1 for balanced signs
and0 otherwise, by Theorem2. The complete graph has2^(r-1) balanced sign
assignments: choose one sign at each vertex modulo common reversal and
put epsilon_ij=tau_i tau_j. Therefore the exact number of ordinary
isomorphism classes in this explicitly defined family is

    F_r = 2^(r-1) 3^(m-r+1)
        + (2^m - 2^(r-1)) 3^(m-r).                     (7)

The first summand consists of classes with Aut=C3 and e=0. The second
consists of rigid classes, each with e=6r. This is an orbit count by linear
phase rank, not a count of all tournaments of order3r.

At r=3, the already inherited count in Session39 §4 is recovered:16
ordinary strong classes:12 with Aut=C3 and4 rigid.
The checker independently re-exhausts all6^3=216 matching presentations and
forms their literal within-cell rotation orbits; it recovers these16
classes and total error72. Thus

    E_s,9 >= 72,
    E_9 = E_s,9 + [z^2]A(z)^2 * 14
        = E_s,9 + 42 >= 114,                            (8)

using a_0=a_1=a_2=1 and the preceding vanishing statements. More generally
(7) gives the explicit family lower bound

    E_s,3r >= 6r (2^m - 2^(r-1)) 3^(m-r), r>=3.

For r=4, the same formula gives720 family classes,504 rigid, and hence
E_s,12>=12096. These are lower bounds, not exact unrestricted strong-error totals. In
particular (8) does not give A000568(9), and no matching-only model is
asserted to exhaust its strong tournaments or its refinement failures.

## 9. An inherited two-prime-cell lower bound

The archived Session40 theorem already counts the two matching-coupled cyclic
regular p-cell family for every odd prime p. Its ordinary class count is

    N_p = ((p-1)! + (p-1)^2)/p.

Exactly p-1 classes have Aut=C_p and the other N_p-p+1 are rigid. Here we
use those proved family facts, not rederive or claim them anew. The cells
have distinct degrees and their matching crossing is equitable, so their
stable cells are exactly the two size-p cells. An affine Aut=C_p class has
e=0, while a rigid class has2p singleton orbits and contributes
2p(p-1). Thus the inherited family gives the new error-series consequence

    E_s,2p >= 2p(p-1)(N_p-p+1),  p>=5 prime.             (9)

For p=5 this is160 at order10; for p=7 it is8568 at order14. These are
selected-family contributions, not exact unrestricted totals. From (4),
E_s,7=14, E_s,8=0, E_s,9>=72, E_s,10>=160, and a_3=2, one obtains

    E_10 >= 6*14 + 2*72 + 160 = 388,                   (10)

because [z^3]A^2=2a_3+2a_1a_2=6 and [z]A^2=2. Positivity of all omitted
strong-error coefficients makes this a rigorous lower bound. It does not
predict A000568(10) without the remaining upper-trace and strong-error data.

## 10. Connection ledger and boundary

Source: actual stable cells, directed internal triangles, and actual
nonuniform crossing matchings. Target: a signed graph and its F3 kernel.
The map preserves the full automorphism permutation action in the one/three
cell domain, hence exact marked-root orbit sizes. It deliberately forgets
matching phase and much of tournament isomorphism type; those phases are
restored in the family count of Section8. Cell sizes plus crossing degrees
alone forget cycle compatibility, as the two nine-vertex controls show.

The SCC map preserves additive refinement error but discards internal
strong-component identity if only the series coefficients are kept.
The finite-order counting improvement needs the full parent atlas and actual
parent automorphism action, not just prior A000568 numbers or group orders.
The positive conclusion is a restricted exact correction; the next open
coordinate is E_s at unrestricted higher orders, including larger regular
cells and more complicated equitable crossings.
