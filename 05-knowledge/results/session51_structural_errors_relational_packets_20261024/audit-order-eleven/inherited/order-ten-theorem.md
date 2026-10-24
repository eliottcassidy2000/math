# Exact structural refinement error at order ten

24 October 2026. **PROVED:** exhaustive strong-order-ten stable-cell classification,
its error sum, and the two-prime-cell matrix formula. **FINITE-EXACT:** independent
small-sector enumeration and direct automorphism controls. The total-order-ten
coefficient additionally uses the separately certified Session50 input
E_s,9 = 716. No full order-ten tournament enumeration or order-ten class count
is used or produced.

## 1. Statement and inheritance

For an ordinary tournament T, start directed color refinement with one color and
repeatedly retain the current color and the outneighbor counts into every current
color. Its eventual cells are called stable cells. For an actual automorphism
vertex orbit O, let C(O) be its stable cell and put

    e(T) = sum_O ( |C(O)| - |O| ).

This is the excess in the pointed refinement trace, not the different symmetry
defect of a root-changing matrix. Write E_s,n for the sum of e over ordinary
strong n-vertex tournament classes, and E_n for the sum over all ordinary classes.
Converse tournaments are not identified.

**Main theorem.** The positive-error strong order-ten classes have exactly one
of the following three stable-cell patterns:

- 5+5: 84 rigid classes, each with error 40; total 3360.
- 7+1+1+1: four classes containing the unique non-vertex-transitive regular-seven
  type, each with error 14; total 56.
- 3+3+3+1: 32 rigid classes, each with error 18; total 576.

Consequently

    E_s,10 = 3360 + 56 + 576 = 3992.

The inherited SCC identity E(z)=A(z)^2 E_s(z), together with
E_s,n=0 for n<=6, E_s,7=14, E_s,8=0, E_s,9=716 and
(a_0,a_1,a_2,a_3)=(1,1,1,2), gives

    E_10 = E_s,10 + 2 E_s,9 + 6 E_s,7 = 5508.

Thus an independently obtained pointed upper trace U_10 must satisfy

    a_10 = (U_10 - 5508)/10.

This packet does not compute U_10 and does not assert a numerical a_10.

Closest inherited mechanisms are the equitable-cell and signed-triangle
holonomy theorems, the regular-seven hostile, the exact order-nine correction,
and the SCC additive-error transfer in the Session50 packet. The two-C5
matching-only orbit count is inherited from Session40, as archived in Session50.
The degree-two C5 crossing sector, all singleton-decorated triangle sectors,
and their exhaustive order-ten combination are the new steps here.
The proof reproduces the needed small mechanisms rather than appealing only to
finite outputs. The finite regular-seven classification is independently rerun.

## 2. Basic restrictions and exhaustive partition list

Every stable cell induces a regular tournament: internal outdegree is constant
by equitability, and counting internal edges makes that degree (m-1)/2.
In particular its size m is odd. For distinct cells of sizes a,b, let r be each
first-cell vertex's outdegree into the second and s each second-cell vertex's
indegree from the first. Double counting gives ar=bs. If gcd(a,b)=1, r is 0 or b,
so the crossing is uniform. Every singleton is uniformly attached to every cell.

The complete odd partitions of ten are

    9+1; 7+3; 7+1+1+1; 5+5; 5+3+1+1; 5+1+1+1+1+1;
    3+3+3+1; 3+3+1+1+1+1; 3+1+1+1+1+1+1+1; 1+...+1.

These are possible size patterns, not a claim that each is realized by stable
cells. The patterns 9+1 and 7+3 have just two cells and a uniform crossing,
so they are not strong. A regular five-vertex tournament is uniquely C5:
fix v, order its two outneighbors a->b and its two inneighbors c->d;
regularity forces b->{c,d}, a->c, d->a and hence the cyclic tournament.
Its full automorphism group is C5 (see Section 4), which is transitive.
Outside 5+5, a five-cell has only coprime-sized other cells; all its internal
rotations extend independently, and the remainder contains at most one triangle.
Those patterns have zero error.

For one- and three-cells, a nonuniform crossing between triangles is a
permutation matrix or its complement. In oriented cyclic coordinates on F3,
its matching is x->epsilon*x+c, where epsilon is +1 or -1. Rotations by k_i,k_j
preserve this crossing exactly when k_j=epsilon*k_i. On the graph of nonuniform
crossings, a balanced component permits one independent C3 rotation, while a
negative cycle forces every rotation in its component to vanish. Every
unbalanced triangle therefore has three fixed vertices and contributes 6;
a balanced triangle contributes zero. This is the full automorphism group when
the cells are actual stable cells, because every automorphism preserves each
such cell individually. With at most two triangles the support has no cycle.
Thus every remaining pattern except 3+3+3+1 has zero error.

Only the three stated patterns remain. This is the completeness argument;
finite computations below never enumerate unrestricted order-ten tournaments.

## 3. Seven-cell substitutions: 56

All crossings in 7+1+1+1 are uniform, so the tournament is a substitution into
a four-vertex quotient Q with one distinguished vertex replaced by a regular
seven-tournament R. Every regular tournament is strong: different SCCs occupy
disjoint degree intervals, contrary to all degrees being equal. Therefore the
substitution is strong exactly when Q is strong.

There is one ordinary strong four-vertex tournament and its automorphism group
is trivial. Here is a direct proof. Every degree is 1 or 2, with two of each.
Let the degree-one vertices be a,b with a->b, and the other vertices c,d.
Then c,d beat a, and b beats exactly one of them, say c. Hence d beats b,
and degree two at d forces c->d. All six arcs are determined. An automorphism
preserves the two degree-pairs and their internal edge, so fixes all four.
There are consequently exactly four rooted strong quotient classes.

Each gives genuine stable cells 7+1+1+1. The seven-cell has degree 3+k, where
k is the root's quotient outdegree, necessarily 1 or 2. A singleton has degree
in 0..2 if it loses to the seven-cell, and in 7..9 if it beats it. Thus the
seven-cell is separated at the degree step. The three singletons cannot be one
stable triangle: that would require identical relations to the seven-cell,
contrary to 1<=k<=2. No stable two-cell exists, so they become separate.
The full automorphism group is exactly Aut(R), independently on R.

The independent finite regular-seven certificate generates all 2640 labelled
regular tournaments and partitions them under all 7! permutations. It obtains
three classes with automorphism orders 7,3,21 and labelled orbit sizes
720,1680,240. The order-3 class has vertex orbits 3,3,1 and error
(7-3)+(7-3)+(7-1)=14; the other two are vertex-transitive and have zero error.
Thus each of four rooted quotients contributes 14 and this sector contributes56.
No substitution is counted twice: its seven-cell is canonical, and an
isomorphism induces a rooted quotient isomorphism and an internal R-isomorphism.

## 4. Two five-cells: 3360

Represent each C5 by x->x+1,x+2 modulo five. Its automorphisms are precisely
translations. Indeed an automorphism preserves v's two ordered outneighbors;
the unique source of their transitive two-tournament is v+1. Hence it preserves
successor and is determined by the image of one vertex.

Let M be the crossing matrix from the first cell to the second. It is a binary
5-by-5 matrix with constant row and column sum d. Uniform d=0 or5 is nonstrong.
The two cell degrees are 2+d and 7-d, always different. Name the lower-degree
cell first. Then d is uniquely 1 or2. Conversely every such matrix produces
a strong tournament with actual stable cells 5+5: different degrees separate
them, each is internally regular, and the crossing is equitable. Strongness
follows from internal strongness and arcs in both directions between the cells.

Isomorphisms preserve the degree-ordered cells and act by translations (a,b)
in C5^2. A stabilizer element on a coordinate axis would make all rows or all
columns equal, and constant row and column degree would then force d=0 or5.
So a nonuniform matrix has no nonidentity axis stabilizer. A subgroup of C5^2
is a vector subspace. It can therefore be trivial or a one-dimensional line
of nonzero slope, but cannot be all of C5^2. The latter case would contain axes.
A nontrivial line projects transitively to each cell, so e=0; a trivial
stabilizer means full rigidity and e=10*(5-1)=40.

For d=1 there are 5!=120 matrices. For d=2, interpret M as a simple 2-regular
bipartite graph on five named row and five named column vertices. Its components
are even cycles using k row and k column vertices, k>=2. There are only the
component-size partitions 5 and 3+2. A k+k bipartite cycle has k!^2/(2k)
labelled possibilities. Thus the exact number is

    R(5,2) = (5!)^2 * ( 1/10 + (1/6)*(1/4) ) = 2040.

Every translation (a,b) with a,b nonzero has five orbits on matrix positions,
one position in each row and each column per orbit. An invariant d-regular
matrix is exactly the union of d of these five orbits, so it has binom(5,d)
choices. There are sixteen such nonidentity translations. Each matrix with
nontrivial stabilizer is fixed by exactly four of them. Therefore exactly
4*binom(5,d) matrices have nontrivial stabilizer, each with orbit size five.
The rest have orbit size25. This gives

    d=1: four nonrigid classes; (120-4*5)/25 = four rigid classes;
    d=2: eight nonrigid classes; (2040-4*10)/25 = eighty rigid classes.

The total error is (4+80)*40=3360. There are96 strong classes in this stable-cell
sector. It is incorrect to count d=3,4 again: reversing the names of the two
cells gives precisely the same ordinary classes in the canonical d=2,1 sectors.
It is equally incorrect to omit d=2: matching-only crossings give just160,
not the exhaustive3360.

## 5. Three triangles plus one singleton: 576

Positive error requires an unbalanced support cycle. On three triangle cells
this forces all three crossings to be nonuniform. Their row sums are1 or2;
orient the majority edge i->j when a vertex in cell i beats two vertices in j.
The majority tournament Q on the cells is either transitive or cyclic.
Let t_i=1 when every vertex of cell i beats the singleton, otherwise t_i=0.
The nine-vertex core is always strong: each triangle is strong and each pair
of cells has arcs in both directions. Adding the singleton preserves strongness
if and only if the t_i are nonconstant.

For a fixed intrinsically named majority/attachment skeleton, there are eight
matching-sign assignments and27 phase assignments. Four sign assignments have
negative product. Their phase-translation map F3^3 -> F3^3 has zero kernel,
hence is bijective: each negative sign assignment is one isomorphism class,
with trivial full automorphism group. Each such class has error18. For a
balanced sign assignment the kernel has size3, so there are three phase orbits;
four balanced signs yield twelve zero-error classes.

It remains essential to prove the skeletons are genuinely distinguished and
that their proposed cells are actual stable cells.

### 5.1 Transitive majority

Order the three cells by increasing degree in the nine-vertex core. These
core degrees are3,4,5, and the decorated degrees are3+t_0,4+t_1,5+t_2.
All eight attachment strings give actual cells3+3+3+1:

- If two triangle degrees tie, they are adjacent in this order and have
  different t values. The singleton has a different degree, and its color
  separates the tied triangles at the next step.
- In a nonconstant attachment, the singleton degree is3 or6. It can tie with
  a triangle only for t=(0,1,1), when it ties with the lowest triangle, or
  t=(0,0,1), when it ties with the highest. In the first case their outcounts
  into the middle triangle are0 and1; in the second they are3 and2.
  Thus the next refinement separates them too.
- With constant attachments the singleton is a source or sink and degrees
  already distinguish all four cells.

The intended four-block partition is equitable, so refinement never splits an
internal triangle. The preceding distinctions prove it is exactly the stable
partition, not merely an equitable refinement chosen by hand.

The transitive majority order intrinsically fixes all cell names. Exactly six
attachment strings are nonconstant; each has four rigid phase classes. Hence
there are24 positive-error classes and their contribution is24*18=432.
Their72 balanced companions have zero error.

### 5.2 Cyclic majority

Each triangle vertex has degree4+t_i. For nonconstant t, the singleton degree
is3 or6, distinct from every triangle degree. The unique triangle with the
minority t value is separated by degree. The two equal-t triangles have
different outcounts1 and2 into that triangle, so refinement separates them.
Thus these too have exactly the intended stable cells3+3+3+1.

The cyclic majority quotient has rotation group C3. Its nonconstant binary
attachments have precisely two orbits, according as sum(t) is1 or2. Each
attachment breaks every nontrivial quotient rotation; the uniquely marked
triangle plus the oriented cycle fixes all three names. Each of the two
skeletons therefore has four rigid phase classes. Their total is8*18=144,
with24 balanced zero-error companions.

For constant t, all nine core vertices have equal degree, and uniform singleton
adjacency never splits them. The actual stable pattern is9+1, not3+3+3+1;
the tournament is also nonstrong. These presentations are a required hostile
control, not extra classes in this sector.

Combining the two majority types gives32 rigid strong classes and total
432+144=576. Isomorphism cannot mix their skeletons because it preserves the
actual stable cells, crossing-majority directions, and singleton relations.

## 6. Completion and the exact guards

The ten odd partitions in Section2 are exhaustive. Their only nonzero strong
sectors were fully classified in Sections3--5, proving E_s,10=3992. Distinct
sectors cannot intersect because actual stable-cell sizes are isomorphism
invariants. Within each sector, the explicitly described group acts on a
presentation set which is exhaustive, and every isomorphism belongs to that
group. Ordinary isomorphism, not converse identification or arbitrary block
relabeling, is used throughout.

The checker independently validates: tournament axioms; actual stable cells;
strong connectivity; exact full vertex automorphism groups by degree-pruned
arc backtracking; literal orbit partition sizes; positive and zero-error
controls; all constant-attachment exclusions; and the normal/-O invariance of
every check. It enumerates2640 labelled regular sevens,2160 canonical-degree
five-by-five crossings, all64 four-vertex quotients, and13824 named
triangle/singleton presentations. It never enumerates all2^45 order-ten
labelled tournaments or builds a full order-ten class atlas.

## 7. General two-prime-cell family formula

For any odd prime p, use two regular carousel tournaments R_p, where each vertex
beats its next(p-1)/2 cyclic successors. Their full automorphism group is C_p:
the outneighborhood is a transitive tournament whose unique source is the
successor, so every automorphism preserves successor. Let R(p,d) be the number
of binary p-by-p matrices with every row and column sum d.

The same degree-order normalization uses1<=d< p/2. The crossing stabilizer is
trivial or a full-slope line in F_p^2. There are(p-1)^2 nonidentity translations
with two nonzero coordinates, each fixing binom(p,d) matrices, and every
nontrivially stabilized matrix is counted p-1 times. Hence exactly
(p-1)*binom(p,d) matrices have nontrivial stabilizer. The class counts are

    rigid:       ( R(p,d) - (p-1)*binom(p,d) ) / p^2;
    transitive:  (p-1)*binom(p,d) / p;
    total:       ( R(p,d) + (p-1)^2*binom(p,d) ) / p^2.

Every rigid class contributes2p(p-1), and every other class zero. Thus this
specified strong family has exact total refinement error

    sum_(1<=d<p/2) 2(p-1)/p * ( R(p,d) - (p-1)*binom(p,d) ).

For p=3 all three-by-three nonuniform matrices have nontrivial stabilizer,
so an individual edge cannot freeze a triangle; a frustrated signed cycle is
needed. For p=5, degree-two matrix edges alone already produce eighty rigid
classes. This distinguishes edge rigidity from loop frustration and provides
a clean source of higher-order error.

At p=5 this family exhausts the5+5 sector because the regular-five tournament
is unique. For p>=7 it is only a selected family: other regular internal
p-tournaments exist. Composite p also has additional subgroup types, so the
prime formula is not asserted there. Neither formula predicts unrestricted
higher-order errors from cell sizes or crossing degrees alone.

## 8. Information contract and scope boundary

Source: actual stable-cell partitions, internal oriented tournaments, and literal
crossing matrices/matchings. Target: class-weighted refinement error. The map is
exact finite quotienting by the permitted internal orientation-preserving
actions, after proving actual-cell recognition. Preserved predicates are
ordinary isomorphism, strongness, orbit sizes and e. Lost data in degree-only
quotients are matrix translation stabilizers, matching signs, and phase orbits.
Those are precisely the sidecars supplied in the proof.

The order-ten structural correction is exact, but its use in a_10 prediction
still needs the independently computed upper trace and an authenticated
order-nine parent atlas. No all-order scalar recurrence, asymptotic speedup,
unmarked blue/black-metagraph sufficiency, or literature priority is claimed.

## 9. Full defect distribution and second moments

The classification proves more than an aggregate sum. The positive-error
strong-order-ten class polynomial is

    4 t^14 + 32 t^18 + 84 t^40.

The pinned Session50 regular-nine reroot certificate gives the regular-sector
histogram: three zero-error classes, three classes with e=18, two with e=36,
and seven with e=72. `distribution_check.py` independently reconstructs this
histogram from the receipt's85 almost-regular order-eight parents, canonicalizing
only eight-vertex deletions. Completeness of that85-parent slice is inherited
from the authenticated Session50 atlas, not inferred merely from reroot closure.
Adding the one seven-cell substitution and four rigid triple classes gives

    strong positive-error order9: t^14 + 7 t^18 + 2 t^36 + 7 t^72.

The smallest defective strong tournament has order7. Therefore below order14,
any tournament has at most one defective strong component. Ordered arbitrary
prefix and suffix components contribute exactly the multiplicative A(z)^2
factor separately for each positive defect value, with no collision of two
defective components. At order10 this gives the full positive-error polynomial

    12 t^14 + 46 t^18 + 4 t^36 + 84 t^40 + 14 t^72.

Thus exactly160 ordinary order-ten classes have nonzero refinement error,
whose total is5508. In particular the structural correction is sparse in
class support, not only numerically small. This is not a claim that all the
remaining classes have no automorphisms; it says their stable cells equal their
automorphism orbits.

More generally let S(z,t)=sum_[C strong] z^|C| t^e(C), with singleton strong,
and let F(z,t)=sum_[T] z^|T| t^e(T), including the empty tournament. SCC
additivity and unique ordered decomposition give the exact formal identity

    F(z,t) = 1/(1-S(z,t)).

Put M_k(z)=sum_[T] e(T)^k z^|T| and
M_s,k(z)=sum_[C strong] e(C)^k z^|C| for k>=1. Differentiating using
t*d/dt and evaluating at t=1 yields

    M_1 = A^2 M_s,1,
    M_2 = A^2 M_s,2 + 2 A^3 M_s,1^2.

These are raw moments; ordinary second derivatives instead give factorial
moments e(e-1), with the analogous consistently interpreted formula. The
cross-component second-moment term first appears at order14. The bivariate
series retains the full defect distribution rather than collapsing it to the
first-moment correction.
