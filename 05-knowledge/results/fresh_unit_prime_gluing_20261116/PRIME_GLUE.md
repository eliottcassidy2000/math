# An obstruction to the proposed exact prime gluing

Let H(T) count directed Hamiltonian paths of a tournament. Write A,B
for induced subtournaments of a proposed host S, W=A∩B, p=H(A),
q=H(B), r=H(W). Set H(empty)=1. Use the positive odd-cycle-pack
expansion of H: each pack contains vertex-disjoint directed odd
cycles and has weight 2 to the number of cycles. The packs contained
in A or B weigh exactly p+q−r. Thus if the proposed gluing satisfies
H(S)=p+q−1, the **entire** weight budget outside those packs is r−1.

**Proposition.** If p,q>1 and W is transitive, such an induced
gluing is possible only when p=q=3. If p,q are primes at least 11,
the gluing also cannot use any overlap with at most three vertices.

## Transitive overlap

Choose any directed triangle in each side. Both triangles exist
because a tournament with no directed triangle is transitive and
has H=1. If the two are disjoint, their joint pack gives forbidden
mixed weight 4. If they share one vertex, their five-vertex union
is strong and has a directed Hamiltonian cycle (Camion's theorem).
This cycle is odd and uses vertices outside each side, giving
forbidden mixed weight 2. They cannot share all three vertices
because W is transitive. Thus every pair shares two vertices of W;
each individual triangle meets W in exactly two, and fixing once
on each side shows all triangles on both sides use the same pair
a,b of vertices in W.

In either input, deleting a or b leaves a transitive tournament.
In the common order of the remaining vertices, a,b each has an
insertion position. The only cyclic possibility is a backward
arc b→a with a nonempty interval of m vertices, each admitting
a→x→b. Outside that interval all vertices are inert to the strong
component. Its cycles close the backward arc by an increasing
subset of the interval, and the odd cycles use odd subsets; hence
H=1+2^m. Every interval vertex is outside W, since it forms a
triangle with a,b. If either interval has at least two vertices,
take three interval vertices using both sides. In the subtournament
on these three vertices take a directed Hamiltonian path and close
it through b→a. This mixed directed 5-cycle exceeds the zero budget.
Thus both intervals have size one and p=q=3. Conversely the
four-vertex bracket with two interior vertices realizes H=5.

## Cyclic triple overlap

It remains, for the second claim, to consider W a directed triangle.
If p,q are primes at least 11, each side has a unique nontrivial
strong component, by the product rule along the ordered strong
components. It contains W and has size ≥5: on four vertices
the only possible odd cycles are triangles, and there are at most
two (each directed triangle uses three arcs and the standard
outdegree count gives a maximum of 2). Thus H≤5 at order four.

Every directed triangle in a strong tournament of size at least
five extends to a strong induced 5-set. Starting at the triangle,
if an outside vertex sees both directions, adjoin it and continue
if still below five vertices. If there is no such vertex, all
outsiders point uniformly toward or away from the set. Strongness
forces a reverse arc from the away side to the toward side;
adjoin its endpoints. If this makes a six-set, undo the first
single-vertex addition: the reverse pair sees the original triangle
uniformly and makes a strong five-set with it.

Take such five-sets A',B' inside the two strong cores. Their
seven-vertex union is strong (paths through W), so its spanning
directed 7-cycle spends all of the mixed budget r−1=2. We show
one more mixed pack is forced. Every outside vertex in each
five-set lies on a directed triangle within that set: otherwise
all its in-neighbors point to all its out-neighbors, making a
return from the out-neighbors impossible. If directed triangles
touching their respective exteriors are disjoint or intersect
once, they give a mixed pack or mixed 5-cycle as above. Otherwise
every such cross pair shares two vertices of W, so the same fixed
pair occurs in every exterior-touching triangle on both sides.
In particular all four exterior vertices lie between these two
shared vertices in the backward-edge bracket sense. Three
exterior vertices using both sides then give a mixed 5-cycle,
again as above. In both cases this is additional to the 7-cycle,
contradicting budget 2.

Optional sharper finite certificate: `overlap3.cpp` enumerates
the 18 free orientations of two five-sets over a fixed shared
triangle, retains the 153664 completions with both sides strong,
and computes the mixed weight by three exact subset DPs.
The minimum is 10, attained by 72 completions. This implies
H(S) ≥ p+q+7 in the cyclic triple case, but is not needed for
the analytic no-gluing proof.

## What must change

This obstructs *induced witness retention*, not arithmetic
coverage. An additive prime transfer must either bring a
genuinely cyclic interface of at least four vertices or alter
the witnesses themselves. An exterior repair next to an inert
common overlap cannot evade the mixed-pack budget by making
the window larger: the obstruction is at least one forced
positive interaction. An attractive constructive inversion is
to program the shared reservoir's pack weight to equal mixed
weight plus one, rather than ask an independent seam to give
a signed cancellation of positive packs.

## Positive exact families: the overlap pays for all interactions

Code a labeled tournament by the bits for pairs (i,j), j increasing
and then i increasing through i<j, starting with bit zero. Bit one
means i→j. The following hosts meet J literally: deleting u
or v gives the two induced summands. The joint deletion gives r.

| order | code | u,v | p,q | H(S) | r |
|---:|---:|:---|:---|---:|---:|
| 6 | 1105 | 4,5 | 13,13 | 25 | 5 |
| 7 | 470924 | 1,4 | 29,23 | 51 | 13 |
| 7 | 393231 | 1,3 | 23,17 | 39 | 9 |
| 7 | 1053296 | 0,6 | 17,19 | 35 | 9 |
| 7 | 1615491 | 2,5 | 31,23 | 53 | 13 |

Thus four common vertices suffice for some prime summands, and five
for the illustrated distinct ones. A blanket J statement beginning
at primes 11 already fails at (11,11), since H=21 is impossible.
Neither obstruction forbids constructing the target by another route.

A physical infinite extension comes from code 393231. Inflate its
common vertices 0 and 2 to transitive modules of sizes s,t≥1
retaining all outside orientations. Put X=2^s, Y=2^t. Then

    H(S)=16XY−12X−3Y+5,
    p   = 8XY− 4X−2Y+3,
    q   = 8XY− 8X− Y+3.

In particular both induced counts lead to exact J for all s,t.
Only the displayed seed and some other small settings have been
checked as prime pairs; no infinitude of prime pairs is asserted.

Proof: after dropping both modules from a Hamiltonian path, keep
the ordered exterior vertices. In each gap place neither module,
one run of either module, or one run of both in their unique
allowed order. Check each used adjacency against the tournament
matrix. The resulting skeleton counts, indexed by numbers
(1,1),(2,1),(1,2),(2,2) of runs in the two modules, are

    S: (39,20,29,16)
    A: (23,12,14,8)
    B: (17, 8,15,8).

No other positive run counts occur. A transitive module of size s
contributes one filling for one run, and 2^s−2 for two runs:
choose any nontrivial two-part ordered set partition of its
vertices and write the two parts in transitive order. Substitution
gives the formulas. The coefficient difference is exactly −1 at
(1,1) and zero otherwise. More generally replace the transitive
modules by arbitrary tournaments M,N. Their one-run responses
are H(M),H(N). In the same skeleton count, all higher-run responses
still cancel (these are ordered spanning path-cover counts of
the modules). The equation becomes H(S)=p+q−H(M)H(N).
This is why a richer shared reservoir changes the offset and
cannot simply be treated as a free programmable source.

The inverse identity

    4H−11 = (16X−3)(4Y−3)

shows exactly what a target would have to supply. Likewise
p−2=(4X−1)(2Y−1), q−2=(8X−1)(Y−1). This is a physical construction,
but the powers-of-two restrictions give only O(log² N) distinct
host outputs up to N, not a universal Goldbach step. To obtain
absolute coverage one needs new port profiles that move the two
factors without changing the target offset, or compensate the
resulting H(M)H(N) term.

Reproduce `python3 check_glue.py` (only standard library) for the
exact host permutations and inflation skeletons. For the optional
sharp overlap-three margin, run
`c++ -O3 overlap3.cpp -o /tmp/overlap3 && /tmp/overlap3`.

One can keep the exact inverse when changing modules. Set u=H(M),
v=H(N), R=P₂(M)+2u and T=P₂(N)+2v, with P₂ counting two nonempty
ordered paths. Then

    4H−11uv=(16R−3u)(4T−3v),
    p−2uv=(4R−u)(2T−v),
    q−2uv=(8R−u)(T−v).

This follows by substituting P₂(M)=R−2u and P₂(N)=T−2v into
the four skeleton terms. For example a singleton M gives
H=29T−19v, a concrete selection line in a physical module
profile. R counts ordered two-path covers where either path
may be empty. For an ordinal sum, restricting each labeled path
to each summand and concatenating back proves R is multiplicative.
Thus (H,R) is a physical coordinatewise multiplicative semigroup.
This does not prove it meets every target line, but it replaces the
one-off embedding problem by an exact algebra of admissible sources.

For the singleton line, a necessary condition is v≡3H mod 29.
Choosing a fixed small representative for each residue cannot suffice.
For fixed v, there are only finitely many lists of nontrivial
strong components of a tournament with count v. They contribute
factors ≥3 to the ordinary count, and each has bounded order:
vertex pancyclicity and the odd-cycle expansion already give
at least 1+2 floor((n−1)/2) for any strong component of order n.
Singleton components do not change v but each doubles T.
Therefore the T-fiber for fixed v is a finite union of exact
dyadic rays. A dense target selector of this form must vary v
without a bound, or change the first module as well.
