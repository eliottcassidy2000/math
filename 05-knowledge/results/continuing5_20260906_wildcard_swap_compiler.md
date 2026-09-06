# Faithful column repairs: small-prime compilers and a linear swap barrier

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The constructive result is a complete seven-column repair class: after one
specified row swap in the canonical fourteen-cycle board, every column
ordering admits a nonincreasing sequence of at most three column
swap-and-affine moves to an actual no-three-in-line board. Nine starting
orbits require a plateau move. Freezing the original row addresses instead
leaves no successful column ordering at all.

The general results explain the limits of this operation. A classical
modular-slope lemma yields an affine-plus-(p-4)-swaps compiler for arbitrary
column permutations. Conversely, a board containing a complete modular line
requires a number of non-affine row/column swaps growing linearly with p
before it can be successful. A completed reciprocal graph supplies a
different native object with exactly one integer triple, together with an
exact factor-based decoder for swapping its zero column. No all-prime
successful construction or optimization running-time theorem is claimed.

## 1. Inheritance, operation and retained predicates

The inherited target is the actual integer determinant predicate in
`{0,...,p-1}^2`, with p an odd prime. A successful board has no three distinct
integer-collinear points. Modular collinearity is used only where explicitly
stated. Every row/column permutation preserves the bipartite incidence graph
and its cycle profile; it need not preserve this geometric predicate.

Closest proved mechanism:
`05-knowledge/results/continuing4_20260906_wildcard_affine_columns.md`, with
its independent audit, gives the exact forbidden-translation interval union
for each affine column orbit. The canonical hostile is its p=5 ten-cycle:
every affine column map fails while four unrestricted column permutations
succeed. The corrected near miss is extending that single non-affine repair
while forgetting the row addresses or the affine realization inside an orbit.
The least-used sidecar is the actual repair word, supplemented by integer
layers within a modular line.

The conditional column theorem
`overnight11_20260906_no3line_rowfreeze.md` and the restart consumer
`overnight12_20260906_no3line_count_restart.md`, both in
`05-knowledge/results/`, require fresh columns uniform conditional on the
history. The purposeful moves below are outside that sampling rule. Nothing
here contradicts their probability bounds or replaces conditional uniformity
by marginal uniformity. The Guy--Kelly primary-paper route remains the
inherited one; its independence heuristic is not a dependency.

Anchor: a faithful within-board compiler. Niche: modular slope products and
permutation cycles. Wildcard: reciprocal conics and retained integer factors.
The concept board, updated after the hostile probes, is:

| Object or operation | Preserved predicate/data | Essential missing coordinate |
|---|---|---|
| Column permutation | Incidence skeleton and row addresses | Actual integer triple set |
| Affine quotient | Existence of a successful affine realization | The selected affine parameter pair |
| One swap followed by affine relabeling | A concrete reachable column word | Completed-move triple count; intermediate swaps can increase it |
| Modular matching triple | Three exact column assignments | Modular zero is not integer zero |
| Unchanged points on a modular line | Native points surviving the repair | Occupied Euclidean line layers |
| Reciprocal zero swap | All untouched conic points | Integer products and the carry in a coordinate sum |

This uses the maintained meta-pattern "Quotient by the actual repair, then
test the native survivor." Its counterindication is explicit: the quotient
does not inherit geometric success as a pointwise invariant.

## 2. The faithful affine-orbit graph and a general compiler

Let H=AGL(1,p), acting on column values by y -> ay+b modulo p. Permutations
compose from right to left. A column assignment sigma is the map taking
each source column to its physical integer column. The vertices are left
cosets H sigma. Every coset has the unique representative

    r(i)=(sigma(i)-sigma(0))/(sigma(1)-sigma(0)) modulo p,

so r(0)=0,r(1)=1 and there are (p-2)! vertices. Join H sigma to H sigma tau
for every source-column transposition tau, omitting loops. This is well
defined. In physical coordinates the swap exchanges sigma(i),sigma(j),
rather than blindly reusing the source labels i,j.

For a fixed board B define

    F_B(H sigma)=min_(h in H) X(h sigma B),

where X counts unordered integer-collinear triples. The inherited interval
decoder computes both this minimum and its attaining affine parameters.
The vertex is successful exactly when F_B=0. The action of H itself does
not preserve X, which is why a minimum and an attaining realization are
retained. A charged move is one column transposition followed by an arbitrary
affine column relabeling. Monotonicity below refers to the endpoints of
these completed moves, not the intermediate bare transposition.

**Exact word-distance law.** If c(w) counts all permutation cycles, including
fixed points, then

    d(H sigma,H pi)=min_(h in H) [p-c(sigma^-1 h pi)].       (1)

A permutation on p letters needs exactly p-c(w) transpositions: each swap
changes the cycle count by one, and each nontrivial cycle of length l can be
built with l-1 swaps. A path reaches the target coset precisely when its word
is sigma^-1 h pi for some h. This proves (1). Interspersed affine maps do
not lower the cost further: conjugating a transposition by an affine map is
still a transposition. Right multiplication by any permutation is a graph
automorphism because it conjugates the complete transposition set. Thus
distances from one vertex determine the diameter.

**General upper bound.** For p>=5 this graph has diameter at most p-4.
The modular ingredient is classical: every point of a permutation graph
over an odd finite field belongs to a modular collinear triple. This is
[Cooper, Collinear Triple Hypergraphs and the Finite Plane Kakeya Problem,
Proposition 4](https://people.math.sc.edu/cooper/collinear.pdf). For clarity,
fix i and suppose all slopes (f(j)-f(i))/(j-i), j!=i, were distinct. They
would be all nonzero field elements, whose product is -1. Their numerator
and denominator products are equal, so their product is also 1, a
contradiction in odd characteristic.

Apply this fact to pi sigma^-1. An affine map g agrees with that permutation
in at least three positions. Therefore sigma^-1 g^-1 pi fixes at least three
positions. If nonidentity, its remaining at most p-3 positions have at least
one cycle, giving transposition length at most p-4. Identity has length zero.
This proves the bound constructively by choosing the matching triple and
then decomposing the residual cycles. It makes no claim that the modular
triple is an integer bad event. The verifier retains the F4 Frobenius graph,
which has no modular collinear triple, as a hostile to dropping oddness.

Consequently, if a fixed-row board has any successful column ordering, then
from every current ordering some affine-plus-at-most-(p-4)-swaps word reaches
one. The premise is essential: this compiler covers an operation space; it
does not create a successful permutation when the space contains none.
The modular lemma is cited rather than claimed as new. The contribution
here is its use with the faithful repair distance and native decoder.

## 3. Exact five- and seven-column operation laws

At p=3 the affine group is S3. At p=5, the general upper bound is one, so the
six-vertex graph is K6. Thus for any board in the 5-by-5 grid, if any column
ordering succeeds, an affine map and at most one column swap suffice from
every ordering. This is stronger than a single witness for the inherited
ten-cycle, but is still a statement confined to five columns.

At p=7 a complete exact coset calculation improves the general upper bound
three to two. The radius profile from every vertex is

    distance 0: 1 vertex; distance 1: 21; distance 2: 98.  (2)

For p>=7, all binom(p,2) transposition neighbors are distinct. Otherwise a
nonidentity affine permutation would be the product of two distinct
transpositions, fixing at least p-4>=3 points; a nonidentity affine map fixes
at most one. The exact p=7 coverage in (2) is certified both by graph BFS and
independently by the minimum cycle-count formula (1). No p>7 diameter formula
is extrapolated from it.

Now retain the canonical fixed-row cycle board

    B_p={(i,i),(i,i+1 modulo p):0<=i<p}.

The p=5 full-column space has four successful permutations. At p=7,
**all 5040 column permutations fail**. The certificate lists all 120 affine
cosets; each has a full forbidden-translation mask for all six multipliers.
Separate literal determinant and primitive-direction counts verify the entire
space. Its best triple count is one. For example the normalized representative
(0,1,2,5,3,6,4) attains one under the affine map (a,b)=(3,5); its 21 neighbor
orbits have minima 2,3,4,5 with multiplicities 8,7,3,3. This is an exact
fixed-row obstruction, not failure of the abstract C14 skeleton.

An actual successful C14 board appears after swapping rows 0 and 5. From
the identity column assignment perform source-column swaps (0,1), then
(1,2), followed by y -> 6y+3 modulo 7. This gives

    rho=(5,1,2,3,4,0,6),
    pi =(2,1,3,0,6,5,4).

The points (rho(i),pi(i)) and (rho(i),pi(i+1)) have no integer collinear
triple. Both row and column degrees are still two and the incidence graph is
still one fourteen-cycle. The row change is therefore a load-bearing address
change, not a new cycle profile. The complete 20328-word bank with at most
one ordinary row transposition, at most one column transposition and arbitrary
affine columns contains no success. This last bank does not allow free affine
row relabelings; it is not a lower bound for every two-shore operation model.

## 4. A complete monotone repair class after the row change

Fix rho=(0 5) as above, retain the C14 incidence pattern, and allow any column
assignment. There are again 120 affine column cosets. Exactly one is
successful. Every vertex has a nonincreasing F_B path to that vertex in at
most three charged moves, with exact shortest nonincreasing-distance profile

    0:1, 1:21, 2:91, 3:7.                              (3)

The finite proof is constructive. The certificate gives every normalized
representative, its minimum and an attaining affine pair, its successor,
the source-column transposition, and the following physical affine map.
Each successor lowers the stored distance by one and does not increase the
literal integer triple count. The final distance-zero board passes the
literal no-three predicate. All 120 records are independently checked against
the interval decoder and actual determinants.

For an arbitrary starting realization, first choose the attaining affine
map of its coset; this cannot increase X. Then follow the certified completed
moves. Nine nonzero vertices have no strictly improving neighbor, so a
universally valid strict-descent rule would stop. Each has a permitted
equal-score first move in the compiler. Seven vertices need three
nonincreasing moves although their unrestricted graph distance is two:
minimum word length and monotone repair length are different observables.

This is a rigorous finite repair class, with the target checked at its actual
endpoint. It asserts neither success for other row addresses nor a general
greedy method. It also explains why a cycle-profile sidecar alone was too
coarse: the unreordered and reordered boards have the same C14 skeleton and
different full-column success sets.

## 5. An all-prime lower bound from unchanged native points

Let B in [p]^2 contain a complete non-axis modular line. Permit arbitrary
affine row and column relabelings at zero cost, interleaved with k ordinary
row/column transpositions in total. Suppose the final board is successful.
Then, for p>=5,

    k >= max(1, ceil((p-2 floor(sqrt(2p)))/2)).           (4)

Proof: row and column operations commute across shores, and affine maps
normalize the transpositions within each shore. Write the final operations
as A_r tau_r and A_c tau_c, where tau_r,tau_c are products of k_r,k_c
transpositions, k_r+k_c=k. Let R,C be their supports. Their total support is
at most 2k. The original complete modular line has one point in each row and
column. At least p-|R|-|C|>=p-2k of its points are untouched by the tau's.
After the final affine maps these points lie on another complete modular
line ell, at their actual standard integer representatives.

The inherited lattice-layer theorem supplies a primitive tangent vector
(u,v) with r=|u|+|v|<=floor(sqrt(2p)). The integer functional v x-u y has one
residue modulo p on ell and at most r occupied integer levels in the box.
Each such level is an actual Euclidean line and can retain at most two
points in a successful board. Hence p-2k<=2r<=2floor(sqrt(2p)). With zero
non-affine swaps the complete modular line survives, and the inherited
equal-parity midpoint obstruction rules out success for p>=5. This proves
the extra lower bound one in (4).

The retained predicate is pointwise, so no independence or sampling law is
used. A sharper direction-specific sidecar keeps each occupied level length
ell_j of the complete final affine-image modular line ell, giving

    |R|+|C| >= p-sum_j min(2,ell_j)
              =sum_j max(0,ell_j-2).                    (5)

These are necessary bounds, not sufficient repair criteria. Moved points can
create additional bad lines, and the supports must actually arise from
permutations. Formula (5) is the exact deletion deficit for this one parallel
line predicate, not an equivalence for full no-three-in-line.

In particular, no fixed number of swaps, even with free affine operations
on both shores, can repair this complete-modular-line family for arbitrarily
large primes. This does not assume that unrestricted permutations eventually
succeed. Equivalently, any successful p-point permutation graph agrees with
any affine field permutation in at most 2floor(sqrt(2p)) positions. Its
distance in number of unequal positions from every affine permutation is
therefore at least p-2floor(sqrt(2p)). The linear upper and lower operation
bounds leave a substantial gap and make different hypotheses.

## 6. A native reciprocal object and an exact zero-swap decoder

The completed reciprocal permutation f(0)=0, f(x)=x^-1 for x!=0 is a
classical modular example discussed in the same
[Cooper paper](https://people.math.sc.edu/cooper/collinear.pdf), Section3.
Here its integer lift gives a useful different starting object.

**Proposition.** For every odd prime p its standard integer graph has exactly
one collinear triple:

    (0,0), (1,1), (p-1,p-1).                            (6)

The nonzero points lie on xy=1 over F_p, which a line meets in at most two
points. Thus a modular collinear triple must contain the origin and a pair
(x,y),(-x,-y), giving exactly (p-1)/2 modular triples. The determinant of their
standard representatives with the origin is p(x-y). It vanishes exactly when
x=y, hence x^2=1; the two choices give the same triple (6). For p>=5 this
object contains no complete modular line and therefore avoids the hypothesis
of the linear barrier, while still having a defect to repair.

Swap the output values 0 and a, where 1<=a<p, and put b=[a^-1]_p. The moved
points are replaced by P=(0,a), Q=(b,0); retain

    U={(x,[x^-1]_p):1<=x<p, x!=b}.

**Exact decoder.** Every integer bad triple is in exactly one of these types:

1. P,R,S with R=(x,y), S=(x',y') in U, where either
   y+y'=a and xy=x'y', or
   y+y'=a+p and x(p-y)=x'(p-y').
2. Q,R,S, where either
   x+x'=b and xy=x'y', or
   x+x'=b+p and y(p-x)=y'(p-x').
3. P,Q,R with R=(x,y) in U and a x+b y=ab.

To prove the first case, modular collinearity through P forces
y+y'=a modulo p: multiply its determinant by the nonzero yy' and use
xy=x'y'=1 modulo p. The two distinct y values have sum either a or a+p.
Substituting each exact sum into the integer determinant yields the displayed
product equalities. The second case is the transposed argument, and the
third is the exact line through the two intercepts. Three points entirely
in U cannot be collinear. These cases are disjoint by which moved points
they contain, proving an iff decoder with O(p^2) pair tests.

The factors are native integers. Replacing their equalities by congruences
would lose the very information needed to decide success. The bounded exact
bank gives successful a values

| p | Successful swaps of output 0 with a |
|---:|---|
| 3 | 1 |
| 5 | 1,4 |
| 7 | 1 |
| 11 | none |
| 13 | 1 |

At p=11 all ten swaps leave at least one actual bad triple. This is the first
odd-prime hostile to universal repair by a zero-column swap, not a claim that
all single transpositions or all column permutations fail. The strongest
survivor is the exact factor-and-carry decoder. The next meaningful question
is to find a provable family with empty factor packets or a specified larger
move that removes them; a larger blind census is not the present result.

## 7. Exact certificate and reproduction

The standalone source has no repository imports. For p=3,5,7 it exhausts
all (p-2)! normalized affine cosets and their p(p-1) parameter pairs, checks
their exact interval counts against literal determinants and independently
against anchored primitive-direction counts, and retains the cycle graph.
It independently compares quotient BFS with permutation cycle lengths.
The complete reordered p=7 class supplies all successor witnesses in (3).
The 20328-word small repair bank, all reciprocal zero swaps at p<=13, and
short-direction/occupied-layer controls for every non-axis modular line at
the nine primes 5 through31 are explicit bounded universes. Every safe subset
of those modular lines is additionally checked at p=5,7. No randomized or
floating computation is used.

Run the source directly, or after filing it in 04-computation:

    python continuing5_20260906_wildcard_swap_compiler.py
    python -O continuing5_20260906_wildcard_swap_compiler.py

It writes the exact JSON beside itself outside the repository, or into
05-knowledge/results in the filed layout. All **73164** gates remain active
under -O. Normal and optimized runs produce byte-identical LF transcripts.

- Source SHA256: `d15877188ab7714390cf9334838151c9b20be39ca553d333530a14ce9cdccad1`.
- Output SHA256: `ede7846755e9eab28d6fea80f3ace8d9f7e94c40c07d0e9d0ef6a0333228d0f3`.
- JSON SHA256: `cb41dab50c3505c48dac732f67013337bb936bddf1f0d762a08e48e3a337a982`.

The exact proofs and finite certificates are separate from
unresolved asymptotic existence, large-prime diameters, and algorithmic cost.

All files were created outside the repository. No scarce theorem IDs,
maintained navigation, Git state or previous frozen artifacts were changed.

Independent [proof and exact referee](continuing5_20260906_wildcard_swap_compiler_audit.md) passes.

Filed checkpoint provenance: the [raw-byte manifest](continuing5_20260906_manifest.json)
pins the final report, source and output. Reviewed candidate-report hashes
above refer to the pre-promotion bytes. Source-location defaults and local
links were made portable where necessary; all emitted outputs were replayed
as raw bytes. The independent audit supplies the stated promotion basis.
