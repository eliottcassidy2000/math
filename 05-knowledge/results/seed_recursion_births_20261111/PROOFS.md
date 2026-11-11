# Stable seed recursion, head completion, and composition boundaries

**11 November 2026. PROVED IN NOTE / FINITE-EXACT, as specified below.
Universal odd-coefficient seeds and Chan–Pak remain OPEN.** No external
novelty claim, independent-person audit, or formal verification is claimed.

## R0. State and inherited inputs

H(B) counts directed Hamiltonian paths of a tournament B. For distinct a,b,
Q_ab counts ordered pairs of nonempty disjoint directed paths spanning B,
the first ending at a and the second starting at b. Put W_ab=Q_ab+Q_ba.
U_b counts such ordered covers with the second starting at b, now allowing
the first path to be empty. For a nonempty selected set O put

    W(O)=sum_{a<b in O} W_ab, U(O)=sum_{b in O} U_b,
    delta(B,O)=H(B)-W(O).

A unit seed has delta=1. Its coefficient is the odd part c of U=c*2^j.
The inherited tail B -> TT_t -> r, with r -> O and B\O -> r, has

    H(tail(B,O,t))=1+2^t U=1+c*2^(t+j), t>=1.               (R0)

This is [THM-5126 / unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law](../../../01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md).
An unselected vertex is stable when replacing it by TT_a preserves delta=1
for every integer a>=1. With P_j(a)=j!S(a,j), the inherited criterion is
that its entire defect coefficient vector be (1,0,...,0), not merely that
its first few replacements work. See E2/E10 of the
[endpoint compiler](../endpoint_coverage_compiler_20261111/PROOFS.md).

The [Session258 core classification](../session258_endpoint_seed_20261111/seed_decision/CORE_CLASSIFICATION.md)
is CITED: every positive unit seed has a unique ordered decomposition
TT_a -> (S,O) -> TT_b with S strong or a singleton, O inside S, and
U_global=2^b U_core. Consequently transitive padding does not supply new
odd coefficients. The inherited complete ordinary-base census through
order seven has 95 internal stable-port presentations and 49 internal
stable-pair presentations. Presentations count all marks on an ordinary
isomorphism representative; they are not automatically marked orbits.

## R1. Full stable sets commute, beyond the 49 pairs

**FINITE-EXACT skeleton classification, with an all-size consequence.**
Group the 95 presentations by (ordinary core,selected set). There are
60 states. Their numbers of stable vertices are distributed as follows:

| Stable vertices in a state | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|
| Marked core presentations | 35 | 18 | 5 | 1 | 1 |

Every full stable set is transitive. Every such set is jointly stable:
replace all its vertices by arbitrary nonempty transitive blocks, with
independently chosen sizes a_1,...,a_m. The resulting marked graph is
still a unit seed. Thus every coordinate increment is a legal recursive
move and these moves commute. Keeping a chosen block decomposition is
part of the construction state.

**Finite certificate proves the infinite assertion.** Collapse maximal
consecutive runs inside each block to its port label. Let r be the number
of ordinary, uninflated vertices. Between ordinary vertices a quotient
word cannot repeat a marked symbol: the marked subtournament is acyclic.
A one-path word therefore visits each port at most r+1 times, and an
ordered two-path cover at most r+2 times. Each word contributes
product_i P_(j_i)(a_i), since an ordered partition of a transitive block
into nonempty paths uniquely specifies its runs. Hence delta and U have
finite coefficient tensors. The exact word enumeration in `cones.py`
finds, for each of the 60 states,

    d_(1,...,1)=1; all other d_j=0.                         (R1)

It also records the full nonnegative U tensor. The bound above is a
proof of exhaustiveness, not an empirical truncation in the block sizes.
The P_j form a triangular, independent basis on positive integers, so
(R1) is also the exact criterion for simultaneous stability in this
acyclic-port model. Literal all-size-two expansions were checked by two
different endpoint-count algorithms. See `cones.json` and `verification.json`.

The U tensor is a finite recursive register. In each coordinate,

    P_j(a+1)=j(P_j(a)+P_(j-1)(a)), P_0(a)=0 for a>=1.

The coordinate operators commute. This provides an executable infinite
family, but no assertion that its odd outputs exhaust the odd integers.
For example the five-port core (n,code,O)=(7,32,65) has

    U=3*2^(a_1+...+a_5), c=3.

Its five degrees of freedom never change c. Counting parameters without
tracking the output and its 2-adic valuation exaggerates coverage.
Here the stored slope tensor factors as
3 product_i(2P_1(a_i)+P_2(a_i)), proving the displayed formula directly.

## R2. No new old stable vertex along any of the 95 single inflations

**FINITE-EXACT all-size classification for these sources.** For each
stable v and every different unselected old w, let d_ij be the two-port
defect tensor. After v is replaced by TT_a, w is stable if and only if

    sum_i d_ij P_i(a)=0 for every j>=2.                     (R2)

The j=1 coefficient is already identically 1. Two marked vertices always
induce an acyclic tournament, so the finite tensor is complete. Expand
P_i(a)=sum_(b=1)^i (-1)^(i-b) binom(i,b) b^a. A nonzero column becomes
f(a)=sum_b e_b b^a. For q the largest base with e_q!=0, choose A with

    |e_q|q^A > sum_(b<q) |e_b|b^A.

The inequality remains strict at every a>=A, since each b/q<1. Testing
all smaller positive a and all columns gives an exact finite root test.
No numerical root approximation is used. `port_transport.json` includes
the chosen column, bound, and complete tensor for all 322 transports.

Result: an originally stable w stays stable for every a; an originally
unstable w never becomes stable for any a. There are 322 literal a=2
controls. This classifies single inflations of the inherited sources.
It does **not** exclude activation after several different inflations,
reselection, edge changes, or changing the strong skeleton. Clones inside
the chosen stable block are themselves stable, by TT_a[TT_b]=TT_(a+b-1).

## R3. Seven head families are closed; two entire coefficient lines

**PROVED by finite exact kernels plus R0.** For every r,s>=1 there is an
explicit tournament with

    H=1+(3*2^r+1)*2^s,
    H=1+(21*2^r+1)*2^s.                                    (R3)

The inherited seeds cover s>=5 and s>=3, respectively. The following
Hamilton kernels fill every missing exponent. In each row use the base
in `head_completion.json`, replace its indicated vertex by TT_a, and set
a=r-start+1. Entries are coefficients of P_1(a),P_2(a),P_3(a).

| Coefficient c | Exponent s | start | Hamilton kernel | Earlier r values, handled separately |
|---|---:|---:|---|---|
| 3*2^r+1 | 1 | 1 | 15,6 | none |
| 3*2^r+1 | 2 | 1 | 29,12 | none |
| 3*2^r+1 | 3 | 3 | 201,96 | H=57,105 |
| 3*2^r+1 | 4 | 2 | 209,96 | H=113 |
| 21*2^r+1 | 1 | 1 | 87,42 | none |
| 21*2^r+1 | 2 | 1 | 173,84 | none |
| 2*3^r+3*2^r-1 | 1 | 1 | 23,18,4 | none |

Use P_1=1, P_2=2^a-2, P_3=3^a-3*2^a+3 to check the formulas.
For example 201+96P_2(r-2)=9+24*2^r. Every witness includes literal
adjacency rows and a port, rather than a scalar recurrence alone.
`verify_new.py` independently recovers each entire coefficient vector
from physical graph counts up to its proved degree bound. The three
prefix graphs and 68 literal line instances also pass independent counts;
42 small instances pass direct path backtracking.

`realize.py` emits the actual tournament for every r,s in (R3). The
graphs have O(r+s) vertices. The third coefficient family now has every
exponent except s=2 covered; the constructor rejects that unresolved head.
None of these proofs constructs seeds for every odd c.

## R4. Exact fixed-kernel decisions and their boundary

The search input is the complete strong H<=375 inventory inherited from
[Session173](../session173_modular_clocks_20261103/counts527/RESULT.md),
with completeness grounded in
[THM-4744 / interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md).
The inventory is an inherited FINITE-EXACT input, not recomputed here.
Its LF-normalized SHA-256 is pinned in `head_completion.json`; 8,525 relevant strong
ports were tested in this search. Busch's strong lower bound excludes
orders >=13 for H<=375. The primary source is
[Busch, 2006](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

**All-order reduction, conditional on the inventory's completeness.**
Let a fixed tournament with a marked vertex have nonconstant Hamilton
kernel h=(h_1,h_2,...). Its strong components are totally ordered and
their Hamilton counts multiply. The marked vertex belongs to a unique
strong component S. If S is a singleton, its transitive inflation leaves
the count constant, contrary to the hypothesis. Otherwise every kernel
coefficient is g times that of (S,v), where g is the product of the
counts of the other components. Thus g divides gcd(h_j), g is itself
attainable, and H(S)=h_1/g. Conversely, every such factor and strong-port
kernel gives a witness by ordinal composition. Exterior factors are
decided recursively from the same complete small-H inventory. Unlimited
transitive padding contributes factor 1 and cannot evade the reduction.

Consequently the following are globally impossible **fixed one-port
transitive-inflation kernels**, at all base orders, given that input:

| Family target at parameter r>=1 | Shifts r=a+q tested and excluded |
|---|---|
| 7+14*2^r | q=0,1,2,3 |
| 3+12*3^r+18*2^r | q=0,1 |
| -3+8*3^r+12*2^r | q=0,1,2 |

These are exactly the shifts whose first value is <=375. They are the
remaining heads for c=7*2^r+3 at s=1, c=6*3^r+9*2^r+1 at s=1, and
c=2*3^r+3*2^r-1 at s=2. The conclusion concerns a uniform model;
it does not say these target integers are missing Hamilton counts.
Larger shifts and other graph operations are outside this exclusion.

## R5. Two concrete reasons naive composition fails

**Exact hostile, not a failed seed-stability claim.** Let
B=[8,33,3,6,15,29] and inflate ports v=2,w=0. Its two-port Hamilton
tensor has entries (i,j,h_ij):

    (1,1,23),(1,2,11),(2,1,18),(2,2,12),(3,1,4),(3,2,4).

With w uninflated, F(a)=4*3^a+6*2^a-1. The coefficient of P_2(b) at w
is G(a)=4*3^a-1. At a=1, G=(F-1)/2, so changing b from 1 to 2 takes
23 to 45=2*23-1. At a=2 it takes 59 to 129, not to 117. The exact
failure is G-(F-1)/2=2*3^a-3*2^a, zero only at a=1. Stable delta pairs
do not imply transport of a different Hamilton-doubling interface.

**Same-clock two-port growth (PROVED, using Camion's theorem).** In a
nontrivial strong tournament choose distinct v,w with v->w. Start a
Hamilton cycle at v, traverse it, and append v,w. The resulting quotient
path visits v,w twice and every other vertex once. Therefore h_(2,2)>=1.
Inflating both ports to TT_a yields H>=P_2(a)^2=(2^a-2)^2. It cannot
equal a prescribed O(3^a) family for every a. This rules out a tempting
same-size second block repair of a ternary head. It says nothing about
asynchronous sizes, varying templates, or differently parametrized models.

## R6. A finite fixed bank cannot cover the Hamilton spectrum

The classical Camion input in R5/R6 has a short proof: take a longest
directed cycle in a strong tournament. If it is not spanning, each
outside vertex must dominate all cycle vertices or be dominated by all
of them; otherwise it inserts between consecutive cycle vertices.
Strong connectivity forces an edge from an outside vertex dominated by
the cycle to an outside vertex dominating the cycle. Inserting these
two vertices into any cycle edge makes a longer cycle, a contradiction.
A strong tournament of order >=3 has an initial directed triangle.

**PROVED for the specified grammar.** Fix finitely many nontrivial strong
marked unit-seed skeletons, with at most m transitive inflation blocks,
and allow any positive sizes for which the result remains a unit seed,
then one final tail. The number of distinct resulting Hamilton counts
<=X is O((log X)^(m+1)). Fixed finitely many exceptional graphs or fixed
Hamilton-inflation head templates do not alter density zero. For the
60 cones above m<=5, giving O((log X)^6).

Proof: Camion supplies a Hamilton cycle through any port v. Cut it at v
and repeat v at the end to get a quotient path with two visits to v.
Hence its single-port expansion has H>=H(S)+P_2(a)>=2^a+1.
In a simultaneous expansion this single-port graph is an induced
subtournament. Hamilton counts are monotone under adding tournament
vertices: insert each new vertex into each old Hamilton path at its
first legal position; deletion recovers the old path, so the map is
injective. Thus H(B)>=2^(a_i)+1 for every coordinate i.

For every selected b, splitting any Hamilton path immediately before b
injects it into U_b. Consequently U(O)>=H(B). Since the final output is
1+2^t U(O), every a_i and t is bounded by log_2 X. Counting parameter
tuples proves the assertion. Exterior source/sink padding either leaves
H unchanged or adds its terminal length to t, by R0's core theorem.

This is a genuine obstruction to proving all odd Hamilton values by
these fixed cones alone. It is **not** a density theorem for their odd
coefficients: removing an unbounded power of two from U can change size
estimates. Unbounded ordinal products, changing templates, selected
vertex births, edge reversals, and other recursion grammars are not
included. To escape this result one must justify an operation outside
the stated bank, not just add more parameters inside one cone.

## R7. The seed birth metagraph has unavoidable new roots

Vertices at level n are marked unit seeds (B,O). An upward edge adds one
tournament vertex, leaving all old orientations and old selection bits
unchanged; the new vertex may be selected or unselected. Its reverse is
a deletion that leaves a nonempty unit seed. This is an intrinsic,
directed, levelled graph, not itself a tournament.

`birth_atoms.py` enumerates all ordinary tournaments through seven, all
nonempty O, and every legal deletion. Dynamic programming then tests
whether repeated legal deletions can reach the singleton seed.

| n | Unit seed presentations | No legal deletion | Can descend to singleton |
|---:|---:|---:|---:|
|1|1|1|1|
|2|2|0|2|
|3|6|0|6|
|4|13|0|13|
|5|28|4|24|
|6|67|7|40|
|7|204|43|63|

The four first nontrivial roots are:

| Adjacency rows | Selected mask | H | U | c |
|---|---:|---:|---:|---:|
| [2,24,3,5,13] | 5 | 11 | 29 | 29 |
| [2,24,3,5,13] | 24 | 11 | 47 | 47 |
| [4,1,18,7,11] | 9 | 9 | 32 | 1 |
| [4,1,18,7,11] | 18 | 9 | 32 | 1 |

The last two have c=1, as does the singleton seed, but no permitted
deletion. Thus a coefficient-only quotient loses a concrete construction
obstruction. The 54 nontrivial root presentations through seven cannot
be born from smaller unit seeds by this one-vertex grammar. This does
not make their coefficients unattainable, nor preclude a multi-vertex
move, temporary nonunit state, or edge change. Whether infinitely many
root types are needed is OPEN.

## R8. Correct symmetry quotient and new sequences

`marked_orbits.py` identifies marks only under actual automorphisms of
the ordinary base. It enumerates all degree-preserving permutations and
checks every pair orientation. Therefore it retains the entire marked
isomorphism type and safely preserves deletion feasibility and all
construction responses, unlike quotienting by H or c alone.

For n=1,...,7 the true marked-orbit sequences are:

    unit seeds:           1,2,4,9,22,59,194;
    deletion roots:       1,0,0,0,4,7,43;
    singleton reachable:  1,2,4,9,18,32,53;
    strong unit seeds:    1,0,1,3,8,24,98.

Strong stable r-subset orbits for r=1,...,5 are:

| n | r=1 | r=2 | r=3 | r=4 | r=5 |
|---:|---:|---:|---:|---:|---:|
|1|0|0|0|0|0|
|2|0|0|0|0|0|
|3|1|0|0|0|0|
|4|2|1|0|0|0|
|5|7|3|1|0|0|
|6|20|10|4|1|0|
|7|63|35|14|5|1|

Thus the 95 port presentations reduce to 93 marked port orbits; the
49 pairs remain 49 orbits. The 60 full stable-set presentations reduce
to 58 marked core types. These are exact refinements of ordinary
tournament enumeration, not additional terms of A000568.

The inherited padding transform also applies to marked orbits. If K_m
is the number of strong marked seed orbits, then A_n=sum_(m<=n)(n-m+1)K_m.
Automorphisms fix the ordered exterior singleton components individually.
For r-subsets of stable vertices, combine the core r-subset orbit counts
with choices of exterior vertices; the same binomial convolution as E10
applies. This avoids repeating padding calculations. Only the core
counts through seven are supplied here; higher core counts remain open
computational work.

## R9. The precise remaining proof obligation

The spectrum statement is [Chan–Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf):
all positive odd Hamilton counts occur except 7 and 21, originating with
bof and Royle. This packet does not settle it. A sufficient next mechanism
would give legal, endpoint-decorated descent for every odd coefficient
outside a stated seed set, with a decreasing rank, and handle every
point below its ray threshold except the two forbidden counts.

There is also an exact way to organize what remains. Every odd N>1
has a unique decomposition N=1+c*2^s with c odd and s>=1. If c=1,
the nearly transitive tournament on s+2 vertices already realizes N.
Otherwise write uniquely c=1+d*2^r with d odd and r>=1. Thus the
remaining spectrum decomposes into the disjoint two-parameter grids

    N(d,r,s)=1+2^s+d*2^(r+s), d odd positive, r,s>=1.       (R9)

The only conjectured omissions are (d,r,s)=(1,1,1) and (1,2,2), giving
7 and 21. R3 closes the **entire grids d=3 and d=21**. Increasing r
implements c->2c-1 inside these certified source families; increasing
s uses the appropriate head or tail constructor. For example the first
grid's coefficients are 7,13,25,49,97,... . This unique parametrization
turns the broad coverage problem into filling all odd-d grids, with the
two exceptions explicit. It is an equivalent indexing, not a proof that
a legal operation can change d or that all the other grids are filled.

The progress here is explicit: commuting source-valid moves, no-activation
certificates, seven completed heads, two complete infinite coefficient
lines, a small exact counterexample to interface transport, a fixed-bank
density obstruction, and a birth-root classification through seven.
The [arithmetic companion](ARITHMETIC.md) supplies analogous exact birth
laws and explains which coordinates must be retained; it supplies no
unproved lift from an arithmetic parent pair to a tournament seed.
