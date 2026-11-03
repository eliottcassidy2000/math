# TM1–TM3: connected tournament memory and sharp boundary budgets

**PROVED all-size construction; FINITE-EXACT seed identities and literal
graph controls.** This extends the incoming S172 strong-memory construction
by replacing its single feedback vertex with a transitive feedback block.
It resolves that packet's stated 2-vertex-strong extension and gives every
vertex-connectivity r. It does not prove an unrestricted spectrum theorem.

## Inherited seed and the extra positive information

Use outgoing-neighbour masks on vertices 0,...,4:

    B0=(16,1,3,7,14), B1=(20,1,2,7,14).

The S171/S172 packets prove and check that these strong, nonisomorphic
tournaments differ only in arc {0,2} and have the common unordered
directed-path-cover polynomial

    P(z)=9z+33z²+31z³+10z⁴+z⁵.

For each seed add one vertex beating only its old vertex 0. Let the result
be Q0 or Q1. The previously known Hamilton counts are 17 and 19. The full
profiles, independently checked here, are

    P_Q0=17z+95z²+138z³+73z⁴+15z⁵+z⁶,
    P_Q1=19z+109z²+150z³+75z⁴+15z⁵+z⁶.

The useful strengthening is coefficientwise positivity:

    P_Q1-P_Q0=2z+14z²+12z³+2z⁴.                       (1)

These finite identities can be proved by exhaustive path/subset recurrence;
the separate audit enumerates permutations and cut sets, canonically
deduplicating unordered covers. The distinction between full-profile
equality before the read and positive full-profile difference afterwards
is what lets the read survive a more highly connected wrapper.

## A modular substitution lemma, including its scope

Replace vertices of any fixed directed graph by nonempty directed graphs,
with every cross-arc between modules uniformly governed by the outer graph.
Then each coefficient of the global path-cover polynomial is a polynomial
with nonnegative integer coefficients in the module path-cover coefficients.

**Proof.** For clarity first label the global path slots. In each slot,
collapse maximal consecutive runs in the same module. This gives a finite
ordered skeleton of module labels, with no adjacent repetitions and all
cross-steps legal in the outer graph. If module i occurs q_i times in this
skeleton, assign to those occurrences an ordered cover of module i by q_i
nonempty directed paths. There are q_i! P_(i,q_i) choices. Concatenation
and maximal-run collapse are inverse maps. For unlabelled global paths,
fix a multiset of nonempty skeleton words. Its stabilizer permutes copies
of identical words. This group acts freely on the occurrence assignments
of any fixed collection of internal covers: moving a global word moves
at least one module occurrence, which carries a distinct nonempty internal
path. Thus the contribution is

    [product_i q_i! / |stabilizer|] product_i P_(i,q_i),

with a nonnegative integer coefficient, by the same free action on the
product of permutation sets. Sum over the finite skeleton multisets.
This proves equality and coefficientwise monotonicity directly, without
assuming that conversion from the labelled-slot basis preserves signs. QED.

In particular equal module profiles remain equal under every *uniform*
modular substitution. Arbitrary external arcs to individual vertices are
not covered. The local read deliberately uses such a nonuniform attachment
inside the selected seed; after it is formed, the wrapper is uniform.

## The construction and exact connectivity

Fix r>=1 and M>=2 ceil((r+1)/5). Choose a split 1<=a<M with
5a>r and 5(M-a)>r. For a bit word h of length M, put

    A=B_h(0) ⊕ ... ⊕ B_h(a-1),
    L=B_h(a) ⊕ ... ⊕ B_h(M-1),
    C=TT_r,

and orient the three groups cyclically A -> L -> C -> A. Here TT_r is a
transitive tournament and ⊕ is the ordinal sum. Call this W_(h,r).
Vertices, bits, block positions, the split and local vertex-0 ports are
retained explicitly. There are 5M+r vertices.

**Connectivity and reconstruction.** Every deletion of fewer than r vertices
leaves all three groups nonempty, hence strong: a directed route can pass
around the three-group cycle regardless of internal orientations. Deleting
all of C leaves A⊕L, which is not strong. Furthermore C is the unique set
of at most r vertices whose deletion destroys strong connectivity. Any
different such set leaves part of C, and cannot empty A or L because both
have more than r vertices. Thus the vertex-connectivity is exactly r.

The unique minimum cut C is intrinsic. Delete it; the ordered strong
components are precisely the original strong seeds. Their nonisomorphic
types recover every bit h(i), so the 2^M wrapped tournaments are pairwise
nonisomorphic even without labels. This uses the whole set C as its canonical
mark, replacing the incoming packet's single strong articulation vertex.

## Memory and one-vertex read theorem

All W_(h,r) have the same complete path-cover profile, by the substitution
lemma and equality of the two seed profiles. On the fixed vertex labels,
their arc distance is exactly the Hamming distance of h and h'. Thus any
binary code of length M embeds with its original labelled distance. No
statistical source model is needed for these statements.

To read h(i), add one vertex to that block with the same Q attachment as
above. Its arcs to every other block and to C are the uniform arcs belonging
to the selected module. This is a fixed one-vertex induced extension using
the named local port; no old arc is changed or port reused. Every queried
tournament still has connectivity r and the same unique cut C.

For fixed r,M,a and address i, the response H depends only on h(i), since
all unqueried module profiles coincide. Equation (1) and the substitution
lemma show that the response for bit 1 exceeds that for bit 0. More precisely
its gap is at least

    Delta_i >= 6*9^(M-1).                              (2)

To see strictness and this bound, restrict to Hamilton paths with exactly
one run in each of A,L,C. There are three allowed orders of these groups.
The Hamilton counts of A and L are products of their seed counts; querying
one block changes its factor from 17 to 19, giving a difference of
3*2*9^(M-1) because H(TT_r)=1. Every other skeleton contributes a
nonnegative difference. It is a disjoint count, so no reserve is reused.

This proves **arbitrarily high vertex-connectivity, common complete profile,
recoverable hidden bits, unchanged arc-distance code, and bounded one-vertex
reads simultaneously**. A scalar Hamilton response has a two-element
alphabet for this particular reader, irrespective of its large numerical gap.
No autonomous tournament update or fixed-cost rewrite of the stored word is
asserted; the read preserves the original tournament as an induced subgraph.

## Exact Hamilton formula for the feedback block

Let a_j=j! P_(A,j) and b_j=j! P_(L,j), taking index zero and out-of-range
indices to be zero. Put c_q=q! S(r,q), where S is the Stirling number of
set partitions. This is the number of ordered q-path covers of TT_r:
each vertex subset has one increasing directed path.

Then

    H(W)=sum_(q=1)^r c_q [a_(q-1)b_(q-1)
             +a_(q-1)b_q+a_q b_(q-1)+3a_q b_q
             +a_(q+1)b_q+a_q b_(q+1)+a_(q+1)b_(q+1)].  (3)

For proof, a collapsed Hamilton word must follow A,L,C,A,L,C,... . Fix
its number q of C-runs. The possible A/L-run counts are precisely the
seven displayed pairs; (q,q) has three starting positions, and each other
pair has one. Assign the ordered internal path covers. This also proves
the formula when some terms vanish because a group has too few vertices.

For r=1 this reduces to the incoming S172 formula using H and F(2).
For r=2, M=2, the four twelve-vertex storage tournaments have

    H=152757.

Each address's one-vertex reader gives 585033 or 645915, a gap of 60882.
The full profiles of all four storage graphs, all single-vertex deletions,
and the unique bad two-vertex deletion are checked on the literal adjacency
matrices. Orders 11–14 (r=1–4) and their query extensions are also checked.

## A sharp boundary for the odd-spectrum programme

The new connectivity does not remove the seed's divisibility obstruction.
For B0 and B1, the ordered nonempty path-cover counts are

    9,66,186,240,120,

all divisible by 3. In any uniform outer substitution containing M such
seed modules, every Hamilton skeleton therefore contributes a multiple of
3^M. Hence

    3^M divides H(W_(h,r)).

This remains an excellent memory carrier but cannot supply the prime-count
part of the {7,21} completeness conjecture. An extension toward that target
must break the uniform-module divisibility mechanism, not merely increase
connectivity or code distance. The same port-level skeleton analysis tells
us both why the memory theorem succeeds and why this spectrum route stops.

## TM2. Arbitrary new vertices can remove at most two module factors each

There is a stronger obstruction, even for nonuniform attachments. Start
with any tournament containing M disjoint seed modules B0/B1, each having
uniform arcs to every old vertex outside that module. Other old vertices
and their tournament may be arbitrary. Add q new vertices with arbitrary
arcs to all old vertices and to each other; preserve every old arc. Then

    3^max(M-2q,0) divides the new Hamilton-path count.       (4)

**Proof.** Collapse each maximal run in an old seed module along a Hamilton
path, retaining every other vertex individually. Call a seed module touched
if one of its runs is immediately next to a new vertex in this skeleton.
Each new vertex has at most two neighbours in the path, hence at most 2q
distinct seed modules are touched. For any untouched seed, all run endpoints
meet old vertices or other old modules with uniform cross-arcs. Its internal
choices are therefore exactly the ordered j-path covers for its number j
of runs, a multiple of 3. Touched modules may have arbitrary endpoint
restrictions; their number of internal choices is still an integer. Fixing
the skeleton factors the count over modules; every term has at least M-2q
untouched seed factors. Summing proves (4). QED.

In particular, to obtain H coprime to 3 from such an old source, at least
ceil(M/2) added vertices are necessary, unless old arcs are changed. Merely
giving one new vertex access to all M named ports does not evade the bound:
an individual Hamilton path uses only its two adjacent incidences.

The factor-two budget is sharp, with a concrete coupled reader. Start with
W_((0,0),1), the eleven-vertex storage source in TM1. Add one vertex that
beats precisely old local vertex 0 in both seed blocks (old labels 0 and 5),
and is beaten by every other old vertex, including the feedback vertex.
The resulting twelve-vertex strong tournament has

    H=30431, which is prime.

Its exact adjacency rows are in `tournament_certificate.json`; subset DP
and a separate nonmemoized enumeration of directed Hamilton paths both give
30431, and trial division through floor(sqrt(30431)) verifies primality.
This changes the source through one actual induced extension using two
distinct named ports. It is not claimed to be the first known realization
of that prime count or a repeatable scalar return rule. The two-block read
now couples bits; it is not the single-address reader of TM1.

For three old seed blocks with word (0,1,0), the analogous one-vertex
attachment beating local vertex 0 in all three gives H=12573525, of exact
3-adic valuation one. Thus the surviving factor in (4) is real, even though
the new vertex's adjacency pattern inspects all three modules.

Sharpness extends to every q>=0 and t>=0 with 2q+t>0. Use an ordinal sum
of q copies of the prime extension and t copies of C3[B0,1,1]. In each
prime extension designate only its last reader vertex as new. The old
source has M=2q+t seed modules, and all are uniform modules before these
q new vertices are inserted. The seven-vertex component C3[B0,1,1] has
H=3*9+66=93. Ordinal-sum Hamilton counts multiply, so the final count is

    30431^q 93^t,       v3(H)=t=M-2q.                   (5)

The final tournament need not be strong when this ordinal sum has multiple
components; the lower bound (4) did not assume strongness. Formula (5)
proves its all-size sharpness for M>=2q. It does not prove every target
count can be reached, nor does it give compatible return sources for the
ternary spectrum grammar.

## TM3. New-vertex readers have bounded Boolean interaction degree

Fix the old uniform-module layout, and let each of its M seed modules
independently be B_(h_i), h_i in {0,1}. Fix a q-vertex labelled extension
rule: every new arc is specified independently of the hidden word h, though
its pattern to individual old vertices may be arbitrary. Let R(h) be the
resulting Hamilton count. Its unique multilinear polynomial over Q obeys

    degree R <= 2q.                                    (6)

**Proof.** In the TM2 skeleton decomposition, only modules directly adjacent
to new vertices can be observed nonuniformly. There are at most 2q. Each
untouched module has the same ordered-cover counts for both bit values,
so that skeleton's contribution depends on at most 2q bits. A function of
t binary bits has a multilinear polynomial of degree at most t, by expanding
in the indicator products h_i and 1-h_i. Sum over skeletons to prove (6).
This argument also works for other fixed binary module pairs with identical
full path-cover profiles; their factor-3 divisibility is not required. QED.

For the coupled two-block reader above the exact truth table is

    00 -> 30431,  01 -> 33397,  10 -> 36325,  11 -> 39603,

so its response is

    R(h0,h1)=30431+5894 h0+2966 h1+312 h0 h1.           (7)

Thus one new vertex really can read two bits in one four-outcome integer
response; the nonzero 312 term makes the degree-two bound sharp. The four
outputs are distinct, and their unequal linear coefficients retain block
order. This does not contradict the inherited one-bit lower bound for the
single-address reader, whose response alphabet had only two outcomes.

For all q, take an ordinal sum of q independent two-block coupled readers.
The response is the product of q copies of (7), with top coefficient
312^q, so degree 2q is attained. These sharpness examples need not be strong
after the ordinal sum. The one-vertex M=3 control has zero third mixed
difference, checked on all eight literal source graphs.

Two limitations are essential:

* Low polynomial degree does not bound the number of bits recoverable by
  arbitrary postprocessing of an exact integer response: even a weighted
  linear function can encode every bit. The theorem bounds interaction
  degree, not general information capacity or arbitrary-decoder query cost.
* If two probability laws on the stored word have identical marginals on
  every <=2q coordinates, they give the same expected R for every such fixed
  reader. This is an immediate expectation consequence of (6). It does
  not equate the pointwise responses or their full probability distributions.

For instance uniform even- and odd-parity words of length M>2q have identical
these marginals and therefore identical mean Hamilton responses to every
fixed q-vertex reader. A reader whose raw count is required to take exactly
one constant on even parity and a different constant on odd parity needs
2q>=M: that two-level parity function has multilinear degree M. A nonlinear
decoder of a richer integer response is not subject to that last conclusion.

This sharp locality statement supplies an exact version of the repo's
recurring distinction between matching low-order averages and preserving a
specific source. It is a theorem about the stated modular tournament family,
not an imported concentration theorem for Collatz or a unit-distance gadget.

Dependencies: [S171 seed/read theorem](../session171_knots_information_20261103/information/NOTE.md),
[S172 one-feedback construction](../session172_memory_dynamics_20261103/spectrum/NOTE.md).
No scarce global theorem ID is reserved or claimed. TM1–TM3 are local to this file.
