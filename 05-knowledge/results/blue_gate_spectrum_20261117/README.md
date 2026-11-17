# The transitive blue gate: exact fans, nonplanar lifts, and spectral families

**Proved.** All Hamiltonian frames at the first blue neighbor have
an explicit run grammar. This gives exact endpoint formulas and
deduplication at every size, and genuine blue and black families
in the 7-channel, with verified prime quotient examples. A
pointed six-vertex certificate lifts a literal K3,3 to every
larger order in the **actual** merged metagraph. A separate
[shifted-divisor analysis](DIVISORS.md) supplies all four odd
quotient residues mod 8 but does not prove absolute coverage.
The spectrum conjecture is still open.

## 1. From the transitive node to the first fan

Throughout a line means complement of the nonpath tiles along
a *fixed directed Hamiltonian frame*, followed by isomorphism
and converse identification. This is not the single-arc-flip
graph or the larger switching observer.

The transitive TT_n has one Hamiltonian path, hence just one
literal line out, blue because its tiling is mirror fixed.
For n≥3 its other endpoint Z_n in the natural order has precisely
the upward adjacent arcs i→i+1; all longer arcs point downward.
It contains the directed spanning cycle and is self-converse by
order reversal. The inherited blue gate and tribonacci counts are
in [THM-4731](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md)
and [THM-787](../../../01-canon/theorems/THM-787-transitivity-flow-on-the-metagraph.md).
The early resemblance H(Z_n)=2^(n−2)+1 stops at n=6.

A Hamiltonian frame in Z_n consists of increasing contiguous
intervals, traversed in descending interval order. The adjacent
blocks must not both have length one: exactly then their join
would be a forbidden downward step of length one. Conversely
every such ordered composition gives a frame. The count GF is

    (1−x²)/(1−x−x²−x³).

The initial counts at n=1,2,3,4,5,6,7,8 are 1,1,3,5,9,17,31,57.
For n≥4 Z_n is rigid: the only degree ties are at the two end
pairs, and the arc within either pair fixes their order.
Thus its frames count literal incident endpoints with no
automorphism division. The unique reversing isomorphism is
the natural reflection, so a frame is blue exactly when the
composition is a palindrome. Its GF, including formal constant
term 1, is

    (1+x)(1+x^4)/(1−x²−x^4−x^6).

For completeness let E(t) count the empty composition plus ones
ending in a part ≥2 and with no adjacent ones, and H(t) count
all such compositions. Then E=(1−t)/(1−t−t²−t³) and
H=(1−t²)/(1−t−t²−t³). The palindrome may have no center,
a singleton center (both need E), or a center of length ≥2.
So its GF is (1+x)E(x²)+x²H(x²)/(1−x).
Counts here are incidences, not necessarily different nodes:
mirror, repeated endpoints and loops can identify them.

## 2. The exact endpoint of each frame

For a composition lambda, partition the natural labels into its
contiguous blocks, in increasing order. Let N(lambda) be TT_n with
two arcs reversed across each boundary between a lower block L
and the next upper block U:

    (last L, first U) and (first L, last U).

Indeed the complement of Z_n first makes all nonconsecutive
arcs upward and adjacent arcs downward. Restoring the frame
arcs restores all adjacent arcs *inside* blocks and introduces
one long descending connector across each boundary. This
is a direct all-order labelled formula for all fan endpoints.

An all-order exact deduplication needs only a small permutation
test, not the whole tiling cube. Start from TT_n's distinct
indegrees. In N(lambda), each interior block's net changes
cancel at its endpoints. In the first block add one at each
endpoint (two if it is a singleton), and subtract at the last
block's endpoints. So at most four vertices were modified,
and a tied bucket has at most one unmodified vertex.
At most eight labels participate in ties. Canonicalize by
indegree and permute only tied buckets (at most 8! trials);
also take the minimum with converse. Distinct normalized
keys are an exact all-order node count for the fan, with
multiplicities and colors obtained from the same compositions.

The one-block frame returns to TT_n. With two parts at least
two, exchanging the adjacent labels at the boundary removes
the short reversal, leaving just the long return n→1. All
these frames land at **the same** strong self-converse F_n,
with H=2^(n−2)+1. Among these n−3 frames the palindrome
exists iff n is even and is unique. With one end part equal
to one, the result is a singleton joined to F_(n−1) or its
converse, in one common merged node with two black incidences.
It is the only nontransitive nonstrong node in this fan.
For three or more parts the original frame is open: its
endpoints across the gaps are not consecutive in Z_n.
Complementing the nonpath endpoint arc closes it; thus the
new tournament is strong by its directed spanning cycle.
For n≥7, H(Z_n)<H(F_n) by the tribonacci recurrence, so
F_n is genuinely distinct from the hub, not just a loop.

**Next fan at F_n.** Apart from the natural order, every
Hamilton frame in F_n uses n→1 and has form
S increasing, n, 1, complement(S) increasing, where S is
any subset of the middle labels. For these frames the next
partner is reverse-transitive except for the ascending edges
along its two strands (S,n) and (1,complement(S)).
The natural frame instead gives Z_n with the long return
reversed. The natural frame is blue; in the subset frames
blue means reflection exchanges S with its complement.
Thus there is one blue incidence when n is odd and
1+2^((n−2)/2) when n is even (for n≥4, F_n also
has trivial automorphism group by its degree pairs).
These are next-step descriptions, not claims that every
endpoint lies at exactly distance three.

For a general N(lambda) there is also an exact marked
next-frame compiler. Choose the ordered list of its descending
arcs used by a path. These are selected from the explicit
2(number of parts−1) reversals, not from all arcs.
Between them the runs must be increasing. Their start/end
vertices are forced by the adjacent descending arcs (leave
the two outer endpoints free); incompatible repeats are
discarded, and a coincident start and end forces a singleton
run. Assign each unused vertex to a run if its label lies
in the corresponding interval, then sort within each run.
Reject a run if an ascending consecutive pair is one of the
reversed arcs. Equivalently, inclusion–exclusion over the
bad (run, reversed pair) events pins both endpoints to that
run and forbids its color on intermediate vertices. Reading
the valid runs in order yields each frame once. This compiles
third-step incidences for every second-step endpoint without
asserting a closed count of their merged images.

## 3. A three-run family with actual 7p witnesses

For lengths (a,m,c) all ≥2, swapping the two disjoint
adjacent boundary pairs puts N into a TT order with only
two interlacing backward arcs C→A and D→B, A<B<C<D.
The numbers of vertices in the three gaps AB, BC, CD are
l=a−2, m, r=c−2. Its exact count is

    H(l,m,r) = 1 + 2^(m+1)(2^l+2^r) + 2^(l+r+1)3^m
                    − [l=0]2^r − [r=0]2^l.

There is a literal cycle-packing proof. One odd cycle using
only the first return can choose an odd subset of l+m+1
interior vertices, giving 2^(l+m) cycles; the second alone
gives 2^(m+r). A cycle using both returns consists of
two disjoint forward strands, A→D and B→C, with odd total
interior count. Each side vertex has two choices, each
middle vertex three choices (unused, first or second).
With M=2^(l+r)3^m its number is
(M−[l=r=0](−1)^m)/2 by parity. A pack of two disjoint odd
cycles uses strands A→C and B→D with both interior counts odd.
Its count is
(M−[l=0]2^r−[r=0]2^l+[l=r=0](−1)^m)/4,
using the two parity tests. The signs cancel exactly upon
adding twice the single-cycle terms and four times the
two-cycle packs in the [odd-cycle read](https://arxiv.org/abs/2307.05569).
In particular the ternary digit isn't a numerical analogy:
it is precisely a middle vertex allowed on either strand
or in neither.

This produces a physically specified family two edges from
the transitive class. The frame is blue iff l=r. For example:

| Run lengths | Path count | Quotient |
|---|---:|---:|
| (3,2,6), black | 721 | prime 103 |
| (6,8,6), blue | 3,375,617 | prime 482,231 |
| (5,14,8), black | 4,900,119,553 | prime 700,017,079 |

Keeping l=1,r=4 and m=2 mod 6 gives infinitely many actual
counts 1+36·2^m+64·3^m divisible by 7 and coprime to 3.
Keeping l=r=4 and m=2 mod 6 gives a separate blue family
with the same divisibility, including the prime quotient
in the table. This does not prove infinitely many primes
on either family. All three-run outputs to X number at most
O((log X)^3), since the formula bounds n by O(log X).
They cannot fill the odd spectrum on their own.

The Catalan pair 8,9 enters exactly at middle length two:
with l,r>0 the count is 1+8(2^l+2^r)+18·2^(l+r).
The 8 is a strand term and 9 is the two middle ternary digits;
this is not using Catalan as a universal selector. A different
hostile: using the literal primitive triple as run lengths
(3k,4k,5k), the whole family has H≡5+5·2^k mod 7, never zero.
This follows from 2³=1 and 3⁴=4 mod 7 in the formula.
The right way to use Pythagorean data is instead to build
two physical **square gaps** in a cyclic cover read, as in
the [verified shifted-square example](DIVISORS.md).

## 4. The planar-shadow obstruction persists at every order

The projective-plane interpretation of the five-label
metagraph is a useful base, not a planar drawing of the
actual simple graph. For n=5 the prior five-label calculation gives a literal K5;
[replay.py](replay.py) regenerates it directly from the 64 tilings.
Here is a new marked lift that gives a literal K3,3 at
**every** n≥6, not just nonplanarity of a switching shadow.

The [replay](replay.py) has three left six-vertex strong
tournaments A_i and three right strong tournaments U_j
with a distinguished end v_j. The six classes are distinct
even after converse. For every i,j an actual Hamilton
frame of A_i has complement U_j with its *last* vertex
at v_j (such paths and pointed identifications are checked
independently). No U_j is isomorphic to Z_6.

At order 6+m, put the left class at A_i followed ordinally
by a transitive tail of length m. Extend every marked wire
by that natural tail. The complement is U_j together with
Z_m, where every tail vertex beats every old vertex
except for the one old arc v_j→s_1 to the start of the tail.
It depends only on the right pointed class. It is strong
for m≥1: U_j reaches v_j and hence the tail path; the tail
reaches U_j. So no right endpoint can collide with a left
one. The left SCC core recovers A_i.

For m≥7 the old vertices of the right tournament are exactly
the six vertices of outdegree <6. Every tail vertex has
outdegree ≥6, and the distinguished v_j is the one old
vertex pointing outward. In the converse, the outdegree <6
vertices are instead the top six tail vertices in the
original natural order; their induced graph is a copy of
Z_6 or its converse. This cannot be one of the U_j.
Hence the three right endpoints are distinct even after
converse. The remaining m=0,1,…,6 are checked by exact
canonical labels. All nine edges lift and stay distinct.
Together with the literal five-label K5, this proves
the simple merged supports are nonplanar at **all n≥5**.

This differs from the prior full switching observer, where
frame eligibility has already been forgotten. The useful
“shadow” here is the marked endpoint v_j: deleting it
from the description loses the recipe for lifting a wire.
Neither this obstruction nor an embedding on another
surface classifies all metagraph paths.

## 5. Finite checks of how the fans spread

[The root producer](rooted.cpp) enumerates the Hamilton
frames of every class, assigns each tiling to its class
by the [tiling dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md),
then glues converse through the mirror and takes
loopless simple supports. It does not conflate line
multiplicity with breadth. [Receipt](rooted.out):

| n | Sizes of consecutive BFS shells from TT_n |
|---:|---|
| 3 | 1, 1 |
| 4 | 1, 1, 1 |
| 5 | 1, 1, 5, 3 |
| 6 | 1, 1, 7, 22, 3 |
| 7 | 1, 1, 15, 149, 89, 14, 3 |
| 8 | 1, 1, 29, 813, 1529, 921, 231, 3 |

All nodes are reached in the displayed orders; this is
not a proof of unrestricted connectivity. For n=8 there
are 437,981 nonloop supports; adding the 235 loop supports
recovers the union cardinality of the independently
reported blue and black supports in
[THM-843](../../../01-canon/theorems/THM-843-n8-metagraph-flow-and-defect-shell.md).
The phase reachability there is a different question from
the ordinary graph distance here.

Run `python replay.py` for the literal certificates. For
BFS, build `g++ -std=c++17 -O3 rooted.cpp -o rooted` and use
`python run_rooted.py /path/to/repo ./rooted`. The complete
class input through eight is the archived catalogue; small
orders through six are regenerated from degree-canonical
one-vertex extensions.
