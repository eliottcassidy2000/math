# Two routes for exact Chan–Pak: induced-gap rigidity and a three-slot descent

The two approaches answer different questions. The first rules out a
broad strategy for exact repair; the second builds a new exact inverse
with arbitrary tournament sources. Neither proves the
[Chan–Pak spectrum conjecture](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).

## 1. A fixed positive repair cannot be an induced extension at large absolute height

Let S be any induced subtournament of T. Allow the empty tournament
with Hamilton count 1 for the statements about packings and factors.
By the [Grinberg–Stanley identity](https://arxiv.org/abs/2307.05569),
H(T) sums 2^j over packings of j disjoint directed odd cycles of
length at least three, with all remaining vertices singleton.
Every packing on S survives in T, so H(T)≥H(S).

If a vertex v outside S belongs to a strong component C of order c≥3,
it lies in a cycle of each length 3 through c. One elementary proof:
start with a triangle through v (otherwise its in-neighborhood beats
its out-neighborhood, contradicting strength). Given a proper cycle
through v, a vertex with mixed contacts to it can be inserted into
a cycle edge. Otherwise outside vertices lie in the two pure
dominance classes. Strength ensures a crossing from the class
beaten by the cycle to the one beating it. Insert that pair in place
of any cycle vertex other than v, increasing length by one. This
is the usual vertex-pancyclicity mechanism.

Each odd-length cycle through v in C gives a distinct new packing
on C, each weighted 2. If H(T)−H(S)=2g>0 and we write the strong
components of T in ordinal order, then at least one component C
contains a deleted vertex on such a cycle. Let b_i be the H-values
of those components and a_i the H-values of their intersections
with S. The products are H(T) and H(S), and a_i≤b_i. For the chosen
component i,

    2g ≥ (b_i−a_i) ∏_(j≠i) b_j,
    b_i−a_i ≥ 2 floor((|C_i|−1)/2).

It follows that |C_i|≤2g+2 and ∏_(j≠i)b_j≤g. In particular

    H(T) ≤ g(2g+2)!.

This treats **any number of new vertices** and arbitrary existing
components. For gap +2, all other components have count 1, and the
chosen strong component has order three or four. Its H is 3 or 5.
Thus the only induced +2 increments are 1→3 and 3→5. This is not
an obstruction to reversing old arcs. It tells us precisely why
a unit-repair compiler must change old contacts, use nonnested
companion graphs, or work with unequal and growing steps. The
crude factorial bound already forces the size of a positive
induced increment to grow unboundedly with the absolute endpoint.

## 2. Invert a genuine three-slot gluing, not just a fixed-width code

The [return-gate compiler](../chanpak_power_three_return_gates_20261118/README.md)
used transitive twins. Replace that twin chain by an **arbitrary**
nonempty tournament B, maintaining all its uniform cross-arcs
to the old core and making all its vertices beat the final mark.
Let P_j(B) be its number of ordered covers by j paths, empty slots
allowed. Write h=H(B). For the exact construction, use one of
the following marked cores with outgoing bit rows and apex outmask:

    PLUS      D=(2,20,9,3,13), q=13
    MINUS     D=(2,4,1,7,15),  q=4
    MINUS17   D=(30,156,8,128,12,671,959,277,575,159), q=31.

All arcs B→D go to the q-marked contacts, the others point back.
The attached apex has the same contacts to D and loses to B.
Before D and B put a transitive chain of width m≥3 that beats
both modules. The apex beats the first two spine vertices and
loses to the others. The resulting three hosts are strong,
and have respectively the following exact Hamilton counts:

    PLUS:      h + 27·2^m P₃(B)
    MINUS:     2^m(18P₃(B)+12P₂(B)−3h) − h
    MINUS17:   17·27·2^(m+1) P₃(B) − h.

Proof: in each marked deletion core the crossing graph has a
single shared return gate. In one slot, the skeleton has at
most two visits to B; in the marked two-slot read W there
are at most three. The same literal skeletons that were filled
by onto assignments to transitive runs can be filled with
ordered nonempty path covers Q_j(B), independently of the
interior of B. This gives the exact marked profiles

    core        H(deletion)       W
    PLUS         12P₂−h            36P₃
    MINUS        2P₂+h              24P₃+16P₂−4h
    MINUS17      90P₂+h             1224P₃.

Specifically P₂=2Q₁+Q₂ and P₃=3Q₁+3Q₂+Q₃. The coefficients
are fixed by transitive inputs of widths 1 and 2 for the
one-slot read, and 1,2,3 for the two-slot read (triangular
substitution). Including the apex as the last twin appends
an ordinal singleton to B, multiplying P₂ by 2 while leaving
h fixed. Thus the signed defect κ−2H(deletion) is +h for PLUS
and −h for the other two. The left-spine read is
3·2^(m−2) W plus this defect, by the previous marked-apex
formula. Each vertex of the old core reaches the apex from
its fixed contacts, and the apex via the spine reaches every
vertex; the arbitrary orientations inside B do not affect
that argument. Literal subset-DP controls in the replay use
nontransitive sources, not just the original transitive case.

For PLUS the inverse at prescribed N is particularly clean.
Since P₃(B)≥3h by using a single nonempty slot, only
positive odd h≤N/649 are candidates when m≥3. Moreover
P₃(B) is a multiple of 3 by its odd-cycle expression (every
monomial has positive degree) and is odd by the same expression
(the empty packing contributes an odd number and all
others even). Hence the candidate h must satisfy

    h≡N mod 648,     m=v₂(N−h),
    P₃(B)=(N−h)/(27·2^m).

For any chosen h there is a finite exact test of the last
profile equation. Remove singleton strong components; the
remaining components are bounded in total order by h.
For example the triangle argument gives
H(C)≥1+2|C|/3 in a nonsingleton strong component, and
multiplying the component counts bounds both the individual
sizes and their number. Enumerate these finite SCC profiles.
Restoring t singleton components multiplies P₃ by 3^t
and does not change H. Check the exact remaining power of
3, and reconstruct only if equality holds. If the normalized
core is empty, require at least one restored singleton to
keep B nonempty. This is a
literal target-decreasing profile selector for its own
language, not a universal assertion that every target has
one of these profiles. Knowing h alone is not enough.
For prime N a shared nontrivial divisor of h and P₃ also
rules out that branch.

Two literal nontransitive controls are:
B=(2,4,9,3), (h,P₂,P₃)=(5,32,117), m=4 gives the prime
50549; B=(28,5,32,86,38,75,23), (101,1312,8823), m=3
gives 1905869=7·272267 with prime quotient. These show the
gluing enters both prime axes, not that it covers them.

## 3. What a repeatable descent would still need

The preceding map on arbitrary sources has a real range
restriction even before any congruence test. Iterating any
of these three reads from a **finite** seed bank, with
adaptive m and arbitrary ordinal singleton padding before
each step, produces only X^(o(1)) possible H-values up to X.

To prove this, for an input of order n, each Hamilton path
splits into three ordered optional paths in binom(n+2,2)
different ways. The results determine the original path
and both cut positions. Thus P₃≥binom(n+2,2)H. Ordinal
padding of b vertices additionally multiplies P₃ by 3^b.
At step i the order without padding has already grown at
least linearly in i, so each of the three displayed maps
multiplies the count by at least 2^(m+b) i². The constant
coefficients above cover the first step as well.
A resulting count ≤X therefore imposes d=O(log X/log log X)
on the number of steps and Σ(m+b)≤log₂ X. For each d,
there are at most a fixed seed factor times
3^d binom(floor(log₂ X)+2d,2d) parameter schedules.
Singleton placements do not give extra profiles: ordinal
padding has the same effect on every P_j, and uniform
module substitution reads only the cover profile.
The count of schedules is exp(o(log X)), as claimed.

The next slot genuinely cannot be reconstructed from the old
three-count state. For PLUS with m=3, P₂ of the output reads
Q_j(B) through j=4 with positive top coefficient. Indeed
the crossing pattern from B's out-neighbors to its in-neighbors
has two independent edges, D0→D1 and apex→first spine vertex,
and every other crossing uses the old gate D1 or the apex.
The [matching-based return theorem](../chanpak_power_three_return_gates_20261118/README.md)
then gives exactly four visits across two slots. The following
two order-15 inputs are ordinal sums of the row-mask factors:

    B1: (8,9,3,100,15,23,55)     (144,5,161,7,78,27,47,122)
    B2: (32,97,19,7,11,28,61)    (112,97,67,7,14,28,56,127)

Their first three P reads agree, namely
(17325,5105184,317017071), but P₄ is 8535609600 and
8540125440 respectively. The small factor reads are replayed
and multiplied. Triangular conversion from P to Q preserves
equality through Q₃ and the difference in Q₄. So the two
outputs have the same H but different P₂, an actual state
collision rather than merely a formal dependence on an extra
variable.

The positive lesson is to change which kind of state can
descend: add a side input or a low-cost independent read,
rather than force every target onto one scalar recursion
that raises the slot cost at each stage. The fixed-gap
theorem and the slot-cost theorem give two different tests
for candidate new moves. They exclude specific recursive
grammars, not Chan–Pak.

Run `python3 verify.py` for three literal module reads,
the two prime-axis controls, and induced-gap controls.
