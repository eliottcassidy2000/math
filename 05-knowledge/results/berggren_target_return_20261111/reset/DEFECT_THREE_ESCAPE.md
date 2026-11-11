# A least-defect return at1997 and an actual+4 escape

**PROVED / exact finite certificates:**1997 has minimum positive root
defect3, attained at an indegree-one root; no branching root of positive
defect at most6 occurs at1997,5993,6013; an infinite defect3 source family;
an all-height exclusion of5993 in every two-triangle deletion chart.
**FINITE-EXACT:** actual1997→5995→ready5995 constructions with retained
physical names. **OPEN:** a uniform repeatable target-preserving return
from this new interface, and Chan–Pak.

## 1. The smallest escape uses one in-neighbour

The complete branching classification in
[defect/PROOF §§4–6](../../uniform_return_metagraph_20261111/defect/PROOF.md)
contains52 transitive-residual rays,8 isolated objects, the two three-port
rays7·2^k−5 and9·2^k−3, and five triangle-residual families

    36·2^m−3, 72·2^m+3, 60·2^m−3, 156·2^m+9, 168·2^m+9.

The first rays have k>=1; the last five have m>=0 and m+1 possible
triangle positions. Exact integer inversion of every formula excludes
1997,5993,6013. The triangle-residual formulas are all divisible by3,
whereas these three targets are not. This exhausts ALL strong roots
with at least two in-neighbours and0<b<=6. The d=1 case is deliberately
outside that theorem.

The same exclusion holds for weak branching sources at these values.
The SCC factorization in [OBSTRUCTION O1](../../ready_return_pentagon_20261103/OBSTRUCTION.md)
forces the root into the first SCC when0<b<H. The suffix count divides
both b and H. For b<=6 its only possible nontrivial odd factors are3,5,
and none of the three targets is divisible by either. Thus the suffix
is transitive and the strong root component has the same H and b.

There IS an indegree-one defect3 source of1997:

    B=(2046,2044,2024,2032,4068,4032,
       3968,3840,3584,1024,2048,527),

with marked root0 and unique in-neighbour11. The five-core compiler
describes it as coreC3=(2,4,1), cut3, flank addresses(A,B)=(3,47).
Direct counting givesH1997 and defect3. Since the all-order inverse
excludes every1-/2-ready1997 representative, this is the minimum
POSITIVE defect over all tournaments and all roots. Defect zero is a
different boundary and is not excluded by that statement.

### An infinite source family with the same interface type

Inflate vertex1 of B to TT_m, m>=1. Its exact count is

    H=1134·2^m−271,

with complete module coefficients(A1,A2,A3,...)=(1997,1134,0,...).
Every source is strong and has the same marked-root defect3. Equivalently,
take left address A=2^k−1, k>=2, and fixed right address47; then
H=567·2^k−271.

The mechanism generalizes: put any fixed core C between transitive flanks,
add an apex, and make the first TWO flank letters O. The first flank
vertex r has only the apex as in-neighbour. After deleting r, every path
starting at the apex must enter the first remaining flank vertex and then
traverse the base; there are exactlyH(C) such paths. The indegree-one
identity therefore gives defectH(C), independent of the other cuts.
With a nonempty right flank ending I the graph is strong. Thus fixed
core counts1,3,5,9 supply corresponding indegree-one interface families.
This is not a branching-defect construction in disguise.

### Actual same-count return from the previously constructed1997

Let S be Q_1 from TARGET_RETURN after its two standard births. It has13
vertices and H1997. Retain

    (0,1,2,3,4,6,7,9,10,11)

and map them to B vertices

    (11,2,5,6,7,8,9,10,0,4).

Delete old5,8,12, then insert B's missing vertices1,3. Every retained
arc agrees. The exact trace is1997→241→619→1997. The pre-excursion marked
root old10 becomes B root0, and its original in-neighbour old0 becomes
B apex11. Ten of thirteen current vertices survive. The return reduces
the needed positive defect from431 to the optimal3.

## 2. Content one still does not give the exact+2 target

The standard two births at B root0 give5997=3·1997+6. Delete only the
old apex11 conceptually, leaving the actual residual order

    C3(0,12,13) → vertex1 → C3(2,3,4) → TT(5,...,10).

The two original directed triangles retain their listed cyclic orders.
The triangle-filtration theorem gives first extension content3 because
they are separated. Therefore changing only the old apex's incidences
cannot attain5993 or5995.

Move singleton gaps past triangles using the named transport from
[GAP_TRANSPORT](../fibres/GAP_TRANSPORT.md). If an internal gap remains,
the same divisibility exclusion applies. If the two triangles become
adjacent, every remaining residual has form

    TT_a → C3 → C3 → TT_b.

The exact binary inverse was run on all64 cuts of the six-vertex core

    C3→C3=(58,60,57,16,32,8)

with arbitrary flank addresses, without a length cutoff. It finds NO
5993 realization, at ANY size. This is a complete all-height exclusion
for an apex whose deletion is a chain of exactly two triangles and
singleton gaps. It permits arbitrarily many extra singleton vertices.
It is not an exclusion of general5993 tournaments, which exist.

For5995 the normalized inverse has exactly18 presentations:

    cut in {25,26,28,41,42,44,49,50,52},
    (left,right)=(4,9) or(9,4).

They all have14 vertices and are strong. Normalized presentations are
being counted; source/sink padding can give larger weak representatives.
The nine cuts are related by independent rotations of the two triangles.

Completeness comes from the inherited positive count formula. Its AB
coefficient is positive, so min(A,B)<=sqrt((N−h)/L); the other coordinate
is solved by the exact binary inverse of d*g(B)+cB+e(B mod2). Uniform
axes use the inherited trailing-one normalization. No finite word-length
search is silently substituted for the all-height inverse.

## 3. A named zero-deletion+4 repair

Start from the14-vertex standard5997 child of B. Keep the old apex at11.
Move residual singleton1 left across the new triangle, then move5 and6
left across both triangles. These are five successive three-arc gap
transports, retaining every name and both triangles' internal arcs.
The final residual order is

    (1,5,6, 0,12,13, 2,3,4, 7,8,9,10).

The complete full-tournament traces, with the old apex still fixed, are

    5997,7727,10001,9421;
    9421,12215,12083,8901;
    8901,11163,11437,8047;
    8047,10477,10333,7671;
    7671,9873,9997,7237.

In that residual order, give the retained old apex cut3273. Only ONE
apex incidence differs, and this last reversal gives7237→5995. The final
graph is strong and has the exact desired3·1997+4. Every one of the14
vertices survives. These macros are not individually H-neutral in the
full graph; their prices are stated, and the peak is12215.

The old root0 and its original in-neighbour11 survive, and11 still beats0.
However0 now has several in-neighbours, so the old one-port interface is
lost: its profile is(H,S,q,F,b)=(5995,339,1897,19759,13764).
Neither orientation of this particular endpoint is1-/2-ready.

There is an explicit same-count return to a1-ready5995 representative
while preserving these SAME root and original in-neighbour names. Use
the one-start parent of order13 and rank578 in the inherited grammar,
and append its ready newborn. Retain current labels

    (0,2,3,7,8,9,10,11,12)

and map to that ready target's labels

    (13,4,5,8,9,10,11,0,3).

Insert the five missing target vertices. The trace is

    5995→107→253→545→1085→2165→5995.

The final marked root13 is old0, and its sole in-neighbour0 is old11.
It now admits the standard+2 step to17987. This is a finite, actual
defect3→general→defect1 return with explicit physical names. It does
not make the whole interface class closed under arbitrary further digits.

The complete ready inverse additionally gives:5995 has two strong1-ready
classes and one strong2-ready class;6013 has no1-ready class and exactly
one strong2-ready class. The parent bound at6013 is15 because Busch's
order16 minimum3125 exceeds(6013+2)/2. The inherited ready2 minima at
target orders15,16 complete the same all-order cutoff used in TARGET_RETURN.

## 4. Uniform lifts tested, with their failure coordinate

For the two5995 flank patterns, inflating any singleton flank position
gives dyadic slopes among2538,2718,3150,5067,5427. None equals3402, the
slope required for3(1134·2^m−271)+4. The triangle rotations preserve
these flank responses, so this covers all18 normalized charts.

The more local partial transport menu was also completely checked:
all14 nonempty subsets of the three crossings in either direction,
followed by every8192 old-apex cut. It finds no5993 output. Its18
5995 outputs are weak and have complete module coefficients
(5995,4544,936,0,...), hence

    H_m=936·3^m+1736·2^m−285.

Their positive third-run channel prevents the proposed uniform target
law. Equality holds only at m1: subtracting3402·2^m−809 gives
936·3^m−1666·2^m+524, zero at1 and strictly positive thereafter by
f(m+1)−2f(m)=936·3^m−524>0.

Thus finite source return, gcd repair, and a matching first count have
all been separated from uniform target preservation. The next useful
operation must change more than singleton gap placement and apex cuts,
or supply a different exact source family.

## Reproduction

Run `low_defect_escape.py`, `two_triangle_inverse.py`, and
`gap_return_probe.py`, with Python `-B`. Matching JSON receipts retain
actual graphs, full maps, every core signature, the all-height inverse
answers, and each stated finite universe. The inherited ready inverse
receipts now include5995 and6013. The positive source counts and all
literal output matrices are independently counted by subset DP.
