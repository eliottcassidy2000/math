# Fixed defect: a finite grammar for every branching stratum

**PROVED IN NOTE, independently reviewed; small automata FINITE-EXACT.**
No all-height source-return or Chan–Pak theorem is claimed.
The new scope is all fixed positive ready defects with at least two inports,
including nontransitive residuals; the older atlas treats two inports with
transitive residual and neutral start/delete response.

## 1. Inheritance and the positive quantity

For an actual strong rooted tournament (T,r), put h=H(T), q=H(T-r),
s=s_T(r), F=E_r(2), and b=F-h. Let I=N^-(r) and d=|I|.
For U=T-r, W_ac(U) counts two nonempty ordered spanning paths, the first
ending at a and the second starting at c. The inherited identities in
[ADAPTIVE section3](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md#3-the-defect-detects-indegree-and-the-legal-inverse)
and [READY section1](../../two_entry_funnels_20261103/READY.md) give

    b = sum_(a in I) s_U(a) + sum_(a!=c in I) W_ac(U).             (1)

The first sum is q-s: prefixing r accounts for exactly those deletion
paths whose first vertex lies outside I. In particular (1) is a sum of
counts of specified objects. A capped algorithm must compute this sum
directly, not subtract two separately clipped large path counts.

The inherited pair-capacity bound is

    W_ac(U)+W_ca(U) >= 2 H(U-{a,c}).                             (2)

It is a bound on the PAIR, not each individual W term. An individual
ordered term may be zero. This distinction is essential for weak deletions.

## 2. A bounded residual kernel for every branching defect

Suppose b is fixed and d>=2. Write R=T-({r} union I). Induced insertion
monotonicity gives H(U-{a,c})>=H(R), including H(empty)=1. Summing(2),

    d(d-1) H(R) <= b,       d(d-1)<=b.                          (3)

Thus d is bounded and H(R)<=floor(b/[d(d-1)]). In the unique ordered SCC
decomposition of R, every nontrivial strong component has at least three
Hamilton paths; their product is H(R). Each such component has order at
most its path count, and there are at most floor(log_3 H(R)) components.
Consequently only finitely many ordered nontrivial strong kernel lists
can occur. An elementary bound on their total vertices is

    K floor(log_3 K),     K=floor(b/[d(d-1)]),                  (4)

with zero nontrivial vertices when K<3. Sharper inherited Busch bounds
can replace this deliberately crude bound. Arbitrarily many singleton
components remain in each gap of the ordered kernel list.

This is the same finite-kernel/variable-gap decomposition used in the
[capped metagraph lane](../metagraph/FAMILIES.md); the important new map is
from a bounded ROOT DEFECT to a bounded PATH COUNT after deleting r and
all its inports. The original h itself can grow without bound.

### A sharper exponential bound from the older capacity atoms

[THM-4167, positive exposure capacities](../../../../01-canon/theorems/THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy.md)
expresses c_ac(U)=W_ac(U)+W_ca(U) as a positive TAGGED odd-cycle-packing
sum. Take one tagged packing on I plus its probing ear and an independent
ordinary packing on R. Their disjoint union is a valid tagged packing
on U plus that ear, and the weights multiply. The map is injective since
restriction recovers both factors. Cross-edges do not affect these cycles.
Hence c_ac(U)>=H(R)c_ac(I), and

    b >= delta >= H(R) [Q_I(2)-2H(I)]
                 >= H(R)(2^d-2).                              (6)

Here Q_I(2) is the two-ordered-bin cover count. The bracket sums
H(I[S])H(I\S) over every nonempty proper subset S, so it is at least
2^d-2. Thus d<=floor(log_2(b+2)) and the residual bound improves to
H(R)<=floor(b/(2^d-2)). The older quadratic bound remains valid.

This degree bound is sharp for every d>=2: take the cyclic composition
r->TT_k->TT_d->r, k>=1, with I=TT_d. The backbone is forced into the
first path of every ordered endpoint-cover pair; q-s=0 and b=2^d-2.
For d>=4 equality in the bracket forces I transitive, since every
three-element subset is proper and a directed triangle has H3 rather
than1. At d3 both TT3 and C3 attain the bracket6. This is a new root
corollary of the inherited tagged-capacity theorem, not a new OCF formula.

## 3. Effective regular languages, with exact isomorphism quotient

Fix d, an ordered residual kernel list, literal representatives of those
cores, all their incidences with the d named inports, and the tournament
on the inports. Treat those finitely many vertices as feedback vertices.
The remaining vertices are the ordered singleton gaps. They form a
transitive backbone, with a finite letter recording their incidences
with the inports and their gap position. Incidences with the residual
cores are forced by SCC order. The root beats the backbone and all
residual cores and loses to every inport.

For each path object in(1), list the feedback vertices in their path
order and, for two paths, choose the split between their lists. Every
backbone vertex belongs to exactly one increasing run between consecutive
feedback vertices or at a path boundary. The inherited
[finite-feedback circuit](../../session176_tournament_information_20261104/spectrum/PROOF.md#4-a-growing-nonuniform-interface-circuit-that-counts-all-orphans)
extends directly to these specified path boundaries: each run is empty,
nonempty with its last vertex pointing to its right boundary, or nonempty
with the wrong last incidence. A first run vertex must satisfy its left
boundary. A final empty internal run requires the feedback-feedback arc.
Forbidden endpoint gaps are kept empty. This is a finite nonnegative
integer weighted automaton counting precisely (1), with no orphan loss.

Cap every register at b+1. Truncation commutes with nonnegative sums and
products, so the exact equality test b is retained. The capped vector
therefore gives a finite deterministic automaton. Strongness is regular
as well. A tournament fails strongness exactly when a nonempty proper
vertex set dominates its complement. In this encoding, its backbone
vertices must form an initial segment. A finite automaton can guess its
subset of the fixed feedback vertices (including r), guess the backbone
cut position, and check every required cross-incidence before and after
that position. It also tracks that both sides are nonempty. Complement
this regular failure language to recognize strongness. This is an effective finite
construction; no uniform polynomial state bound in b is asserted.

The residual SCC list is intrinsic after removing the intrinsic set I.
Any rooted isomorphism preserves its order and every singleton position.
For each fixed list the finite group

    S_d times product_i Aut(C_i)

acts on the inport template AND contact letters. There is no permutation
of equal core occurrences at different SCC positions. This gives exactly
the rooted-isomorphism quotient of the representations, not just a
presentation count. For each group element, fixed templates and fixed
letters impose regular restrictions. Burnside's finite average therefore
shows that the ordinary rooted-class generating function, by vertex
order, is rational for every fixed b and d>=2. The finite union over d
remains rational. This is a classification algorithm for each fixed
branching defect, not a fixed finite bank as b increases.

The d=1 case has no bound(3). It retains the inherited endpoint-budget
funnel grammar and is not included in this rationality conclusion.
Nor does the theorem say that a ternary output returns to a low-defect
stratum: the inherited H35 obstruction already refutes that general step.

## 4. Complete first branching strata, including the missing nonneutral roots

For b<=5, (3) forces d=2 and H(R)<=2, hence R is transitive. Order the
two inports a->c. Encode a backbone vertex by the pair of bits saying
whether it beats a,c. The code in `two_inport.py` uses integers
0=00,1=10,2=01,3=11 (low bit is the first coordinate).
The source is strong exactly when the word is nonempty, its last letter
is nonzero, and at least one letter has first bit1. These conditions
ensure that the root reaches a, then c, and that every backbone vertex
can reach an inport; every inport returns directly to the root.

Four disjoint two-run skeletons count(1): the two port orders in a path
starting at an inport, and the two ordered path-bin endpoint choices.
Each has nine raw empty/good/bad run states. Exhaust all reachable
vectors capped at b+1, then retain exactly states that can reach an
accepting strong word. Removing the listed self-loops leaves a DAG.
Every remaining path is retained as an exact regular expression in the
receipt; these disjoint paths certify ALL word lengths.

The resulting ordinary rooted-class series are

    b=2: (x^4+2x^5+2x^6)/(1-x),
    b=3: 0,
    b=4: x^4+x^5+(6x^6+10x^7)/(1-x),
    b=5: 2x^5.                                                (5)

The b2 series is the inherited ready pentagon. The infinite part at b4
is exactly the inherited sixteen neutral tails. The two extra b4 words
are `10` and `01,10`, of orders4 and5. At b5 the complete list is
`10,01` and `01,11`, both of order5. These are nonneutral branching roots;
they were outside the old q=s atlas, not counterexamples to it.
Rooted isomorphisms fix the ordered two inports and the unique residual
chain, so different words are different rooted classes and all are rigid.

At b6 the transitive-residual d2 chart has four additional isolated words
and the inherited31 neutral tails. Its series is

    2x^5+2x^6+(x^6+12x^7+18x^8)/(1-x).

The isolated words, in integer notation, are `11`,`13`,`221`,`231`.
Here `11` means two consecutive letters of value1, not one letter11.
This chart is one of three disjoint sectors in the complete b6 result below.

## 5. Complete three-inport minimum-defect families

For d3,b6, (3) forces H(R)=1 and q-s=0. Since each of the three pair
bounds is at least2, equality forces H(R+v)=1 for each inport v. Each
such one-port extension of the transitive residual is transitive, so
every coordinate of the contact word is of the form1*0*. Moreover
q-s=0 forces the first SCC of U to avoid I. It is the first residual
singleton, so the word starts111. Strongness requires a nonzero last
letter. These are proved necessary restrictions, not empirical filters.

The generic positive run compiler in `multi_inport.py` exhausts this
restricted language for both inport tournaments TT3 and C3. For each,
its complete capped DFA has98 reachable states and two productive states;
the latter are exactly the initial state and the accepting111 self-loop.
Thus the only words are111^k, k>=1. This certifies a complete infinite
classification from a finite exact automaton, rather than extrapolating
the tested word lengths. The two rooted-class families have series

    2x^5/(1-x),

and their counts are7*2^k-5 (TT3 inports) and9*2^k-3 (C3 inports).
Indeed for an inport tournament I, uniform111^k gives
H=(Q_I(2)-H(I))2^k+3H(I)-Q_I(2); use(H,Q(2))=(1,8),(3,12).
The rooted automorphism groups are respectively trivial and C3, so the
two families must not be counted as rigid labelled presentations.

There is now an elementary proof for this sector and its all-degree
generalization in [MINIMAL_DEFECT](../metagraph/MINIMAL_DEFECT.md):
for every d>=3 the minimum b=2^d-2 consists exactly of
r->TT_k->I->r, with I=TT3 or C3 for d3 and I=TT_d for d>=4.
A mixed-contact residual vertex, or a zero-before-one pair in its chain,
creates a new positive tagged five-cycle and contradicts equality in(6).
This proof has been independently reviewed; the automaton is an additional
finite certificate of the d3 case.

## 6. The residual triangle and the complete defect-six classification

The only remaining possibility at b6 is d2,H(R)=3. The residual is
uniquely TT_a->C3->TT_b, a,b>=0. The two ordered inports are fixed, and
the contacts of C3 have24 orbits under its cyclic automorphism group.
An exact transition for this residual block sums its12 assignments to
two ordered path bins, respecting each bin's actual path order. This
extends the same positive run compiler without treating the triangle
as a transitive word.

Equality b6 forces q-s=0. If the first residual gap is nonempty, its
first vertex must beat both inports; if it is empty, all three vertices
of the first residual SCC must beat both. In either case the source SCC
of U avoids I, which preserves q-s=0 for every later residual prefix.
Before the triangle, tagged-capacity positivity and its independent
odd-cycle packing force delta_final>=3 delta_prefix. After it, deletion
monotonicity forces delta_final>=delta_prefix. Therefore the prefix
budgets2 and6 used by `cycle_residual.py` are necessary. The empty
starting state is treated separately: its value includes q-s and is
not pruned as though it were delta.

The complete capped automaton has34 reachable states and7 productive
states. Its language is exactly

    11^a C3_11 11^b w,       a,b>=0,
    w in {empty,10,01,(10,01),(01,10)}.                      (7)

The triangle's three contact letters are all11. The five terminal
tails are precisely the old ready pentagon. Thus the residual triangle
inflates that pentagon to FIVE two-parameter families. The root fixes
the two ordered inports, and the intrinsic SCC order fixes both gaps;
the remaining rooted automorphism group is C3. Each presentation in(7)
is one rooted isomorphism class. Their series is

    x^6(1+2x+2x^2)/(1-x)^2.                                 (8)

The universe of all strong rooted tournaments with branching defect6
is now exhausted: equation(6) permits only d2 with H(R)=1 or3, or d3
with H(R)=1. Summing the three sectors gives

    B_6(x)=(4x^5-2x^6+11x^7+10x^8-18x^9)/(1-x)^2.           (9)

Its counts from orders5 onward are4,6,19,42,47,52,..., and exactly
5n+2 for every n>=8. This is42 disjoint described families: four
isolated objects,31 transitive-residual rays, two three-inport rays,
and five triangle-residual surfaces. This number is an output of the
classification, not an assumed structural correspondence.

The mechanism is informative beyond the count: the first extra SCC
opens a second independent gap, changing eventual constant rooted
counts into linear growth. At higher fixed defects the same finite
kernel/gap grammar applies, but neither this particular polynomial nor
this list of cores is asserted to persist.

### Exact Hamilton values and the information lost by gap collapse

For a ready tail w, let W consist of its two inports and its tail vertices.
Write E for the ordered two-bin cover count of W whose first path ends
at an inport, q=H(W), and C for the sum of the two opposed endpoint covers.
The inherited variable-prefix decomposition in
[Session181](../../session181_multiport_defect_20261104/synthesis/NOTE.md)
gives, for any nonempty prefix tournament P with uniform contacts,

    H=E Q_P(2)+(q-C)H(P),   s_r=q_r=q H(P),   defect=C H(P).

The five tails have exact(E,q,C,q-C) data

|tail|E|q|C|q-C|
|---|---:|---:|---:|---:|
|empty|3|1|2|-1|
|10|6|3|2|1|
|01|5|1|2|-1|
|10,01|13|5|2|3|
|01,10|14|5|2|3|

For P=TT_a->C3->TT_b, let D=2^(a+b). Then H(P)=3 and Q_P(2)=12D.
The five triangle-inflated Hamilton values are respectively

    36D-3, 72D+3, 60D-3, 156D+9, 168D+9.                    (10)

This identifies every size and scalar value in these families. For a fixed
tail and m=a+b, there are exactly m+1 different rooted classes with the
same full displayed response. The scalar data forget the triangle's
position; the two ordered gaps restore it. These values are all divisible
by3, so NONE equals n^2+1, since a square is0 or1 modulo3. This excludes
the whole sector from A002496 without assuming anything about primes.
`pentagon_inflation.py` computes the finite interface coefficients directly
and independently checks50 literal triangle-gap instances. The general
formula follows from the path-bin decomposition, not from interpolation.

## Reproduction and scope

Run `python -B defect/two_inport.py` from the packet directory.
All341 literal words of lengths0..4 (including weak hostiles) are checked
by independent subset endpoint counting against the positive automaton.
The full capped state spaces at budgets2,3,4,5,6 have respectively
736,1217,1972,2765,3681 reachable states. These are not claimed minimal.
The computed finite automata certify their infinite languages; the
general rationality theorem rests on the construction in section3.
Independent implementations in [the metagraph audit](../metagraph/audit_defect.py)
and [the arithmetic audit](../arithmetic/audit_defect.py) review sections1–4
and check literal path words, positive terms and the exact accepted-language
series independently. Section5 has1170 literal DP controls, an exact
weighted-bisimulation/capped-DFA receipt, and the independently reviewed
all-degree proof in MINIMAL_DEFECT. Run `python -B defect/cycle_residual.py`
for section6. Its1368 literal controls are separately repeated with a
different graph gauge and first-vertex DP in
[audit_cycle_residual.py](../arithmetic/audit_cycle_residual.py), with24
additional larger random controls and independent generating-series checks.
Both the metagraph and arithmetic reviewers audited the pruning argument.
