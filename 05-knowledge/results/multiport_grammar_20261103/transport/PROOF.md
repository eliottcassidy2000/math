# Several anonymous endpoint colours, named exits, and exact changing-defect transport

**Status: PROVED, after independent analytic/source/replay audit.**
This is an unbounded exact evaluator for a specified construction language.
It proves neither reusable two-digit return closure nor global Hamilton-path
spectrum coverage. Native Hamilton-complement legality is a separate issue.

## 1. Inheritance and the exact increment

The inherited
[one-bulk forest grammar](../../multiport_defect_grammar_20261103/transport/PROOF.md)
already retains arbitrarily many path components under births, named flips,
forgetting and converse. The newer
[Session182 named birth law](../../session182_physical_returns_20261104/defect/PROOF.md)
already carries all physical joint responses and inner defects when the full
named endpoint matrix is retained. We do not claim either mechanism anew.
The Session183
[terminal-context theorem](../../session183_two_branch_frontier_20261104/mining/PROOF.md)
also already identifies the port debt of its two fixed local recipes.

The extension here is to several anonymous endpoint populations. A later
vertex may see the populations in different directions, without restoring
all their individual names. This requires paired endpoint colours, not just
a component count or separate start/end totals. In contrast to the one-bulk
case, one birth can consume two anonymous components, and the resulting
operator has genuine mixed second derivatives.

The canonical hostile is the inherited two-vertex failure of anonymous
reactivation. Here its two vertices are kept in separate colours, making
the formerly forbidden cut an authorized operation and giving different
correct outputs. The corrected near miss is counting the same component
once as incoming and again as outgoing: that creates a cycle, not a path.
The least-used sidecar is the ordered pair of anonymous endpoint colours.
Our live board is: endpoint pairing; anonymous colours; physical ports;
falling multiplicities; evolving root defect; native/return guards.

The source is an actual named tournament with a colouring of forgotten
vertices. Its target is a multivariate cover polynomial. The map retains
all observations and operations proved below, but loses individual hidden
names, internal positions and nonuniform future incidence inside a colour.
The full construction word or literal adjacency remains a separate source
sidecar. Splitting an existing colour or resurrecting a forgotten name is
not an operation on this quotient.

## 2. The finite-dimensional polynomial state at every finite boundary

Fix a finite set P of active physical names and a palette of q>=1 colours.
Every other vertex has exactly one colour. A spanning directed path cover
is an **unordered** collection of nonempty vertex-disjoint directed paths
covering the tournament; singleton paths are allowed.

For each component retain its paired (start,end), using its physical name
if an endpoint is in P and its colour otherwise. A component with at least
one active endpoint belongs to an unordered anchored signature kappa.
A component with anonymous start colour a and end colour b contributes a
factor z_ab. Define

    C_T[kappa](z) = sum_covers_with_kappa product_(a,b) z_ab^(k_ab).   (1)

There are q² commuting variables. An active vertex cannot occur as an
endpoint in two different components. The pair (p,p) represents its
singleton component; an active name absent from the signature lies in a
component's interior. These are the same physical names throughout a word.

The number g_(p,q) of formal anchored signatures satisfies

    g_(0,q)=1, g_(1,q)=2+2q,
    g_(p,q)=(2+2q)g_(p-1,q)+2(p-1)g_(p-2,q) for p>=2.                         (2)

Indeed the least active name is absent, a singleton, paired with one
anonymous colour in either direction, or paired with one of the other
p-1 active names in either direction. The q=1 specialization recovers
the inherited recurrence. At graph order n, every polynomial has total
anonymous-component degree at most n, so the complete state has at most

    g_(p,q) binom(n+q²,q²)

integer coefficients. This is an upper bound on stored coefficients,
not a rank, a reachability claim, or a polynomial-time algorithm for
arbitrary tournaments. With fixed p,q it is polynomial in n, but the
allowed future cuts are restricted accordingly. Both p and q may grow;
each state and each operation still have finite exact semantics.

The empty graph initializes (1) to the scalar one at the empty signature.
For nonempty T, H(T) is read as described in Section 6. Physical order,
active names, colouring and orientation sidecars are retained separately.

## 3. Birth with a separate incidence bit for every colour

Let x be a fresh active name. Specify a subset of active vertices beating
x and a subset I of colours beating x. Every vertex in a colour in I
beats x; x beats every vertex of the other colours. This restriction is
on each colour separately, not on their union.

For every old cover select either no incoming component or one whose end
beats x, and either no outgoing component or one whose start loses to x.
If both are selected they must be **different components**. Join them
through x. When neither is selected, add the singleton (x,x). The new
endpoints are the incoming component's start (or x if absent) and the
outgoing component's end (or x if absent).

An anchored component is selected by its actual pair in kappa. Each is
unique because an active endpoint cannot be reused. Selecting an anonymous
(a,b)-component acts by partial_ab; it contributes k_ab possibilities.
If the two selected anonymous components have the same type, the correct
factor is k_ab(k_ab-1), exactly partial_ab². If their types differ, their
factor is k_ab k_cd, exactly partial_ab partial_cd. Resulting anonymous
endpoints multiply by their new variable; anchored endpoints change kappa.

In particular the part consuming two anonymous components is precisely

    sum_(b in I, c notin I) sum_(a,d)
        z_ad partial_ab partial_cd.                                  (3)

The derivative acts first; multiplication by z_ad acts afterwards. Its
indices describe an incoming (a,b)-path and outgoing (c,d)-path. This is
an operator on each anchored coordinate, leaving kappa unchanged. The
other terms use at most one derivative and the deterministic changes of
anchored pairs just specified. Selecting the same anchored component twice
is rejected, even when both its endpoints satisfy the cut.

**Birth theorem.** This rule is a positive integral polynomial differential
operator of order at most two, and it produces the complete state (1) of
the actual birth. It is valid for arbitrary component degree and any finite
q,p, not just the presently observed one- and two-component layers.

Proof. Delete x and its incident path arcs from a new cover. If x was
internal this splits its path into two different old components; there is
no shortcut between its former neighbours. If x was an endpoint it leaves
one component, and if singleton it leaves none. This uniquely recovers the
old cover and the selected ordered components. Conversely the selection
rule has the proper indegree/outdegree bounds and creates no cycle exactly
when the two components differ. Colour-uniform incidence makes the number
of valid selections depend only on (1), giving the displayed falling
multiplicities. The two constructions are inverse.

For q=1 there cannot be two selected anonymous components: the one colour
would have to lie both in I and outside I. Thus (3) vanishes and the old
first-derivative grammar is recovered exactly. For q>=2 it is necessary.

## 4. Names, colour maps, converse and flips

**Forget a port.** To forget active p into an existing colour a, replace
that name by colour a in its endpoint pairs. Each pair that becomes fully
anonymous multiplies the corresponding z_ab. Coefficients with the same
new signature are added. The actual vertex remains in the graph.

**Map colours.** For any function phi from the old palette to a new finite
palette, substitute z_ab -> z_(phi(a),phi(b)) and apply phi to anonymous
ends of anchored pairs. This exact pushforward includes colour permutations,
merges and adjoining unused new colours. It does not split an occupied old
colour. A fresh name can later be forgotten into a newly added colour.
Future incidence must respect the resulting coarser or larger palette.

**Converse.** Reverse all anchored pairs and substitute z_ab -> z_ba.
This is exact reversal of every path in every cover.

**Named joins and flips.** For distinct active names a,b, let J_ab join
the component ending at a to the different component starting at b.
If either endpoint is missing or they belong to the same component, the
result is zero. If both resulting endpoints are anonymous, multiply their
variable; otherwise retain their new anchored pair. Then reversing an
actual arc a->b has exact state update

    C_(after)=(I+J_ba-J_ab) C_(before).                              (4)

The inherited cutting proof remains valid with endpoint colours. A cover
using the removed arc has a unique cut preimage; a new cover using the
added arc has a unique cut preimage in the old graph. Active names make
those arcs definite physical edges. The active induced orientation must
be retained to supply the sign and to update the next reversal.

The same-forest argument proves named joins commute and
J_ab²=J_ab J_ba=J_ba J_ab=0. Distinct flips can therefore be multiplied
with their mixed terms retained. Birth and forgetting need not commute
with these operators. Signed intermediate algebraic expressions are not
asserted realizable; a whole authorized flip on an actual state has
nonnegative final coefficients by the bijection above.

Consequently any finite word of these typed operations has exact state
semantics, with no uniform limit on word length or component degree.
This is the claimed unbounded grammar. No accepted-return theorem is
hidden in this compositional statement.

### Bounded-interface evaluation is polynomial, with the input contract retained

Fix upper bounds p and q throughout a supplied typed word of length L,
and let n be its final physical order. Starting from the empty graph, a
fixed explicitly initialized seed, or a **given exact initial profile**,
the whole word and all readouts above can be evaluated in time polynomial
in n and L (using indexed names, with the external name dictionary
counted as part of the supplied input). This is an evaluation statement, not an algorithm for finding
a successful return word or initializing arbitrary tournaments cheaply.

There are at most g_(p,q) binom(n+q²,q²) stored coefficients. Per monomial,
a birth has at most (1+p+q²)² ordered component choices; joins, forgetting,
colour maps and converse have only a bounded number of local changes for
fixed p,q. Each operation therefore takes polynomially many coefficient
updates. A tournament on N>=1 vertices has at most N! 2^(N-1) unordered
path covers, since permutations with cut sets already enumerate all ordered
covers and hence overcount the unordered ones. Every canonical nonnegative
coefficient has O(N log(N+1)) bits. A named join selects at most one ordered
pair from any physical cover, so the total absolute coefficient sum in the
signed flip calculation is at most three times the old cover count.
Birth's positive partial sums are bounded by its final actual cover count.
Thus arithmetic bit lengths remain polynomial as well.

The literal adjacency and current colour/name assignment require O(n²)
bits plus the finite word data. Updating them and testing the SCC guard
by direct reachability is also polynomial. Palette changes must stay within
the fixed bound q; unbounded q or p is permitted by the general theorem but
not by this polynomial-complexity corollary. Forgetting a port does not
delete a vertex, and no unproved hidden-name resurrection is used.

## 5. Two cheap hostile tests and the exact loss boundary

Take no active vertices and give vertex0 colour red, vertex1 colour blue.
The actual tournaments 0->1 and 1->0 have respective states

    z_red,blue + z_red,red z_blue,blue,
    z_blue,red + z_red,red z_blue,blue.                            (5)

Merging both colours gives the identical one-bulk polynomial z+z².
The allowed coloured birth red->x->blue produces a transitive triple
from the first and a directed triangle from the second, so the actual
Hamilton counts are 1 and 3. Separate colours pay this distinction; the
merged state cannot supply the cut. Applying an arbitrary splitting map
back to the merged colour would incorrectly identify these two sources.
One must retain the colours beforehand, or recompute a richer observer
from the literal graph, rather than call such recomputation transport.

There is also a same-component hostile inside (5). The single path
blue->red has both an eligible incoming end and an eligible outgoing
start for this cut. Selecting it in both roles would make the directed
cycle blue->red->x->blue. Its monomial has

    partial_blue,red²(z_blue,red)=0,

not the product 1*1. Selecting two distinct such components, when present,
uses the nonzero factor k(k-1). Thus simply multiplying eligibility totals
without the diagonal exclusion is false.

Paired colours themselves cannot be replaced by separate colour marginals:
the two cover signatures {(a,b),(c,d)} and {(a,d),(c,b)} have the same
start and end marginals. The colour cut I={a,b}, with {c,d} outgoing,
has only b among eligible ends and only c among eligible starts. It joins
two different components in the first, but would
select the same (c,b) component in the second. These can be actual covers
of the same transitive four-vertex tournament ordered a,c,b,d.

## 6. Named endpoint observers and every retained inner defect

For a monomial/signature put d=sum k_ab+|kappa|. For nonempty T:

* H is the total coefficient with d=1.
* F_r=E_r(T,2) is the total coefficient with d=1 or 2 and a component
  ending at the retained physical name r. At d=1 the second slot is empty;
  at d=2 the r-ending component is uniquely the first slot.
* J_T(r,v) is the total coefficient with d=2, one component ending at r
  and the **other** component starting at retained v.
* The actual inner defect is b_r=F_r-H, for every currently retained r.

These are exact integer readouts after every authorized operation. The
rooted bad-join identity gives the positive interpretation

    b_r=sum_(v->r in T) J_T(r,v)>=0.                              (6)

If all incoming hidden vertices form whole colours, the corresponding
summed colour-joint readouts evaluate (6) directly. If a colour has mixed
incidences to an already existing r, colour aggregation need not identify
that summand; the difference F_r-H still evaluates b_r exactly. We do not
assume historical root incidences become uniform when colours are merged.
A forgotten root loses its own named readout, precisely as specified by
the grammar. Choosing it again requires retaining or rebuilding a richer
observer, not a silent scalar reset.

Current H, F_r and J values alone are not a sufficient iterative state:
operators can lower component number, bringing arbitrarily high layers
into those observations. All degrees in (1) remain stored. The inherited
one-bulk cyclic-macro coefficient `(k+1)!k!` already demonstrates that a
fixed component cutoff cannot be justified for unrestricted future words.
Our result extends its allowed incidence language, not that truncation.

## 7. Actual strongness and an infinite changing-defect control

Strongness belongs to the retained actual graph sidecar. For a nonempty
old tournament, let C_first,C_last be its first and last SCCs in the
unique transitive order. The new graph is strong exactly when x beats
some vertex of C_first and some vertex of C_last beats x. Necessity:
otherwise the first component is unreachable from x, or the last cannot
reach x. Sufficiency: from x reach the first component and then every
old component, and from any old vertex reach the last and then x. This
also handles one old SCC; a singleton cannot satisfy both requirements.
A birth from the empty graph yields the strong singleton separately.
This guard is exact, but is not a Hamilton-complement native-path guard.

For a positive arbitrarily long word, start with cycle 0->1->2->0.
Keep root0 active, colour1 red and2 blue. At every step birth x with
all red vertices beating x, x beating all blue vertices and root0, and
forget x into red or blue alternately. Both populations stay nonempty,
so every birth has an old in-neighbour and out-neighbour and preserves
strongness. No old name or edge is deleted from the physical source.

After k births the graph is C3({0},{1},TT_(k+1)): 0 beats only1,
1 beats the remaining transitive block, and that block beats0. The
coloured birth inserts its new vertex between the existing red tail
and blue head of that block. Let n=k+1. A Hamilton path either places
0,1 consecutively between an arbitrary subset of the chain and its
complement (2^n choices), or is the single order 1,chain,0. Thus

    H=2^(k+1)+1, F_0=2^(k+2)-1, b_0=2^(k+1)-2.                  (7)

For F_0, deleting0 leaves a transitive order with1 first and n other
vertices beating0. Its ith vertex, i=1,...,n, has two-slot endpoint
count 2^i. The root-peeling identity therefore gives
`F_0=1+sum_(i=1)^n 2^i=2^(n+1)-1`.
This is a known simple tournament family used as an all-length transport
control, not a new spectrum family. Its changing defect is retained
exactly even though the same old name0 is used throughout.

## 8. Closure and coverage remain different proof obligations

The grammar evaluates a supplied typed word; it does not supply a word
for either desired digit from every output. Neither polynomial coordinate
size at a fixed boundary nor the infinite family (7) implies a physical
return class. Actual adjacency, named exits, source acceptance, strongness,
native guard paths and return conditions remain separate interfaces.

Even a future closed +2/+4 tree from finitely many sources would not prove
full spectrum coverage. The
[Session154 adaptive-branch obstruction](../../session154_gadgets_20261103/closure/NOTE.md)
applies without a finite-state assumption: at most two children with each
H multiplied by at least3 give only O(X^(log_3 2)) constructed values below
X, missing infinitely many primes and 7p values. The present multicolour
extension changes state semantics, not that counting theorem's hypotheses.

## 9. Exact verification contract

The independent direct evaluator enumerates actual spanning path covers
by choosing the unique component containing the least remaining vertex,
then sums directed path endpoint counts obtained by literal permutations.
It does not use the new birth or flip operators to construct the target
profile. The finite bank includes every labelled tournament through
order4, every assignment of each old vertex as active/red/blue, every
active-and-colour-uniform birth, every active flip, every forget target,
colour merge, palette permutation and converse. Endpoint subset DP
independently checks all named one/two-slot observers and strongness guards.
The two-colour, same-component and split hostiles remain explicit. A longer
literal word checks the whole profile at small orders and direct endpoint
DP plus (7) at larger orders. No finite result is extrapolated into a
missing two-branch acceptance theorem.

Both normal and optimized runs pass **216,070 always-active exact gates**
and produce identical **474-byte** stdout. The complete finite universe
has 76 labelled graphs including the empty graph, 5,422 coloured/active
assignments, 67,732 typed births, 3,530 active flips, 14,282 forget operations,
and 5,422 merges and converses. The source checks all the named observers
against independent endpoint subset DP. The longer word has 16 observed
births and ends with H=131073, F_0=262143 and b_0=131070, using 2,236
nonzero profile terms. Its complete profile is independently enumerated
through graph order6 and its endpoint readouts independently DP-checked
through graph order12; the displayed formula (7) is proved for every length.

```bash
python3 04-computation/multiport_grammar_20261103/transport_colours.py
python3 -O 04-computation/multiport_grammar_20261103/transport_colours.py
```

Frozen source: **13,609 bytes**, SHA256
`860928f457f1227d017d677f1412b2e9da017e2604f8ee7d4a6d513f618ee77a`.
Frozen output: [verification.out](verification.out), SHA256
`d415cc1fc86727fb986a7a094c5701f99629372d750deaf8c325d36782b34d7e`.
Semantic receipt SHA256:
`c9601f087a73eb440bf18e940e31b0f9ad2e25dfdb255344a412f6388b8afc2a`.

The independent [audit](AUDIT.md) accepts the all-order proof and exact
implementation. Root promoted this status after that audit; the mathematical
body and executable source were unchanged by the promotion.
