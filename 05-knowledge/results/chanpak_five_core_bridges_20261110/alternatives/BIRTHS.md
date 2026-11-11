# Actual ternary branches at73, and an all-order obstruction at141

**Status: exact finite certificates plus the proved strong-ear reduction.
Chan–Pak remains OPEN.** This extends NOTE's constructive attempts: a
ready-root obstruction does not forbid arbitrary two births, but actual
source retention still fails at the next member of the same explicit family.

## 1. A positive small branch and its complete two-birth universe

Let T_a(k), T_b(k) be the two literal sources from NOTE: free prefix `11^k`
and tail `01,11,10` or `01,01,10`. Both have H=`34*2^k+5`.
All old vertices and old arcs remain fixed in this section. A new outmask O
means the new vertex beats precisely the old vertices in O.

At k1, H73. Exhausting all `2^(2*7+1)=32768` ordered two-vertex attachments
to each source gives:

|Actual source|New total221=3H+2|New total223=3H+4|
|---|---:|---:|
|T_a(1)|4 labelled attachment pairs|4 pairs|
|T_b(1)|2 pairs|8 pairs|

The receipt `two_birth_k1.json` retains every cut pair and actual child rows.
The fourfold families for T_a and all T_b223 cases are one-birth realizations
with an extra neutral singleton. The T_b221 cases are genuinely two-birth:
first outmask0, then outmask8 (or the opposite labelling of the new pair).
Deleting the first new vertex leaves H133, and deleting the second leaves
H73; neither deletion has count221. This is a concrete positive realization
despite NOTE's proof that no supertournament of the source can be ready.
The missing ready interface and the existence of a general prescribed-count
extension are different questions.

The endpoint calculation used to exhaust cuts is exact. For a parent R,
let s_i,e_i be its start/end Hamilton counts and let Q_ij count ordered
spanning covers by two nonempty paths, first ending at i, second starting
at j. If x_i indicates that the new vertex beats i, then

    H(new)=sum_i s_i x_i + sum_i e_i(1-x_i)
           +sum_(i,j) Q_ij(1-x_i)x_j.

Deleting the new vertex from a Hamilton path partitions it into the
starting-at-new, ending-at-new, and interior cases. This proves the formula;
subset endpoint DP computes each Q_ij without a sampling assumption.

## 2. Why the uniform-prefix templates do not lift

Replace the one initial11 vertex by a transitive block of size k, keeping
each new vertex's old comparison with that entire block. If A_r counts
valid orders of the other vertices with r separated placeholders for this
block, then the exact response is

    H(k)=sum_(r>=1) A_r r! S(k,r),

where S(k,r) is a Stirling number of the second kind. Indeed the k ordered
block vertices must be assigned surjectively to r labelled runs; within
each run their order is forced. The finite positive placeholder DP in
`two_birth.py` computes A_r directly from the actual graph. It is a complete
response, not a fitted exponential curve.

The successful base templates have exactly the following five responses:

|Base case|All-k response|Desired all-k response|
|---|---|---|
|T_a,221|44*3^k+51*2^k-13|102*2^k+17|
|T_a,223|103*2^k+17|102*2^k+19|
|T_b,221|104*2^k+13|102*2^k+17|
|T_b,223, first group|111*2^k+1|102*2^k+19|
|T_b,223, second group|26*3^k+79*2^k-13|102*2^k+19|

Every response equals its target at k1 and strictly exceeds it for every
k>=2. For the linear cases subtract the target directly. For a3^k-b2^k-c
difference, its next value minus twice its present value is a3^k+c>0,
and its initial difference is0. Thus none of the finitely many successful
base templates can be the desired uniform-prefix recurrence, even at one
later k. This refutes that precise lift, not every variable-prefix operator.

The next finite test includes all nonuniform prefix choices as well. At k2,
H141, all `2^(2*8+1)=131072` ordered two-vertex extensions per source give:

|Actual source|425=3H+2|427=3H+4|
|---|---:|---:|
|T_a(2)|0|0|
|T_b(2)|0|4|

The four positive427 cases are a one-birth outmask94 with neutral padding.
`two_birth_k2.json` records them, including both possible one-prefix-vertex
inflation responses. Thus a movable split of the two prefix vertices cannot
repair the missing425 branch at this scale. This is an actual source failure,
not a failure inferred from a scalar recurrence or from ready-root absence.

## 3. Strong ears make the141 failure an all-order obstruction

**Lemma (relative strong enlargement).** If S is a proper induced strong
subtournament of a strong U, it has an induced strong enlargement inside U
by either one or two vertices. If there is no one-vertex enlargement, every
outside vertex either beats all S or loses to all S. Strongness of U gives
an edge from a vertex x losing to S to a vertex y beating S. Then
`S -> x -> y -> S` is the required two-vertex enlargement.

For that two-ear, write Q_S(2) for the ordered two-slot path-cover count,
including empty slots. Its exact count is

    H(S[x,y])=H(S)+Q_S(2).                           (1)

If y precedes x in the Hamilton path, all S lies between them, giving H(S).
If x precedes y, they are consecutive, and the old vertices before x and
after y form an arbitrary ordered two-slot cover, giving Q_S(2).

Cutting each Hamilton path of S at any of its n+1 positions injects the
pairs `(path,cut)` into those covers. Hence

    H(S[x,y]) >= (n+2)H(S).                          (2)

The maps are injective because concatenating the two pieces recovers the
old Hamilton path and its cut. This argument also explains why a two-ear
can be far more expensive than an arbitrary two-birth realization: the
relative-strong assumption forces a very specific cyclic attachment.

**Proposition.** Neither T_a(2) nor T_b(2) is an induced subtournament of
any tournament with H425. Moreover T_a(2) is not induced in any tournament
with H427, whereas T_b(2) has the explicit one-vertex427 extension above.

Proof: each source has n8 and H141, so (2) gives two-ear lower bound1410,
already above427. Every later strong enlargement has at least8 vertices
and H>=141 by ordinary Hamilton insertion monotonicity. Therefore any
strong supertournament of total count at most427 can be reached from the
actual source by a chain of strong **one-vertex** extensions: a missing
one-vertex step would force a forbidden two-ear inside the target.

The complete finite strong frontiers under H<=427 are:

|Source|One new vertex|Two new vertices|Three new vertices|
|---|---:|---:|---:|
|T_a(2)|25 labelled states|4 states, H381 or393|empty|
|T_b(2)|28 labelled states|4 states, H381 or411|empty|

There is no isomorphism filtering. The four states in each second layer
include different new-vertex label orders. `frontier.json` retains every
state, its exact cuts from its parent, the count histogram, and each full
empty-child certificate. Every retained count is independently checked by
direct subset Hamilton DP. No425 occurs;427 occurs only in the second
source's first layer. The empty third frontier proves closure, rather than
merely stopping at a chosen order. The ear lemma promotes this exact finite
closure to all strong supertournaments at the stated bound.

Finally, in an arbitrary tournament containing the strong source, the source
lies inside one strong component C. Hamilton counts multiply along the
ordered components, so H(C) divides the total and H(C)>=141. Every proper
positive divisor of425 is at most85, and every proper divisor of427 at
most61. Thus at either requested total the containing component must have
the whole count, already excluded or supplied by the strong-frontier result.

This is not a new missing spectrum value:425 is explicitly attained by the
ordinal sum of S4, S4, and a one-apex source `10^4`, with count5*5*17.
It is an all-order obstruction to preserving either of these two H141
sources while realizing the next desired ternary count. The source must be
replaced, a retained old arc must change, or old vertices must be deleted.

## 4. Scope, checks, and next constructive obligation

`python -B two_birth.py --prefix 1` and `--prefix 2` reproduce the two exact
attachment receipts. Each positive hit is directly recounted, and each
module response is checked against literal one-copy and two-copy graphs.
`python -B frontier.py` reproduces `frontier.json`; `--fail` must exit
nonzero. The frontier tests use no inherited enumeration, no unrecorded
symmetry filter, and no inference from residues. They prune only using the
proved monotonic bound and strong connectivity appropriate to the ear lemma.

The constructive next question is whether another H141 source has a425
extension, and whether an actual count-neutral replacement can transport
one of these sources to it while retaining the next required interface.
The explicit deletion returns in NOTE supply genuine graph operations, but
their counts change; scalar all-target coverage does not fill this gap.

The history lane has now supplied the positive comparison: a different
one-apex H141 source has an explicit two-birth trace `141 -> 215 -> 425`.
Its exact rows and cuts are in
[`../history/alternate141_search.json`](../history/alternate141_search.json).
That source belongs to the current five-core grammar. The all-order negative
and positive examples therefore concern different actual realizers of the
same count, not a missing target or a missing grammar value. Transport to
that alternate marked source remains the constructive requirement.

The cover lane independently reproduced all35,840 strong-frontier cuts
using an even-directed-path cut formula and a first-vertex DP, with no
imports from this implementation. Its126 literal controls and every retained
row/path/count agree; see
[`../cover/frontier_audit.json`](../cover/frontier_audit.json).
Normal and optimized executions of all three retained receipts agree
byte-for-byte. The independent audit is in addition to, rather than a
replacement for, the strong-ear and component-divisor arguments above.
