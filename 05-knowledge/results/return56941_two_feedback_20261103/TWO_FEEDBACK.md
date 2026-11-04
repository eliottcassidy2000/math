# A positive run expansion and the1182 feedback boundary

## F1. Uniform blocks have positive ordered-run coordinates

**PROVED.** Fix a tournament template, one vertex v, and its full external
incidence. Inflate v to a transitive block TT_m, m>=1. There are nonnegative
integers c_r, determined by the template with its ordered visits to v, such
that

    H_m = sum_(r>=1) c_r r! S(m,r),                         (F1)

where S is the Stirling number of the second kind. The sum is finite, with
r at most one plus the number of external vertices. Delete the block from
a Hamilton word, retaining each ordered run position and every external
vertex. Each resulting skeleton with r run positions accepts exactly the
surjections from the m labelled clones to those positions. Inside each
run the transitive order is forced. Restriction and insertion are inverse,
so the coefficients are actual nonnegative skeleton counts.

If deleting q distinguished feedback vertices leaves a transitive base
and the growing block is an **interval** of that base, then r<=q+1.
Deleting those q vertices splits a path into at most q+1 increasing base
pieces, each meeting the interval in at most one run. For q2 this gives

    H_m=c1+c2(2^m-2)+c3(3^m-3*2^m+3),
    c1=H1, c2=(H2-H1)/2, c3=(H3-3H2+2H1)/6.             (F2)

Thus a pure constant-plus-2^m family forces c3=0. A third-run response
cannot cancel by a negative realization coefficient. This is a small
observer for this particular growth operation; actual complements need
not preserve the transitive-base hypothesis or its block interval.

## F2. A literal second port removes1180 units of the excess

Inherit the repaired source from
[RC2--RC4](../return18979_ports_20261103/PROOFS.md). Its base is TT_m then
the11-chain Q with names32,22,2,1,6,8,4,13,5,14,7, plus apex p=21.
The apex beats every clone and beats Q_i exactly when bit_i(1243)=0.
Its count is

    H_m=1243*2^m-1061.

This is the **repaired6325 clone source**, not a free presentation of18979.
The historical edge repair has a count debt, retained in the inherited
source receipt. The one-apex target forced by the dyadic coefficient has
constant-1999, whereas3H+2 requires-3181: the excess is1182.

Birth q beaten by every clone and by Q_1=22, and beating all other exterior
vertices, including p. This is cut code5 in `feedback.cpp`: bit0 controls
all clones; bits1..11 control Q; bit12 controls p. The target is strong for
every m because the old source is strong and q has both an in-neighbour
and an out-neighbour. Direct run coefficients are

    (c1,c2,c3)=(4275,3729,0),
    H(T_m+q)=3729*2^m-3183=3H_m.                          (F3)

This reduces the original1182 excess to a deficit2. `audit_feedback.cpp`
independently constructs the matrices and checks m1..5. Formula(F2)
extends the first three values to all m. It is an actual infinite family
with named incidences, not a numerological equality.

`feedback_ports.json` identifies the m1 matrix literally with the inherited
repaired source by its old names. For m1..7 it retains31 as the first clone,
the spare6, apex21 and the new feedback port; every ordered two-bin pair
response and all five Q4/current F49/F343 charts are recomputed. Later clone
addresses are distinct. The full named generator remains the carrier; the
scalar3H equality is not used to infer any of those port responses.
Their all-size graded transport and the growing third-mode correction are
derived in [TRANSPORT T1--T4](TRANSPORT.md).

Every second-port cut was also tested with up to three individually split
prefix clones. At the correct dyadic coefficient and c3=0, the closest
constant is always target minus2 (`feedback.json`). That bounded test alone
does not exclude arbitrary split-incidence words or a third feedback port.

## F3. Complete two-feedback exclusion on the same clone clock

**PROVED reduction + FINITE-EXACT exhaustive certificate.** There is no
fixed two-feedback template with a uniform transitive interval of size m,
arbitrary fixed exterior chain on both sides of that interval, and strong
realizations satisfying either

    H_m=3729*2^m-3181,  or H_m=3729*2^m-3179              (F4)

for every m>=1. The conclusion also holds for eventual equality on every
consecutive m: (F2) extends such an identity backwards to m1, and strongness
of a nonempty transitive inflation is equivalent to strongness of its
quotient. **The block has size m on the source clock.** This statement does
not silently cover a decreased block clock m-s, multiple growing blocks,
or an exterior depending on m.

Let t be the fixed exterior-chain size. At m1 the total order is t+3 and
the required counts are4277 or4279. Busch's sharp strong minimum at17 is
5625, so t<=13. Encode every base vertex by

    c=[p->v]+2[q->v] in{0,1,2,3},

and retain the orientation of p-q. Place the uniform block in any of the
t+1 base gaps. Apex exchange takes bulk letter2 to1; global converse and
base reversal take0 to3. Hence bulk letters1,3, every exterior word,
every gap and both p-q orientations exhaust the required objects. These
are symmetries of a counting proof, not free operations in a source program.

For a fixed order of the two feedback vertices, delete them from a word
and count its three increasing base gaps. Each of the first two gaps is
empty, nonempty with an unsuitable last boundary, or nonempty with a
suitable last boundary. The last gap needs only its empty/nonempty bit.
This gives18 states per apex order,36 states total. First-entry guards and
last-exit guards enforce the actual incidences. Empty middle gaps require
the correct p-q edge. The implementation evaluates m1,m2,m3 together.

Appending a base vertex cannot decrease H: insert it at the first old
vertex it beats, or at the end if none; deleting the inserted name recovers
the old path, proving injectivity. Therefore a lane is permanently excluded
as soon as one of its three counts exceeds its target. The same argument
justifies pruning a prefix before the block is inserted.

`search_words.py` checks process completion and records:

    all_positions visits 97648782
    pruned              62397763
    prefixes            10833829
    hits                       0.

This proves the finite head of the reduction. An independent literal
subset-DP audit checks38,477 cases: every exterior word through length4,
every block position, **all four** block letters, both apex orientations,
m1..3;120 seeded longer controls; and the positive(F3) family through m5.
Counts fit the declared unsigned storage: each apex order has at most
3^(t+m) assignments and t+m<=16 in the exhaustive head. The separate
readout uses64-bit sums. No floating-point comparison or unproved pruning
criterion occurs.

## F4. The square-sum transfer is about paired runs and distinct names

The inherited square-sum classification in
[SQUARES](../return18979_ports_20261103/SQUARES.md) is complete for its
declared anchor/connector problem, not for every square-sum graph. Its
source-to-target map is useful here: deleting q connectors from a path
leaves at most q+1 ordered old path pieces; reconstruction must retain each
piece's paired endpoints and require all connector names to be distinct.
An arbitrary connected connector graph is not enough; its three-leaf
hostile is connected but nontraceable.

The tournament analogue deletes named feedback vertices, retaining the
same ordered pieces plus directed endpoint incidence. A transitive base
then compresses each piece to its boundary conditions, giving the36-state
compiler above. With a general strong base, the pieces have independent
internal path counts and can revisit a module; their joint endpoint tensor
must replace this small compiler. This is exactly where the familiar
two-endpoint-only carrier fails. No orientation is imposed on square-sum
nonedges, and a Hamilton path is not used to declare the whole graph a
tournament.

The next useful split is therefore **not** just another scalar offset:
permit a controlled nonzero c3, record the actual ordered three-run
incidences, and ask which legal source transformation returns those data.
That is an OPEN source grammar problem. (F3) supplies a positive control;
(F4) supplies a complete uniform-block hostile.
