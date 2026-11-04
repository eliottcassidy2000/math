# Exact dimension of the uniform tournament interface

**Status: PROVED (R1--R4 below), with finite checks described separately.**
This sharpens the realized-rank lower bound in
[Session176](../../session176_tournament_information_20261104/information/PROOF.md)
by recovering the pure-parity mechanism of
[THM-3181, tournament half-grid reciprocity](../../../../01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md).
No historical novelty, general circuit lower bound, spectrum completeness,
or nonuniform-attachment sufficiency is asserted.

## R1. A simple realized basis

Fix n>=1. Let q_T(k) count ordered covers of the n labelled vertices by k
nonempty directed paths. Put

    Q_T(t) = sum_{k=1}^n q_T(k) binom(t,k).

For integer t>=0 this counts vertex colourings by t ordered colours, with
one directed Hamilton word chosen on every nonempty colour class. Its degree
is n, its leading coefficient is one, and Q_T(0)=0. The leading coefficient
follows from q_T(n)=n!. THM-3181 proves
Q_T(-t)=(-1)^n Q_T(t). Thus all profiles lie in

    V_n = span_Q {t^n, t^(n-2), ..., t^epsilon},
    epsilon=1 for odd n, epsilon=2 for even n.

This space has dimension d=ceil(n/2). Ordered sums of tournaments multiply
Q: in each colour, a path is the concatenation of its paths in successive
summands. In particular adjoining a universal source or sink multiplies Q
by t. This uses the actual ordered-sum operation, not formal realizability.

Let B_m, m>=3, have vertices 0,...,m-1 in transitive order, except for the
reversed edge m-1 -> 0. Then

    Q_(B_m)(t) = t^m + t[(t+1)^(m-2) - (t-1)^(m-2)].       (R1)

Proof: a colouring with differently coloured endpoints has the same unique
monochromatic paths as the transitive tournament. If the endpoints have
one common colour and exactly j interior vertices have that colour, the
same holds when j=0. When j>=1 the induced tournament has 1+2^j Hamilton
paths: its increasing path, plus every path using the reversed edge,
specified uniquely by assigning the j increasing interior vertices to the
segment before m-1 or after 0. There are 2^j such assignments. Summing the
increment over colourings gives

    t sum_{j=1}^{m-2} binom(m-2,j) 2^j (t-1)^(m-2-j),

which is (R1). Equality at all positive integers gives a polynomial identity.
The j=0 case must remain separate: a reversed two-vertex tournament still
has one path.

Take W_0 to be the transitive n-tournament. For 1<=j<=d-1 take W_j to be
B_(2j+1) followed by n-2j-1 singleton summands. Their polynomials have
distinct least degrees n-2j: coefficient one for j=0 and coefficient two
for j>=1. Hence they are linearly independent, by successively examining
the least power in a purported linear dependence. All are actual
n-tournaments, proving

    span_Q {Q_T: |T|=n} = V_n,     dimension = ceil(n/2).   (R2)

Every Q_T is monic. The d witnesses are also affinely independent, so the
affine hull is exactly the monic hyperplane of V_n, of dimension d-1.
The binomial/power change of basis is invertible, so both statements also
hold for the ordered-cover vectors q_T.

## R2. Exact number of linear profile coordinates

A linear state and linear decoder defined on the rational profile span
require at least d scalar coordinates; d suffice. On the fixed-order
affine hull, an affine state and affine decoder require d-1 coordinates;
d-1 suffice. These are rank statements. They are NOT lower bounds on bits,
nonlinear encodings of the finite realized set, arithmetic registers,
adaptive oracle queries, or general Boolean circuits.

The distinction matters already for n=1,2: the affine hull is a point and
known order alone determines the profile. At n=3,4 the scalar H is enough
with known order. At n=5,6 one further independent coordinate is necessary
for affine recovery over the whole hull.

## R3. Attaining the bound with actual tournament probes

Use Session176's uniform contexts P_r: a transitive exterior
a_1,b_1,...,a_r,b_r, with T->a_i and b_i->T. Their actual Hamilton counts
h_r satisfy

    h_r = sum_{k=1}^{r+1} c_r(k) q_T(k),  c_r(r+1)=r!.

This triangular skeleton identity and its proof are inherited, not new.
It recovers q_T(1),...,q_T(s) from h_0,...,h_(s-1). Therefore it recovers
Q_T(1),...,Q_T(s) by the binomial definition.

For s=d these observations recover every vector in V_n by interpolation
at the d distinct positive squares 1,4,...,d^2 after dividing Q(t) by
t^epsilon. This proves linear recovery from P_0,...,P_(d-1).

For known n, monicity fixes the top coefficient. For s=d-1 subtract t^n
before interpolation; the remaining polynomial after division by
t^epsilon has degree at most d-2 in t^2. Thus

    P_0,...,P_(d-2) recover the entire actual fixed-n profile.       (R3)

There are no probes when d=1. When d>=2 the largest exterior has 2d-4
vertices. This attains the affine-coordinate bound, not a universal
minimal-size-context or finite-set-query theorem. A larger single context
may encode a finite family injectively without providing affine recovery.

## R4. When the recently compiled gates finish reading a fixed-order input

The independently developed
[gate-flow packet, CF2 and its finite-depth observer theorem](../../circuit_flow_tournaments_20261103/PROOFS.md)
provides actual uniform context operators L,U. On nonempty inputs every
Hamilton readout after at most r such gates lies in the span of
Q(1),...,Q(r+1). The readouts at the successive prefixes of the single
word L^r already span those evaluations: at positive j the highest shift
of L is j(j+1)(2j+1) Q(j+1), which never vanishes. At j=1 the lower shift
vanishes; U can also read Q(0), which is identically zero here. Thus no
negative evaluation enters the Hamilton readouts.

On the actual fixed-order profile span, the precise rank is therefore

    linear rank = min(r+1, d),
    affine variation rank = min(r+1, d-1),   d=ceil(n/2).    (R4)

Indeed evaluation of the allowed powers at distinct positive integers is
a Vandermonde matrix in their squares, with nonzero row factors. For
affine variations omit the fixed leading monomial t^n. R1 shows these
spaces are spanned by actual tournaments, so the ranks are attained there.

For d>=2, reading H at every prefix of L^(d-2) therefore recovers the full
profile when n is known, by affine interpolation. For d=1 order alone
suffices. These depths are sharp for affine recovery over the whole hull
using readouts of bounded gate depth. This does not say a single terminal
readout suffices, nor give a nonlinear encoding/query lower bound.

This is a consequence of the incoming gate theorem plus the sharp realized
rank proved here. The operators and their arbitrary-input depth theorem
are credited to that packet; their composition or spectrum results are
not reproved or extended here.

## Inheritance, tests, and boundaries

Closest proved mechanisms: THM-3181 pure parity and half-grid recovery;
Session176 exact uniform-context equivalence and triangular actual probes.
Corrected near miss: its floor(n/3)+1 construction was valid but not sharp;
three-cycles alone consume too many vertices for the smallest power sought.
The underused sidecar is the older reciprocal path-colour polynomial.

The canonical hostile remains the two H=15 five-vertex tournaments in
Session176, with q(2)=78 and 90. Their cyclic wrappers have 123 and 135
paths, confirming that the extra coordinate is an actual contextual
distinction. Uniform profiles still lose labelled orientation: opposite
two-vertex tournaments have identical Q but one nonuniform vertex can
distinguish them. No statement here changes that boundary.

The mechanism compares all three relevant spaces: the n-dimensional formal
cover-vector space, the d-dimensional realized linear span, and the
(d-1)-dimensional fixed-order affine hull. No arbitrary point of the hull
is asserted to be a realizable tournament profile.

`check.py` exhausts all 1,099 labelled tournaments of orders 1 through 5;
uses subset path and ordered-cover dynamic programs; checks the rank,
affine rank, parity, actual probe counts, and reconstruction; and checks
the explicit basis through n=11. It also checks all gate words through
depth six on every fixed-order power span and affine variation space
through n=12, including saturation. Failure checks stay active under `-O`.
The separate permutation-based audit is in `../rank_audit/`.
