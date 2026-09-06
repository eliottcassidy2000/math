# An exact arithmetic decoder for the actual rank-eleven 11+2 equality branch

**Status: PROVED ELEMENTARILY ON THE STATED DOMAIN + INDEPENDENTLY AUDITED;
FINITE-EXACT controls.**
The theorem decides the equality W=V_dec for an individual row in the stated
actual graph and physical-box domain. Its proof is elementary algebra and
the signed-box lemma. It does not classify every LRC row or prove LRC(14).
The eleventh and twelfth results remain separately frozen. No theorem ID or
external priority claim is requested.

The [independent referee](overnight13_20260906_lrc_entry_decoder_audit.md)
accepts the full proof and implementation without repair. Its separate
complement-semigroup membership routine and literal rational-span engine
pass 31,754 gates, including all 363 retained support certificates and
additional actual-graph hostiles. Normal and optimized outputs agree.

## 1. Inheritance and the missing orientation

The immediate proved mechanism is
[the audited twelfth signed-box theorem](overnight12_20260906_lrc_gcd_semigroup.md),
especially Lemma 2 and its
[independent proof audit](overnight12_20260906_lrc_gcd_semigroup_audit.md).
For two physical coordinates with primitive ratio a<b<=Q, that lemma decides
every bounded relation with a third coordinate by testing only the smallest
possible coefficient of the third coordinate. It uses the exact central
radius `Q(a+b)-(a-1)(b-1)` and the fact that this radius exceeds half the
whole signed support.

The exact entry object is inherited from **THM-3818**, Sections 6.4--6.5:
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
The full relation span W and the actual all-scale decoder graph must be
retained. **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`,
supplies the old finite-box pigeonhole mechanism; its required instance is
reproved in Section 3. MISTAKE-486/490 and the eleventh audit retain the
repaired all-scale atlas and the physical-box boundary.

The twelfth selected collection consisted of 110 supports with two core
labels and one label from the two-component. A complete crossing relation
can also use one core label and both labels of the two-component. Those
eleven omitted supports are the recovered sidecar. The canonical cheap
hostile in Section 5 passes all 110 old tests and has exactly one crossing
in the omitted orientation. This identifies the failed implication in an
attempt to turn the selected collection into an entry equivalence.

The board is: actual decoder components; componentwise weighted kernels;
internal pair heights; minimal distinguished coefficients; both support
orientations; and the finite physical box. The map sends each mixed support
to a two-generator coefficient box while keeping the distinguished label,
its physical coordinate, and its cleared coefficient. It preserves bounded
integer relation existence. Passing to a selected orientation loses an
entire family of actual crossings; restoring all 121 supports repairs that
loss. This is a source-to-target equivalence, not a phase or moment transfer.

## 2. The declared domain and the complete criterion

Fix `Q=91^6=567869252041`. Let n=(n_1,...,n_13) be thirteen distinct positive
integers, with

    gcd(n_1,...,n_13)=1,          sum_i n_i<=Q^2.       (1)

Construct the actual THM-3818 all-scale decoder graph. For a pair of physical
coordinates, reduce to coprime p<q. They are adjacent exactly when p+q is one
of the inherited admissible sums at most 356: each prime divisor is 2 modulo
3 and each prime exponent is at most two. This reproduces the 5,855-entry
atlas. Assume this actual graph has exactly two connected components I,J,
of sizes eleven and two. A proposed split supplied without the actual graph
does not satisfy this domain assumption.

Define W to be the rational span of all integer relation rows z satisfying

    sum_i z_i n_i=0,          |supp z|<=3,          max_i |z_i|<=Q.

Define V_dec to be the rational span of all primitive two-coordinate
relation rows on actual decoder edges. The graph gives dim V_dec=11.
For any pair i,j let

    H(i,j)=max(n_i,n_j)/gcd(n_i,n_j).

For a pair i,j in one component and distinguished label k in the other,
order the first two physical coordinates as A<B and put Y=n_k. Set

    D=gcd(A,B),       a=A/D, b=B/D,
    delta=gcd(D,Y),   c=D/delta, x=Y/delta,             (2)
    r0=a^(-1)x mod b in [0,b-1],
    s0=(x-a r0)/b,
    lower=max(ceil((-Q-r0)/b),ceil((s0-Q)/a)),
    upper=min(floor((Q-r0)/b),floor((s0+Q)/a)).          (3)

The distinguished label in this arithmetic test may be a core label or a
pair-component label. It must not always be interpreted as a physical
outside label. This is the coordinate change that lets the same lemma
handle both orientations.

Call the support positive when

    c<=Q and lower<=upper.                            (4)

**Theorem.** Within the domain (1) and the actual 11+2 graph assumption,
W=V_dec if and only if both of the following hold:

* all 55 pairs within I have H(i,j)<=Q;
* every one of the 121 mixed triple supports is negative under (2)--(4).

The 121 supports consist of `2*binom(11,2)=110` supports with two labels
in I and one in J, and eleven supports with one label in I and both labels
in J. The one pair internal to J automatically has height at most 355<Q
because it is an actual decoder edge. Equations (2)--(4) are applied only
after the internal height gates, so b<=Q is always available.

Two-coordinate crossings are included: a coefficient of the third label
may be zero. There are no nonzero one-coordinate relations because every
physical coordinate is positive. Thus the mixed triples, with these zero
coefficients allowed, cover every crossing of support at most three.

## 3. Proof of both directions

First identify the decoder span exactly:

    V_dec={z in Q^13:
             sum_(i in I) n_i z_i=0,
             sum_(j in J) n_j z_j=0}.                 (5)

To see this, multiply coordinates by n_i, an invertible diagonal map over
the rationals. Each primitive decoder edge relation becomes a nonzero
multiple of e_i-e_j. The incidence differences of a connected graph span
the zero-sum subspace of that component. Applying the inverse diagonal
map proves (5). Also V_dec is contained in W: every decoder edge has
primitive coefficient height at most 355<Q.

For necessity of the internal pair-height gates, suppose an internal pair
i,j has H(i,j)>Q, and select any k in the other component. There are
(Q+1)^3 coefficient triples in {0,...,Q}^3, but their weighted sums lie
among at most Q^3+1 integers, since

    Q(n_i+n_j+n_k)<=Q sum_l n_l<=Q^3.

Pigeonhole gives a nonzero relation on these three labels with coefficients
of magnitude at most Q. Its coefficient at k cannot vanish: otherwise it
would be a nonzero two-coordinate relation at i,j, whose primitive height
already exceeds Q. Its distinguished partial sum is therefore nonzero,
so the relation belongs to W and lies outside (5). This contradicts
W=V_dec. This necessity uses the physical box explicitly.

Now all internal pairs have b<=Q. The twelfth minimal-coefficient lemma
applies to each mixed triple: a bounded relation

    C Y=r A+s B,       1<=C<=Q,       |r|,|s|<=Q

exists if and only if (4) holds. A global sign change makes any nonzero
distinguished coefficient positive, so neither sign is lost. Such a
relation has a nonzero partial sum on the component containing Y, hence
lies outside (5). Therefore equality entry makes all 121 supports negative.

Conversely, assume all the stated gates and all 121 negative tests. Every
bounded relation whose support meets both components is captured by one of
these tests, including the support-two case. Hence no such relation exists.
Every bounded relation remaining entirely in one component lies in (5).
Thus W is contained in V_dec. The opposite containment was proved from
actual decoder edges, so W=V_dec. This proves the iff.

The internal-height rejection need not construct the pigeonhole relation:
the displayed pair, its exact height, the physical box, and the proof above
are a complete obstruction certificate. A positive support test does return
a literal bounded relation: choose an integer h between lower and upper,
put r=r0+hb, s=s0-ha, and use coefficient row (-r,-s,c) on its selected
labels. The source checks its exact weighted sum and coefficient height.

## 4. The exact decoder and its computational scope

[The standalone implementation](../../04-computation/overnight13_20260906_lrc_entry_decoder.py)
first checks the row and box domain, then constructs all 78 possible graph
edges and its components. A different graph type or an out-of-box row is
reported as out of the declared domain, not as a negative equality verdict.
Within the domain it records the internal heights, and then every support
certificate in both orientations. It returns equality exactly when the
theorem says so.

At most 121 mixed-support tests are needed for a fixed thirteen-label row.
Each is a fixed list of gcd, exact division, modular inverse, and integer
rounding operations. The number of arithmetic calls is independent of the
speed height. Their bit cost depends on the integer sizes; this is not a
constant bit-time claim. The routine neither enumerates Q-sized coefficient
boxes nor traverses physical phase grids. Its small coefficient enumeration
appears only in the independent verifier controls.

The retained JSON lists all 121 certificates for every named in-domain
control: labels and orientation, normalized ratio, physical pair gcd,
clearing divisor, minimal distinguished coefficient, normalized target,
and both interval endpoints. Positive tests include the full thirteen-entry
relation row. Negative tests retain whether the coefficient gate or the
interval gate fails. Thus a reviewer can replay the entry decision without
trusting a final Boolean or a sampled rank.

The theorem is an efficient pointwise entry decoder for the given branch.
The number of possible speed rows inside the inherited finite box remains
enormous. No feasible exhaustion of that box, safe-phase classification, or
LRC(14) closure is inferred.

## 5. Both essential hostiles and nonvacuity

Let

    U=(1,4,6,8,10,12,14,15,16,18,22).

The two positive entry controls are the frozen canonical row
`U union 2^45(1,3)` and the unitless row from the twelfth theorem,

    (2,3,4,5,6,10,12,14,15,20,30) union (60Q+1)(1,3).

Both have the actual connected 11+2 graph, lie in the physical box, and
have no positive support among all 121. The independently proved dominance
condition g>2Q max V excludes crossings in both orientations, providing a
second route to W=V_dec for these controls. Their safe phases are inherited
from the twelfth report; safety is not needed for the present entry test.

**Hostile to retaining only 110 supports.** Set

    t=3Q+1,             n=tU union (1,3).             (6)

The source verifies the actual 11+2 graph, physical box, distinctness and
primitivity. All core internal heights are at most 22. Every two-core/
one-pair test is negative. Indeed its minimal distinguished coefficient is

    tD/gcd(tD,s) >= t > Q,           s in {1,3},

because 3 does not divide t. However the omitted orientation has the
literal bounded crossing

    t*1 - 1 - Q*3 = 0.                                (7)

Its coefficients are 1,-1,-Q and it meets both actual components. Among
the eleven opposite-orientation tests it is the unique positive one in
this row: all other core entries are at least 4t>4Q, beyond the signed
support of the physical pair (1,3). Thus the old 110-negative condition
is not an entry equivalence, even inside the precise domain.

**Hostile to dropping the physical box.** Let

    V=(1,355,355^2,...,355^10),
    K=355^10,             g=2QK+1,
    n=V union g(1,3).                                  (8)

The graph still has exactly the two actual components. Consecutive core
coordinates have primitive ratio (1,355), an atlas pair since 356=4*89
and both prime factors are 2 modulo 3 with allowed exponents. These
bounded-height edge rows connect the whole core. Cross edges are absent:
g is coprime to every core power and every reduced cross ratio has sum
larger than 356. The pair itself is an atlas edge.

The inequality g>2QK excludes every bounded mixed relation, so the actual
decoder span equals W. Nevertheless the internal pair (1,K) has primitive
height greater than Q. There is no contradiction: this row has sum greater
than Q^2. The implementation reports it outside the domain and does not
turn its failed height gate into a false equality verdict. This is an
actual decoder/equality example outside the box, not a claim that graph
connectivity alone implies the inherited internal-height bound.

## 6. Complete finite controls and reproduction

The source imports no mathematical producer. In its independent small
engine, every support of size three is examined, the first two coefficients
are enumerated in [-Qtoy,Qtoy], and the third is solved by exact integer
division. Zero coefficients retain every support-two relation. Duplicate
rows are removed as literal integer vectors. This engine uses neither the
minimal-coefficient lemma nor its interval calculation.

The complete bounded universe is all primitive five-element subsets of
{1,...,9} with sum at most Qtoy^2 for Qtoy=4,...,7: **323 physical rows**.
All ten choices of a three-label component are tested for each row, giving
**3,230 splits** against **71,528 distinct bounded relation rows** counted
over the full row/height universe. The direct comparison checks whether
every relation has zero weighted partial sum on the proposed component.
This is the componentwise kernel equality criterion, with no assertion
that the toy split is an actual THM-3818 graph.

Two explicit positive toy controls also retain literal complete relation
banks: at Qtoy=30, components (1,3) and (181,543); at Qtoy=40, components
(1,3,4) and (321,963). Their internal pair relations span both component
kernels and every bounded relation remains internal. The real thirteen-
label controls then independently verify actual atlas typing, all 121
certificates, both orientation counts, and both hostiles above.

    python -B 04-computation/overnight13_20260906_lrc_entry_decoder.py --write-certificates 04-computation/overnight13_20260906_lrc_entry_decoder_certificates.json
    python -B -O 04-computation/overnight13_20260906_lrc_entry_decoder.py

The [full JSON](../../04-computation/overnight13_20260906_lrc_entry_decoder_certificates.json)
is a consequence certificate, not just a gate count. Both runs pass
**25,286 always-active gates** with byte-identical LF stdout. The JSON also
uses LF. Frozen SHA256 identities:

    source: 9b8fe6804a037a0f93a29396940e510691ee3b2765be1eeca1271ed28e5b7b6c
    output: 8932a324412f5c0b7ce8cfc243146af374454a2239aa6059030f00af261bd57d
    JSON:   af177ea9703390a8169d9466d04a5e4d8bb53f66833cfbfd223c94a0006cb847

No Git, shared navigation, or earlier frozen producer was changed.

**Filing:** root integrated these independently audited artifacts in the thirteenth
checkpoint. Reproduction commands are relative to the repository root.
Transcript hashes refer to filed LF bytes; Windows CRLF captures were normalized.
