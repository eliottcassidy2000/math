# One physical C2 context encodes every bounded source response

Session 199. **PROVED, conditional only on the named reviewed source and graft
inputs:** a single explicitly constructed strong C2 graft context distinguishes
all numerical source responses with cover count at most B. The context uses
O(log B) vertices, and its Hamilton-path count is an invertible base expansion.
Three constant-size contexts already distinguish all physical source responses
without any bound. This is contextual tomography, **not odd-target coverage**.

## 1. Historical inheritance and the bridge

The unexpected starting point was the knot connected-sum work, not another
modular covering argument. We examined the surrounding chain:

* `01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md`,
  Sections 1-2: an operation-ready continuation profile, rather than one raw
  observable, is the natural congruence carrier.
* `01-canon/theorems/THM-2191-catalytic-localization-of-the-gordian-metric.md`,
  Section 1 and its response-ideal discussion: common translation and
  nonexpansivity make optimal-context loci upper ideals.
* `01-canon/theorems/THM-2281-common-optimal-context-for-finite-catalytic-families.md`,
  Sections 1-1a and 4: one common context realizes all requested pairwise
  localized distances in a finite packet; arbitrary infinite packets are not
  covered.
* `01-canon/theorems/THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete.md`,
  Section 2: actual forced-leaf contexts recover each selected matching-kernel
  entry. The important ingredient is a physical isolating context, not a formal
  linear functional.

The source citations and their imported scope were checked in
`05-knowledge/reference/CORE-PAPERS-COMPOSITIONAL-RELATIONS.md`. We use no new
external knot or matching theorem: the source results motivate the construction
below. In particular, neither knot prime decomposition nor perfect-matching
leaf attachment is asserted to act on tournaments.

A direct transfer fails: the source letters do not commute, their positive
readouts expand rather than contract, and a tournament cannot have the sparse
forced leaves used to isolate a matching edge. We replace those operations by
an explicit one-entry tournament-source word and an exact radix decoder.

The closest inherited target mechanism is the already proved injective signature
in Session 197 `SOURCE_SIGNATURE.md`. Its rank-three N quotient is **not** a new
result here. The canonical hostile is a fixed small partner that gives equal
path counts for two different physical source responses. The least-used sidecar
is the distinction between a graph-isomorphism class and its four numerical
source responses.

Concept board: physical context selection; source-response observability;
positive radix digits; bounded common context; exact target membership.

## 2. Physical domain and inherited laws

A nontrivial source is a rooted tournament (V,c), |V|>=2, with exactly one
Hamilton path starting at c. Its response is

    v=(D,E,U,d)=(H(V),e_c(V),R_c(V),H(V-c)).

Here R counts ordered two-path covers with a nonempty c-start reserved path
and an optional unreserved path. H(empty)=1. We inherit

    d<=D<=U,  E<=d,  D-E>=d,

and the complete source grammar with seed t=(1,0,2,1):

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

O appends a universal sink, Z appends a vertex beating only the root and keeps
the old root, and J makes that same birth the new root. Words are chronological.
All letters are actual graph operations.

The inherited signature is

    A=D+2d,  B=E+d,  C=U-d.                         (1)

It determines the physical response, not the graph isomorphism class. The exact
physical inverse from `dependencies/SOURCE_SIGNATURE.md`, Section 2, is

    if B is odd:
      d=B when A>=3B, and d=A-B-C otherwise;
    if B is even:
      d=B/2 when 2A<=3B+2C, and d=A-B-C otherwise;
    D=A-2d, E=B-d, U=C+d.                            (2)

This formula is justified only on the physical image. It is not a realization
test for arbitrary integer triples.

For two nontrivial sources the actual strong C2 graft N has

    N(v,w)=(A-B)(A'-B')+AC'+CA'+3CC'.                (3)

The full graph construction and exactly-two-cycle law are in the copied
`BILINEAR_C2.md`, Sections 2-4. Every probe used below is a source, so its N
readout is an actual tournament Hamilton-path count.

## 3. Three small physical partners separate all response profiles

Let

    w0=t, w1=J(t), w2=J(J(t)).

Their respective orders are 2,3,4 and signatures are

    (3,1,1), (5,2,2), (11,5,3).

Write ri=N(v,wi). Then (3) gives

    r0=3A-2B+6C,
    r1=5A-3B+11C,
    r2=9A-6B+20C.

The coefficient matrix has determinant 2. Its explicit inverse is

    A=3r0+2r1-2r2,
    B=(-r0+6r1-3r2)/2,
    C=(-3r0+r2)/2.                                  (4)

Consequently, equality of these three actual C2 path counts for two physical
sources is equivalent to equality of their full response quadruples. It is
also equivalent to equality against every N partner: the reverse implication
follows by specializing the partners, and the forward implication follows
from (3) and the inherited physical injectivity (2).

The same response equality is preserved under every subsequent source word,
and under every inherited C2/C3 graft, because those constructions depend only
on the complete response. This does not say that arbitrary tournament contexts
or graph isomorphism are determined.

### The displayed bank is irredundant

No one of the three displayed probes can be discarded. The following are
physical collisions for each possible pair; all words start at t:

    probes  left word / response       right word / response      common values
    0,1     OJ / (5,1,5,1)             ZJ / (5,2,6,3)             (41,73)
    0,2     OOJJOOZ / (89,9,135,9)     ZOJOZOZ / (93,3,129,3)     (1041,3375)
    1,2     OJOZ / (15,1,19,1)         ZJOO / (5,0,24,3)          (277,501)

Each row follows by applying the source laws and evaluating (3); the verifier
also reconstructs and counts all twelve corresponding graft graphs directly.
The omitted probe separates each pair, as (4) requires. This is irredundancy
of this specified bank, not a proof that every possible two-context bank fails.

## 4. A single bounded common context

Fix an integer cap Bcap>=2. Choose

    k=ceil(log2(14 Bcap)),  x=2^k,
    wx=t . O^k J J O^k.                             (5)

Here k>=1 is essential: O^0 is the identity, not the map which sets E=0.
The partner has 4+2k vertices. Direct application of the physical laws gives

    wx=(1+4x, 0, 2x+4x^2, 1+2x),
    (A',B',C')=(3+8x,1+2x,4x^2-1).                 (6)

Substituting (6) into (3) yields

    N(v,wx)=q2 x^2+q1 x+q0,                         (7)
    q0=A-2B=D-2E,
    q1=6(A-B)+8C=6D-6E+8U-2d,
    q2=4(A+3C)=4(D+3U-d).

For every physical source with U<=Bcap, these are nonnegative integer digits.
More precisely,

    1<=q0<=Bcap,
    0<=q1<=14 Bcap-2,
    0<q2<=16 Bcap-4.                                (8)

For q0, D-E>=d>=E gives nonnegativity, while D is odd, so q0 is positive.
For q1, A-B=D-E+d>=0 and C=U-d>=0 prove nonnegativity, and dropping -6E,
using D,U<=Bcap and d>=1 proves the bound. The q2 bounds follow similarly.

Since x>=14Bcap, the two lower coefficients q0,q1 are strictly below x. Thus
there are no carries out of either lower position, and from n=N(v,wx) one
recovers

    q0=n mod x,
    q1=floor(n/x) mod x,
    q2=floor(n/x^2).                                (9)

No upper digit bound is needed for the last formula. Put L=q2/4 and S=q1/2.
Then

    C=3L+3q0-2S,
    A=-8L-9q0+6S,
    B=-4L-5q0+3S.                                  (10)

Equations (2) and (10) recover the full physical response exactly. In particular:

> **Common-context theorem.** For each Bcap>=2, the one actual strong C2 graft
> context with partner (5) maps distinct physical source-response quadruples
> of U<=Bcap to distinct Hamilton-path counts. The N graft adds 5+2k vertices
> to the input source and has exactly two Hamilton cycles.

This is one common context for every source in the cap, not a separate context
chosen for each pair. The bound is effective. Since x<28Bcap, (8)-(9) also give

    N(v,wx)<16 Bcap x^2<12544 Bcap^3.               (11)

No sharpness or optimality of these constants is claimed. The construction
uses O(log Bcap) additional vertices and an O(log Bcap)-bit scalar readout.

## 5. Exact failure boundaries and connection contract

A small fixed context need not separate sources. At x=2, the actual wx has
response (9,0,20,5). The different physical sources

    t . O Z J = (11,4,12,5),
    t . Z O Z = (9,1,11,1)

both give N(v,wx)=931. The cap condition fails. This is an actual physical
collision, not an unrestricted-cone example. We make no assertion that the
constant 14 is a necessary threshold or that every fixed context has a collision.

Conversely, valid-looking radix digits do not create a physical source. The
formal tuple (1,0,3,1) yields (q0,q1,q2)=(1,28,36), but is not a source:
D=1 forces a transitive source with U a power of two. Thus the number
36x^2+28x+1 is a formal readout of this tuple, not a graph construction from
that tuple. We do not claim this number is absent from other tournament families.

The exact bridge is:

    source: finite common-context/physical-probe mechanisms in the older
            knot and matching work;
    target: the physical start-one source-response set with a cover cap;
    map:    v -> the actual C2 graft N(v, t.O^k.J.J.O^k);
    retained: full numerical response and every response-determined
              source-word/C2/C3 continuation;
    lost: graph-isomorphism type, unrestricted tournament contexts;
    sidecar: the cap, the physical source grammar, and the physical
             signature inverse;
    hostile: the 931 collision and the nonphysical (1,0,3,1) tuple.

The resulting probe is an injection into a sparse positive scalar image.
Neither (7) nor decoding (9) says that a prescribed odd target has physical
source digits. A common context solves simultaneous **distinguishability**,
not the existential coverage problem. Exact odd coverage outside 7 and 21
remains open.

## 6. Validation

Run `python3 verify.py` and `python3 -O verify.py` in this directory. The script
uses explicit exception checks, not assertions. It verifies the identity on a
four-vector basis against the original four-coordinate N expression; the
source response and the three-probe inverse; the exact bank of 126,910
nontrivial source profiles with U<=5,000; all 40 source words through length
three by literal subset path counting; 26 literal C2 probe graphs, twelve literal irredundancy grafts, two literal
small-context collision grafts; and a full
no-carry selector at Bcap=2 with 17 vertices. The literal counter separately
counts Hamilton cycles, always obtaining two for these grafts. It includes the
physical collision and the fake-source control above. Frozen dependency hashes
are in both receipts. These finite checks supplement the all-order proof.
