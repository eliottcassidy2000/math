# One selected complement recovers integral three-run weights through order seven

2026-10-29, 10:00 UTC continuation. **PROVED elementary reconstruction and
orbit formulas; FINITE-EXACT complete strong-source census through seven.**
Same-agent review with alternative enumeration and integer certificates.
This concerns arbitrary weights on a fixed tournament's ordered three-port
space. It does not classify all strong tournaments, identify every physically
realizable attachment, or close the full unnamed metagraph carrier.

## K1. The operation, source and integral certificate

For a tournament T on n>=3 named vertices, let D=n(n-1)(n-2). The matrix
M_T has ordered triples of distinct vertices as rows and columns; its entries
count spanning covers by three ordered, nonempty directed paths, including
singletons. This is exactly THREE_RUN's native matrix. A state v is an
arbitrary column of D port weights. The response is M_T v.

For a legal Hamilton path P retain the vertex alignment and form the actual
relation U=T^op xor E(P). The two observations use the SAME weight vector:

    O_(T,P) v = (M_T v, M_U v).                         (K1)

This is not a product of independent port marginals. Global converse would
give M_T^T; the actual path complement gives a generally different M_U.
Targets need not remain strong. Their full relations, including any ordered
strong-component decomposition, are retained in the calculation.

**Integral reconstruction criterion (PROVED).** If full-column-size minors
d_1,...,d_s of O have gcd1, O has an integer left inverse. Consequently it is
injective over every Z/mZ, and over any commutative coefficient ring.

Indeed, let S_j select the rows of the corresponding square matrix B_j.
Choose integers u_j with sum u_j det(B_j)=1. Then

    L = sum_j u_j adj(B_j) S_j,       LO=I.              (K2)

This is the usual determinantal-ideal/Bezout mechanism, not a new abstract
linear-algebra theorem. It supplies an unusually small certificate for the
actual path operators: the receipt uses at most five minors per source.
The largest determinant in the complete census has102 decimal digits.
No prime factorization or guess that only primes2 and3 matter is used.

The closest proved mechanism is [THREE_RUN T4--T5](THREE_RUN.md).
The hostile is the81/27 pair; a second hostile below defeats arbitrary
choice of a complement path. The corrected near miss is confusing full
rational rank with recovery at every modulus. The less-used sidecars are
the July [free Hamilton-frame lemma](../../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md)
and [retained-path memory N1--N3](../charged_recurrence_memory_20261028/COMPLEMENT.md).
The latter already proves that repeating a fixed complement merely alternates
two states: a bad path cannot repair itself by repetition.

Concept board: actual relation / common port alignment / determinantal
ideal / unnamed path orbit / source symmetry / precision. The anchor is
native complement closure; the niche is integer recovery without factoring;
the wildcard is an exact Frobenius term in the unnamed state count.

## K2. Complete finite theorem, including the sharp depth distinction

**FINITE-EXACT.** For every strong tournament with3<=n<=7, either M_T already
has an integer inverse or there is a Hamilton path P for which O_(T,P) has
an integer left inverse. Thus a single suitably selected complement suffices
for every modulus simultaneously. The receipt records the path and coprime
integer minors for each of the396 source classes.

| n | Strong classes | Modular primes | D | M_T injective over Q | M_T invertible over Z |
|---:|---:|---:|---:|---:|---:|
|3|1|1|6|1|1|
|4|1|0|24|1|1|
|5|6|3|60|3|3|
|6|35|15|120|15|7|
|7|353|197|210|176|16|

The last column is exact: every excluded source has a recorded replay test
at a prime where its rank drops. Thus the minimum complement depth for
recovery at all moduli is zero on the last column and one on all other
sources. Rational depth zero is the preceding column. These are two
different optimization questions.

The source bank is inherited from the audited
[strong-atom atlas](../strong_atom_continuation_vitali_20261026/README.md),
but its exhaustiveness is checked again without assuming its class counts.
The verifier tests strong connectivity, every possible proper module,
automorphisms and pairwise nonisomorphism using degree-ordered canonical
forms. It sums n!/|Aut(T)| and compares with the independent recursion

    s_1=1,
    s_n=2^binom(n,2)-sum_(k=1)^(n-1) binom(n,k)s_k 2^binom(n-k,2). (K3)

The proof of K3 uses the unique first strong component. The matched labelled
counts at n=3,...,7 are2,24,544,22320,1677488. This closes the universe check;
the numerical coincidence of a matrix rank with a class count plays no role.

Source rank certificates contain integer null vectors and matching modular
minors. Both independent cover enumerations from THREE_RUN are rerun on
every matrix. All legal paths of every rationally deficient source through
six are tested; at seven a successful path is supplied for each deficient
source, without claiming an exhaustive classification of every path.

The chosen integral witnesses sometimes leave the strong locus: this occurs
for16 of28 depth-one sources at six and168 of337 at seven. Those are counts
for the stored choices, not impossibility results for strong targets.
The native metagraph must accommodate these actual target relations.

## K3. The cheapest hostile paths already occur at five

For cyclic T5 with rows[6,12,24,17,3], all15 path choices give rational
rank60, as proved previously. Their full integral classifications differ:

    10 paths: Smith(O)=1^60;
     5 paths: Smith(O)=1^59,2.                           (K4)

The five parity-defective paths form one Aut(T5) orbit. For example
P=(0,2,4,1,3) gives target[28,9,18,4,10]. At precision2^k its observer
has image size2^(60k-1), for k>=1, and kernel size2. An extra binary
distinction is needed for this fixed path, at every positive precision.
Choosing one of the other two path orbits removes that loss integrally.

A nonunit Smith factor does not create ambiguity when exact integer outputs
are available: an injective integer matrix can be inverted on its image
using rational coefficients and the image guard. The missing bit here is
about finite-precision observations and an integer-linear inverse valid
over every coefficient ring. These claims must not be interchanged.

Now let B5 have rows[16,5,9,3,14]. Its module{1,2,3} is a directed triangle;
the quotient is C3, with the other two vertices singleton blocks. Its15
paths divide into

    9 paths: Smith(O)=1^60;
    3 paths: Smith(O)=1^59,3;
    3 paths: rational rank58.                           (K5)

For the bad path P=(1,2,3,0,4), the target is[22,28,24,17,0]. The common
kernel has dimension two. Unlike T5's four run-symmetric boundary vectors,
this plane carries the standard two-dimensional S3 representation on run
names: its character is2,0,-1 on identity, transpositions, and three-cycles.
Both basis vectors are fixed by every VERTEX automorphism of B5.

In particular the loss survives restricting to vertex-invariant weights.
One kernel vector has six entries+1 and six entries-1. Its positive and
negative parts are distinct nonnegative weight laws of equal mass6, fixed
by Aut(B5), with identical source and bad-complement responses. The exact
port lists are in complement_details.json. Switching to the good path
(1,0,4,2,3) changes their response at start triple(0,1,3) from0/0 to0/2.
The three bad choices form one path orbit and all miss this same plane.

This also limits a phase-only repair. The hidden plane is trivial under
vertex transport, so specifying another vertex rotation phase
does not determine its two amplitudes. A joint run observation or an extra
realization equation is needed. An LCA rule that reconstructs relative
vertex frames alone cannot recover this particular invariant weight loss.
This is a scoped obstruction, not a rejection of LCA incidence recovery
when its additional equations are actually available.

This is the smallest rational arbitrary-path failure among strong sources.
For modular-prime sources, the first such failure is at six: rows
[2,4,9,19,39,15] with P=(0,1,2,3,4,5) give stacked rank112 of120.
Another path restores all120 coordinates integrally. The complete deficient
six-vertex path census contains16 rank-deficient choices across six sources.
None of these failures refutes the existential selected-path result K2.

## K4. An unnamed path-orbit selector, without averaging away correlations

Fix T and write G=Aut(T), acting on its port coordinates. A vertex-unnamed
weight state is an ORBIT [v] under G. It is not the same thing as a vector
fixed by G. The distinction is essential even in characteristic dividing |G|.

Suppose P satisfies K1's integral criterion. Its G-orbit is a legal unnamed
selector. For each Q in that orbit, express the source and complement
response vectors in the vertex order supplied by Q, retaining them as one
PAIRED record. Keep the multiset of these records, forgetting the names of Q.

**Orbit reconstruction lemma (PROVED).** This multiset determines [v] over
every coefficient ring. In a fixed reference P-frame it is precisely the
multiset of F(gv), g in G, where F is the injective stacked operator with
its output coordinates reordered into that frame. Equality of two multisets
gives F(gv)=F(hw) for some g,h, hence v and w lie in the same orbit. The
converse follows by permuting the multiset. No division by |G| is performed.

The Hamilton-frame freeness lemma makes each path orbit a torsor. This does
not select one vertex-equivariant physical path. One can select the orbit
intrinsically, for example by the least framed source code among successful
orbits. The two response vectors must retain their common frame and pairing.
An arithmetic average or separately canonicalized endpoint vectors are
different observers; the lemma does not grant them reconstruction.

Consequently K2 also establishes a single successful UNNAMED PATH-ORBIT
probe at every strong source through seven for this specified weighted,
path-framed observer. This is stronger than testing invariant weight laws,
but still distinct from HYP-9416's free class-valued reset protocol.
The full tournament relation is supplied, not reconstructed from M3 alone.

## K5. Exact minimum counts and the surviving Frobenius term

For weights in Z/mZ, the minimum set carrier for the above observations on
a fixed successful source is its orbit space. If c(g) counts cycles of g
on the D ordered port triples, classical Burnside counting gives

    K_T(m)=|G|^(-1) sum_(g in G) m^c(g).                 (K6)

Injectivity in K4 proves minimality for this observer, rather than merely
using an orbit count as a proposed carrier size. At m=p^k, this is a finite
sum of exponentials p^(c k), so it obeys the recurrence with characteristic
polynomial product_(distinct c)(X-p^c). The actual cycle indices give

    K_T5(m)=(m^60+4m^12)/5,
    K_B5(m)=(m^60+2m^20)/3,
    K_Paley7(m)=(m^210+14m^70+6m^30)/21.                 (K7)

For example T5 has K_(k+2)=(p^60+p^12)K_(k+1)-p^72 K_k.
These formulas count arbitrary weighted port states; they are not A000568,
nor counts of realizable attachments with a fixed cover mass.

The field interpretation is literal and inherited from
[Session69's1--2--4 bridge](../session69_prime_ports_power_aliquot_20261026/session69-root/CONNECTIONS.md).
Use F8=F2[alpha]/(alpha^3+alpha+1) and identify exponent i modulo7 with
alpha^i. Paley arcs are j-i in{1,2,4}; Frobenius squaring sends i to2i.
Together with field multiplication it gives the actual21-element group.
Its six nonidentity translations have30 port cycles, and its fourteen
order-three elements have70. Thus the term14m^70 retains the Frobenius-fixed
weight laws. Quotienting by translations alone gives

    K_C7(m)=(m^210+6m^30)/7,
    K_Paley7(m)=(K_C7(m)+2m^70)/3.                      (K8)

This is the precise correction to naive division by three. It is a transfer
of a known field action to the newly certified weighted observer, not a new
identification of Paley7. The F64/F256 exceptional-isogeny phase results in
FROBENIUS remain separate operators with their own12/16 clocks.

The twelve-vertex81/27 control has these port-cycle histograms:

    G81: {148:36,440:18,444:8,520:12,776:6,1320:1},
    G27: {444:8,520:12,776:6,1320:1}.                    (K9)

The actual labelled G27 is a subgroup of G81, checked directly. K6 therefore
also counts their ambient unnamed weight orbits, exposing the symmetry lost
by cyclic triples. No minimality for their three-run observations is claimed:
no full1320-by1320 response or integral complement certificate was computed.
The original two-run mass9011898 and modulo-three ranks24/22 are replayed.

## K6. Larger hostile probes, reproducibility and the next questions

The same integral certificate succeeds on three specified order-eight
sources C3[C3,C3,T2], C3[S6,1,1], C3[T5,T2,1], where S6 is the prime
source in K3, and on the recursive order-nine source C3[C3,C3,C3]. The last
has automorphism group81. These four are **FINITE-EXACT examples**, not
complete size-eight or size-nine censuses. The stored successful paths and
coprime minors are replayed in complement_scout.json. They do not settle
the earlier nine-vertex reversal-marking obstruction, a different language.

In the nine-vertex example the source rows are
[58,60,57,464,480,456,135,263,71]. Its native ranks modulo2,3,1009 are
432,255,504. The standard path(0,1,2,3,4,5,6,7,8) raises these to502,504,504:
two binary directions remain invisible. Moving the path to
(0,1,3,4,5,6,7,8,2) gives coprime integer minors, hence full recovery
at every modulus. This is an actual moving-selector result on weighted
ports of the recursive nine-vertex atom, with its scope explicitly fixed.

Reproduction, from the repository root:

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/complement_verify.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/complement_verify.py
```

The generation options --build-ranks and --build-units belong to
complement_kernel.py; normal verification does not rewrite receipts.
The verifier uses exact integer determinants and null vectors, two native
cover constructions, a separate orbit-mass exhaustiveness argument, full
five-vertex Smith forms, and explicit F8 arithmetic. This is a same-agent
audit, not an external review or formal proof-assistant verification.

**OPEN selected-path conjecture:** for every strong T of order n>=3, some
legal P makes the three-run operator O_(T,P) integrally left-invertible.
The census and larger probes support it; they do not prove it. Its strong-
target variant allows zero moves for sources already integrally recoverable;
it is a separate all-size question. A counterexample may force more
path orbits, more complement depth, or a different run arity.

The general composition frontier still needs physically realizable joint
port laws and complete complemented relations, with quotient walks and
their shared alignment. Frobenius and LCA can recover some coordinates
only where the relevant equations and common sources have been supplied.
Our minimum for arbitrary weighted states is not automatically a lower
bound on a smaller realization-constrained attachment family.

Incoming0b1eb8940d/c2943c341d supply full Paley7 Smith groups and the
clock-preserving pairing completion law. Those results are inherited, not
recomputed as discoveries. Meta-patterns used: distinguish representations
from carriers; compose witness fibres before taking scalar observations;
test the actual integral target rather than an intermediate rank.

**11:00 UTC continuation:** [RUN_RECURRENCE R1--R7](RUN_RECURRENCE.md)
certifies the strong-target variant through seven. Zero moves are essential:
all Hamilton complements of C3 are transitive, but its source response is
already invertible. The inherited March interval family supplies new all-size
ordered-cover recurrences and a joint-port query compiler. Its actual strong
complements can leave that family. The general selected-path conjecture and
realization-constrained closure remain open.

**12:00 UTC continuation:** [DEFECT_TRANSFER D6](DEFECT_TRANSFER.md) now
computes the complete twelve-vertex 81/27 three-run response matrices:
equal mass 80614764, ranks 430/435 over F3, maximum entries 396/540.
This supersedes the earlier computational boundary about uncomputed full
matrices; three-run minimality and integral complement recovery there are
still unclaimed. D2--D4 also give exact query recurrences for bounded backward
cutwidth, retaining component pairing and the literal source relation.
