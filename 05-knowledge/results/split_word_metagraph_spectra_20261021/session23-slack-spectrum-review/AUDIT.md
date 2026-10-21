# Independent audit: slack-three split-word path spectra

Status: PASS, with the inherited exact symmetry-slack-three classification
stated as a dependency. No required mathematical repair was found in the
producer's path-orbit formulas, generating functions, converse action, local
incidence bijection, exact degree support, or the added balanced-cut/matching bijection. This review made no repository
or Git changes. All review files are outside the repository.

## 1. Proof findings

### Exact path partition and quotient

Let x -> W -> y -> x be the inherited core. Since y's only out-neighbor
outside W is x and in fact y has no out-neighbors in W, avoiding y -> x
forces y to be last. Similarly x must be first. These paths are exactly
x, P_W, y, yielding one orbit and 3^t literal paths.

Every other path has the unique form P_L,y,x,P_R. Here L and R partition
W, and both induced path sequences visit the ordered SCC factors in their
inherited order. The two empty-side cases are allowed. Thus the product
sum is exact, not an upper bound. K1 contributes two local choices; C3
contributes twelve literal choices, arranged in four free rotation orbits
indexed by the cardinality j=0,1,2,3 of its left subset. Using the FULL
inherited group C3^t, independent local actions give a genuine bijection:
one direct orbit plus the product of the local split-size alphabets.
This proves H=3^t(1+2^(a+2t)) and h=1+2^(a+2t).

The excluded W=K1 boundary is essential: the core is C3, its endpoints
are no longer fixed, and h is 1 rather than 3. The empty boundary is also
outside the inherited normal form.

### Fixed-spectrum and padding counts

For e=a+2t, there are e-t letters, t of which are triangles, and vertex
order e+t+2. Hence the coefficient is binom(e-t,t), not binom(e,t).
The resulting Fibonacci total and weighted recurrence are exact.
Canonical uniqueness of the inherited middle word and of the SCC word
justify ordinary class counting without a hidden orbit divisor.

Arbitrary neutral words on the two sides preserve h and supply p(z)^-2,
p=1-z-z^3. Mirrored neutral words supply p(z^2)^-1 in the self-converse
channel. Their exponent is vertex order, not e: padding changes raw H
and |Aut| in tandem but does not change h. This is an important marker
separation in the bivariate formula.

### Converse, blue paths, and the literal square

The endpoint swap sends W to its reversed atom word. On the full split
orbit alphabet it sends j_i to s_(l+1-i)-j_(l+1-i). This conclusion does
not count fixed literal paths of an arbitrary selected anti-automorphism;
it acts on the actual path-orbit carrier. A central odd-size atom prevents
any fixed split state. An even-length palindrome has 2^(e/2) such states.
The direct orbit is fixed in both cases. The stated b formulas follow.

For W=V reverse(V), complementing and reversing the second split half
identifies the full split carrier with X times X, and reflection becomes
coordinate swap. Its diagonal is precisely the blue split states and
its free two-orbits are exactly the two-element subsets of X. This is
an equivariant bijection, not numerology.

The graph interpreted here is the fundamental paired-black graph. The
set of two-element subsets of X labels LOCAL INCIDENCES at one tournament
node; it does not make X the vertices of an external metagraph component.
Loops, neighbors and multiplicities require the nonpath-complement target.
The finite hostile in Section 3 below directly demonstrates this loss.

### Complete degree support and exact denominator consequence

The three disjoint positive-degree families within symmetry slack three are:

- Odd-core self-converse channel: k=2^(e-1), e>=2, series O_e/q.
- Even-core self-converse channel: k=binom(2^r,2), r>=1, series E_(2r)/q.
- Non-self-converse channel: k=1+2^e, e>=2, series (P_e/p^2-Q_e/q)/2.

Here q=p(z^2). These exhaust the layer because every class is either
self-converse with one of the two core parities, or non-self-converse.
Mirrored neutral padding has even vertex order and preserves core parity.
Each family is nonempty: use an odd central atom for the first, 2r
singleton letters for the second, and asymmetric outer padding for the
third. The degree sets are disjoint, including the isolated r=1 value 1.
The displayed examples D_1=z^4/q, D_2=z^5/q, D_4=z^5/q, and
D_6=(z^6+z^8)/q check exactly.

The proposed exact recurrence refinement is valid. The positive core
polynomials do not vanish at the relevant positive root. Irreducibility
of p and q makes the self-converse denominator exactly q. Coprimality
prevents the p^2 and q poles cancelling across the non-self-converse
subtraction, so its denominator is exactly p^2 q. Minimal eventual
recurrence orders are therefore 6 and 12 respectively, only for these
slack-three slices, not for all nodes at the indicated degrees.

The producer's self-contained irreducibility argument for q is sound:
f=z^6+z^2-1 reduces modulo 2 to the square of the irreducible cubic
z^3+z+1, forcing any proper monic integral factor to have degree 3.
Modulo 3 it factors into irreducible factors of degrees 2 and 4,
(z^2+1)(z^4-z^2-1), which forbids a degree-three factor. The quartic has
no roots; a hypothetical factorization into two monic quadratics has
nonzero constants 1 and 2, and its cubic and linear coefficients force
both linear coefficients to vanish. Its quadratic coefficient would then
be 0, contradicting -1. This proves the needed irreducibility without
an inherited computational claim.

### Balanced cuts, Fibonacci products, and the matching bridge

The added class-level square is valid and different from the path-level
square. For a fixed half exponent r, strict positivity of the letter
exponents makes a boundary at exponent r unique. Words admitting it are
exactly U reverse(V) with U,V of exponent r. Those missing it have exactly
one C3 letter crossing from exponent r-1 to r+1, hence are exactly
L C3 reverse(R) with L,R of exponent r-1. Consequently

    P_(2r)(z)=z^2 R_r(z)^2+z^5 R_(r-1)(z)^2.

Both descriptions are unique; reversal swaps the two word coordinates.
Their diagonal class counts give the two self-converse parity channels.
The class coordinate count F_(r+1) must not be confused with the split-path
coordinate count 2^r at one even self-converse core.

The cumulative balanced-core interpretation of
1,2,6,15,40,104,273,... is correct only with the separately adjoined
actual empty tournament. It is not the invalid empty-middle C3 expression.
The old Fibonacci identity and path-matching interpretation are explicitly
credited to metagraph_spectrum_torus_session16_20261020.md, Section 5.
The new map is an object-level realization, not a new Fibonacci identity.

For a pair of tilings of lengths r and r+1, repeatedly remove an initial
dimer from the second member and swap. An initial monomer terminates at
equal-length words of exponent e; the zero case is the unit. An inverse
from U,V at stopping exponent e is explicit: begin with (U,monomer V),
then apply (A,B) -> (B,dimer A) exactly r-e times. Thus no omitted binary
history or multiplicity is concealed. The endpoint r is part of the fixed
source universe, while e is readable from the output core. Vertex order
and Hamilton spectrum are not preserved because e varies; no matching-move,
face-flip, torus, or metagraph adjacency equivalence follows.

The added uniform sampler is also valid for its declared measure. The
numbers of continuations beginning K1 and C3 are F_e and F_(e-1), with
sum F_(e+1); recursively choosing by these ratios gives uniform ordinary
strong types at fixed e. Conditional on vertex order, the triangle count
and total number of letters are fixed, so uniform position choice is exact.
It does not sample labeled objects or converse-merged types uniformly.

### Inheritance and scope

The raw all-singleton case is not new: THM-284-linear-coefficient-formula.md
proves the single-backedge count H=1+2^(skip-1). The Structural Analysis
section of THM-115-h21-permanent-gap.md records the same historical formula.
Neither its raw-H statement nor its whole historical narrative supplies a
classification of all h. The producer properly credits this inherited part
and restricts new explicit polynomials to symmetry slack three.

No unrestricted strong fixed-h atom census, global degree-support theorem,
external wiring theorem, or synchronized-triangle product rule is proved.
The full C3^t hypothesis is indispensable and cannot be carried to the
known diagonal-C3 synchronization example at symmetry slack six.

## 2. Independent finite controls

Reproduce from the parent workspace:

    python3 session23-slack-spectrum-review/check_literal_ports.py --receipt session23-slack-spectrum-review/normal.json
    python3 -O session23-slack-spectrum-review/check_literal_ports.py --receipt session23-slack-spectrum-review/optimized.json
    cmp session23-slack-spectrum-review/normal.json session23-slack-spectrum-review/optimized.json

The transcripts normal.txt and optimized.txt also byte-match. The checker
imports no candidate or repository verification implementation, uses no
third-party package, and uses explicit exceptions rather than assertions.

Its principal independent method is different from the producer's quotient:
for each literal Hamilton path, write the entire source tournament in the
path's own vertex coordinates and retain the nonpath-edge bit mask. Two
paths have the same mask if and only if an actual source automorphism takes
one to the other. Distinct masks therefore quotient the paths WITHOUT
assuming or computing the automorphism group. Reflection of those masks
is the literal staircase involution; its fixed masks are blue.
Full actual automorphism and anti-isomorphism sets are then computed
separately as a cross-check, not used to obtain the path-mask quotient.

Exact finite universe and totals:

- Every allowed middle word of vertex weight at most 9, plus W=(3,3,3,3).
- 58 cores, maximum order 14.
- 49,076 literal Hamilton paths, 6,250 distinct path-coordinate orbit masks.
- 66 actual padded tournaments, total order at most 9.
- All fixed-e word coefficients for 2<=e<=20, including both parity channels.
- Every classified (left word, core word, right word) signature through
  total order 14: 3,967 ordinary signatures and 2,036 converse-merged
  signatures. Direct converse quotient matches all 106 nonzero
  degree/order coefficients of the stated rational degree channels.
- 69 literal local black-incidence target receipts for W=11,1111,33,1331,3113,
  using actual-vertex nonpath complementation and independent canonical
  relabeling of targets, with converse also allowed for merged destinations.
- The weighted balanced-cut/straddling partition through half exponent 10.
- All 20,736 unequal tiling pairs across r=0,...,10, including 12,816
  pairs at r=10, each decoded and reconstructed by the inverse; exact
  equality with cumulative balanced-core image sets at every r.
- 98,169 explicit checks, all passed in both normal and optimized mode.

These checks validate the scoped formulas and expose local/global losses.
They are not an exhaustive census of all tournaments of the stated orders.

## 3. Concrete hostile: local K4 indexing does not determine neighbors

Both W=1111 and W=33 have e=4, h=17, b=5, and six local black incidences,
canonically indexed by the edges of a four-element set X. Yet literal
nonpath complementation gives:

- W=1111, order 6: five distinct merged destination nodes, with incidence
  multiplicities 1,1,1,1,2.
- W=33, order 8: six distinct merged destination nodes, all multiplicity 1.

Neither tested source has a black self-loop. This last observation is
FINITE-EXACT for these sources only, not a no-loop theorem.

There is a fully replayable collision already at W=1111. Use the source
vertices 0,...,5, with endpoints 0,5 and middle singleton vertices 1,...,4.
Its adjacency rows are (30,60,56,48,32,1), where bit j records an arc to j.
The distinct local reflection pairs of split states

    {1101,0100} and {1001,0110}

have representatives P=(1,2,4,5,0,3) and P'=(1,4,5,0,2,3). Their source
path-coordinate masks are 91 and 525; the complement masks are 932 and
498. The respective target adjacency rows in those path coordinates are

    (18,20,57,51,32,3),
    (10,60,57,16,33,9).

The permutation (5,0,1,2,3,4), read as old vertex -> new vertex, is an
ordinary isomorphism between these two targets. Thus the collision does
not even need converse merging. Full per-port literal paths, source masks,
target masks, target adjacency, and merged canonical keys are retained in
normal.json and optimized.json under literal_port_receipts.

### Stronger same-order and same-scalar wiring hostile

The parent-requested bounded comparison W=1331 versus W=3113 also
finishes exactly with score-pruned canonicalization, not full 10!
permutation enumeration. Both cores have

    (n,H,|Aut|,h,b)=(10,585,9,65,9),

and hence the same 28 local black incidences. Each has 28 distinct actual
converse-merged nonpath-complement destinations, all with multiplicity one,
and no black-loop incidence. Thus their multiplicity profiles alone do not
separate them. But their ACTUAL DESTINATION SETS ARE DISJOINT: intersection
zero between the two 28-element sets.

This is a stronger loss-of-wiring witness at the same order and the same
full displayed scalar signature, not merely between different graph sizes.
Both local incidence carriers have the abstract square-minus-diagonal
form for an eight-element set, while every external destination differs.
The normal and optimized receipts include every literal source/path/target
for these two cores and compare their complete canonical merged keys.
The same-order result, like the smaller multiplicity collision, is
FINITE-EXACT and does not assert a universal disjointness or no-loop law.

## 4. Publication status

The proof, scope boundaries, provenance attribution, and exact finite
controls pass. The producer may incorporate the finite wiring hostile
with its finite status preserved. This review has not changed or reserved
any theorem ID, edited the repository, or performed Git operations.
