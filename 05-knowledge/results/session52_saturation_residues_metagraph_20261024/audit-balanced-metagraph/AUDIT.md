# Independent audit of balanced fixed-path ports and the regular-cycle lift

24 October 2026. **PROVED audit of stated identities; FINITE-EXACT independent controls.**
No tournament-class atlas, target canonicalization census, external count sequence,
or producer arithmetic engine enters the independent scalar computation.
The frozen producer receipt is identified in `comparison-normal.json`.

## Verdict and actual contribution

The all-size algebra is valid and the numerical claims through even n=12 pass
independent coefficient computations. The ordinary score polynomial itself is
not new: THM-554, `tiling-score-partition-function`, already proves exactly this
product in the reverse path gauge. THM-555, `tiling-cycle-moment-hierarchy`, and
its archived critic output already contain balanced totals 5, 157, 51949 and
204019829 through n=10. This audit found and reported that predecessor, and the
producer added the precise inherited/reproduced labels and pinned sources.
A prior engine is therefore not to be described as absent or newly invented.

The useful refinements here are the complete endpoint/color split, exact
balanced-target state truncation, n=12 census, and actual C-equivariant
cycle-marked regular extension. The latter gives Eulerian orientations on
K_(n+1) minus its Hamilton cycle, not merely a numerical correspondence.
The class census 13333 ordinary balanced n=10 classes and 229 self-converse
ones is a separately declared Session51 input. It is not computed by these
coefficient engines, and no A000568 value or full wiring is inferred.

## Independent methods and finite results

`audit.py` eliminates degree variables from highest to lowest, with all fixed
path arcs already included in the exponent state. This is a separately written
coefficient recurrence with a different state convention and endpoint handling
from the producer's first-to-last past-wins recurrence. Its cache sizes at n=10
and n=12 are 17338 and 396960, respectively. The blue polynomial is independently
rebuilt by constructing the actual free-edge orbits of reflected converse on
all n vertices and then projecting each of the two orientation monomials onto
the first n/2 coordinates. The producer's symbolic half-polynomial is not used.

The endpoint quadruples below are ordered (LL,LH,HL,HH):

- n=10: (44967565, 46961451, 67123248, 44967565).
  Total 204019829; internal 67123248; blue internal 22904;
  blue boundary 18795; total blue 41699.
- n=12: (2301770367447, 2371435845579, 3213825974560, 2301770367447).
  Total 10188802555033; internal 3213825974560; blue internal 6803648;
  blue boundary 5765829; total blue 12569477.

Binary-order literal full-matrix enumeration exhausts 1, 8, 1024 and 2097152
fixed-path words at n=2,4,6,8. It rebuilds outdegrees from matrices, rather than
using the producer's Gray-code incremental score calculation. On every balanced
word it verifies literal C, R, their commutation, the endpoint transport law,
internal iff HL, blue boundary iff LH, the regular completion, closure of the
actual specified Hamilton cycle and equality of every arc under transported C.
Literal summaries and deterministic universe digests are in `normal.json`.

## All-size proof checks

For n=2p the known transport law is

    d_C(i) = 2p-1-d(i) + [i=0] - [i=n-1].

Interior high/low statuses interchange. The first coordinate remains allowed
exactly when high and the last exactly when low. Thus HL is necessary and
sufficient, not just necessary. Adding z with z→high and low→z is the unique
regular completion. The fixed path closes through z exactly for HL. On C,
all old nonpath arcs reverse; the two z cycle arcs remain fixed; the other z
arcs reverse as their old vertices exchange high/low. Hence the map intertwines
C with reversal of every off-cycle arc. Removing the fixed cycle gives one
outdegree p-1 at every vertex of K_(2p+1) minus C_(2p+1), exactly Eulerian
orientations. Restoring the cycle and deleting z is the inverse.

For a regular parent R and deleted-vertex orbit O represented by z,
Aut(R-z)=Aut(R)_z, and every directed Hamilton cycle has one reading starting
at z. Therefore the internal card degree is

    c(R)/|Aut(R)_z| = |O| c(R)/|Aut(R)|.

`parent_cards.py` independently verifies actual permutations, automorphisms,
cycles and card canonical forms at n=2,4,6. At n=6 the three regular parents
have (Aut, cycles, vertex-orbit sizes, internal card degrees):

- (7,17,[7],[17]);
- (21,24,[7],[8]);
- (3,15,[3,3,1],[15,15,5]).

Their degrees total 60. These ordinary card identities do not say that C
preserves an unpointed parent isomorphism class.

## Rotation fixed sets and quotient guardrails

For a rotation with gcd d>1 and orbit length L, each vertex orbit contributes
(L-1)/2 internal outgoing edges at every vertex independent of choices, with
2^[d(L-1)/2] orientations. Each pair of distinct vertex orbits gives L matching
edge-orbits, except adjacent quotient-cycle pairs give L-1. Thus the producer's
coefficient formula has residual target ((d-1)L-2)/2. The d=3 cycle has three
edges, not a doubled cycle; d=1 must be handled separately because all deleted
cycle edges lie in one vertex orbit. It gives 2^((n-2)/2). Burnside and the
prime n+1 congruence follow with the stated quantifiers.

Independent actual directed edge-orbit factors verify all proper rotation
fixed sets needed at n≤12, including f9(1)=8 and f9(3)=80. Additional finite
controls, not a computation of the n=14 scalar census, give
f15(1)=64, f15(3)=22144 and f15(5)=5327104. Different steps with the same gcd
agree. These are saved in `extra-rotation-controls.json`.

Native internal edges are C-orbits, D/2. Merging class vertices does not remove
or identify lines. If one deliberately further quotients lines by R, CR has
no fixed word for n≥4 because (0,n-1) is a reflection-fixed free coordinate
and C negates its bit. The number of optional R edge-orbits is therefore
(D+B)/4, where B is the internal blue port count. At n=6 these are 18 rather
than 30 native C-lines; at n=10 they are 16786538 rather than 33561624.
Both constructions are legitimate, but they are different, and the packet
uses the native one. At n=2 C is the identity, so free-pairing formulas do not
apply.

## Hostile controls and reproducibility

`hostile.py` rejects nine plausible errors: full converse substituted for C;
reversed LH internal criterion; inverted endpoint correction signs; reversed
regular-completion z-arcs; omission of interior z-arc reversal in the lift;
free C pairing at n=2; odd internal count; odd blue-internal count; and replacing
native C-lines by optional R-edge-orbits. All checks use explicit exceptions
and survive Python optimization. Three additional packet controls reject engine-byte
tampering and forged n=10 internal or blue-boundary counts even after their
checksums have been refreshed. These are in `hostile-packet-normal.json`.

Normal and optimized independent scalar, card and hostile receipts are
byte-identical. A real filesystem copy of the independent checker reproduces
both independent scalar and card receipts. The producer's manifest, pinned
sources, full normal and optimized run, and a real-copy replay are separately
verified by the completion receipt. No symbolic link stands in for a copy.

Reproduction, from the audit directory:

    python3 audit.py --output normal.json
    python3 -O audit.py --output optimized.json
    python3 parent_cards.py --output parent-cards.json
    python3 -O parent_cards.py --output parent-cards-optimized.json
    python3 rotation_controls.py --output extra-rotation-controls.json
    python3 hostile.py --output hostile-normal.json
    python3 -O hostile.py --output hostile-optimized.json
    python3 compare.py --producer producer --output comparison-normal.json
    python3 -O compare.py --producer producer --output comparison-optimized.json

The full producer replay and manifest pins are retained in `producer/`.
An actual-copy replay alone is not an independent mathematical method; that
role is served by the separately written coefficient, matrix and orbit checks.
