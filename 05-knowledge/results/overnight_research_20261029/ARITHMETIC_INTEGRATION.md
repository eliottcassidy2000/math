# Fifteen kernels, high ranks and the retained arithmetic coordinates

Status: **CITED / inherited exact certificates / elementary derivations**,
with open arithmetic transfers stated at their point of use. We credit the
ICARM Elliptic Curve Rank Leaderboard, its contributors and NSF Grant
DMS2425401 for the public records.

## Fifteen has three different mechanisms

The Ogg set is `{2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}`. Its property
concerns rationality of **all** supersingular j-invariants in characteristic
p, and the corresponding genus-zero modular curves. It is not the set of
supersingular primes of one fixed elliptic curve. The primary reference is
[Aricheta, Theorems1.1 and1.3](https://arxiv.org/abs/1809.07421), already
recovered in the repository's October22 and Session32 work.

For two symplectic planes over F2, there are `(2+1)(4+1)=15` Lagrangians.
A single symplectic plane has three. The concrete genus-two bridge uses
six branch points: even subsets modulo complement give a four-dimensional
symplectic space, and its Lagrangians correspond to perfect matchings of
those six points. This is the direct Richelot-kernel bridge recovered and
audited in
[incoming Session97](../session97_overnight_arithmetic_carriers_20261029/supersingular/RICHELOT_FROBENIUS_AND_FANO.md).
Frobenius must fix the unordered kernel; rationality of j does not ensure
that. At11 the incoming packet exhibits supersingular products with15
versus3 rational kernels. Equal cardinalities do not identify Ogg primes
with the kernel set.

A third fifteen is recovered from August: for curve273, fifteen good
primes give thirty independent local mod-two observations. Each
`E(F_p)/2E(F_p)` has dimension at most2, so at least15 primes are needed
for this form of rank30 certificate. The inherited full finite-group
checker realizes this lower bound, was replayed here, and its thirty
ordered points exactly match the current public curve. This is an
observation-capacity statement, not a supersingular-prime identification.
Its bank is `{43,61,101,211,223,241,263,271,283,311,313,457,521,569,571}`.

## An explicit quadratic selector turns K6 into K3,3

In the even-subset quotient V6, choose an unordered partition A,A^c into
two sets of size3. Define

    q_A([S]) = |S|/2 + |S intersection A| mod2.

Complementing S changes this by `3+|A|=6=0 mod2`, so it is well-defined.
Its polar form is `|S intersection T| mod2`. A nonzero class has a unique
two-element representative; q_A vanishes precisely on pairs crossing
the chosen partition. Therefore the Lagrangians on which q_A vanishes
are exactly the perfect matchings of K3,3, of which there are six.

This gives a specific **missing selector** behind the symplectic15 versus
split-quadratic6 distinction. There are ten possible3+3 partitions.
If instead A has size1, the same formula gives the other quadratic type:
its five nonzero isotropic vectors are the five edges through that
vertex, and no three form a matching. These six choices and the ten
split choices give all sixteen quadratic refinements.

Thus the earlier Petersen/K5/K3,3 and Hamming15 investigations have a
real contact through the K6 matching model and a quadratic refinement.
An unmarked F1 count alone loses that refinement. The preceding packet
already proves that symplectic and split-quadratic counting polynomials
agree at q=1 but have different first derivatives.

The whole fifteen-kernel geometry has symmetry Sp4(2)=S6 of order720,
which has no element of order7. This obstructs identifying a nontrivial
order-seven Ogg/Monster action with an automorphism of this marked
geometry. It does not obstruct every conceivable arithmetic relation.
The incoming genus-three extension explains where7 first enters the
branch-point model: eight points give135 Lagrangians, of which105 are
matchings and30 are extended-Hamming structures. Frobenius type(1,7)
fixes two nonmatching kernels, with cyclic Fano normalizers of order21.
Its example over F5 preserves the actual field and branch-point action.

## The same integers test different admissibility rules

The incoming all-order
[cyclic-kernel formula](../session97_overnight_arithmetic_carriers_20261029/supersingular/CYCLIC_KERNEL_RECURRENCE.md)
uses an odd cycle length n, one fixed branch point and the even-subset
symplectic module. Its number of Frobenius-stable Lagrangians is zero
if any prime p dividing n has even ord_p(2). Otherwise it is

    R(n)=2^[(1/2) sum_(d|n,d>1) phi(d)/ord_d(2)].

The proof keeps every reciprocal pair of simple binary factors; a
self-reciprocal factor obstructs the required half-dimensional isotropic
space. This is a classical cyclic-code mechanism with explicit marked
arithmetic realization, not a consequence of an unmarked count.

Consequently R(25)=R(27)=R(57)=0: ord_5(2)=4 and ord_3(2)=2.
In contrast our weighted integer joins realize H25,H27,H57, while the
fixed native simultaneous-C3 cuts retain H9 and H27 only. These are three
different operations and types. In particular H27 surviving the native
clock is compatible with the length27 cyclic module having no invariant
Lagrangian. A proposed bridge must specify which action it transports.

Among the fifteen Ogg primes the admissible cycle lengths are exactly
{7,23,31,47,71}; the ground characteristic is a separate parameter.
Lengths7 and23 yield the marked extended Hamming and Golay codes.
Lengths5,11,13 are obstructed by orders4,10,12. Thus12 appears here as
a multiplicative-order obstruction, whereas the unit clock below uses
12 to preserve or erase a phase. The predicates are not interchangeable.

## Rank 21 and the current high rank records

The dated public snapshot has6,127 selected records, including563 with
reported lower bound at least22 and66 at least30. It includes23 records
at31 and one at32. These are selected leaderboard entries, not an
unbiased height sample. The stored receipt retains the source hash,
rank histogram and factual high-rank metadata.

The [record4239](https://elliptic-rank.icarm.cloud/curve/4239) reports32
independent points and identifies a specialization of a rank17 K3
fibration. Its exact rank32 upper claim is conditional on GRH and BSD.
The incoming Session97 local-character certificates independently prove
the lower bounds30,31,32 for records273,302,4239; their verifier was
replayed here and rejected twelve hostile mutations. Point signs are
correctly invisible to these mod-two observations. The subsequent
[phase2 packet](../session97_overnight_arithmetic_carriers_20261029/ranks/phase2/README.md)
certifies all66 listed point banks:64 over F2 and two over F3. The latter
two had incomplete binary observations, not proved point dependence.
Changing the observation prime retains their curve and points while
repairing the certificate. This full-bank verification is credited to
that packet; after integration we also replayed its complete coverage
checker, passing140,334 gates for2,005 listed point directions.

The later [joint-halving continuation](../session97_overnight_arithmetic_carriers_20261029/halving/phase2/README.md)
now resolves the two binary gaps by actual rational halves. Each named
rank30 subgroup has a 2-saturated enlargement of index8192 and quotient
(C2)^13. The pinned unchanged verifier was replayed here from exact Git
blobs, passing125,459 checks and rejecting six hostiles. The thirteen
directions are genuine global joint divisibility for those source banks;
that conclusion required exact lifts and a final full-rank local bank.
It did not follow from the original finite-character kernel alone.

This is an arithmetic counterpart to retaining content in the
[matched-clock carrier](EISENSTEIN_CLOCK_JOINS.md), but not an identification
of its single gcd coordinate with thirteen elliptic directions. Their
source groups, maps, output predicates and saturation questions differ.

The [PPVW paper](https://math.mit.edu/~poonen/papers/bounded-ranks.pdf)
predicts finitely many curves of rank greater than21, rather than none.
The balance is `5/6-(r-1)/24`, so21 is the vanishing-exponent threshold.
This is a theorem about their random-matrix model and a heuristic for
curves. Finite higher-rank records do not refute it. A proved nonisotrivial
family with generic rank at least22 would have a different implication.

[FAMILY_HEIGHT](FAMILY_HEIGHT.md) now supplies a concrete family-level
input absent from the raw leaderboard: the rank32 record fibration has
j-map degree24 and Theta(X^(1/12)) distinct curves up to invariant height
X. Neither its generic rank nor a conditional rank-jump distribution is
deduced from that count. The number21 in tournament questions still has
a different predicate; the true bridges are local characters, pairings,
Frobenius, and retained source maps.

## The fruit equation has an arithmetic lifting test

The corrected prior equation `sum A/(B+C)=4` uses
`y^2=x^3+109x^2+224x`, with exact excluded pair sums and positive chamber.
The September note's rank1 certificate, torsion6, and correction to the
79/80/81-digit triple remain the inherited starting point. Rank and
positivity must not be confused with a generator/saturation certificate.

The incoming
[fraction-cubic packet](../session97_overnight_arithmetic_carriers_20261029/fractions/README.md)
adds a genuine3-isogeny phi and a named point H with `3H=phi(G)`.
Its cyclic subgroup has the all-order rule

    nH lifts rationally through phi iff 3 divides n.

The positive example27H lifts to9G. Conversely35H has three positive
real roots in the quotient description but no Q2 lift: the residual
cubic is rootless. This is a local obstruction, not an everywhere-local
Sha class. It is a much stronger arithmetic connection than matching
the numerals27 or21 across unrelated structures. The next Sha-directed
step must find and certify local solubility at every place before a
global lifting failure can serve that role.

## The 12 and 13 annotation remains source specific

The [retained clock proof](../weighted_crt_clocks_20261028/CLOCKS.md)
connects Zagier's level-two Fricke expression to the actual unit
`u=(5+sqrt29)/2`, with `g58^24=u^12`. At13 the ring splits and
`u^12=+1`; at11 it is inert and `u^12=-1`. The sampled trace has13
states per fiber modulo13, while adding the primitive clock recovers
all169 states. This is the precise11+1=13-1=12 mechanism. The separate
Zagier Apéry recurrence with parameter11 or12 retains an accessory
coefficient and is not the same modular object.

The analogy with the new port theorem is now exact at the level of
representation design: a restricted family of probes can lose a
coordinate which a newly admitted operation reads. Neither analogy
supplies a missing arithmetic Selmer map.
