# Session126: prime registers, divisor balance, and retained symmetry

**Status:** audited component proofs and exact finite computations; classical
inputs credited. No claim of historical novelty, a five-family equivalence,
or a solution of a global tournament-spectrum conjecture.

The user's ordered fourteen fractions are Conway's PRIMEGAME. The first ten
primes are its register labels, not ten interchangeable metagraph vertices.
This session makes the requested connections explicit and tests where they
stop preserving structure.

## 1. The exact connections

- [Prime-register dynamics](primegame/fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md):
  four data counters and six control primes, with first-enabled priority.
  The first ten output exponents are2,3,5,7,11,13,17,19,23,29.
- [Actual radius-one CA](primegame/fractran_ca/ACTUAL_RADIUS_ONE_CA_REPORT.md):
  an explicit finite-alphabet local rule simulates the program on every
  positive integer, including spectator primes. Equality holds at READY
  observations, not at every microstep. The construction is not minimal.
- [Complete balance classification](balance_dynamics/README.md): all538
  balanced-input types modulo outside-prime renaming;39 preserve balance
  for one step. All orbits that stay balanced forever are explicit two-cycles.
  The seed2 orbit has exactly16 balanced visits, last60 at update70.
- [Six-cell CA and the metagraph](metagraph_ca/README.md): Rule150 is the
  unique interacting reflection-symmetric linear nearest-neighbor binary
  rule commuting with complement. In the specified perimeter coordinates,
  7 maps to26,17 to63, and55 through22 to37. These are the three nonzero
  fixed patterns, not the three pure-blue merged nodes.
- [Five families](five_families/README.md) and [the value60](five_families/A5_FUS_MEETING_POINT.md):
  constructibility, Fermat products, five convex regular solids, A5's ten
  Sylow3 subgroups, and the square-sum onset set are related by explicit
  maps or separated by exact counterexamples.
- [Inherited definitions](inheritance/FRONTIER.md) recover the actual
  metagraph and F,U,S conventions before applying any dynamics.

## The direct origin of {7,17,55} in the fraction list

Exactly three program fractions have a unit numerator or denominator:
1/17,1/7,55/1. Selecting their nonunit sides recovers {7,17,55} directly
from the fixed rack. The first removes a17 control token, the second removes
a7 data token, and the third creates one5 data token and one11 control token.
This is a canonical selection rule for this program, not an arbitrary set
of labels. It does not by itself identify these three instructions with the
three balanced exponent shapes or with three pure-blue tournament fibres.

## 2. Why the six-cell model is more than a numerical coincidence

Let T cyclically shift the six perimeter bits over F2. Write A=T+T^-1
and B=I+A. Since T^6=I, A^3=A and B^3=B^2. Thus B^2 is a projection
onto the two-dimensional space of2-periodic patterns. Its four fixed masks
are0,26,37,63; each has16 preimages. The three requested masks select its
three nonzero basins, but many other triples do too.

The genuine inherited reflection R fixes two coordinates and swaps the other
four in pairs. It commutes with B; complement does too. However the ten-node
merged observer forgets too much: the pure-blue fibre {11,46,53} has B-images
37,0,17 in different merged nodes. Therefore no deterministic B-action exists
on those ten nodes. The minimum exact refinement for B, complement and
reflection has40 states, precisely the R-orbits. This is a proved repair of
the quotient, rather than an analogy between equally sized sets.

The original three pure-blue fibres are {0}, {16}, {11,46,53}. A pure-blue
node means every tiling in its fibre is reflection-fixed. It does not mean
one of three fixed vertices of an involution on a ten-point set. Indeed an
involution on ten points cannot have exactly three fixed points, by parity.
The literal masks7,17,55 are not representatives of those three fibres in
this gauge; arbitrary labels from another drawing remain a separate question.

## 3. A complete arithmetic dynamical result

For N>=2, F counts proper nontrivial divisors, U their prime elements, and S
their squarefree elements. The inherited iff is

    F=U+S  <=>  N=p, p^3, or p^2*q*r,

with distinct primes in the last shape. U is a subset of S: the equation is
not a disjoint set decomposition. Every balanced number has total exponent
count Omega in {1,3,4}.

Outside primes are spectators. Factoring N=A*M with A supported on the ten
active primes gives T(N)=T(A)*M. The balanced shapes require at most three
outside labels, so538 symbolic types cover all positive magnitudes. The
proof classifies every balance-preserving rule, not just finitely many small
integers. The only orbits staying balanced forever are

    11 <-> 13,
    11*a^2*b <-> 13*a^2*b,

where a,b are distinct primes in {2,5} or outside the active ten-prime set.
Every other balanced start exits balance within three updates. Returning
later is not excluded by this restricted-graph theorem.

For seed2, a separate divisor-loop argument proves a permanent Omega>=5
tail after update145. Checking the preceding finite window then proves the
last balanced state is60 at update70. The same argument shows N_t>=32 for
every t>=70. With Gerbicz's cited square-sum cycle existence theorem, G_(N_t)
is Hamiltonian exactly outside times {0,1,19,20,69}. The five negative cases
have direct degree-one/zero obstructions. This is an observer corollary of
a known all-order graph theorem, not a new proof of that theorem. At CA
READY time547 the observed integer is60; intermediate microsteps have no
such asserted arithmetic or graph meaning.

## 4. What the fives and tens actually preserve

The five given orders {15,16,17,23,25} are the first five nontrivial square-sum
path orders, not a complete five-element exception family. Five convex
Platonic solids is a complete classical classification. Five Fermat primes
means five currently known examples, not a proved complete list.

The five inputs split differently under different predicates:

- Constructible regular polygons:15,16,17.
- PRIMEGAME seeds joining the seed2 tail:15,16,23.
- Divisor-balanced inputs:17,23.

In the one-sided XOR model, the polynomial (1+x)^t modulo2, evaluated as a
binary integer, is the product of Fermat numbers selected by the binary bits
of t. This is Rule60, or a parity-compressed Rule90 row. It is not the raw
six-cell periodic rule above. The values15,17,51,85 are exact Fermat products
appearing in the program, but prime-label renaming preserves the machine
while destroying that literal arithmetic naming. Therefore these products
alone do not explain its prime-generation algorithm.

The least mixed balanced input60 is also the order of the icosahedral
rotation group A5. Its ten proper nontrivial divisor values differ from the
ten Sylow3 subgroups: the latter all map to the one order value3. The A5
action on those subgroups is nontrivial; conjugation on numerical order
values is trivial. The ten-edge K5 / ten-axis correspondence is equivariant,
but it cannot be transferred through the order-forgetting map or to ten
converse-merged tournament nodes without retaining additional data.

## 5. Research board and stopping boundaries

Anchor: recover the exact prime program and the metagraph before attempting
an identification. Niche: classify balance-preserving dynamics for all prime
supports. Wildcard: test the requested masks under an explicitly chosen
symmetry-compatible CA, then derive the minimal observer that actually works.

The six live concepts were prime registers, divisor exponent profiles,
marked-path tilings, quotient observers, Fermat/XOR arithmetic, and A5 geometry.
The decisive hostile controls were priority order, the unit exit, spectator
primes, a pure-blue fibre with nonconstant CA image, and nontrivial versus
trivial group actions on two ten-element carriers.

The retained-data rule common to these lanes is precise: a proposed observer
must give the same future observation on every member of each fibre. When
it fails, refine by future experiments rather than equating matching counts.
The40-state refinement is one completed example. Allowing every named
coordinate flip requires all64 tilings, even though the initial observer
has ten values. This is an interface-information result, not a general
Boolean-circuit lower bound.

Open directions left explicit: a meaningful alternative encoding of the
user's7/17/55 labels; a natural geometric action compatible with more of the
program; and extensions beyond the stated balanced subgraph. No claim here
settles general FRACTRAN orbit behavior, the completeness of Fermat primes,
LRC(14), or the global tournament Hamiltonian-path spectrum.

## Sources and concurrent work

The identification uses Conway's own1987 chapter, printedp4 Theorem1 and
pp9–10 proof: [primary scan](https://gwern.net/doc/cs/computable/1987-conway.pdf).
It lists the user's exact tail15/2,1/7,55. A different published version
uses15/14,15/2,55; the step maps already differ at14 (105 versus15).
All proofs and step numbers here refer to the exact ordered user/1987 list.

The [inheritance pins](inheritance/SOURCE_PINS.json) credit prior metagraph,
divisor-balance and marked exponent-box mechanisms. Concurrent S70 work at
adfb2645207756de07a639e8632ac8264241026e independently investigated the two
published lists and a state-free phase-entry law. Its preregistration is
not used as a proved dependency here. In particular 'phase entry' cannot
be replaced by 'every configuration with no control prime': the seed2
orbit's next state15 already supplies a hostile test of that wording.
