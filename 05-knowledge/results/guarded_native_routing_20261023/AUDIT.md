# Proof and exact-computation audit

23 October 2026. **PASS in the declared scopes.** This is a same-agent proof
review with separately constructed exact matrix, flow and port controls.
It is not an independent-subagent audit. No new theorem ID is reserved or
silently promoted, and no general Collatz or unrestricted wiring conclusion
is inferred from the packet.

## 1. What is proved, and what is supplied

PROOF.md proves polynomial congestion and bounded-cost native realization
for the specific pinned cyclic-cell binary-word carrier. Its primitive
moves were already authenticated by the incoming coupled-phase theorem.
The additional residue guard is explicitly a guard on the integer number
of ones. It is not claimed to equal the original electrical phase.

The compiler keeps the source host and actual owner labels, actual endpoint
cell pins, every integer bit increment, augmented phase, source weights,
common-phase kernel lifts when present, converse, and paid event/work
clocks. Unpinned common phases are retained as actual lifts; they are not
divided out. Extra internal pins give a product only when the guard really
factors into interval guards. An indexed arbitrary owner pin, a correlated
guard, or a finer target mark needs substitution into the accepted family
before any product or connector conclusion is used.

The native host is nonempty for the stated odd q and odd number of cells;
the finite programs reconstruct literal Hamiltonian paths. Empty residue
classes are excluded, singleton fibres are treated separately, and the
short-word range r<m has its own bound. The host need not have prime q:
q=9 is an explicit positive control. A modulus-six count guard and a
modulus-eight count guard on ternary cells retain augmented moduli six and
24, respectively. Their names do not authenticate electrical guards.

## 2. Proof checks that matter for the quantitative conclusion

The repair length m-1 is enough for every directed pair of augmented phases
because binary increments can realize every residue in that many steps.
The moving splice never mixes individual owners between cells: each whole
cell still contains each of its q actual owners exactly once. The same
construction on arbitrary Hamiltonian path prefixes could repeat or omit
owners; no such unguarded splice is claimed.

The complementary path is itself accepted and reconstructs both original
endpoints when the two bounded repair words are retained. This is an
injection for each stage. A non-idle phase edge belongs to at most m-1
moving stages; the two ends are counted separately. Repeated appearances
are counted as occurrences rather than erased. This justifies the
congestion estimate, which is stronger than the inherited connectivity
and diameter statements.

Local sorting only occurs within m differences. For nonzero block weight
residue the integer weight is unique; for zero residue the only choices
are the two constant blocks. A directed primitive edge cannot repeat in
one expansion: the two sorting legs have opposite inversion directions.
The number of possible containing windows is at most m, and the number
of endpoint completion pairs within one window is at most C_m^2. These
are deliberately generous bounds; they do not assume an optimal route.

The stationary-demand compiler uses a CONSTANT length majorant D, rather
than multiplying an empirically averaged length by four. The first and
second leg marginals are uniform all-pair demands by row sums and
stationarity, even for a nonreversible return kernel. Reversed paths have
the same native capacities. Whole-route length is at most 2D, so its
length-weighted congestion is at most 4A. For positive nonuniform weights,
the exact Metropolis conductance and variance comparisons give the stated
factor b/a. No depth-independent weight ratio is presumed.

Event length, changed-cell work, changed-owner work, and original excursion
length are four distinct retained quantities. Changed-owner work is q
times changed-cell work. The compiler proves paid dilation, not equal
duration. One atomic seam flip is an authenticated coupled move, not a
claim that every partial one-cell rotation is accepted or an old elementary
metagraph edge.

## 3. Exact test universe and independent constructions

All checks use integers or rational numbers and explicit exceptions, so
Python optimization cannot disable them. Frozen normal outputs are replayed
byte-for-byte in both normal and optimized modes by reproduce.py.

| Program | Declared finite universe | Checks |
|---|---|---:|
| routing.py | 49 native families, 65,373 endpoint pairs; 48 complete pair universes and four selected pairs in the modulus-24 case | 10,054,411 |
| audit.py | independently constructed bit-word transition matrices and weighted demands, four small families | 13,976 |
| electrical.py | all 1,866 actual paths of four specified order-eight sources, all nonempty endpoint/phase and target-refined fibres | 5,411 |
| ports.py | all 120 permutations of the five-state fibre; three selected 22-state permutations and two exact nonreversible kernel decompositions | 10,239 |
| carry.py | cyclic cells M=3,9,15,27,81; four complete pair universes and three selected pairs at M=81, with all 2187 paths there checked | 635,471 |

The routing receipt checks every native state even in its selected-pair
modulus-24 case. It records exact congestion only for COMPLETE pair
universes. Transfer counts are checked independently against binomial
residue sums, with rank/unrank controls at r=128 and r=1024. Deep transfer
checks do not stand in for a full deep route census.

audit.py constructs moves from integer bit masks, separately from the
producer's phase-word proposals. Exact rational Schur elimination checks
the Poincare matrix inequality and the PSD lazy chain, both uniformly and
with nonconstant weights. The producer routes are consumed as certificates
and tested against that independent matrix. The weighted flow audit
enumerates all x,y,z triples for the 5-, 22- and 7-state controls and checks
both weighted marginals and the full microstep length. It performs another
64,076 producer-certificate checks, reported separately rather than silently
added to its own count. These finite matrix checks support, but do not
replace, the all-r proof.

ports.py independently constructs the ternary bit-move graph. Its exact
augmenting-matching decomposition includes a nonreversible kernel, so the
proof is not inadvertently limited to symmetric return laws. Every edge
exchange preserves bijective occupancy at each clock, and the entire
terminal permutation is checked. A scarce-target demand fails the matching
test. PORT_PERMUTATIONS.md proves the general unit-capacity statement;
its cost depends on |X| and is not advertised as polynomial in r.

electrical.py uses the committed native atlas and exact potential as
inputs. It does not independently prove their full enumeration. It rebuilds
the four source matrices, enumerates every Hamiltonian path, authenticates
both merged owners and converse using the original forward-bit convention,
and retains actual outside-clock tuples in every move. The saved minimum
window widths include lower-width disconnected components. Full-window
refresh is a trivial upper bound, so the packet does not misinterpret it
as all-size locality.

## 4. Incoming work changes the next obligation

Commit a201515762's [native-move audit](../native_guard_moves_farey_20261023/farey/AUDIT_NATIVE_MOVES.md)
authenticates strong-four transpositions and two selected-port rectangles.
Its exact original electrical census overlaps ours; that overlap is credited
in ELECTRICAL.md. The different costs are not conflated: 105 fibres require
six CHANGED OWNERS in the incoming menu, whereas 216 require a CONTIGUOUS
window of width six in ours. The incoming higher-dyadic debts and scarcity
of distinct actual paths defeat amplitude doubling of the two rectangles.

The new permutation compiler answers a narrower ownership question: a
uniform-stationary return on a COMPLETE accepted fibre can have distinct
terminals on every trial and can be implemented by legal exchanges. It
keeps the occupied multiset fixed. It neither supplies missing selected
ports nor cancels a nonzero current at another target. That distinction is
the map from the new routing result to the incoming obstruction, with its
lost information and required capacity sidecar explicit.

Commit 96981b1561's [Newton compiler audit](../session44_guarded_native_moves_20261023/audit-newton/AUDIT.md)
authenticates all-ring table coordinates and the cyclic-wrap relation
D^3=-3TD, including the mod-nine setting where the earlier invertible
Fourier division cannot be used. The underlying guard table can therefore
be retained exactly in either setting, but this provides no connector
between accepted actual states. Its Newton-index-span hostile also warns
against identifying coefficient indices with a minimal observer.

The original potential metadata had a LF/CRLF hash discrepancy. Both exact
views were identified, and 96981b1561 repaired that incoming metadata.
Our script pins the mathematical LF bytes, explicitly verifies the old
CRLF view, and accepts only those documented hashes. No mathematical claim
was changed. This is not a reason to weaken other raw source hashes.

### Final incoming extension: paid native carry and universal target guards

The later commits 3b3727544c and 45c5d81756 supplied stronger authenticated
objects during this session. We read the coupled-move theorem, universal
electrical-guard theorem, native lift-holonomy theorem and second review,
spectral/continuation theorem, and endpoint-observer supplement.

The all-character electrical result identifies merged targets at every
fixed source in the complete order-eight atlas. It upgrades the meaning of
our target-fixed finite test: the gap bound 1/6 applies while retaining the
whole original modulo-eight harmonic bank. The one-character test still
has only its stated weaker guard. General source/target reconstruction and
connectivity do not follow from this finite equivalence.

The native coupled-packet theorem supplies a genuine three-slot assignment
cycle with heterogeneous eligibility. We added its abstract eligibility
pattern as a negative control for the permutation compiler. Pair exchanges
are valid in our common accepted family; they are not asserted to connect
those differing slot families. The incoming joint-lift obstruction
delta-D(r) in B/D(K), and its retained kernel multiplicity, remain the
correct first test when shared owner values are not yet fixed.

CARRY_GAP.md is the positive quantitative response to the incoming native
holonomy example. The literal rotating-cell chart is exactly
Z/M x Z/(M/3), not a ternary quotient. Its two authenticated generators
change separate product coordinates. This proves the full character
spectrum and the minimum of the two cyclic gap terms. The all-size family
retains all class-based electrical guards because the full source/partner
isomorphisms are supplied, not because a numerical guard was guessed.
The claims concern the complete orbit of the DECLARED rotation subgroup,
not an unproved classification of all source automorphisms or all paths.

The normal chain uses probability 1/8 for each signed generator and is
PSD by at least one-half holding. Its exact gap is
(1-cos(2pi/M))/4, with rational lower certificate 2/M^2. The second factor
is absent at M=3; its nominal moves are identities and retain their clock
cost. For general positive generator rates the slower kernel term is kept
explicit, so vanishing rates cannot be hidden by a full-rank statement.
Shortest signed displacements on the two odd cycles give the diameter and
translation-invariant directed load used for the paid return bound.

carry.py independently constructs the actual source and literal partner
matrices, checks every path in each declared chart, and tests the two source
automorphisms on every pair of owners. Partner and converse squares are
checked for all states through M=27 and four selected states at M=81.
All pairs are routed in the first four controls; only three pairs are
routed at M=81. The spectrum is an all-size character proof, not a claim
that these finite route checks establish a spectral identity by experiment.

Incoming spectral compression and Jacobsthal work also confirms the clock
boundary: a compact address need not have native execution depth comparable
to its length. Our transfer sampler retains the actual integer position,
and our event bounds are paid in that clock. We do not transfer a rational
two-state gap to an algebraic higher-order spectrum.

The final concurrent commit 4c4f68ba69 arrived during publication. Its
stable-padding theorem and rectangle-bank theorem were read before the
second clean rebase. Endpoint-fixed partner transports now extend the carry
family to every n>=8, with gap at least 8/(n-2)^2; this is a proof corollary
of the incoming closure, not a newly enumerated padding census. The incoming
packet separately supplies constant-size two-slot cross-orbit moves at
every n>=8. Their existence is already proved; what remains is uniform
coverage/connectors, not finding any cross-orbit move at all.

The final rectangle bank strengthens the earlier obstruction: the three
relevant source/target charge differences have exact order four in the
original modulo-eight charge quotient, so paired black rectangles need
even amplitude to preserve the entire bank. The actual scarce sides have
capacity one. This already excludes a nonzero allowed move at modulus
eight, without waiting for a modulus-256 debt. We retain that stronger
incoming result. The change to the universal-guard document itself is a
provenance clarification about repository input paths, not a theorem
correction. Our frozen earlier proof pin remains explicit.

## 5. Hostiles and precise remaining questions

- A guard accepting only 0^r and 1^r has constant automaton width but needs
  r(q-1)/q changed phases for a connecting move. Complete-state connectors,
  not just small encodings, are the necessary extra hypothesis.
- A one-slot ternary repair cannot connect phase 1 to phase 0. The m-1
  repair is supported by its residue arithmetic, not a guessed short path.
- For arbitrary owner words, a prefix/suffix splice may duplicate owners.
  A native connector certificate must include those resources and ports.
- For selected-port wiring, aggregate mass can exceed actual distinct
  witness capacity. Full-fibre permutations do not create more witnesses.
- Weights must remain positive and their ratio must be paid. Weighted
  mass is not automatically a collection of authenticated unit-copy ports.
- A native Collatz correspondence and universal entry estimate remain
  OPEN; this packet proves no universal Collatz convergence statement.

The most direct next target is a full electrical-state connector certificate
with bounded legal support, retaining indexed owners, target marks and
higher-current debt. A second target is a polynomial-depth distinct-port
schedule for the cyclic-cell carrier. They are different questions, now
separated from the already proved polynomial single-demand routing bound.
