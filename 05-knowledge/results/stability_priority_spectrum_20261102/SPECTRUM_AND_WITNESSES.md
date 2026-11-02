# Which path and arithmetic information survives a quotient

Status: **PROVED** for the elementary propositions and exact calculations
identified below; **CITED/INHERITED** for named earlier theorems;
**OPEN** for general Hamiltonian-value realization and the remaining
weighted orders. This synthesis distinguishes a numerical coincidence
from a map preserving the requested mathematical predicate.

## Three different spectra

The square-sum **order spectrum** asks which vertex counts admit a path.
The tournament **value spectrum** asks which integers occur as H(T).
The merged metagraph asks about classes and operations on presentations.
They can interact, but neither their inputs nor their values have the same
type. In particular the square-sum exceptions15,16,17,23 are vertex
counts, whereas the tournament omissions7,21 are path counts.

The earliest recovered mechanism is [THM-002 / odd-cycle collection formula](../../../01-canon/theorems/THM-002-ocf.md),
reviewed in March2026 and attributed to the earlier Grinberg--Stanley and
Irving--Omar work. For a binary tournament,

    H(T)=sum_k 2^k alpha_k(Omega(T)),

where Omega's vertices are literal directed odd cycles modulo rotation,
and conflicts mean shared original vertices. This is also the May27
[THM-326 / universal identity](../../../01-canon/theorems/THM-326-staircase-hp-equals-ip.md).
It keeps multiplicity, not just the supports of cycles.

The old five-vertex witness is especially informative: mask40 under
lexicographic edge-bit convention has four triangles and three directed
five-cycles. All seven cycles conflict, giving H=1+2*7=15. Collapsing the
three five-cycles to their common support gives the false value11.
`verify_carriers.py` reconstructs this exact witness. Thus a genuine5/11
connection is an identifiable quotient error, not an unexplained integer
pattern. The correction already existed in [THM-209](../../../01-canon/theorems/THM-209-independence-polynomial.md)
and Session66; it is recovered here, not claimed new.

At five binary vertices the odd-cycle sum is spectral. At six, disjoint
cycle information first becomes available and the native H=13/17
cospectral pair separates it from closed-walk moments. The repaired
[THM-1780 / binary spectral boundary](../../../01-canon/theorems/THM-1780-H-leaves-the-spectral-ladder-at-n6.md)
states precisely this reconstruction failure. The June15
[THM-505 / non-spectral defect](../../../01-canon/theorems/THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md)
already sought the missing simple-cycle/overlap coordinates. Its historical
dimension guesses must be read with its later corrections. No complexity
separation follows merely from a spectral collision.

## Five-vertex stability needs native incidence

The concurrent [Session139 weighted zeta boundary](../session139_incidence_moduli_and_stability_20261102/riemann_five_eleven/WEIGHTED_ZETA_BOUNDARY.md)
sharpens the distinction: weighted tournaments first lose H from adjacency
zeta at order5, because the determinant mixes a five-cycle with a disjoint
digon/triangle. The restoring scalar is

    M23=sum_edges e q_e * tau(A[e complement]).

Its exact rational witness has equal adjacency moments but Delta differing
by1/3456, and M23 differing by1/6912. Our native replay reproduces both.
This fixed-order obstruction is distinct from the binary order-six one.

The sharp five-vertex inequality can be copied across all five-subsets S:

    sum_|S|=5 Delta(A[S])
      >=(7/8) binom(n-2,3)V+(3/4) binom(n-3,2)tau.

This is an exact incidence count: each edge and triangle has the displayed
multiplicity. But it bounds a new observable, the sum of local gaps.
It does not identify that sum with the full gap. For an ordinal sum of
R=floor(n/2) fair pairs and at most one singleton,

    sum_|S|=5 Delta(A[S])
      =binom(n-2,3)V-(n-4)binom(R,2)/16.                  (1)

Proof: on each five-set H=1 and Delta=V_S-M22_S. Each pair variance appears
binom(n-2,3) times; each pair of disjoint fair edges, of product1/16,
appears n-4 times. At n=6 the edge-normalized local sum is21/32 while
the global gap is9/16. At n=11 they are115/96 and781/1024. Both reject
the proposed direct local-to-global lower bound. Equation(1) and the
values are verified by separate native determinants and path sums.

This explains the useful connection to the quasi-Riemann work. The
[October5 preprint](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex)
uses a specified sixth-power replication identity to extract one target
from a mean-square bound. That is an analogy of proof structure: copies
need a proved relation to the target and an error budget. Its claimed
11/12 zero-free half-plane is **CITED, not audited here**. It supplies no
weighted tournament inequality or RH consequence. Our successful
[n>=238 proof](STABILITY_238.md) retains the native binary-rounding baseline
and disjoint group incidence instead of assuming a local-to-global map.

## Every modulus still does not settle the integer spectrum

The inherited [Session136 local-closure theorem](../session136_full_tournament_spectrum_20261102/local_closure/PROOF.md)
shows that strong one-apex H-values already fill every odd residue target
compatibly at all finite moduli. Yet13 is absent from that restricted
integer spectrum. The unrestricted omissions7,21 also survive; see
[THM-1370 / all-order omissions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
whose correction credits earlier external work. Thus no condition testing
membership in a finite collection of allowed H residue sets can classify
the exact global integer spectrum. Simultaneously quantifying over all
moduli without a bounded-resource or coherence condition does not fix
this: every one of the named odd targets lies in the profinite closure.

There is a stronger operational arithmetic object. Concurrent
[Session139 all-modulus context theory](../session139_incidence_moduli_and_stability_20261102/all_moduli/PROOF.md)
uses the path-cover profile p_T(c) and the factorial transform

    Psi_T(x)=sum_c (c-1)! p_T(c) (x)_c.

Two tournaments have equal H responses mod M in every homogeneous
one-port context iff these functions agree on Z/M. This retains future
responses, rather than only the present integer. At odd primes p the
ordinary profile F_T(x)=sum p_T(c)(x)_c suffices; at composite moduli this
replacement can fail. Its literal L4/L6 modulo16 hostile is retained.

The actual H=5 and H=11 modules have profiles x^4+4x^2 and
x^5+8x^3+2x. Their supports on F5 are disjoint. Their ordinal join is
therefore the zero context-function modulo5, not merely a tournament
with H=55. This is an operational5/11 mechanism with an explicit product,
and is credited to that concurrent proof. Its complete prime-state results
at3,5,7,11 do not establish complete realization at every prime power.

The new [FRACTRAN clock formula](FRACTRAN_CLOCKS.md) is a different
all-modulus result: division increases required input precision at exactly
the primes dividing active denominators. Compatible CRT input charts share
that selector. In tournament profiles factorial divisibility can remove
response coordinates; in arithmetic programs division can demand more
coordinates. Both are primewise laws, but their operations and directions
of information loss differ.

During final integration, the concurrent [guarded-horizon theorem](../stability_priority_moduli_20261102/FRACTRAN_PRECISION.md)
supplied the general exact finite-clock modulus and independently recovered
the same one-step formula. It minimizes over realizable priority words,
improving the sufficient iterated clock. Its reciprocal family has h+2
trace classes and returns periodically; our [1/2,2/3] family has
1+h+floor((h+1)^2/4) trace classes and decreases to halt. Neither family
supplies a generic all-time finite carrier.

## The 81/27 metagraph control is stronger than a count collision

Before merging converse classes, a tournament's ordered Hamilton-path
presentations modulo isomorphism number H(T)/|Aut(T)|. The automorphism
action on an ordered spanning word is free: fixing it fixes every vertex.
This is the presentation-fibre mechanism of
[THM-4729 / closed and open tilings](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md).
Literal nonpath-arc complementation pairs the closed and open presentation
states; quotient incidence must still be retained.

The native twelve-vertex pair has common labelled cyclic-triple set of
size67 and common H=462915, but automorphism orders81 and27. Thus its
ordered-path presentation fibres have sizes5715 and17145. This exact
threefold discrepancy survives both of those scalar summaries. It blocks
a metagraph carrier based on H and cyclic triples alone, including any
claim that such a carrier reconstructs A000568 without additional data.
We replayed the literal automorphism action, not a guessed group order.

Odd H, parity of vertex order, path-profile parity F(-x)=(-1)^n F(x),
and complement symmetry are separate observations. The new exact
even/odd stability recurrence describes the fair-pair extremal locus;
it is not an identification of those four parity notions.

## A precise tournament abstraction for square-sum paths

**Proposition.** For any finite labelled simple graph G, its oriented
Hamiltonian paths are in bijection with labelled transitive tournament
orders whose consecutive cover pairs are edges of G. The map takes a
path word to its total order; its inverse reads that unique total order.
For at least two vertices, reversal identifies the two orientations of
each undirected path.

This creates the requested abstraction while preserving the acceptance
test. Its vertices are the original graph vertices; the binary relation
is precedence in the candidate order; there are no ties; the orientation
gauge is the chosen direction of reading. The required sidecar is the
original adjacency predicate, and for square sums the integer labels.
Every graph, including the empty graph, has transitive tournament orders.
What starts at25 is acceptance of at least one marked order for each
square-sum graph Q_n, not the existence of a tournament.

The [Gerbicz reconstruction and exact threshold record](../ninja_pairs_square_sum_number_system_20261031.md)
credits the all-order square-sum construction. Our independent native DFS
recovers exactly1,15,16,17,23,25 among orders1..25, with explicit positive
witnesses and exhaustive failures in that finite universe. The all-order
tail n>=25 remains inherited, not inferred from this finite check.

An arbitrary finite FRACTRAN execution prefix likewise has a transitive
order on its **time indices**. This loses whether machine configurations
repeat, whether each fraction was legal, and whether it had priority.
Recovering execution requires configuration labels and the first-applicable
transition guard. Calling the chronological order Hamiltonian proves none
of those properties. For the specific [1/2,2/3] program we instead prove
termination by the native integer rank v2(N)+2v3(N).

## The earlier 12 and 13 annotation

The [quadratic unit clock theorem](../cursor_cm_phase_20261028/TRACE_CARRIER.md)
is the closest old arithmetic mechanism. For u^2=5u+1 modulo13,
u^12=1. On all169 ring states, the sampled trace has13 classes; adding
Tr(ux) restores169 because the discriminant29 is invertible modulo13.
Our checker replays every state, so this is not inferred just from the
numerals12 and13.

[Zagier, section8](https://people.mpim-bonn.mpg.de/zagier/files/tex/TracesSingModuli/fulltext.pdf)
gives j_2^*=t+2^12/t+24 for t=(eta(z)/eta(2z))^24. Forgetting t in favor
of this invariant generically identifies two reciprocal choices. That
two-choice Fricke quotient is distinct from the13-element trace fibres
above. The same section's Heegner coordinates retain beta mod2N with
beta^2=-d mod4N before the stated quotient; which labels may be forgotten
depends on the observed invariant. This is a precise source for the
clock/compatible-chart analogy, not a theorem transporting tournament
stability to modular forms. :codex-annotation{index="1"}

## Transfer ledger and live next tests

| Source | Map and preserved query | Lost data and required sidecar | Decisive control |
|---|---|---|---|
| Weighted tournament | Independent rounding plus retained matching bounds Delta | Disjoint group allocation and binary baseline must remain | Fair pairs; negative -BR term |
| Literal odd cycles | Independence polynomial at2 computes binary H | Individual overlaps needed for stronger operations | Five-vertex15 versus collapsed11 |
| Tournament | Psi mod M preserves homogeneous-context H responses | Nonhomogeneous ports, exact order and integer magnitude | L4/L6 at16;81/27 presentations |
| Ordered fractions | Least input modulus preserves one-step branch/output | Later precision and priority |30/42 charts must agree mod6 |
| Square-sum graph | Marked total orders encode Hamiltonian witnesses | Integer labels, adjacency and distinct visits | Repeated1,8 walk; empty graph |

The next substantial targets are: use actual incidence distributions for
weighted orders10..237; classify realizable factorial-transform functions
at prime powers with compatible source witnesses; and seek bounded-resource
integer realization beyond profinite spectrum density. None is replaced
by a numerical coincidence, a chronological order, or a scalar average.
