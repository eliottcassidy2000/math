# Scale credits, finite-horizon coverage, and tournament compression

Research checkpoint, 2026-10-21. **PROVED** scoped bounds, representation
lemmas and the three-presentation classification; **FINITE-EXACT** supplied
certificates, enumeration and timings. Universal 153N coverage, Collatz
convergence and complete metagraph adjacency remain **OPEN**.

## What moved

1. A uniform, cancellation-free native energy bound now covers every fixed
   arrival horizon at sufficiently large scales. An explicit all-depth
   sufficient region is `N >= 2^(8d+12)`; sharper six-band certificates give
   depth one for every `N>=1024`. This changes a finite census into an
   infinite scale result, with the depth restriction retained.
2. Two new moving-parent credit families have linear-sized supply. Their
   signs are valid at every cut, even for parents that have not arrived:
   the contribution is then zero. At N8192, 342 such pairs plus 38 early
   frontier pairs complete the remaining 159 intervals without changing
   the prior 16,983 completions or counting any credit twice.
3. Every one of those 380 base sources now has a checked paid child and
   ROOT suffix. The 374 supplied parent words generate infinite guarded
   three-scale rays, keeping sources, mod15/mod35 and separate native/odd
   clocks. There are 360 first-step and 20 second-step payments.
4. Strong tournaments with exactly three Hamiltonian-path automorphism
   orbits are now classified in every order: only the cyclic R5 occurs.
   All ordinary three-presentation types are ordered K1/C3 stacks with
   exactly one R5. Their generating function is `z^5/(1-z-z^3)^2`.
5. Exact integer normalization and packed future responses accelerate the
   inherited A000568 computation without changing its mathematical states.
   On this host, bound 200 took 7.94s versus 50.30s; pair-free cores through
   100 took 2.01s versus 15.96s. These are single-run local timings, not
   asymptotic complexity claims.

## Inheritance and working board

Anchor: source-owned Collatz coverage. Niche: fast exact tournament
enumeration. Wildcard: small presentation fibres as a route to metagraph
classification. LRC(14) remains the repository anchor but no artificial
transfer to it is claimed in this packet.

The session inherited origin/main through `5ce22d4ecb75`, including the
audited black-deficit center bounds, current pure-blue classification,
balanced choice-network quotients, and conductor/fibre sidecars. A second
fetch before filing found no newer incoming commit. The
[preceding joint-flow packet](residual_flow_fibre_clocks_20261020k.md)
is the closest source-owned mechanism: N4096 was completed there, while
N8192 retained 159 sibling-bank intervals. We do not claim those older
results as new.

| Live concept | Invariant / operation | Missing information / hostile | New consequence |
|---|---|---|---|
| Guarded inverse rays | `s -> 64s+35` or `64s+49`, clock +6 | Phase rotates; fixed parents supply only boundedly many annulus addresses | Three phases restore, but fixed-parent rays cannot give linear supply |
| Moving-parent credit | One exact source pair and birth cut per owner | Negative eventual supply may arrive too late | Linear families plus an explicit early supplement at8192 |
| Inverse-word clock cost | Clock increases in steps of6 for a fixed branch | Long-horizon multiplicity remains uncontrolled | Polynomial-in-log-scale count at fixed depth and an explicit logarithmic-depth cone |
| ROOT word compilation | Two literal odd steps reconstruct a primitive edge | First discarded valuation must authenticate | 380 grounded base owners and all guarded lifted rays |
| Burnside history denominator | Odd part of N! clears every raw history | Aggregated state is not a class count; excess is not a cycle | Exact integer updates with lossless packed responses |
| Hamilton-word score classes | Cycle-cut orbits and intrinsic shore words | Cycle-cut index need not be odd | Four distinct words exclude index2, completing tc3 classification |

These are connections of mechanisms, not identifications of different
objects. In particular, the Collatz branching alphabet and tournament
response keys have different semantics; both become safe to compress only
after the discarded coordinates are proved irrelevant or retained as a
sidecar. No tournament is imposed on tied or nonbinary data.

## Source credits: a useful positive result and its exact obstruction

The [scale supply packet](scale_credit_supply_20261021a/README.md) proves:

* `ceil(3N/2)<=u<2N`, `u mod72 in {9,27,63}` admits
  `s=(16u+3)/3`, a primitive edge of clock4 and interaction -3.
* `N<=u<ceil(3N/2)`, `u mod72 in {17,35,71}` admits
  `s=(64u+19)/3`, a primitive edge of clock6 and interaction -6.

Each has exactly N/48 parents when144 dividesN. The literal parent and
source, first arrival, terminal, mod15/mod35 and clock are retained; exclude
a new source already in the stopped frontier. Before the parent arrives,
its child cannot have arrived, so the pair contributes zero. After the
child arrives the sign is the stated negative value. No conjecture about
the parent is needed for this **sign** identity.

Timeliness is different. At8192 there are171 parents in each family and
513 units of eventual -3 credit, but only8 units at cut1. Forty-six residual
intervals meet that cut. The exact cut obstruction leaves38 after the
optimal linear-bank allocation. Thirty-eight distinct early frontier pairs
remove it, completing all159 residual intervals. The -6 owners are split
into exactly two declared units; assignments never spend either unit twice.

For N>63665, every linear-family parent lies beyond the fixed frontier, so
none supplies depth-one credit. Conversely, a fixed627-parent bank has at
most1254 inverse addresses in any64-annulus. Each address has three-phase
zero-sum interaction under repeated scale lift; a sign seen at one scale
cannot be asserted at all three. These facts stop a tempting but invalid
uniform extrapolation of the8192 completion.

The [child packet](credit_children_20261021a/README.md) discharges ROOT
grounding for these actual owners. It does not discharge a uniform arrival
deadline or turn a finite parent bank into an all-source family.

## The residual is bounded on an infinite region without spending credits

Let A be arrived new-source norm, M unarrived new-source norm, and
I=R'+J' the full interaction after ownership reallocation. Then
`Delta E=A+I`, `A+M=153N`, with every negative pair appearing once in J'.
For `2^K<=N<2^(K+1)` the
[finite-horizon proof](finite_horizon_budget_20261021a/README.md) gives

\[
B(K,d)=\sum_{r=0}^d2^r{\lfloor(K+8+r)/6\rfloor+r\choose r},
\qquad R'\le1881B(K,d)^2-A-J'.
\]

When `1881B²<=153*2^K`, this is `R'<=M-J'`, exactly sufficient for the
original153N budget. It is not the stronger `R'<=M`. The proof overcounts
guarded inverse words by their integer clock cost, then uses the Q-norm
triangle inequality. It never relies on arrival of an unobserved source.

The explicit region `N>=2^(8d+12)` handles an unbounded logarithmic horizon
as N grows. The sharper fixed-depth thresholds begin `(d,K)=(1,10),(2,16),
(3,22),(4,28),(5,34),(6,41)`. For each fixed d this covers an infinite tail
of scales. The unresolved region is large d relative to logN. Replacing
these quantifiers with a bound at every finite d for one fixed N would be
an unsupported step.

## Tournament classification and what the count forgets

The [three-presentation packet](metagraph_classification_20261021a/README.md)
combines audited atom/center classifications with a new four-word
obstruction. For a step-two-invariant Hamilton cycle on2m vertices with m
odd, the two rotation orbits have different full scores. Thus every
automorphism preserves the two shores. Alternating words and the two
concatenations of their internal Hamilton cycles give four distinct
Hamilton-word orbits. This handles the index-two case that the earlier
rotation argument could not discard.

Together with the cited Hamiltonian bypass for order>=7 and the current
pure-blue/deficit-two classifications, this proves strong tc3 iffR5. The
ordered SCC word then reconstructs the complete classified stratum. The
packet keeps all three literal Hamilton words, endpoint pairs and target
adjacencies; node counts alone do not recover them.

The full degree-three node generating function in the fundamental merged
black graph now has only one remaining center polynomial:

\[
\frac{z^5A(z)^2-z^5A(z^2)}2+A(z^2)R_6(z)
  +z^{10}A(z^2)^2(1+z+z^3),\qquad A(z)=\frac1{1-z-z^3}.
\]

Here R6 counts strong self-converse deficit-six centers; the inherited
order bound makes it finite of degree<=4096. Its coefficients are still
unknown. Adjacency, ports, bridges and black-core connectivity remain
separate obligations. HYP-9352 now records the scoped promotion, the
previously proved pendant law, and the refutation of its quasi-polynomial
expectation by exponentially many degree-zero nodes.

The [enumeration packet](tournament_enumeration_20261021a/README.md) is an
orthogonal improvement: multiplying raw Burnside histories by odd(N!)
makes every update integral, and that is the least uniform scalar for
those individual histories. Packed future-gcd responses are lossless with
stage/degree/bound/mode retained. The resulting2,994,814 transition
contributions through200 are exactly the old ones. Thus this improves the
computation of A000568 and pair-free-core coefficients, but does not
classify the individual objects those sums count.

## Audits and next decisive obligations

All five packets have normal/optimized immutable replay and independent
proof/API review. The finite controls include actual sources, deliberately
bad guards, wrong clock/word inputs, labelled tiny tournament universes,
independent partitions/orientations, and the index-two hostile. Before
promotion, peer review repaired parent-word authentication and a duplicated
internal-cycle sidecar; the known finite claims were replayed after repair.

| Packet | Body checks | Scope |
|---|---:|---|
| Scale supply |157146|source-owned finite completion and all-scale guarded formulas|
| Finite horizons |20285|all-scale fixed-depth proof and forward-prefix controls|
| Paid children |545725|380 owners,374 parent words, guarded rays|
| Tournament enumeration |3243955|all modes, independent controls, exact coefficients|
| Three presentations |21964|all-order proof with inherited inputs and literal finite receipts|

The strongest next Collatz obligation is now an arrival-sensitive estimate
for the moving-parent count or for the unresolved long-horizon interaction,
with each literal pair owned once. More negative eventual mass alone will
not settle it. The direct clock bound can be sharpened by retaining source
bands, terminal size, and compatible residues in the composition count;
it must still handle the depth/scale quantifiers honestly.

For the metagraph, classify strong deficit-six centers and transport their
actual ports through the ordered padding, then test connectivity claims.
For enumeration, seek a summation that removes future-response distinctions
mathematically; another repacking cannot reduce the state count. These are
precise remaining targets, rather than claims of universal coverage or a
complete metagraph theory.
