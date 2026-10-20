# Defect transitions, black rerouting and rooted Petersen reductions

2026-10-20. **PROVED:** the scoped structural results below.
**FINITE-EXACT:** the declared computational controls.
**OPEN:** a smaller overlap-preserving counting quotient, general black-core
connectivity, and an unavoidable support-preserving reduction for arbitrary
snarks. No Collatz or LRC consequence is claimed.

## The results

### Exact addition of one defect

[The budget transition](budget.md) proves
\[
 S_{k+1}=\bigcup_{s\in S_k}R_k(s),\qquad S_0=\{[K_1;1]\},
\]
where states are capped weighted clone-core orbits. Each branch retains
two actual endpoints, their offsets inside clone chains, the target-cap
refinement of every gap, and the final canonical clone mergers.
Reversing a quotient arc is not the operation: it may reverse many actual
arcs. A single actual arc changes the canonical core order by at most four,
and the bound is sharp already at order five.

This is an unrestricted structural state recurrence, with no next-layer
core census as input. It gives all fixed-budget counting series and hence
A000568:
\[
 A_k(x)=1+\sum_{s\in S_k}\frac{x^{m(s)}}{(1-x)^{h(s)}},
 \qquad
 a_n=[x^n]A_{\lfloor n(n-1)/4\rfloor}(x).
\]
It still stores individual capped weighted-core orbits. The unresolved
compression is the overlap pattern of their successor sets; a constant
division factor or a count of heavy coordinates cannot replace it.

The first two transitions recover all 20 and 2,691 states independently
compiled in the preceding packet. There are also 991 direct comparisons
with literal full-tournament arc reversals through order six.

### A grading and the third universal coefficient

[The deficiency theorem](strata.md) gives
\[
 2k+1-h=2(k-r)+(2r-t)+u,
\]
for a good canonical block order: \(r\) backward arcs, \(t\) incident
endpoints, and \(u\) empty or light endpoint gaps. At deficiency \(D\),
the nonmatching and nonunit defect motifs occupy at most \(3D\)
endpoints; all remaining defect edges form an ordered unit-weight matching.
The endpoint weights are at most \(\lfloor D/2\rfloor+1\).

This places paths and claws in a precise hierarchy. Both four-vertex
three-edge trees enter at collision excess two; a triangle enters at three.
The isolated vertex in triangle-plus-isolate is gap data, so graph
complementation does not preserve the grading.

The complete deficiency-two grammar is implemented at budgets one through
four, with independent complete-compiler equality at budgets one and two.
[The local collision classification](strata_audit.md) then counts that
stratum at every budget and proves, for fixed \(k\ge2\),
\[
 a_{n,\le k}=
 \frac{n^{2k}-3k n^{2k-1}
       +\tfrac12 k(17-9k)n^{2k-2}}{2^k k!}
 +O_k(n^{2k-3}).
\]
This counts tournaments of feedback at most \(k\), not unrestricted
A000568 at fixed \(k\). The new third coefficient follows from a proved
canonical ownership rule, rather than extrapolation from four values.

### Local black detours and the remaining blue quotient

[The rerouting proof](blue.md) connects the three distinct mixed neighbors
of every trivalent pure-blue node through a common black neighbor. Thus
none of its three blue edges is a bridge, even in the simple support graph.
Together with the inherited pure-blue classification, every non-leaf
bridge must have two mixed endpoints.

A separate boundary-preserving construction supplies black detours for
blue edges at padded mixed-center stacks. The six-vertex tournament
\(1\to R_4\to1\) defeats the proposed universal fixed-boundary shortcut.
An internal shared-edge cut repairs it. The family
\(R_4\to Z\to R_4^{op}\), with pure-blue center \(Z\), supplies the stronger
obstruction to methods fixing both original endpoints.
Neither this failure nor an abstract rewiring is a genuine non-leaf bridge
of the tournament metagraph. General black-core connectivity remains open.

### Petersen reductions beyond flower strips

[The rooted-pole proof](petersen.md) adds a graph-minor certificate to
the inherited four-port support theorem. In a colourable four-pole with
support exactly
\[
 \{(a,a,b,b):a,b\in\{0,1,2\}\},
\]
a two-colour Kempe subgraph contains vertex-disjoint paths connecting the
prescribed port pairs. A crossed path would allow a forbidden boundary
word. Deleting the other edges and suppressing these paths gives the
rooted two-wire minor.

Consequently adjacent-vertex deletion poles of suitable uncolourable cubic
donors can be inserted into two specified edges and then removed while
preserving extendability in every exterior and an explicit minor model.
Petersen insertions include both authenticated Blanuša graphs and may be
iterated. A connected colourable cubic three-pole similarly reduces to one
cubic vertex. This realizes the repaired seven-port junction with actual
graphs and ordered ports, rather than only a parity-compatible formula.

## Inheritance and concept board

The session began from the [fixed-feedback compiler](../defect_boundary_20261020/README.md),
[future-response compression](../tournament_future_boundary_20261020/README.md),
[pure-blue atom proof](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md),
and [boundary-support synthesis](../boundary_support_20261020_synthesis.md).
These were recovered from the live shared Git surface before derivation.
The incoming [endpoint transport theorem](../merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md)
arrived during the session. Its exact row-current constraint rejects our
abstract bridge rewiring; the actual black-detour certificates remain valid.
The corrected near miss was the shared-endpoint equality junction, which
violates cubic seven-pole parity; the surviving repair uses one cubic vertex.

| Lane | Object and operation | Preserved target | Lost data / required sidecar | Decisive test |
|---|---|---|---|---|
| Anchor | Capped cores; one actual arc flip | Exact next-budget state set | Actual endpoint offsets and refined gaps | Independent endpoint compiler |
| Anchor | Defect graphs; fixed deficiency | Bounded exceptional motif shapes | Global endpoint order, gap weights, duplicate ownership | Normalize all deficiency-two templates |
| Niche | Blue lines; path rotation | Explicit black detours | Hamiltonian path and boundary triple | Six-vertex fixed-triple hostile |
| Niche | Cubic poles; substitution | Boundary support plus rooted minor | Ordered ports and disjoint path certificate | Kempe-path audit and Blanuša banks |
| Wildcard | Response quotients across these lanes | Only the explicitly retained observable | Support differs from counts and object transport | Uniform-factor and abstract-wiring hostiles |

The least-used sidecar brought forward is the incidence structure behind
the matching dictionary: shared endpoints change the defect grading and
also change the cubic junction type. The map preserves marked endpoint
incidence; it does not identify tournament feedback with Tait colourability.
The new rooted-pole argument supplies a genuine operation on the cubic
side without asserting that the tournament operation transports across it.

The common useful method is to retain the boundary needed by the next
operation. Its limits are equally concrete: the budget cap theorem is
uniform over all permitted lifts, whereas the marked successor relation is
existential over lifts; Boolean pole support does not retain multiplicities;
and local black detours do not alone prove global connectivity.

## Reproduction and unresolved targets

See [VALIDATION.md](VALIDATION.md) for the checked universes, controls,
commands and proof audits. The scripts use explicit failures under normal
and optimized Python. These are research proofs and exact finite checks,
not Lean formalizations.

The next strongest targets are:

1. Extend the proved deficiency-two ownership rule to deficiency three,
   while seeking a smaller quotient that retains general successor overlaps.
2. A replacement for fixed-boundary black rerouting at short outer stacks
   and at strong mixed nodes.
3. An unavoidable rooted support-preserving pole in a sufficiently broad
   irreducible snark class. The present reductions certify their generated
   families; they do not establish that every snark admits such a move.
