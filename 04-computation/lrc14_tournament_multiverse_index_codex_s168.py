"""Build the LRC tournament-technique multiverse index.

This is a synthesis script, not a proof.  It makes the index reproducible:
the same card data emits the navigation annex and a compact Tournament
Analysis fingerprint over technique families.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Iterable


DIMS = [
    "predicate",
    "exact_scale",
    "topology",
    "endpoint",
    "arithmetic",
    "harmonic",
    "metagraph",
    "formal",
    "residual",
    "transfer",
]

CATEGORY_ORDER = [
    "proof_sheaf_zipper",
    "tournament_metagraph",
    "packet_families",
    "topology_boundary_discrepancy",
    "harmonic_fourier_certificate",
    "arithmetic_sieve_series",
    "analytic_sieve_smoothing",
    "state_lift_geometry_cross",
    "computation_formal_forum",
]


@dataclass(frozen=True)
class Card:
    handle: str
    category: str
    name: str
    vertices: str
    observable: str
    preserves: str
    guardrail: str
    pull: str
    anchors: str
    dims: tuple[str, ...]


def c(
    handle: str,
    category: str,
    name: str,
    vertices: str,
    observable: str,
    preserves: str,
    guardrail: str,
    pull: str,
    anchors: str,
    dims: Iterable[str],
) -> Card:
    dims_tuple = tuple(dims)
    bad = [d for d in dims_tuple if d not in DIMS]
    if bad:
        raise ValueError(f"{handle} has unknown dimensions {bad}")
    return Card(
        handle,
        category,
        name,
        vertices,
        observable,
        preserves,
        guardrail,
        pull,
        anchors,
        dims_tuple,
    )


CARDS = [
    c("LTM-001", "proof_sheaf_zipper", "Controlled-kernel quotient ledger", "quotient fibers / forgotten coordinates", "fiber-constancy, reconstruction, dual annihilation, residual emission", "the LRC predicate through each compression", "scalar quotient with unnamed kernel", "Attach a kernel ledger to every new reduction.", "HYP-2990, LTI-037, LTT-001", ("predicate", "formal", "residual", "transfer")),
    c("LTM-002", "proof_sheaf_zipper", "Cocycle-sheaf exactness matrix", "C0 packet labels, C1 emitted cocycles, C2 incompatibilities", "boundary map rank and named cohomology classes", "all emitted proof debt at C1", "cochain metaphor without finite matrix", "Build the C1 matrix on HYP-2963 packets.", "HYP-3004, LTI-108", ("predicate", "topology", "harmonic", "formal", "residual")),
    c("LTM-003", "proof_sheaf_zipper", "Zipper tooth DAG", "proof teeth / handoff arrows", "which tooth reconstructs, stops, descends, or emits debt", "local-to-global proof assembly", "tooth order treated as theorem", "Give every row a tooth path and stopping reason.", "HYP-2993, HYP-2990", ("predicate", "metagraph", "formal", "transfer")),
    c("LTM-004", "proof_sheaf_zipper", "Residual sector naming", "F7, THM-572, harmonic/state-lift sectors", "whether a failure has a typed residual address", "non-AP/GW failure information", "anonymous unknown bucket", "Define F7 as finite state-lift or harmonic cohomology.", "HYP-2987, THM-572", ("predicate", "residual", "formal")),
    c("LTM-005", "proof_sheaf_zipper", "Proof-obligation tournament", "certificates, obligations, proof carriers", "retention vector dominance", "which proof route is stronger for a packet", "runner tournament used where proof-carrier tournament is needed", "Add carrier tournaments before scalar rankers.", "AGENTS.md, HYP-2987", ("predicate", "metagraph", "formal", "transfer")),
    c("LTM-006", "proof_sheaf_zipper", "Robbins bridge checklist", "certificate pieces / assembly bridges", "whether deleting a field disconnects the proof manifest", "certificate reconstructibility", "floating Fejer or numeric certificates", "Add bridge fields to Fejer manifests.", "HYP-2981, LTT-044", ("predicate", "formal", "transfer")),
    c("LTM-007", "proof_sheaf_zipper", "Certificate handoff braid", "q-witness, Fejer, Ramanujan, endpoint, moment, smoothing carriers", "handoff compatibility and proof-arrow closure", "multi-clock packet proof data", "one carrier declared universal", "Prove the six HYP-2987 atlas arrows.", "HYP-2987, LTI-036", ("predicate", "harmonic", "endpoint", "transfer")),
    c("LTM-008", "proof_sheaf_zipper", "Meta-proof debt ledger", "hypotheses / theorem obligations / forum deltas", "which claim discharges which named debt", "research-route dependencies", "new idea not linked to a debt", "Add debt tags to future forum posts.", "CONCEPT-MAP, poke-forum", ("formal", "metagraph", "transfer")),
    c("LTM-077", "proof_sheaf_zipper", "Curried packet functional tower", "packet evaluators / quotient fibers / lane functions", "lost-coordinate derivative under partial evaluation", "proof-route data before scalarization", "row to scalar evaluator", "Add curried_call_signature and lost_coordinate_function to packet manifests.", "HYP-3002, LTI-152, LTT-059", ("predicate", "exact_scale", "formal", "transfer")),
    c("LTM-009", "tournament_metagraph", "Source-marked A000568 fiber", "observer-source marked tournaments", "source cone plus deleted-root class", "loneliness predicate at a moment", "unmarked isomorphism class", "Attach endpoint owners to source fibers.", "HYP-2486, THM-381", ("predicate", "metagraph", "endpoint")),
    c("LTM-010", "tournament_metagraph", "Phase-swept tournament spectrum", "isomorphism-class path over phase", "class measure, binding scale, edge flips", "magnitude-aware tournament shadow", "single apex class", "Compute spectra for AP, GW, K33, petal rows.", "HYP-2928, LTI-018", ("predicate", "exact_scale", "metagraph")),
    c("LTM-011", "tournament_metagraph", "Fixed Hamiltonian path tiling", "tie path plus above/below tiles", "tile interval and path address", "local switch support", "hidden dependence on chosen tie path", "Declare tie Hamiltonian path in scripts.", "LTT-006, THM-354", ("metagraph", "topology", "formal")),
    c("LTM-012", "tournament_metagraph", "SC-NS spine-ribs-sea metagraph", "tournament isomorphism classes / flip graph", "class adjacency and complement behavior", "deformation neighborhood", "visual color slogans without chain data", "Place LRC packet classes on the metagraph.", "metagraph parity, A000568", ("metagraph", "formal", "transfer")),
    c("LTM-013", "tournament_metagraph", "Edge-flip stress graph", "scalar rankers / row pairs", "order flips under one changed coordinate", "unsafe ranker detection", "ranker accepted before flip count", "Store flip counts beside every scalar feature.", "HYP-2984, LTI-101", ("predicate", "formal", "metagraph")),
    c("LTM-014", "tournament_metagraph", "Good-cut SCC support", "fixed-path cut intervals / SCC condensation", "n minus SCC count and support intervals", "strong-connectivity support", "bucket 1 or raw connectivity shortcuts", "Add good-cut profiles to packet spectra.", "THM-354, LTT-009", ("metagraph", "topology", "residual")),
    c("LTM-015", "tournament_metagraph", "OCF activity coimage", "odd-cycle compatibility graph", "activity sectors and coimage maps", "cycle-packing structure beyond H", "H count without observer labels", "Compare OCF activity at q=14,27,41.", "THM-002, HYP-2618", ("harmonic", "metagraph", "formal")),
    c("LTM-016", "tournament_metagraph", "Path-homology boundary tournament", "directed paths / boundary operators", "which cycles survive boundary quotient", "cycle-space homology signal", "homology rank as certificate", "Run path-homology on spectrum representatives.", "path homology notes, HYP-2995", ("topology", "metagraph", "residual")),
    c("LTM-017", "tournament_metagraph", "Walsh W-polynomial channels", "arc monomials / Fourier channels", "even/odd channel support", "channel-specific tournament data", "single scalar evaluation", "Pair W channels with Haar zeta packets.", "THM-080, THM-081", ("harmonic", "metagraph", "formal")),
    c("LTM-018", "tournament_metagraph", "Metatournament of proof routes", "claims, proof routes, failed quotients", "which route defeats another on retained data", "research-route conflict structure", "linear chronology mistaken for proof order", "Add meta-edge reasons to route conflicts.", "repo-as-tournament, LTI-020", ("metagraph", "formal", "transfer")),
    c("LTM-019", "packet_families", "Labelled packet classifier schema", "primitive rows with route labels", "q-witness/AP-GW/petal/K33/covering/Fejer/unknown", "family route", "classifier label without theorem obligation", "Extend with Haar, Ramanujan, spectrum fields.", "HYP-2963, LTT-040", ("predicate", "exact_scale", "formal")),
    c("LTM-020", "packet_families", "Family-sporadic decision tree", "F0-F7 families and sporadic rows", "route split plus terminal residual", "counterexample classification", "F7 undefined", "Make F7 concrete and finite-checkable.", "HYP-2961, HYP-2956", ("predicate", "formal", "residual")),
    c("LTM-021", "packet_families", "Skeleton gate equality atom", "qdiv=14 boundary skeletons", "AP/GW skeleton versus impostor shell", "boundary equality class", "AP-like profile confused with AP/GW", "Use skeleton gate before K33 claims.", "HYP-2960, HYP-2951", ("predicate", "endpoint", "topology")),
    c("LTM-022", "packet_families", "Exact Farey node", "M=p/q, qdiv, excess e=14p-q", "exact scale and unit-excess branch", "front-end route scheduling", "product/summand value alone", "Require M, q, e in every packet record.", "HYP-2931, HYP-2984", ("predicate", "exact_scale", "arithmetic")),
    c("LTM-023", "packet_families", "NORK pinch template", "pinch sites / apex apertures", "local pinch or boundary debt", "few-apex local obstruction", "pinch without global route", "Add NORK tags to HYP-2963 records.", "HYP-2966, HYP-2964", ("topology", "endpoint", "residual")),
    c("LTM-024", "packet_families", "Few-apex packet lift", "small lifted packet changes", "which apex changes and which proof route survives", "finite-core lift data", "wide analytic estimate before lift address", "Audit few-apex residuals for strict open mass.", "HYP-2968", ("predicate", "exact_scale", "residual")),
    c("LTM-025", "packet_families", "Boundary-gap bridge graph", "safe interval endpoints / owner transitions", "positive bridge or zero-length taut vertex", "strict witness versus boundary debt", "raw divergence", "Classify zero bridge packets by moment/state-lift.", "HYP-2965, HYP-2975", ("topology", "endpoint", "predicate")),
    c("LTM-026", "packet_families", "Source-spectrum pullback", "Farey-indexed rooted packet spectra", "source-local spectrum before quotient", "global proof data over source fibers", "apex spectrum without source packet", "Attach source-spectrum IDs to F6/F7.", "HYP-2954, HYP-2928", ("predicate", "exact_scale", "metagraph")),
    c("LTM-027", "topology_boundary_discrepancy", "Regular-open Haar-Baire split", "finite safe/danger arc unions", "strict-open mass versus closed boundary atom", "topological witness type", "boundary-only equality counted as open", "Use exact interval unions in every zero-open claim.", "HYP-2948, HYP-2949", ("predicate", "topology", "formal")),
    c("LTM-028", "topology_boundary_discrepancy", "Endpoint-owner skeleton", "active endpoint-owner pairs", "owner transfer and pair-sum pattern", "AP/GW boundary skeleton", "Haar mass without owner labels", "Prove zero-open skeleton rigidity.", "HYP-2951, HYP-2970", ("endpoint", "topology", "predicate")),
    c("LTM-029", "topology_boundary_discrepancy", "Tope-cocircuit wall complex", "endpoint arrangement topes / cocircuits", "open cell or boundary cocircuit", "wall-combinatorial certificate", "oriented-matroid words without packet labels", "Search minimal forbidden wall packet.", "HYP-2986", ("topology", "endpoint", "residual")),
    c("LTM-030", "topology_boundary_discrepancy", "Haar zipper zeta switch", "2 x 2 fixed-margin packet tables", "zeta=T00-T01-T10+T11", "mixed local cocycle", "row/column margins alone", "Compute zeta signatures on packet bank.", "HYP-2991, LTI-106", ("topology", "harmonic", "formal")),
    c("LTM-031", "topology_boundary_discrepancy", "Haar rectangle product atlas", "dyadic rectangles and product classes", "orthogonal zero, atom, owner strip, cross handoff, nested descent", "2D discrepancy handoff type", "product class without endpoint/period labels", "Test the no-nonzero-coefficient lemma.", "HYP-2992, LTI-047", ("topology", "harmonic", "transfer")),
    c("LTM-032", "topology_boundary_discrepancy", "Colored discrepancy reservoir", "mod-14 color grids / mixed switches", "independent color-compatible Haar switches", "resonance-aware discrepancy", "component count K as scalar", "Replace K by bounded switch count.", "HYP-2595, HYP-2594", ("topology", "harmonic", "arithmetic")),
    c("LTM-033", "topology_boundary_discrepancy", "Boundary-moment multichart", "exact-period denominator charts", "multi-chart moment and cap-cell pressure", "covering residual discharge", "one covered chart as obstruction", "Couple moment charts to Haar coefficients.", "HYP-2969", ("topology", "arithmetic", "endpoint")),
    c("LTM-034", "topology_boundary_discrepancy", "Taut bridge current", "positive bridges and taut zero vertices", "bridge length plus endpoint current", "AP/GW equality as current atom", "strict and taut cases merged", "Compare taut currents with cocircuit atoms.", "HYP-2975, HYP-2970", ("endpoint", "topology", "predicate")),
    c("LTM-035", "topology_boundary_discrepancy", "Kernel homotopy support radius", "open components / boundary defects", "support radius and smoothing admissibility", "stable open certificate", "kernel deformation losing labels", "Add radius preconditions to smoothing dispatcher.", "HYP-2984", ("topology", "harmonic", "formal")),
    c("LTM-036", "topology_boundary_discrepancy", "Two-dimensional Haar product rule", "Haar rectangles / staircase Walsh masks", "coordinatewise product and xor tiling", "multiplicative descent under labels", "tiling shadow treated as theorem", "Instantiate product classes on HYP-2963.", "HYP-2989, THM-351", ("topology", "harmonic", "metagraph")),
    c("LTM-037", "harmonic_fourier_certificate", "Fejer interval certificate", "packet-keyed rational trig intervals", "degree, center, atom formula, sign interval", "formal positivity certificate", "floating negative value", "Move Fejer manifests to interval backend.", "HYP-2981", ("harmonic", "formal", "predicate")),
    c("LTM-038", "harmonic_fourier_certificate", "Toeplitz PSD dual", "Fejer/Toeplitz matrices and dual vectors", "negative pairing or PSD violation", "dual certificate", "matrix scalar without packet key", "Formalize familywise PSD certificates.", "HYP-2974", ("harmonic", "formal", "predicate")),
    c("LTM-039", "harmonic_fourier_certificate", "Ramanujan exact-period projector", "exact-period denominator modes", "prime-power period labels", "period-aware harmonic packet", "squarefree divisor quotient", "Build endpoint-aware projectors.", "HYP-2979, HYP-2978", ("harmonic", "arithmetic", "exact_scale")),
    c("LTM-040", "harmonic_fourier_certificate", "Resonance-channel Fourier split", "congruence-compatible frequencies", "which Fourier channel sees color resonance", "off-scalar resonance evidence", "low-frequency shadow alone", "Merge resonance signs with Haar product signs.", "HYP-2595, HYP-2867", ("harmonic", "arithmetic", "topology")),
    c("LTM-041", "harmonic_fourier_certificate", "Krawtchouk shadow audit", "orthogonal basis coordinates", "which packet info a basis erases", "anti-scalar basis guardrail", "basis coefficient as proof", "Use as loss detector for new transforms.", "Krawtchouk reflections", ("harmonic", "formal", "residual")),
    c("LTM-042", "harmonic_fourier_certificate", "Riesz product witness route", "positive trigonometric products", "weighted test function witness", "candidate lonely-point certificate", "unlabelled global product", "Revisit with Fejer packet backend.", "T827, Riesz reflections", ("harmonic", "predicate", "formal")),
    c("LTM-043", "harmonic_fourier_certificate", "Spectral shadow dual", "low-frequency shadows and dual ranks", "cheap pre-split plus lost-label report", "screening certificate candidate", "spectral shadow as final proof", "Compare with Ramanujan projectors.", "lrc14_spectral_shadow_dual", ("harmonic", "formal", "exact_scale")),
    c("LTM-044", "harmonic_fourier_certificate", "Danger-count moment dual", "danger counts and moment functionals", "moment inequality with packet owner labels", "count-moment dual certificate", "raw count moments", "Attach moment duals to packet signatures.", "HYP-2971, HYP-2970", ("harmonic", "endpoint", "predicate")),
    c("LTM-045", "arithmetic_sieve_series", "Mobius product packet ledger", "squarefree, unit, product-support fibers", "mu, unit, prime-power side labels", "multiplicative packet data", "mu or mu^2 scalar without exact-period repair", "Run fiber-mixing on any mu feature.", "HYP-2899, HYP-2982", ("arithmetic", "exact_scale", "formal")),
    c("LTM-046", "arithmetic_sieve_series", "mu(n)/n cancellation gauge", "weighted divisor packets", "partial-sum cancellation with packet labels", "signed arithmetic cancellation clue", "global Mertens-style scalar jump", "Attach cancellation ranges to qdiv packets.", "analytic-sieve prompts", ("arithmetic", "harmonic", "exact_scale")),
    c("LTM-047", "arithmetic_sieve_series", "mu^2(n)/phi(n) density guardrail", "squarefree totient-density packets", "Euler-factor density plus lost prime-power report", "large-sieve weight precondition", "squarefree blindness", "Record which exact periods survive the weight.", "HYP-2982, LTI-086", ("arithmetic", "exact_scale", "formal")),
    c("LTM-048", "arithmetic_sieve_series", "Singular-series packet factor", "major-arc arithmetic factors", "local density by labelled packet class", "circle-method middle layer", "singular series detached from endpoints", "State packet labels in every major/minor split.", "HYP-2982, LTI-094", ("arithmetic", "harmonic", "exact_scale")),
    c("LTM-049", "arithmetic_sieve_series", "Farey mediant obstruction monoid", "Kuratowski/Wagner and packet obstruction counts", "componentwise mediant as disjoint-union density", "additive obstruction composition", "edge-density numerator/denominator swapped or scalarized", "Test AP/GW obstruction mediants under packet union.", "Kuratowski prompt, HYP-2932", ("arithmetic", "exact_scale", "transfer")),
    c("LTM-050", "arithmetic_sieve_series", "Totient curvature CRT residual", "Euler factor and CRT packets", "curvature from prime-power recombination", "why phi recurrences bend", "denominator-only quotient", "Feed curvature into exact-period ledger.", "HYP-2900, HYP-2899", ("arithmetic", "exact_scale", "residual")),
    c("LTM-051", "arithmetic_sieve_series", "Pisano Fibonacci band ladder", "modular recurrence fibers", "period band and quotient fiber", "sequence periodicity as packet mark", "next-term chasing", "Use bands only with qdiv/Farey labels.", "Pisano reflections, HYP-2523", ("arithmetic", "metagraph", "exact_scale")),
    c("LTM-052", "arithmetic_sieve_series", "Zeckendorf boundary normal form", "Zeckendorf packet decompositions", "canonical summand carry", "boundary decomposition clue", "normal form without endpoint meaning", "Try on AP/GW taut currents.", "HYP-1902", ("arithmetic", "endpoint", "transfer")),
    c("LTM-078", "arithmetic_sieve_series", "Summand multiplicand Farey fiber merge", "Farey packets with additive antidiagonal and product hyperbola fibers", "which operation graph preserves the witness", "typed p+q and p*q proof currencies", "sequence shadow collapsed to one scalar", "Add summand_fiber_id and multiplicand_factor_fiber to packet records.", "HYP-3003, LTI-153", ("arithmetic", "exact_scale", "predicate", "transfer")),
    c("LTM-053", "arithmetic_sieve_series", "Euler Glaisher Witt parity transform", "partition parity channels", "doubled or forgotten boundary labels", "GF(2) channel discipline", "visible color parity slogans", "Test exact packet fibers for parity channel flips.", "HYP-2685, metagraph parity", ("arithmetic", "metagraph", "formal")),
    c("LTM-054", "arithmetic_sieve_series", "p-adic residual tree", "prime-power carries and CRT recombination", "carry owner, shell height, recombined address", "local-prime residual split", "squarefree or denominator-only quotient", "Build tree for q=14,27,41,63,84,168.", "THM-568, HYP-2929", ("arithmetic", "endpoint", "residual")),
    c("LTM-055", "analytic_sieve_smoothing", "Major-minor arc labelled split", "arc regimes and kernels", "which packet label enters each regime", "circle-method architecture", "arc estimate with lost endpoint/state labels", "Split late walls into labelled regimes.", "HYP-2982, HYP-2983", ("arithmetic", "harmonic", "formal")),
    c("LTM-056", "analytic_sieve_smoothing", "Large-sieve controlled precondition", "Selberg/large-sieve weights", "upper-bound weight plus exact-period survivor set", "screening estimate", "sieve weight as final equality carrier", "Log surviving exact-period fibers.", "HYP-2982, LTI-068", ("arithmetic", "harmonic", "formal")),
    c("LTM-057", "analytic_sieve_smoothing", "Upper-bound quadratic sieve analogue", "quadratic residue buckets / pair tensions", "which pair channels remain after upper bound", "second-order sieve prefilter", "quadratic average losing endpoint owners", "Compare pair tensions with C27/K33 handoffs.", "sieve prompts, HYP-2982", ("arithmetic", "harmonic", "endpoint")),
    c("LTM-058", "analytic_sieve_smoothing", "Kaczynski approach-class boundary", "boundary functions and approach classes", "which limiting approach preserves a packet", "exceptional-set typing", "exceptional set as anonymous null set", "Define approach class for true-wide residuals.", "HYP-2983, HYP-2679", ("topology", "harmonic", "residual")),
    c("LTM-059", "analytic_sieve_smoothing", "Admissible smoothing dispatcher", "smoothing policies / proof clocks", "policy preserves certificate or emits boundary defect", "kernel route selection", "global smoothing before packet labels", "Prove dispatcher lemma over HYP-2963.", "HYP-2985, HYP-2984", ("topology", "harmonic", "formal")),
    c("LTM-060", "analytic_sieve_smoothing", "Saddle/explicit-formula debt", "saddle points, parabolic cylinder approximants, explicit formulas", "analytic asymptotic plus residual bound label", "late analytic proof clock", "complex estimate detached from packet manifest", "Make every saddle term name its packet residual.", "Goldbach/Kaczynski prompts", ("harmonic", "formal", "residual")),
    c("LTM-061", "analytic_sieve_smoothing", "Exponential-sum core checksum", "finite exponential sums over packet labels", "phase cancellation with endpoint and period tags", "analytic core invariant", "continuous smoothing without discrete checksum", "Pair sums with Ramanujan projectors and Fejer atoms.", "exponential-sum prompts", ("harmonic", "arithmetic", "formal")),
    c("LTM-062", "state_lift_geometry_cross", "K33 state-lift residual", "K33 incidence and finite state lifts", "hidden lift or forbidden state construction", "p>=3 unit-excess route", "positive row without state address", "Define lift data emitted by failed carriers.", "HYP-2908, THM-572", ("residual", "arithmetic", "transfer")),
    c("LTM-063", "state_lift_geometry_cross", "C27 unital shell transfer", "antipodal C27 shells and branch-local unital charts", "hole/double owner and pair completion", "p=2 petal/two-block route", "global unital claim without branch chart", "Formalize p=2 transfer rigidity.", "HYP-2937, HYP-2942", ("endpoint", "arithmetic", "transfer")),
    c("LTM-064", "state_lift_geometry_cross", "Graph-minor Kuratowski obstruction", "minor cores / Kpq incidence walls", "minor containment plus exact packet label", "state-lift obstruction analogy", "edge-density-only graph analogy", "Define packet-category K33 containment.", "HYP-2932, HYP-2945", ("residual", "metagraph", "transfer")),
    c("LTM-065", "state_lift_geometry_cross", "Octahedral Hodge current", "divergence, curl, harmonic currents", "local conservation and harmonic residue", "F7 current language", "current not attached to wall/packet", "Express C27/K33 handoffs as current packets.", "octahedral current threads", ("topology", "residual", "transfer")),
    c("LTM-066", "state_lift_geometry_cross", "Unit-distance endpoint-ear recursion", "endpoint masks / deletion ears", "endpoint-universal deletion and spine extension", "cross-domain zipper model", "unit-distance geometry treated as LRC proof", "Name LRC packet spine-flop events.", "HYP-2620, THM-408", ("endpoint", "metagraph", "transfer")),
    c("LTM-067", "state_lift_geometry_cross", "Irreducibility convolution lift", "factor allocations / coefficient lifts", "no hidden convolution factor", "hidden-lift guardrail", "single witness value as global factor proof", "Borrow no-lift language for zero-open packets.", "HYP-2452, HYP-2450", ("residual", "formal", "transfer")),
    c("LTM-068", "state_lift_geometry_cross", "Matroid circuit-tope transfer", "topes, cocircuits, directed cycles", "cycle/cocircuit compatibility", "oriented-matroid wall proof data", "matroid rank shadow", "Build circuit polytope for wall packets.", "HYP-2986, circuit polytope notes", ("topology", "metagraph", "residual")),
    c("LTM-069", "state_lift_geometry_cross", "Covering pressure cap-cell residual", "cap cells, pressure fronts, covering packets", "late covering pressure and moment chart", "F6 covering route", "one chart or scalar pressure only", "Tie cap pressure to Fejer degrees.", "HYP-2969, THM-398", ("topology", "endpoint", "residual")),
    c("LTM-070", "computation_formal_forum", "Exact rational interval engine", "rational endpoints and finite unions", "exact safe mass and components", "zero-open claims", "float topology", "Use exact endpoints for every audit.", "HYP-2948, HYP-2975", ("formal", "topology", "predicate")),
    c("LTM-071", "computation_formal_forum", "Stale-quotient fiber-mixing test", "feature fibers with mixed routes", "same scalar, different LRC route witness", "quotient unsafety proof", "negative result deleted from index", "Add fiber-mixing tables to scalar proposals.", "HYP-2978, LTI-100", ("formal", "predicate", "residual")),
    c("LTM-072", "computation_formal_forum", "Tournament fingerprint payload", "any finite carrier tournament", "score histogram, SCCs, 3-cycles, H paths, tie path", "comparable exploratory output", "fingerprint presented as proof", "Use shared helper in future scripts.", "AGENTS.md, LTT-045", ("formal", "metagraph", "predicate")),
    c("LTM-073", "computation_formal_forum", "Canonical class cache", "small-n tournament representatives", "canonical edge word under relabeling", "isomorphism-class comparison", "raw labelled bitstring", "Cache achievable classes for n<=8 projections.", "A000568, tournament iso work", ("formal", "metagraph", "exact_scale")),
    c("LTM-074", "computation_formal_forum", "Lean local core", "formal arithmetic and interval lemmas", "checked local theorem fragment", "proof kernel", "unformalized key arithmetic step", "Formalize shell collapse or state-lift alternative.", "HYP-2929, LRCApexShell.lean", ("formal", "arithmetic", "predicate")),
    c("LTM-075", "computation_formal_forum", "Forum delta protocol", "POKE posts / index entries", "new technique, guardrail, and pull hook visibility", "agent coordination", "private insight not posted", "Post a delta for every new reusable carrier.", "poke-forum, LTI-105", ("formal", "transfer", "metagraph")),
    c("LTM-076", "computation_formal_forum", "Namespace checkpoint repair", "HYP/T/LTI/LTT/LTM IDs", "unique ID and incoming-signal integration", "shared live research surface", "duplicate IDs silently inherited", "Reserve honest stubs and checkpoint after claims.", "AGENTS.md, CONCURRENT-SESSIONS", ("formal", "transfer", "residual")),
]


def category_vectors() -> dict[str, Counter]:
    vectors: dict[str, Counter] = {cat: Counter() for cat in CATEGORY_ORDER}
    for card in CARDS:
        vectors[card.category].update(card.dims)
    return vectors


def rotated_priority(i: int, j: int) -> list[str]:
    shift = (3 * i + 5 * j + 1) % len(DIMS)
    return DIMS[shift:] + DIMS[:shift]


def orient(vectors: dict[str, Counter], a: str, b: str) -> str:
    ia = CATEGORY_ORDER.index(a)
    ib = CATEGORY_ORDER.index(b)
    for dim in rotated_priority(ia, ib):
        if vectors[a][dim] != vectors[b][dim]:
            return a if vectors[a][dim] > vectors[b][dim] else b
    return a if ia < ib else b


def adjacency(vertices: list[str]) -> dict[str, set[str]]:
    vectors = category_vectors()
    adj = {v: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner = orient(vectors, a, b)
        loser = b if winner == a else a
        adj[winner].add(loser)
    return adj


def score_hist(adj: dict[str, set[str]]) -> dict[int, int]:
    return dict(sorted(Counter(len(v) for v in adj.values()).items()))


def directed_3cycles(adj: dict[str, set[str]]) -> int:
    total = 0
    for a, b, c0 in combinations(adj, 3):
        if b in adj[a] and c0 in adj[b] and a in adj[c0]:
            total += 1
        if c0 in adj[a] and b in adj[c0] and a in adj[b]:
            total += 1
    return total


def scc_sizes(adj: dict[str, set[str]]) -> list[int]:
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    low: dict[str, int] = {}
    sizes: list[int] = []

    def strong(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adj[v]:
            if w not in indices:
                strong(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack.remove(w)
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for v in adj:
        if v not in indices:
            strong(v)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: dict[str, set[str]]) -> int:
    vertices = list(adj)
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    out_masks = [0] * n
    for v in vertices:
        out_masks[idx[v]] = sum(1 << idx[w] for w in adj[v])

    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            avail = out_masks[last] & ~mask
            while avail:
                bit = avail & -avail
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += count
                avail -= bit
    return sum(dp[-1])


def edge_flips_against_tie_path(adj: dict[str, set[str]]) -> int:
    flips = 0
    for i, a in enumerate(CATEGORY_ORDER):
        for b in CATEGORY_ORDER[i + 1 :]:
            if a not in adj or b not in adj:
                continue
            if a not in adj[b]:
                continue
            flips += 1
    return flips


def canonical_edge_word(adj: dict[str, set[str]], vertices: list[str]) -> str:
    n = len(vertices)
    best = None
    for perm in permutations(range(n)):
        bits = []
        for i in range(n):
            for j in range(i + 1, n):
                a = vertices[perm[i]]
                b = vertices[perm[j]]
                bits.append("1" if b in adj[a] else "0")
        word = "".join(bits)
        if best is None or word < best:
            best = word
    assert best is not None
    return best


def emit_markdown() -> str:
    lines = [
        "# LRC Technique Multiverse Index",
        "",
        "**Purpose:** a large pull-from / contribute-to annex for LRC14 technique",
        "reuse.  `00-navigation/LRC-TECHNIQUE-INDEX.md` owns compact `LTI-*`",
        "handles, and `00-navigation/LRC-TOURNAMENT-TECHNIQUE-INDEX.md` owns",
        "`LTT-*` tournament-specific cards.  This annex owns `LTM-*` cards:",
        "creative cross-links that future agents can instantiate, refute, merge,",
        "or promote into the main indexes.",
        "",
        "## Exact Tournament-Analysis Setup",
        "",
        "Vertices are technique families, not runners.  Inside a family, a later",
        "agent may choose runners, gaps, endpoint owners, sections, wall events,",
        "residues, cover arcs, Fourier modes, Haar rectangles, matroid circuits,",
        "Fejer atom banks, state-lift obligations, or proof obligations as the",
        "secondary vertex set.",
        "",
        "Pairwise observable: for two technique families, compare how much of the",
        "following LRC predicate payload each retains:",
        "",
        "```text",
        ", ".join(DIMS),
        "```",
        "",
        "Binary gauge: rotate the priority order by the pair of families, take the",
        "first dimension where their retained-card counts differ, and orient toward",
        "the family retaining more of that dimension.  If all dimensions tie, use the",
        "declared Hamiltonian tie path:",
        "",
        "```text",
        " > ".join(CATEGORY_ORDER),
        "```",
        "",
        "This intentionally prevents one scalar score from deciding every edge.",
        "",
        "## Contribution Protocol",
        "",
        "When adding a card, keep these fields: handle, family, technique, vertex",
        "candidates, observable/gauge, preserved predicate, guardrail, next pull,",
        "and anchors.  If a card fails, do not delete it; add the failure as a",
        "quotient guardrail.",
        "",
        "## Multiverse Cards",
        "",
        "| Handle | Family | Technique | Vertex candidates | Observable / gauge | Preserves | Guardrail | Next pull | Anchors |",
        "|--------|--------|-----------|-------------------|--------------------|-----------|-----------|-----------|---------|",
    ]
    for card in CARDS:
        lines.append(
            f"| {card.handle} | {card.category} | {card.name} | {card.vertices} | "
            f"{card.observable} | {card.preserves} | {card.guardrail} | "
            f"{card.pull} | {card.anchors} |"
        )

    lines.extend(
        [
            "",
            "## Immediate Pull List",
            "",
            "1. Promote any `LTM-*` card that gets a theorem-facing computation into the",
            "   main `LTI-*` index.",
            "2. Add packet examples for cards that currently have analogy evidence only.",
            "3. For every new scalar, run LTM-071 before trusting it.",
            "4. For every new tournament, record LTM-072 fingerprints and the tie path.",
            "5. For every analytic estimate, attach LTM-022, LTM-039, and LTM-055 labels.",
            "6. For every zero-open residual, force a choice among LTM-021, LTM-031,",
            "   LTM-062, LTM-065, and LTM-069.",
            "",
            "## Challenged Assumption",
            "",
            "The failed assumption is that tournament vertices should be runners, or that",
            "a continuous proof route must first be scalarized before it becomes a",
            "tournament.  The useful abstraction is stricter: any finite set of proof",
            "obligations with a pairwise observable becomes a tournament once a gauge is",
            "declared.  Continuous data enters through cutoffs, signs, ranks, endpoint",
            "ownership, exact-period labels, or packet-class changes.  The theorem risk",
            "is always the same: what did the gauge forget?",
        ]
    )
    return "\n".join(lines)


def emit_report() -> str:
    vertices = CATEGORY_ORDER[:]
    adj = adjacency(vertices)
    summit_vertices = vertices[:8]
    summit_adj = {v: {w for w in adj[v] if w in summit_vertices} for v in summit_vertices}
    counts_by_category = Counter(card.category for card in CARDS)
    dim_coverage = Counter()
    for card in CARDS:
        dim_coverage.update(card.dims)

    lines = [
        "S168 LRC tournament-technique multiverse index",
        f"cards={len(CARDS)}",
        f"categories={len(vertices)}",
        "category_counts=" + repr(dict(counts_by_category)),
        "dimension_coverage=" + repr(dict(dim_coverage)),
        "",
        "Tournament Analysis",
        "vertices=technique families, not runners",
        "pairwise_observable=rotating first-difference over retained LRC predicate dimensions",
        "dimensions=" + ", ".join(DIMS),
        "tie_hamiltonian_path=" + " > ".join(vertices),
        "score_hist=" + repr(score_hist(adj)),
        f"directed_3cycles={directed_3cycles(adj)}",
        "scc_sizes=" + repr(scc_sizes(adj)),
        f"edge_flips_against_tie_path={edge_flips_against_tie_path(adj)}",
        f"hamiltonian_path_count={hamiltonian_path_count(adj)}",
        "",
        "Summit-8 isomorphism projection",
        "n=8 vertices=" + ", ".join(summit_vertices),
        "canonical_edge_word=" + canonical_edge_word(summit_adj, summit_vertices),
        "score_hist=" + repr(score_hist(summit_adj)),
        f"directed_3cycles={directed_3cycles(summit_adj)}",
        "scc_sizes=" + repr(scc_sizes(summit_adj)),
        f"hamiltonian_path_count={hamiltonian_path_count(summit_adj)}",
        "",
        "Top next pulls",
        "1. C1 emitted-cocycle matrix on HYP-2963 packet banks.",
        "2. Packet schema with exact M/qdiv, Haar class, Ramanujan projector, Fejer manifest, endpoint owner, state-lift debt.",
        "3. Source-marked tournament spectra for AP, GW, K33, petal, covering, and hard Fejer rows.",
        "4. Fiber-mixing tests for mu, mu^2/phi, divisor, spectrum, and any proposed scalar.",
        "5. Formal F7 residual sector as state lift, harmonic current, or named cohomology class.",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--emit-markdown", action="store_true")
    args = parser.parse_args()

    if args.emit_markdown:
        print(emit_markdown())
    else:
        print(emit_report())
