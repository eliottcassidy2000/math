# Proof scope, independent review, and reproduction

**PROVED IN NOTES:** the arbitrary-module wrapped transfer, protected
observer congruences, all-parameter primitive host family, finite inverse
and sparse-language bound, odd-wrapper divisor/release, and finite physical
identity/separation principles. **FINITE-EXACT:** declared graph searches,
prime certificates, residue-state return and small hostile controls.
**OPEN:** unbounded determinant compensation, arbitrary-target entry,
global marked primitivity after every outer wrapper, and Chan–Pak.
Ordinary Collatz is not solved or identified with these graph moves.

## Sources, priority and boundaries

The starting snapshot is087740676b20. The previous determinant equations
are in target_collatz_transitions/POLARIZED; the NR stabilizer law for
strong parents is inherited from adaptive_metagraph_reformulation/ANCESTRY.
The direct NR proof here adds the explicit nonstrong indegree≥1 boundary
and source-root exception, then applies it to the new exits. The C3 wrapper
transform is inherited THM-3134, and its run expansion is THM-3121.

The A/B neutral pair and its full polynomial are inherited from Session310;
rigidity and external-mark congruence are inherited from Session311. New
work protects a mark in a different leaf, proves its primitive arithmetic
guard, and supplies induced host growth. No claim is made that full Q
preserves selective marks inside the changed leaf. Global asymmetry is
proved separately from marked content.

The cyclic-triangle invariant is inherited from the corrected THM-4803
file cited in history/MARKED_HYPERGRAPH. Its extension to arbitrary
protected-port birth programs is proved by the explicit new-triple update.
Historical reflections supplied pointers; the proof dependencies are the
scoped canon or exact-result files linked in the notes. Neighboring LRC,
six-uniform, planar-Jacobian and port papers are not treated as theorems
about tournament target coverage.

Incoming work through7888efd8f9 is separated in INCOMING.md. No scarce
theorem ID was reserved. There is no new canonical claim of a solved prize.

## Exact universes and controls

| Artifact | Universe / positive result | Hostile or boundary retained |
|---|---|---|
| pair_growth.py | All256 single-vertex extensions of each specified order8 module;242/239 profiles; zero compensated pairs | Only one birth on each side; no arbitrary-growth impossibility |
| pair_growth2.cpp | All131072 two-vertex extensions per side;56614/58244 profiles; all response-pair solutions via exact modular inverse | No strongness filter or sampled attachment bank; counts response pairs, not all labelled pair multiplicities |
| pair_check.py | All12 resulting module profiles/strongness; two independent10! permutation counts; six26-to28-vertex formulas; four trial-division prime certificates; named24-to28 induced embeddings | Large output H uses universal theorem, not28-vertex DP;285 complete small NR automorphism controls include indegree0 |
| boundary_check.py | Exact W(A),W(A²),W(A^191) values; context3→content1 release; all extension masks actually nonconstant | Third slot removes an apparent457 divisor; old induced subset is not still a uniform module |
| ports/port_return.py | Complete cyclic-itinerary convolution from≤8-vertex profiles; literal small controls; finite elementary prime/order/character certificate; all-parameter gcd proof | Arithmetic guard sharp; global C3 and point stabilizer differ; primitive right-module phase does not imply primitive outerΠ |
| clock/refill.py | Literal named two-flip cycles, old induced embeddings, finite run signatures, prime/rigid7-site quotient, actual target inverse | Positive wraps and mod3 trap; exact inverse includes its complete size bound; no scalar-address injectivity claim |
| clock/entry_probe.py | Every A95-seed address for six determinant outputs, allowing every initial left size;380 states; zero entries | Not an arbitrary-seed exclusion or a global nonattainment result |
| clock/phase_certificate.py | Exact finite vector return at indices2 and338 modulo3,5,49 | Proves all future period336 for this orbit; word(0,4+336t,0) has exact v7=1; no prime-cofactor assertion |
| history/hypergraph_probe.py | All labelled tournaments of orders2–5 for fixed-star and indegree-one fibres;1098 birth identities;260 short multiport nodes | Same-H/star H9 rearm changes a cyclic triple; disjoint neutral-flip hostile first occurs at6 in the stated edge-pair universe |
| history/transitive_probes.py | Stirling physical basis through cap6; tensor reconstruction; six Π basis probes; sharp simultaneous9/3/5 control | Needs certified run caps; does not apply to unrestricted cyclic variable modules or interior marks |

The C++ counts fit signed64-bit arithmetic: at order10, H≤10! and
P2≤11·10!, while every product used in the modular-inverse search remains
below2^63. The Python witness verifier uses arbitrary-precision integers.
The C++ exhaustive search was run once after compilation with−O3, and its
search completeness and bounds were independently reviewed. It is not
described as an independent rerun of the whole universe.

## Independent reviews

The ports proof was reviewed independently by the history agent through
its exact modulus certificate and by the clock agent for the six-round
reservoir phases. The clock construction, named arcs and rigidity were
reviewed by the history agent; the root checked address recovery, the
finite inverse, weighted-word bound and all residue guards. The history
congruence and simultaneous physical separator were reviewed by the clock
agent and the root. The root determinant universe, modular pairing,
literal embedding and NR boundary were reviewed by the ports agent.

The initial history replay scripts used Python assertions. Before this
packet was integrated those21 assertions were replaced by explicit
require/RuntimeError checks. Normal and optimized modes both reject a
deliberate false guard. This is checker hardening, not a change to a
mathematical result. Other Python guards are likewise explicit.

## Reproduction

From the repository root, first reproduce the complete two-birth search:

```bash
c++ -O3 -std=c++17 05-knowledge/results/symmetry_port_return_20261113/pair_growth2.cpp -o 05-knowledge/results/symmetry_port_return_20261113/pair_growth2.bin
05-knowledge/results/symmetry_port_return_20261113/pair_growth2.bin 05-knowledge/results/symmetry_port_return_20261113/pair_growth2.json
```

Then run the nine Python scripts in dependency order:

```bash
python3 -B 05-knowledge/results/symmetry_port_return_20261113/pair_growth.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/pair_check.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/boundary_check.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/ports/port_return.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/refill.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/entry_probe.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/phase_certificate.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/history/hypergraph_probe.py
python3 -B 05-knowledge/results/symmetry_port_return_20261113/history/transitive_probes.py
```

The compiled executable is intentionally excluded from the committed
packet. Receipts and source hashes are retained. VALIDATION.json records
normal/optimized receipt agreement and repository checks; SHA256SUMS covers
the packet, including the separate history manifest. Relative source links
are checked against Git objects as well as working files because the
session worktree is sparse.
