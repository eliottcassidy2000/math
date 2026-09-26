# Crossroads crossing session: reproduction and validity audit

**PROVED arguments / independently audited / FINITE-EXACT computations,
2026-09-26.** Current conclusions are in the [research board](crossroads_crossing_20260926_board.md).
Collatz and the full 1/3--2/3 conjecture are not resolved here. Previously
attached manuscripts remain external, unverified proof claims; their
terminology is not a substitute for checking their hypotheses and measure.

## Reproduction

From the repository root:

```text
python -X utf8 -B 04-computation/experiments/crossroads_crossing_20260926_audit.py
python agents/check_docs.py
```

The [replay harness](../../04-computation/experiments/crossroads_crossing_20260926_audit.py)
runs all eight lanes normally and with Python `-O`. Checks use explicit
exceptions rather than removable assertions. Both outputs must equal
the retained output. The [manifest](crossroads_crossing_20260926_manifest.json)
pins eight scripts, eight outputs and four inherited dependency files,
normalizing UTF-8 line endings to LF. Hashes establish reproducibility,
not the truth of a proof. Updating a retained output requires inspecting
the discrepancy; `--refresh` refuses to silently overwrite a mismatch.

## Universes and decisive controls

| Lane | Declared scope and controls | Result |
|---|---|---|
| [colour](crossroads_crossing_20260926_colour.md) | All35 supplied symbols; candidate rules and100,000 arithmetic carry cases | First34 near-fit, symbol35 disagreement; equal-length factors refute binary balance; marked-unit repair |
| [potential](crossroads_crossing_20260926_potential.md) | 50,000 flux checks;665,811 height-language checks;152 all-growth realizations, L<=16;100 resets | Exact crossing identity, rational finite-height certificate, and hostile growth/reset blocks |
| [arithmetic](crossroads_crossing_20260926_arithmetic.md) | 792 primitive triples,80 square-hypotenuse targets; r_2 through2000; r_4 through300;372 Mobius and1520 denominator checks | Fibre bijection, signed dynamics, support distinctions and denominator controls |
| [arithmetic run](crossroads_crossing_20260926_arithmetic_run.md) | 100,000 odd starts below200,000;199,999 DFA inputs;1024 fixed-u edges;2499 terminating pairs <=9999 with10,000-step cap | Exact coalescence criterion, eight-state predicate, alternating pair boundary, stopping-time identities |
| [dyadic](crossroads_crossing_20260926_dyadic.md) | All32,766 binary words of lengths1--14; actual odd paths from starts<=300 | Unwrapped affine identity and residue uniqueness; analytic seam, derivative and Fourier corrections separately proved |
| [resource](crossroads_crossing_20260926_resource.md) | Five prime sets, eight forms,25 rational-shadow blocks, up to3620 steps and7620-bit starts | Exact parity, core avoidance, CRT precision, affine endpoint and endpoint valuation equality |
| [independent phase audit](crossroads_crossing_20260926_phase_audit.md) | 72 perturbed languages L<=24;324 all-growth words;16 polynomial shadows, nine forms, four prime sets | Independently implemented controls importing none of the scripts under audit |
| [integration](crossroads_crossing_20260926_integration.md) | Three finite poset cells and one29-step same-band segment | Repairs perfect-balance overclaim and carry-free depth; explicitly retains first-return scope boundary |

These are finite universes, not a search for a counterexample over all
integers. The convergence checks in the run lane verify the stated
stopping-time relation on terminating examples; its general theorem is
conditional on termination. No density is silently promoted to frequency
along one orbit.

## Proof audits and dependencies

[THM-4507](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md)
was reserved and checkpointed before promotion. Its proof received two
independent reviews: rational-word legality, denominator units, polynomial
root avoidance, every precision loss, positivity and finite-core avoidance,
and the accelerated and maximal-rise reset extensions. The selected
negative rational cycle and feature vector stay fixed while positive
starting integers vary with requested block length. This quantifier is
essential. Fibre-dependent boundedness suffices because one fibre is fixed.

The independent review of the board also separated the fixed-count
two-chain parity poset from the proposed excursion forest. The0-or1
source coefficient uses a bijection with binary words; an arbitrary forest
does not inherit it. This is now explicit in the research board.

THM-4476 supplies reciprocal summability only for the hypothetical
nonperiodic positive orbit. THM-4483 supplies the inherited expanding-cycle
obstruction; this session extends its finite-feature scope rather than
claiming the cycle mechanism anew. THM-4503 supplies the source-law warning.
Incoming THM-4506 closes a worst-case landing-multiplicity route, leaving
an orbit-coupled average as a distinct open target.

Analytic literature checks use the primary references linked in the lane
notes: DLMF for Bernoulli Fourier series, theta coefficients and Chebyshev
identities; Romik for Pythagorean dynamics. Kontorovich--Miller's source-
population limits are not the fixed-orbit conditional argument here. The
January2026 near-conjugacy paper was checked for scope, not independently
verified in full. No literature-priority claim is made for the elementary
identities or obstruction theorem.

## Corrections and stopping boundaries

Demonstrated errors were repaired in the incoming crossing note, its
synthesis and navigation routes, and recorded in
[MISTAKES](../../01-canon/MISTAKES.md). The integration note records the
finite balance constants, entropy distinction, logarithmic inequality,
run invariant, carry-sensitive depth, graph-root exception, authoritative
colour mismatch, and modularly invisible draft carry error.

The positive outcomes are exact predicates, identities and obstructions.
The remaining proof attempt requires unbounded information inside a fixed
finite-feature fibre and a justified bound on costs across actual
excursions. Merely choosing a different return map is insufficient:
THM-4507 already covers the maximal-rise reset map.
