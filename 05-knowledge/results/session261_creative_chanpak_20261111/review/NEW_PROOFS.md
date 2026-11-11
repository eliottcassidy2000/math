# Focused independent review of Session261 new proofs

**PASS with one harmless endpoint-lemma boundary clarification.** Review scope is the four named new proof interfaces below. Inherited OCF, Redei, Camion and older constructions were accepted as named inputs. This is proof inspection, not a new numerical census or an audit of the inherited corpus.

## Universal quantitative defect: EXTENSION4.md

The main theorem and exact nonnegative fibre formula are correct under the stated inherited OCF input.

- Inclusion-exclusion on undominated vertices gives the stated signed pair sum. Grouping by X union Y is legitimate: no crossing edge is exactly the condition that X is a union of induced connected components. Each component contributes 1+(-1)^size, giving the positive identity even for empty graphs and graphs with isolated vertices.
- The first and last strong components have the required colours. Camion gives the needed one-sided endpoint paths in distinct components; in a single nonempty strong component, a Hamilton cycle crosses from B to A and cutting there gives both endpoints simultaneously.
- Every admissible nonempty component in the application has order at least two, so no bipartite isolated vertex occurs. Absence of intercomponent graph edges gives every required B-to-A join. Fixed component paths and their distinct block orders really yield c! distinct literal paths; W is recovered from the used vertices. No bijection with dominating sets is assumed.
- The factorial comparison is exact at c=1,2 and strict thereafter. Summation counts every path on its unique vertex set, with no cross-fibre collision. Closing through z is invertible once the marked cycle is specified. Complementary OCF factors include the empty complement with H=1.
- The equality criterion follows term by term from nonnegative brackets. The induced matching lower bound is therefore valid. Empty U and empty A or B cause no problem: both path and nonempty admissible sums vanish.

**Clarification needed:** the standalone endpoint lemma should say “nonempty tournament R.” An empty R satisfies its neighbour conditions vacuously but has no path beginning in A and ending in B. Every use in the theorem is already nonempty, so this does not affect any main conclusion. Also, “Equivalently” before (2) is best changed to “Consequently”: (2) is the advertised consequence, whereas (1) retains the additional intermediate statistic P.

## Packing grammars

The resource matching theorem correctly separates upper bound from simultaneous witness realization. Fresh spine vertices and disjoint core supports remove collisions; arbitrary inter-witness edges do not destroy the selected triangles. The feedback-cover converse works with t witness copies: a core deletion set of size below t leaves some resource support untouched, and fewer than t spine deletions cannot hit all its copies. Core triples and t=0 are handled explicitly.

The forbidden-subsequence bound uses at most two spine vertices per selected triangle. Monotone prefix-recognizer state changes bound exceptional positions; a reachable accepted self-loop alphabet really is freely repeatable and hence subject to the resource theorem. The full-size alphabet classification correctly forces one varying coordinate per triangle block and a single global fixed-bit cut. The clone argument uses disjoint external segments and distinct clones. Finite numerical classifications and receipts were not independently rerun.

## Integral return gates

The increasing lattice closure terminates by the Noetherian property of Z^D and yields exactly the integer span of the orbit, not a reachability certificate. Integer column membership of the stacked defect and constant constraints is the appropriate necessary condition. Rational rows may be cleared of denominators without changing the zero constraints. The distinction between integral feasibility and physical reachability is maintained. Numerical ranks and twelve-pair exclusions were not independently rerun; their certification depends on the supplied exact computations.

## Spectral normal form

The local common flags and diagonal labels give simultaneous triangularization in the tensor incidence algebra, not merely an arbitrary total-order triangularization. Products retain letter-independent diagonal j^L. Both the initial safe annihilator and the sharper argument are valid.

For the sharper bound, the spectral-projector columns have diagonal one and support only on comparable later positions, so they give an invertible incidence-algebra basis change. The transformed block preserves label eigenspaces. Each strict same-label step raises total grade by at least two, and the stated grade ranges give maximum chain lengths k+1 at label zero and k-j+2 at positive label j. These yield the advertised annihilator, eventual polynomial degrees and transient bound, including k=0. “Jordan chain steps” should be understood as maximum block size (number of chain positions); the displayed exponents use the correct convention.

No claim reviewed supplies arbitrary target-count coverage; the OPEN Chan–Pak boundary is preserved throughout.

## Follow-up: hereditary equality and effective finite basis

**PASS** for `bridges/EQUALITY_SKELETON.md` and `bridges/FINITE_EQUALITY_BASIS.md`, by focused proof inspection. The earlier EXTENSION4 wording fixes are reported applied in checkpoint2.

Old W-fibres retain p, admissibility and component count under rooted induced extension; their complementary H increases by the inherited unused-vertex OCF inclusion. All new fibres are nonnegative. This proves the required rooted monotonicity without an unjustified edge-change monotonicity claim.

The order-five formula, wrong-nesting P4 paths, P5 exclusion and connected bipartite chain-graph deduction are correct. In particular the exclusive-neighbour/common-neighbour P5 has no extra cross-shore chords. Twin-class triangles give three literal paths through their common opposite neighbour; the forced strict nesting directions exclude mixed triangles. The induced 3K2 and feedback-pair obstructions follow from individual nonnegative fibres. The supplied isolated-vertex counterexample correctly prevents promoting the necessary skeleton to sufficiency.

For the finite-basis reduction, a fixed edge ab disjoint from a U-triangle has p({a,b})=1 and complement H=3 on the six-vertex rooted witness, producing slack at least 4. Otherwise deleting the chosen edge leaves a transitive spine with eight attachment letters. This reduction is exhaustive for positive slack because empty G has zero slack.

The four domination flags are sufficient: a prior spine A can be dominated by a later spine B, but never conversely across the transitive order; selected core membership must remain fixed alongside the evolving coverage flags, as represented by the four initial core choices. The Hamilton bin assignments uniquely represent paths. The resulting signed finite representation need not be positive internally; nonnegative output from the theorem is exactly what identifies a nonzero avoiding word as a positive witness.

Reachable-row linear-span saturation decides identically-zero series and can retain literal word witnesses. Greedy single-letter minimization gives true subsequence minimality by upward closure. Every later synthesized word avoids every earlier basis word, even after minimization. An infinite synthesis would therefore contradict the supplied minimal-bad-sequence proof of the finite-alphabet subsequence lemma. Final finite rooted containment pruning yields the minimal rooted forbidden types, not merely a word basis tied to one chosen core. No explicit basis-size or order bound is silently inferred.

The numerical automaton receipt and full basis synthesis were not rerun. `constructors/TWO_PORT_KERNEL.md` was absent at this follow-up inspection, so no assessment of that proof is included yet.

## Follow-up: bilinear two-port kernel

**PASS** for the logical all-size extension and arithmetic boundaries of `constructors/BILINEAR_KERNEL.md`, conditional only on its cited previously proved all-chain formula and literal quotient identification. The file became available after the preceding inspection.

Compressed maximal runs give finitely many nonnegative skeleton coefficients: once each run is replaced by its port label, internal ordered spanning path covers independently reconstruct the original path. On a transitive module C_r(TT_a)=r! S(a,r), positive for a>=r and asymptotic to r^a for fixed r. Hence any positive coefficient with r>=3 contradicts the all-chain O(2^a) formula after fixing b>=s; nonnegativity prevents cancellation. The symmetric argument excludes s>=3. The four functions 1, (2^a-2), (2^b-2), and their product are independent, so coefficient identification is valid for all modules, not merely tested sizes.

Expanding x=37C_2+72H confirms 37 times the claimed bilinear law equals xx'+329HH'. The cut injection counts (n+1)H distinct ordered two-slot covers with empty slots allowed: concatenation and cut position recover the original path. Thus x>=(37n+35)H, and the positive-definite norm gives finitely many positive integer candidates (x,h) and then a rigorous finite order bound. The integrality test and actual tournament enumeration retain the physical profile constraint.

A nonnegative row offset forces n<=3; the four small isomorphism types and their stated profiles yield exactly the listed offsets. Thus the unique unit-offset row is genuinely all-order. The residue statement is correct because 37 is 2 modulo 7 and its inverse 4 is itself a quadratic residue. The chain increment formulas and exchanged sidecar agree with the all-chain expression. No numerical coefficient receipt was independently rerun, and no general coverage result is implied.
