# Audit and reproducibility boundary

2026-10-25. All claims below have explicit domains. The main proof files
and independent reviews were frozen before the root's paired replay.
The reproduction scripts use explicit exceptions; Python -O does not
disable their gates. The manifest normalizes CRLF to LF and pins both
the current packet and its immediate proof/code inputs. It is a source
integrity check, not a substitute for mathematical review.

## Proof and control ledger

| Claim | Independent audit | Exact finite universe / hostile |
|---|---|---|
| h=3, 1->3, every n>=12 | [h3 review](h3/INDEPENDENT_REVIEW.md) rederives the affine degree cocycle and induced-orbit degree sum | Producer91 full affine cells n12..100 plus120,160; formulas through500; independent Gaussian h3/h4 check n12..40; genuine blue-word exceptions at n8 and n10 |
| h=4, 3->1, every n>=12 | [h4 review](h4/INDEPENDENT_REVIEW.md) checks each of three unequal arc witnesses | Producer109 full affine cells n12..120 and327 literal models; formulas through1000; extra port5 makes the new first join illegal |
| Closed endpoint-four, n>=15 | [Inheritance audit](PORT_FRONTIER.md) checks the unfiltered partition into base cases and the residual h/chord branches | Inherited24 complete affine systems at15..26 replayed with the independent Gaussian solver; all5040 permutations verify the n7 target has9 blue paths and no black path |
| Integral adaptive gauge | [Counting review](counting/INDEPENDENT_REVIEW.md) independently derives exact odd division, J>=0, b'<=b+k and destination-gauge merging | Ordinary/core0..26,2642 intermediate states; unmerged core oracle0..22; full core44 with719 intermediate states;374 transfer squares,430 negative ungauged channels; independent105-cycle ceiling and mask-deletion hostile |
| Pair-free and interpolation parity | [Parity review](parity/INDEPENDENT_REVIEW.md) checks canonical weight orbits separately from cut orbits | Eight odd actions and an even hostile; full tournament classes0..6; signed interpolation0..20; independent weighted/cut orbit agreement0..5,53 weight-parity checks through total order9; all C5-equivariant maps to two points |
| Fixed binary precision and low-bit formulas | [Precision review](precision/INDEPENDENT_REVIEW.md) checks quadratic cycle bounds, equality cases, multiplicity patterns and modular units | Producer31583 odd partitions0..50,306 modular comparisons K1..6; independent integer-response counts0..60,183 comparisons;1326 sharpness cases,950 triples,609 four-cycle cases, totient controls through500, budgetsK1..128 |

Producer h3/h4 checkers reuse the inherited affine-source constructor.
The independent h3 review reconstructs its Gaussian equations without
importing those producers and covers both new ledgers. The port-frontier
control deliberately reuses the *independent* predecessor boundary
solver and separately checks the order-seven hostile literally.

The integral counting checker compares its states with a rational
f-response computation and its totals with unmerged arithmetic oracles.
The precision review uses that integer response route as a separate
count oracle from both full rational partition summation and modular
bounded-cycle evaluation. The parity review's literal weighted quotient
route does not import the producer. These independent paths are stated
at their actual scope; no finite sample supplies an infinite theorem.

## Types, quotients and retained information

The port result concerns a *marked path* in a literal labelled source.
Blackness is tested at each actual intermediate tournament. Every word
is a permutation, every join is directed, and the path-edge XOR has the
declared literal target. Clone configurations are not contracted away.
No assertion equates a local repair theorem with global connectivity.

The adaptive state is an arithmetic aggregate of decorated permutation
cycles. It is not a tournament, a port carrier or an orbit of a metagraph.
Its moment quotient retains every allowed future aggregate at a
decreasing horizon. The new product lattice contains the actual moment
lattice and can be larger. The proof is stronger than integrality on
reachable histories: all its transfer matrices are integral on that
product lattice. It does not establish positivity there.

Canonical clone inflation uses full automorphism orbits of positive
integer weights on a pair-free quotient. Dropping these orbit sizes is
valid only for the stated parity observable; the first integer failure
is explicitly checked at order4. Constant subset toggles descend exactly
when their supports are invariant. A modular depth of n therefore does
not add n operations. The positive correction B still fails the C5
equivariant quotient test.

The fixed-precision filter concerns rational Burnside summands in the
2-adic integers. Its denominator is split before reduction; only an odd
unit is inverted. Discarded terms individually vanish at the declared
precision, so no cancellation heuristic is needed. Order2 is excluded
from the even baseline and handled explicitly. Increasing precision
requires restoring terms; the full integer complexity is not bounded
by the fixed-K statement.

## Status and inheritance

The all-odd-action law was supplied by the user and was independently
proved in the incoming session63 packet. It is credited as recovered
input. The April odd-order valuation theorem is also recovered input.
The current proofs add the explicitly stated consequences and bounds;
bounded searches do not establish historical priority.

The n7 pure-blue obstruction remains **REFUTED** evidence against the
old unfiltered all-order detour claim. The present n>=15 theorem has a
different quantified domain. Missing lower orders, other endpoint
layers, open-path transport, a positive state recurrence, and a smaller
unrestricted structural A000568 recurrence remain **OPEN**. No famous
open problem is declared solved by analogy.

Run replay.py to reproduce all saved paired controls. Run verify.py
read-only to check source pins, paired PASS receipts, Python syntax and
local links. Only explicit --refresh writes either receipts or pins.
