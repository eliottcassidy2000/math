# What the older enumeration work contributes now

2026-10-28. This is a targeted recovery across March, April, June and October,
plus the July/August connections routed by the October inheritance audit.
It is not a claim to have read every historical session. Embedded dates and
messages identify provenance; current proofs and reproduced computations
decide mathematical status. Historical timings are not current benchmarks.

The [incoming integration](INCOMING.md) credits the independent terminal fold,
faster degree-bucket counter and representative-only root selector received
before publication. Our isolated fold is not presented as the newest positive
scalar frontier.

## Mechanisms recovered

| era and source | retained mechanism | consequence for this continuation |
|---|---|---|
| March7/8, [a000568_enum.py](../../../04-computation/a000568_enum.py) | enumerate odd partitions; integer common denominator | legitimate exact baseline; no claim that every memoization improves it |
| March8, [GMP v4](../../../04-computation/a000568_gmp_enum_v4.c) | thread partition traversal, bucket by exponent, postpone huge shifts | batch operations sharing a coefficient before doing large arithmetic |
| April4, targeted MSG-522 in the kind-pasteur inbox | historical GMP/CRT comparison reports partition traversal and modular-inverse overhead | fewer bits per arithmetic operation may still mean more total work; no same-host GMP speed claim here |
| April5, [pruning proposal](../../../04-computation/a000568_dp_pruned_s25.py) and [asymptotic series](../../../04-computation/a000568_series_exact_s26.py) | tempting relative-error shortcut, with reversed exponent heuristic | corrected below; exactness needs absolute omitted-mass control |
| June2, [exact Burnside repair](../../../04-computation/a000568_exact_burnside.py) | exact Fraction sum, explicit wrong float values at14/15 | independent current check through40; gives correction lineage |
| June15, [THM-514 / unit-tail kernel](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md) | isolate the soluble1-tail; next prime tail sees divisibility counts | our last-stage polynomial is this inherited kernel in the integral DP's coordinates |
| June15, [THM-516 / phi reframe](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md) | gcd exponent is a totient-weighted divisibility quadratic form | nonlinear arithmetic dependence is retained by future insertion responses, not an ordinary graph cycle analogy |
| October20, [future-response packet](../tournament_future_boundary_20261020/README.md) | merge histories by degree and all remaining gcd responses | the main reduction in the prefix is inherited; its saved Fraction computation independently checks our coefficients through200 |
| October21, [degree-local integer kernel](../degree_local_tournament_enumeration_20261021/session22-enumeration/README.md) | packed responses, O(d)-scaled integral transfers | production baseline reused here; neither packing nor the local gauge is new |
| October25, [adaptive moments](../unprotected_ports_cancellation_20261025/counting/README.md) | only finitely many moments can be asked by a bounded future | contract the terminal operator into P_d(2^k), with a measured gain; the general moment principle is inherited |
| October26/28, [unnamed carrier](../unnamed_carrier_dissection_20261026/README.md), [native Paley](../session72_even_graph_unifying_carriers_20261026/session72-cross-domain/native-paley/README.md), [previous continuation](../frobenius_lca_metagraph_20261028/README.md) | reset/complement quotient, failed pair gluing, phase trees, reset trace | compute the smaller singleton-menu closure; make the pair failure visible in a ten-vertex Hamilton count |

The old March message MSG-069 records a large speed difference between two
historical algorithms. The April MSG-522 records a large advantage for its
then-current GMP implementation. Those are useful leads, not measurements
repeated here. gmpy2 is unavailable in this environment; the present report
compares only the current standard-library kernels on this host.

## Repair: exact integers need absolute control

For n>=3, the partition (3,1^(n-3)) alone contributes

    c_n = 2^[(n-3)(n-2)/2+1] / (3(n-3)!).

At n=4, c_4=4/3. For n>=4, c_(n+1)/c_n=2^(n-2)/(n-2)>1. All Burnside terms
are nonnegative. Therefore the error from keeping only the identity term is
greater than one at every n>=4, including arbitrarily large n. It may be a
tiny *fraction* of the answer while remaining a huge absolute integer error.
This does not rule out a rigorously bounded, growing truncation scheme; it
rules out the historical implication that relative asymptotics provide exact
rounding for free. A bound for one omitted term also does not bound the sum.

The old pruning prose put the exponent extremes backwards: many1-cycles
maximize the identity exponent binomial(n,2); one odd n-cycle has exponent
(n-1)/2. Its prose also mixed 1/z_lambda weights with a second n! division.
Warning headers now preserve these proposals while routing readers to the
repair. The June exact script already diagnosed float precision loss at14/15;
this continuation closes the unwarned historical-source gap, not priority on
that earlier numerical correction.

## What was improved and measured

Three alternating runs per bound in one Python3.13 process; every full
coefficient vector agreed. Times include the whole kernel, exclude JSON I/O,
and are medians. [Raw timing samples](benchmark.json) are kept separately from
mathematical equality receipts.

| family / bound | inherited local DP seconds | terminal polynomial seconds | ratio old/new |
|---|---:|---:|---:|
| ordinary100 | 0.04466 | 0.03786 | 1.180 |
| ordinary200 | 2.09915 | 2.00009 | 1.050 |
| ordinary250 | 11.38186 | 10.94782 | 1.040 |
| pair-free core50 | 0.03390 | 0.01413 | 2.399 |

The two-stage {1,3} contraction is a negative performance control: at100 its
median0.05535s is slower than the0.04466s baseline despite fewer prefix
transitions. It increases terminal work. Its outputs remain exactly correct.

The extension through300 used1,683,206 peak prefix states and44,675,329 prefix
transitions. It took54.21s in a separate, unpaired run, so that time is not
used as a controlled speed ratio. The inherited baseline was also run
through300 and agreed at every coefficient. a(300) has12,887 decimal digits.
This is a verified local extension, not a world record or evidence of a
polynomial-time algorithm. At larger bounds the unchanged prefix dominates;
another terminal tweak cannot remove that bottleneck.

The constructive route yielded a larger improvement. Using the inherited
root-switch block theorem to skip known rows regenerated the exact same
order-nine bank in47.72s, computing191536 rows instead of1716608. A matched
order-eight control changes only whether already assigned rows are skipped;
both modes include parent setup and exclude the old generator's extra
Hamilton/prime audits. Three alternating samples have medians7.4959s versus
1.4112s, a5.31x ratio, with all ordered representatives identical. See
[matched samples](benchmark-roots.json). This is a tested implementation
improvement derived from an inherited theorem, not a new equicount identity.

## Connections that preserve a specific target

| source -> target | actual map / preserved predicate | loss and required sidecar | decisive next test |
|---|---|---|---|
| native atom -> group-frame tree | choose a base; retain relative automorphisms on tree edges; preserves all joint frame orbits | arbitrary ports need compatible stabilizer witnesses; complement changes the atom | P7 marker counts and81/27 transporter, both reproduced |
| history -> scalar enumerator | sum weights by future response, then evaluate terminal polynomials | loses individual realizations and path connectivity | compare every coefficient and keep the slower two-stage control |
| rooted graceful monomials -> signed coefficient tensor | push formal terms through the incidence substitution | cancellations depend on signs and shared variables, not total degree | inherited July28 P4 coefficient-floor hostile; no graceful existence consequence |
| Petersen/K5 contraction -> boundary path tensor | sum over each cluster's traversal possibilities and glue matching ends | one/two runs, cyclic order and global connectivity survive as boundary data | inherited August boundary compiler; test connectedness after contraction rather than just degree |
| square-sum/square-class selection -> linear syndrome | parity of prime valuations determines the square predicate | Fano dependence need not be graphic; a linear syndrome does not retain order or positivity | recover the seven-prime Hamming bank; add a signed endpoint equation before claiming a path construction |
| machine configurations -> bounded continuation classes | identify configurations with identical chosen future outputs up to horizon M | discarded tape or alignment bits can be read later | N4 is an explicit finite transition system; make no unbounded decidability inference |

The July/August records and their corrected limits are linked in the inspected
[October inheritance audit](../tournament_future_boundary_20261020/inheritance.md).
These are mechanisms to transfer, not isomorphisms between the entire problems.
In particular an equality of counts is weaker than an intertwining of allowed
operations. Royle et al.'s [equicount theorem](https://arxiv.org/abs/2204.01947)
concerns **automorphism-sign-even** graphs, not graphs with all vertex degrees
even. It transports the scalar A000568 value but does not supply a port map.

## Proposed next problems — OPEN

1. For the precise singleton/reset/complement menu N4, is the infinite-horizon
   carrier discrete on all X_n? It is now finite-exact through6. Find the least
   separating extension horizon h(n), or a pair surviving every context. The
   inherited strong two-port hinge decoder uses a different product and does
   not automatically settle this literal path-presentation question.
2. Which strong atoms admit small bases whose frame transport remains usable
   across complement partners? Require an explicit translator for path orbits
   and multiplicities. P7's full size-seven closure is a mandatory test.
3. Choose the contraction boundary adaptively using a cost estimate involving
   coefficient bit lengths and number of polynomial groups, not state count
   alone. Can divisibility-profile factorization shrink the middle-stage
   bottleneck without enlarging the final arithmetic? Compare against the
   current GMP kernel before making broader speed claims.
4. Can a terminal observable of the actual black/blue metagraph be contracted directly,
   with certified low tensor rank, to recover the reset trace faster than
   generating the graph? Preserve the81/27 reset masses and blue fixed-point
   flag; do not identify this question with already available scalar Burnside.
   Root-switch block skipping is a positive constructive example for its
   particular operator, not yet a contraction of the full black/blue graph.

The unifying proposal is a compiler with two explicit choices: retain witness
transport when future operations can inspect a realization; contract weighted
responses when only a terminal sum is requested. The proofs above specify when
each choice is sound. No universal compression theorem, Collatz resolution,
graceful-tree proof, or LRC(14) advance follows without the missing hypotheses.
