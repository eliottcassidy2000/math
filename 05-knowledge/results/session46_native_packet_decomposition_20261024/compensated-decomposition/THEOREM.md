# Genuine charged compensation lowers native atomic support

Status: FINITE-EXACT native paths, complete pinned domains, unit-capacity resources, direct image witnesses, and optimal bottleneck bounds. The positive routes permit later reoccupation of a released resource; they do NOT satisfy strict temporal nonreturn. This packet concerns compatible partial observer packets. It does not rewire the fully occupied canonical atlas or construct its unchanged remainder.

## 1. Main result and the resource/coordination trade-off

The Session 45 sixteen-source circuit remains primitive in its original complete pinned family: the only full-dyadic accepted states are its two endpoints, so its atomic changed-slot budget is sixteen. A genuine additional source option changes this answer.

| Circuit | Added charged source owners | Selected slots | Maximum changed slots per move | Maximum changed owner clocks | Number of moves | Total changed owner clocks |
|---|---|---:|---:|---:|---:|---:|
| Original 16 | none | 16 | 16 | 32 | 1 | 32 |
| Original 16 | 1394 | 17 | 8 | 16 | 3 | 44 |
| Original 16 | 3015 | 17 | 6 | 13 | 4 | 50 |
| Original 16 | 1394, 2905 | 18 | 6 | 12 | 4 | 48 |
| New 18, original cost 43 | 1394 | 19 | 6 | 17 | 4 | 59 |
| New 18, original cost 37 | 2071 | 19 | 6 | 13 | 4 | 49 |

The table displays the initially certified routes, not necessarily the least-work or least-peak-clock routes. Exact event/work optima and improved attaining paths are proved in [variants/THEOREM.md](variants/THEOREM.md): the single-1394 route improves from44 to40 total clocks, the cost43 extension from59 to55, and the cost37 extension trades four events/49 clocks/peak13 against five events/57 clocks/peak12.

Every added source returns to its initial literal choice. All moves retain the full original nodewise dyadic harmonic bank and the selected ordinary per-slot X phase modulo 8. The last row has its own packet and image file; it is included by the `--include-eighteen37` verification option. Exact receipt values, rather than this summary, are authoritative for all clock totals.

The maximum changed-slot budget six is optimal for all three endpoint pairs even if every one of the supplied twenty-eight sources is available as a charged compensator, with arbitrary initial buffer choices and arbitrarily many moves. An integer invariant of every possible four-slot move proves this stronger lower bound.

For the original circuit, one additional source is the least possible reservoir size for any reduction below sixteen. In the twenty-source cost-two supply, one buffer has sharp budget eight and two buffers are necessary and sufficient for budget six. For the original sixteen-source circuit, across all twelve possible single-buffer owners in the twenty-eight-source supply, the minimum achievable maximum owner-clock cost is thirteen. Two cost-two buffers achieve twelve, which is optimal even with all twenty-eight slots available: some move must change at least six slots and each changed permutation changes at least two clocks. Thus spare charged capacity lowers atomic coordination, but generally increases total routing churn.

## 2. Native carrier and exact binary state convention

The targets are U=1965 and V=2234. A bit zero chooses the actual `before` path to U; a bit one chooses the actual `after` path to V. The original sixteen source order is

    1425,1457,1477,1480,1638,1691,1974,2837,
    176,1549,2071,2176,2185,2581,2945,3000.

Its initial mask is 65280 and its final mask is 255. Buffer bits are appended. For the main certificates they initially equal zero and return to zero. Source 3015 has owner-clock cost three. Sources 1394 and 2905 have cost two.

All twenty-eight canonical source matrices are reconstructed from the original atlas. Their source automorphism groups are trivial; hence these particular source kernels have no nontrivial hidden lift choices. This is not a claim about arbitrary native sources. Literal matrices, physical owner bijections, paths, and all four actual complement/converse representatives are retained in `packet.json`.

Only the first owner, last owner, and both owners of one indexed internal arc are shared between slots. The remaining owner symbols are explicitly source-private. Every unchanged indexed owner clock is pinned. The matrices in different source frames need not agree on undeclared induced arcs, and no common-host tournament amalgam is claimed.

Completing all common fixed clocks gives fifty-seven actual Hamiltonian paths. Fifty-six belong to the U/V pair, precisely two for each supplied source. The only additional path belongs to source 1923 and has target 1559. Since the initial margin at 1559 is zero, no margin-preserving state can select it. Thus the twenty-eight binary U/V domains are complete after all physical target margins are imposed, not artificially truncated menus.

All fifty-seven normalized F/R resource orbits are distinct. Each has four literal closed ports; their union has 228 ports. Each supplied route uses at most one option of each source at a time. An option may be released and later restored; that is temporal routing, not two simultaneous uses of one resource. The union of alternatives needed by the seventeen-slot routes has 34 resources/136 closed ports, with 17 resources/68 ports per state. The eighteen-slot route has 36/144 across its union and 18/72 per state. Each nineteen-slot extension has 38/152 across its union and 19/76 per state.

## 3. The exact row-current law

The effective paired-charge group is

    Q = (Z/2)^4 plus Z/128.

Its authentic source labels are four Y parities together with X modulo 128. The common targets have equal Q label, and X_U=X_V=146. For a state b write

    c(b) = (sum_s b_s, sum_s b_s ell(s)) in Z plus Q.

Two supplied states have equal target margins and the same owned full-dyadic current signature exactly when their c-values agree. The first coordinate is an integer, not a modular margin.

For a move b to b', put u=sum_s (b'_s-b_s)e_s and v=e_U-e_V. Its actual signed edge change removes 2(b'_s-b_s) from sU and adds it to sV. Because sum u_s=0, its full physical-row signed Laplacian is

    Delta L = 2(u v^T + v u^T).

No source/target copies are split: all contributions are accumulated in the original physical row. Here the source set is disjoint from the two targets, simplifying but not weakening the row test. Each source row is a signed copy of 2v; the target rows are signed copies of 2u. The certificates retain direct grounded integer vectors satisfying

    L z = 2u + 256 w,
    L z_target = 2v + 256 w_target,

including every coordinate of the correction vectors w. Independent multiplication over the full original edge export checks every equation. These equations annihilate every original harmonic character modulo 256, not merely a guessed basis. The separately authenticated original exponent-256 saturation extends the result to every genuine dyadic level. Odd-primary characters and all ordinary per-port characters are outside the claim.

The chosen integer X often has literal target debt ±256 or ±512 during a move. This is not a failure of the genuine higher-dyadic bank: the appropriate higher harmonic lifts scale. Evaluating one unchanged integer X at indefinitely increasing moduli is a different test.

The checker also reconstructs every signed edge change, verifies the complete displayed row identity, and directly evaluates every physical row against authentic generators. The five generators organize the exhaustive state search; the full integer image witnesses independently certify each performed move.

## 4. Explicit four-move routes for the sixteen-source circuit

For the single buffer 3015, the complete guarded fiber has exactly six states:

    255, 65280, 85840, 92372, 100534, 114918.

One shortest six-slot route is

    65280 -> 85840 -> 92372 -> 100534 -> 255.

The positive and negative source triples at its successive steps are:

1. +(1638,1974,3015), -(2185,2581,3000).
2. +(1477,2581,2837), -(176,1549,2071).
3. +(1457,1691,3000), -(1974,2581,2945).
4. +(1425,1480,1974), -(2176,3000,3015).

The buffer changes at steps one and four and is restored. At the intermediate states, projecting away the buffer does not yield a legal state of the original circuit: the original sixteen slots have both margin and charge debt, compensated by the extra source. This is precisely why the old two-state projection obstruction does not apply. Neutral resource workspace cannot produce these intermediate signatures.

For the two cost-two buffers 1394 and 2905, the complete guarded fiber has twenty-five states. A four-step route is

    65280 -> 260104 -> 210993 -> 86071 -> 255.

Its triples are:

1. +(1394,1480,2905), -(176,1549,2071).
2. +(1425,1638,1691), -(1480,2945,3000).
3. +(1457,1477,2945), -(2176,2581,2905).
4. +(1480,1974,2837), -(1394,2185,2945).

Every step has six changed slots and twelve changed owner clocks. Both buffers are restored. The single cost-two buffer 1394 instead has a twelve-state fiber and route

    65280 -> 73058 -> 66471 -> 255,

with changed-slot sizes 8,8,6. Exhausting each of the four outside cost-two owners 1394,1921,2905,3369, with both possible initial choices, proves that every such single-buffer menu has exact bottleneck eight. This establishes the stated minimum two-buffer reservoir for six-slot routing within the cost-two supply.

## 5. A complete small-move obstruction, even with all charged supply

Distinct source labels rule out a nontrivial two-slot move. Fixed target margins force every changed-slot count to be even. Any four-slot move must therefore replace two selected sources by two unselected sources with the same Q sum.

The checker enumerates all 378 unordered source pairs and all their equal-Q collisions. There are exactly twenty-seven possible four-slot trades on the complete supplied twenty-eight-owner universe. Every one preserves the INTEGER statistic

    I(b) = b_1638 + b_1974.

The only four-slot trades involving these owners are

    (1638,1841) <-> (1974,2176),
    (1638,2176) <-> (1974,2945).

Both preserve I. The original sixteen-source endpoints have I=0 and I=2. The original cost-43 eighteen-source endpoints change I from zero to one after unchanged outside coordinates are canceled. The cost-37 eighteen-source endpoints change I by minus one. Restored buffers cannot change any of these endpoint differences. Consequently no sequence of moves changing at most four slots can connect any of these endpoint pairs, even with all twenty-eight slots flexible. Odd support sizes do not offer an escape, so the lower bound is six.

This invariant is stronger and much shorter than checking finitely many routes. An optional independent reachability search over all 4096 initial settings of all twelve original-circuit buffers finds no four-slot route, but that computation is not needed for the proof.

The result concerns these endpoint pairs, not every fiber of the full supply. Equality of global charge invariants and generation of a difference lattice do not in themselves imply connectivity under binary capacity constraints.

## 6. The two eighteen-source extensions

The inherited primitive cost-43 support has positive sources

    1549,1691,1974,2087,2185,2553,2837,2945,3000

and negative sources

    1457,1477,1841,1921,2176,2581,2835,3015,3360.

Append buffer 1394 at bit eighteen. The complete nineteen-slot fiber has nineteen states. Its certified route is

    261632 -> 437825 -> 274217 -> 272671 -> 511.

All four moves change six slots; their owner-clock costs are 13,17,14,15. The zero-buffer domain has only its two original endpoints, so one buffer is necessary and sufficient.

The distinct cost-37 primitive support has positive sources

    1394,1425,1457,1480,2176,2185,2581,2835,2945

and negative sources

    1477,1549,1691,1921,1974,2837,2905,3000,3369.

Append buffer 2071 at bit eighteen. Its complete guarded fiber has thirty-eight states, with route

    261632 -> 351784 -> 479592 -> 332150 -> 511.

Again each move changes six slots and the buffer returns to zero. The route has its own literal source order and four direct image witnesses. The common critical-pair invariant proves global six-slot optimality for this endpoint pair as well.

## 7. Two different meanings of “without reusing ports”

The positive statements use standard unit-capacity injection: no state double-occupies a normalized resource, no transition counts a removed resource twice, and no hidden lift supplies a second copy. An option may be released and later reacquired. The buffer restoration and some repeated original-slot flips genuinely use this permission. The certificate records those events; calling the resource a buffer does not erase their cost.

Under the stricter rule “once a normalized resource is released, it may never be occupied again during this route,” NONE of the positive compensated routes qualifies. In fact a short conditional obstruction proves that no restored-buffer decomposition exists in these exact complete pinned menus:

1. A binary auxiliary slot that changes and later returns must reacquire its initial resource, violating strict temporal nonreturn. Hence every restored auxiliary slot is constant.
2. A binary original slot can change at most once. Any second change reacquires its initial resource. Thus the original changes are conformal: every intermediate state changes only a subset of the endpoint changes.
3. The primitive original sixteen- and eighteen-slot domains have only their two guarded states. No nonempty proper subset of changes is legal. Therefore the first nontrivial move must change all sixteen or eighteen original slots at once.

This proof retains all supplied buffers and is independent of their initial orientations. It relies on the exact pinned binary domains and fixed target margins; it makes no claim after releasing pins or supplying genuinely new path alternatives. So the result is a precise dichotomy, not an unqualified resolution of a request forbidding all temporal port reuse: charged compensation lowers simultaneous atomicity only under returnable unit-capacity resources. Strict temporal nonreturn preserves the primitive atom.

## 8. Availability and the unchanged remainder are real hypotheses

A positive state or route is not a proof of a global completion. To install one of these partial-packet routes, its initial selected resources must be available, and a fixed unchanged remainder must avoid the entire union of resources used anywhere along the route. That same remainder must satisfy the declared shared-owner compatibility constraints and supply any current required to balance the enlarged initial packet. Each intermediate packet preserves that initial packet's owned current signature, so an already compatible remainder remains compatible throughout.

No such complete remainder is constructed here. In particular, the full canonical all-path atlas already occupies the alternative resources; one cannot merely add these compensators to that atlas and call the result a native rewiring. The claims are exact selected-packet alternatives and conditional exchanges against an explicitly compatible, disjoint unchanged remainder.

This is the all-size conceptual rule supported by the examples: authenticated charged resources enlarge the actual fixed-margin charge fiber; legal small moves must connect its chosen endpoints while returning added coordinates. Resource availability, exact native domains, and fiber connectivity are separate obligations. Counting labels, introducing neutral storage, or finding a lattice identity cannot discharge them automatically.
