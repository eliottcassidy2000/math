# Exact native receiver thresholds, including strict resource non-reacquisition

Status: **FINITE-EXACT, original order-eight native packet, frozen exterior.**

## 1. Scope and answer

Use exactly the nineteen-slot receiver from the frozen native-repair packet:

- Sixteen independent cheap sources
  `1394,1477,1480,1549,1638,1691,1921,1974,2176,2185,2581,2837,2905,2945,3000,3369`.
- The actual target-class equality `1457 = 1841 = 3015`, so this triple is either its three original U records or its three original V records.
- Every other core and helper full record is frozen at an accepted exterior assignment. Its population and charge constants are subtracted, not discarded.

The full native inventory remains 82 occupied units / 328 inverse-converse-closed words, within 164 reserved units / 656 words. All original paths, indexed/block pins, kernel lifts, source/K/converse partners, owner-4 shared-tree restrictions, original dyadic guards and capacity-one resource reservations remain in force. The 652-owner source-frame host is not a new global Hamiltonian observer.

**Theorem.** Across every pair of accepted endpoints in every exact receiver `(grade, Q-charge)` fibre:

1. The sharp uniform ordinary connectivity budget is **16 changed slots and 34 literal owner-clock occurrences per atomic event**.
2. The sharp uniform strict no-resource-reacquisition budget is **the same simultaneous 16 slots and 34 clocks per event**.
3. Both assertions already have their sharp lower witness in grades 8..11, and hold for all grades 0..19.

The bounds are event-local bottlenecks, not total route cost. Zero-event routes suffice for identical endpoints. A single accepted atomic replacement releases changed alternatives and acquires their authenticated partners simultaneously; a different microscopic transition contract would require a separate analysis.

This does **not** establish a 34-clock theorem when exterior helpers/outside records are free, when additional relational restrictions are imposed, for unrestricted helper-assisted connectivity, or for all-size native wiring. In particular, adding arbitrary pins can delete the intermediate records used here. The abstract interval criterion below still applies to the new accepted domain, but its primitive moves must be recomputed.

## 2. Complete classification, not branch-boundary minima

There are exactly `2^17 = 131072` actual assignments, `23444` nonempty fibres, and `474896` unordered distinct endpoint pairs. The largest fibre has 22 records. The grades 8..11 contain 334504 endpoint pairs. The all-charge complete grades remain exactly 8,9,10,11.

The exact uniform thresholds by grade are:

| Receiver grade | Ordinary slots/clocks | Strict slots/clocks |
|---|---:|---:|
| 0,1,18,19 | 0/0 | 0/0 |
| 2,17 | 4/8 | 4/8 |
| 3,16 | 6/12 | 6/12 |
| 4,15 | 8/18 | 8/18 |
| 5,14 | 10/22 | 10/22 |
| 6,13 | 12/26 | 12/26 |
| 7,12 | 14/30 | 14/30 |
| 8,11 | 16/34 | 16/34 |
| 9,10 | 14/30 | 16/34 |

Each displayed pair is a simultaneously attainable bound on every event of the same route. It is not formed by optimizing the two coordinates on different paths. The checker computes the **separate exact ordinary and strict optimum for every endpoint pair**, not just the uniform maximum or minimum cross-branch distance. Their deterministic streaming digest is frozen; the certificate retains all primitive trades, complete lower-witness fibres and exact per-fibre summaries without a bulky all-pairs payload.

Encode the sixteen independent bits in increasing source order and the triple bit as bit 16. For a changed mask let `b` be the changed independent-bit count and `t` the changed triple bit. Its actual event costs are

    slots = b + 3t,    clocks = 2b + 8t = 2 slots + 2t.

Equal grade forces the slot count to be even. Thus all native event vectors are totally ordered coordinatewise: at even support `s` the only possible costs are `2s` and `2s+2`, both below any cost at support `s+2`. Minimizing the clock bottleneck therefore returns the unique nondominated joint slot/clock bottleneck. This special cost-chain argument is explicit; it is not assumed for arbitrary weighted receivers.

## 3. Strict routes and the primitive interval criterion

Let `D` be any finite accepted domain of binary full-record choices with disjoint reservations. For endpoints `x,y`, define their accepted interval

    I_D(x,y) = {z in D : every coordinate on which x and y agree is fixed at that value}.

A strict no-reacquisition route is exactly a path whose changed coordinates move once, toward their endpoint values. Necessity: changing a binary coordinate twice returns to its initially occupied full record and reacquires its reserved resource. A coordinate whose endpoints agree cannot change at all. Sufficiency: each changed coordinate acquires only its other authenticated full record once, and disjoint source reservations prevent acquiring another slot's previously occupied resource.

Call a distinct endpoint pair **primitive** if `I_D(x,y) = {x,y}`. This is the binary full-record analogue of conformal/Graver primitive moves, not a new general algebraic principle.

**Primitive uniform criterion.** For any rectangular per-event budget `(S,C)`, every ordered endpoint pair in `D` admits a strict route if and only if every primitive pair itself has direct event cost at most `(S,C)`.

Proof. A strict route between primitive endpoints has no permitted interior state, so must use their direct event. Conversely, split any nonprimitive pair at a proper accepted interval state `z`. Both subpairs have fewer changed coordinates, and their changed-coordinate sets partition the original set. Inductively route `x` to `z` and `z` to `y`; concatenation changes each endpoint-different coordinate exactly once and never changes another coordinate. This terminates in primitive direct events, all within the same joint budget. The argument uses all accepted full records, so no projected compatibility is silently supplied. QED.

Consequently the strict uniform rectangular threshold is the pair of coordinate maxima over primitive moves. Those maxima are simultaneous because the same primitive-event collection supplies the routes. In a general weighted domain, a particular endpoint pair may have a smaller optimum and may have a nontrivial Pareto frontier; the uniform criterion does not identify its optimum. Here the cost chain and full interval dynamic programming determine every individual joint optimum separately.

After removing coordinates at which both endpoints are V and identifying reversal, there are exactly **2543 distinct primitive balanced trades**, occurring in **463840 primitive endpoint pairs**. Their complete cost histogram is:

| Slots | Clocks | Distinct primitive trades |
|---:|---:|---:|
| 4 | 8 | 2 |
| 6 | 12 | 48 |
| 8 | 16 | 240 |
| 8 | 18 | 10 |
| 10 | 20 | 516 |
| 10 | 22 | 97 |
| 12 | 24 | 313 |
| 12 | 26 | 453 |
| 14 | 28 | 36 |
| 14 | 30 | 632 |
| 16 | 34 | 196 |

All 196 maximal moves change the triple and thirteen independent slots. There is no primitive move above 16/34. In particular every 18-slot endpoint difference, and every 16-slot triple-preserving difference, admits a proper accepted interval split. The table and criterion give a structural explanation of the strict upper bound; exact enumeration supplies the native finite fact that no larger primitive survives.

Stripping common V coordinates does not change a move's primitivity: all intermediate choices keep those coordinates fixed, and adding/removing their common grade and charge constants is a bijection of the relevant intervals. All primitive trades and their literal joint costs are listed in `certificate.json`.

### Exact primitive-context generating functions

Let a normalized primitive trade have actual support `s`, independent-bit support `b`, and triple-change flag `t`. Its two disjoint sides each have grade `s/2`. The unused independent bits can be common U or common V freely; if the triple is unchanged, its common U/V value is also free. Thus its exact endpoint-pair completion polynomial by grade is

    z^(s/2) (1+z)^(16-b) (1+z^3)^(1-t).

This is an identity for **primitive endpoint pairs**, not for all endpoint pairs. Every context gives exactly one unordered endpoint pair; removing the common V bits recovers the normalized primitive uniquely. Summing the formula over all 2543 trades reproduces every grade and every joint cost bin of the 463840 primitive endpoint pairs. The checker compares the formula with the independent full-fibre pair enumeration.

There is an additional support-additivity explanation at the critical cost. Equal-grade endpoints in nineteen slots differ on at most eighteen slots. Any nonprimitive pair recursively splits into at least two primitive events with disjoint changed supports. Every primitive has at least four changed slots, so each event in such a decomposition has at most `18-4=14` slots and hence at most 30 clocks. Therefore a strict endpoint optimum of 34 occurs **if and only if the pair itself is a 16-slot/34-clock primitive**. There are no extra decomposable critical pairs.

The 196 maximal normalized trades flip the triple and thirteen independent bits, leaving exactly three independent common bits. Consequently the complete strict-critical endpoint-pair polynomial is

    196 z^8 (1+z)^3.

The exact counts are 196,588,588,196 in grades 8,9,10,11 respectively, totaling 1568. By contrast only two ordinary endpoint pairs in the whole receiver require 34 clocks: the grade-8/charge-53 two-state fibre and its grade-11 complement. All these counts are checked separately from the structural argument.

## 4. Sharp lower controls and ordinary/strict separation

### 4.1 The inherited two-state control

At grade 8 and residual charge 53, the whole receiver fibre consists of exactly these two V-source sets:

    {1477,1638,1691,1921,1974,2176,2905,3000}
    {1394,1457,1549,1841,2185,2581,2837,3015}

They are disjoint eight-subsets. Their only connecting event changes 16 slots and costs 34 clocks. Hence no smaller slot budget, and no smaller clock budget, can work uniformly even ordinarily. The direct replacement is strict. This is a frozen-exterior lower bound only.

### 4.2 Sharp ordinary grade-9 threshold

At grade 9, charge 651, the complete fibre in compact mask notation is

    38386, 88077, 92768.

The three pair costs are respectively

    38386 <-> 88077: 16/34,
    38386 <-> 92768: 14/30,
    88077 <-> 92768: 8/16.

Below 14 slots or below 30 clocks, mask 38386 is isolated. At 14/30 the two smaller edges connect the fibre. Full enumeration proves every other grade-9 fibre also connects at 14/30. The complementary grade-10 result follows from the literal isometry and is independently counted.

### 4.3 A much cheaper ordinary route that must reacquire resources

At grade 9, charge 3, take endpoint masks `31171` and `75320`. The complete fibre has eleven states:

    7907,20452,27198,31171,39819,56652,61224,66460,70748,75320,99971.

Their accepted endpoint interval contains only the two endpoints. Thus their strict optimum is the direct **16/34** event. Their ordinary optimum is only **8/18**, attained by

    31171 -> 7907 -> 99971 -> 75320,

whose event costs are `6/12`, `8/18`, `8/16`. This literal route changes source 1691 three times and sources 2945 and 3369 twice, hence reacquires departed alternative resources. Any route under 16/34 between these primitive endpoints must violate strict endpoint-monotonicity, although it need not use this particular detour. The frozen certificate contains the complete native full-record arrivals and changed-clock/source lists for both routes.

This witness separates an **individual endpoint optimum** (ordinary 8/18, strict 16/34) from the **uniform grade-9 optimum** (ordinary 14/30, strict 16/34). Neither should be substituted for the other.

For every other nonsingleton grade, the certificate likewise retains a complete ordinary lower-witness fibre, components below each coordinate threshold, and a primitive strict lower-witness interval.

## 5. Native authenticity and exact physical guards

The new packet includes the frozen native-repair packet byte-for-byte. Its original independent verifier is re-executed: all 82 canonical menus and full resource-kernel lifts, native target classification, source/K/converse alternatives, original indexed/block pins, owner-4 shared tree, all 164 disjoint normalized alternatives, and black capacity are authenticated from the original native carrier, not from charge labels alone.

The new constructor independently rejoins the original electrical labels and literal path differences, checks all source frames and alternatives, and reconstructs the actual receiver variables. Every independent base bit has literal cost 2; the triple has costs 2+3+3. Native target values U and V agree in each original potential, so every switched source row is unchanged. The target-row change is computed directly as `2 sum(f(source)-f(U))`, with its opposite at V. All six original generating potentials and their moduli are accumulated for **every one of the 131072 states**, and checked constant on each exact `(grade,Q)` fibre. The frozen exterior contributes a common constant. This verifies actual target diagonals and the paired-black factor 2, including the radical coordinate not explicitly named in Q.

The original complete dyadic filtration audit is also bundled and re-executed. Its established saturation equivalence is inherited, rather than inferred from generator spot checks. The guard scope remains original aggregate nodewise dyadic currents at fixed margins, not individual marked-port phases, odd-primary currents or exact integer-X currents.

Every route vertex is a fully authenticated assignment; every route edge is an accepted simultaneous full-record replacement. No extra source or auxiliary kernel unit is introduced. Exterior records never move. Fixed core/helper populations remain separately exact. The literal grade-8 and grade-9 controls replay all 82 records, capacity-one reservations and full physical generator rows at every arrival.

## 6. Exact algorithms and portable evidence

`check.py` uses only the Python standard library and explicit exceptions, never correctness-critical assertions.

1. Enumerate all 17-bit assignments and compare every fibre multiplicity with an independent weighted generating-function recurrence (the triple has grade 3, not grade 1).
2. For each fibre, form the complete paired-cost graph. Kruskal trees give every ordinary minimax endpoint bottleneck.
3. Compute every strict endpoint bottleneck by exact dynamic programming over all accepted proper interval splits, increasing changed support. Reconstruct every selected strict route and verify disjoint target-directed changes and exact endpoint restoration.
4. Independently check that the largest strict endpoint optimum equals the largest primitive cost in every fibre.
5. Digest all pair optima with their actual endpoint masks, retain all 23444 exact fibre summaries in plaintext `fibres.tsv`, and all 2543 primitive trades, grade histograms, lower cuts and literal physical witnesses in `certificate.json`.
6. Check complement under `(grade,q) -> (19-grade,2032-q)` in `Q = Z/128 x F2^4`. This is a counting/isometric identity, not an executed physical move.

The pair digest streams compact JSON rows `[grade,q,mask_i,mask_j,[ordinary_slots,clocks],[strict_slots,clocks]]` followed by newline, in grade/charge order and then ascending masks with `i<j`. Its SHA-256 is `6e40c8385f1a7e842f703e6247b31f223aa32cad522fe3310b9bfef75f7ea2b4`. An independently implemented full endpoint-monotone DAG and ordinary minimax calculation obtained the same digest. Clock values uniquely recover the simultaneous slot optimum as `2 floor(clocks/4)` in this specific receiver; the frozen digest nevertheless includes both coordinates explicitly. Every check recomputes all endpoint data rather than trusting the summary.

`reproduce.py` checks the closed manifest, executes the new checker, the frozen native-repair checker, the original complete native authentication and the original complete dyadic audit, each in normal Python and `-O`, requiring byte-identical output for each pair. `verify_manifest.py --relocate` copies exactly the manifest-listed files into a fresh directory and runs that full suite there. No network, Git, original repository, parent directory or absolute workspace path is required.
