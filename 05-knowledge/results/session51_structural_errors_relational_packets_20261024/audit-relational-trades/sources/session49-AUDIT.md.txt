# Independent audit: exact receiver route thresholds

Status: **PASS, finite-exact native threshold classification and proved finite-binary structural criterion.** This is a new Session49 audit, not a modification of the earlier frozen native-hole or repair-theory audits.

## Scope and preserved data

The receiver consists of the sixteen original cost-two sources

1394, 1477, 1480, 1549, 1638, 1691, 1921, 1974, 2176, 2185, 2581, 2837, 2905, 2945, 3000, 3369,

and the three linked sources 1457, 1841, 3015, whose target bits must all agree. Every exterior full record, including every helper, is fixed. A state retains the original literal U/V path record at each source. No new port, source, kernel lift, population allowance or external helper transition is introduced.

The compact code uses the listed sixteen source bits in increasing-source order and bit16 for the triple. This is seventeen independent bits but **nineteen occupied source slots**. The grade of a compact mask is its low-sixteen-bit population plus three times bit16. The charge is in C128 × C2^4, encoded as16X+Y with the inherited Y1,...,Y4 coordinates.

The independent audit reconstructs labels directly from the original electrical potential arrays, and costs from disagreement of the original full U/V paths. The triple charge is1263 and its cost is2+3+3=8 clocks. For a pair differing at b base bits and triple indicator t, its exact event cost is

    slots = b+3t,       clocks = 2b+8t = 2 slots+2t.

The original complete native authentication and complete dyadic-filtration audit are bundled unchanged and re-executed. They retain all source/K/converse records, full kernels, indexed/block pins, literal owner4 shared-tree relation, disjoint normalized alternatives, original paired-black factor2 and target diagonals. The independent threshold script additionally checks exact original potential-row agreement for every assignment in each grade/charge fibre, not just charge equality by definition.

## Exact finite answer

There are131,072 accepted receiver assignments,23,444 nonempty grade/charge fibres, at most22 states in a fibre, and474,896 unordered endpoint pairs. Empty charge fibres make no connectivity assertion. Singleton fibres need zero events.

The exact simultaneous uniform thresholds over all fibres at each grade are:

| Grades | Ordinary slots/clocks | Strict slots/clocks |
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

Thus the uniform threshold over all grades is16 slots and34 clocks for either ordinary or strict routing; within grades9 and10, forbidding resource reacquisition raises both uniform thresholds. “Uniform” is essential: this does not claim every individual fibre requires16/34.

At fixed grade every event changes an even number of slots, and its clock cost is either twice that slot count or two more. Consequently all actual event-cost pairs are a single coordinatewise chain. The reported simultaneous thresholds are not formed by combining incompatible minimizers of two independent objectives. The audit checks this chain across all474,896 pairs and provides exact joint path/cut certificates.

There are2,543 distinct balanced interval-indecomposable trades after removing common V bits, including196 with cost16/34. There can be many embeddings of one trade into different fibres; the number2,543 does not count these embeddings.

## Independent algorithms and agreement

No threshold producer code is imported or executed. The main script implements four independent computations:

1. Dense Prim spanning trees give ordinary uniform minimax thresholds per fibre.
2. Minimax Floyd closure gives the ordinary optimum for every endpoint pair and agrees with each Prim result.
3. A full source-relative subset-poset DAG computes strict optima from every source to every target. Edges may jump between any accepted states whose switched-coordinate sets increase.
4. A separate endpoint-interval recursion computes strict optima by splitting at each proper accepted interval state. It agrees with the DAG for every pair, and its indecomposable maximum agrees with each fibre's strict threshold.

Each ordinary grade maximum has explicit slot-below and clock-below disconnected-component certificates. Each strict grade maximum has an interval-indecomposable lower witness. AUDIT_RESULT.json records all grade histograms and witness fibres; FIBRE_SUMMARY.tsv retains the complete23,444-row exact fibre summary. The independently derived list of all2,543 primitive trades and all grade endpoint-optimum histograms are also compared to the frozen producer certificate, after computation. The canonical complete endpoint-optimum stream has SHA256

    6e40c8385f1a7e842f703e6247b31f223aa32cad522fe3310b9bfef75f7ea2b4

This independently matches the producer's stream. The stream orders grade, charge, ascending compact endpoint codes and serializes one compact JSON row [grade,charge,start,finish,[ordinary slots,clocks],[strict slots,clocks]] plus a newline for every unordered pair.

## Exact native obstruction: grade9, charge3

The complete fibre has codes

    7907,20452,27198,31171,39819,56652,61224,66460,70748,75320,99971.

For start31171 and finish75320, the exact ordinary optimum is8 slots/18 clocks. A realizing path is

    31171 → 7907 → 99971 → 75320,

with event costs(6,12),(8,18),(8,16). Source1691 switches three times, and sources2945 and3369 each switch twice. These are actual released-resource reacquisitions in the authenticated binary menus, not merely repeated projected labels.

Below either8 slots or18 clocks, the start component is exactly

    {7907,20452,27198,31171,39819,56652,61224},

which excludes75320. The accepted endpoint interval for31171,75320 contains only those two endpoints. Hence every strict route must be the direct16-slot/34-clock event. The audit independently checks the literal path, both below-threshold cuts, and the complete accepted endpoint interval.

The grade8/charge53 fibre remains an independent sharp control: it has exactly two accepted states, both grade8, and their only event costs16/34. Its complete literal V-source lists are retained in AUDIT_RESULT.json.

## General finite-binary theorem, including joint budgets

Let Ω be any nonempty finite accepted subset of {0,1}^n representing full literal records. Each coordinate has two distinct nonempty resource alternatives, with all alternatives disjoint across coordinates and across the two choices. Unchanged resources remain continuously occupied during an atomic event. An event may replace any subset of coordinates simultaneously, provided its arrival is in Ω. Use any number of nonnegative additive per-coordinate event costs, including actual slot and clock costs.

For x,y∈Ω define

    IΩ(x,y) = {z∈Ω : z_i=x_i=y_i wherever x_i=y_i}.

Call a distinct pair interval-indecomposable when IΩ(x,y)={x,y}. Let M_j be the maximum j-th event cost among all such pairs, with maximum0 for a singleton domain.

**Claim.** A coordinatewise event budget B connects every pair of Ω without reacquiring a released resource if and only if B_j≥M_j for every cost coordinate j. These bounds are jointly sufficient.

**Proof.** In a binary coordinate, a second switch during the route reoccupies the alternative released by the first switch. Thus a strict route switches each endpoint-differing coordinate exactly once and every endpoint-fixed coordinate zero times. Conversely, such a monotone route has no reacquisition, since distinct source alternatives have disjoint reservations. Every strict intermediate therefore lies in IΩ(x,y).

If x,y are interval-indecomposable, a nontrivial strict route has no possible intermediate; it must use the direct x→y event. Necessity follows in every cost coordinate simultaneously.

For sufficiency, induct on the number of differing coordinates. An indecomposable pair uses its direct event, within B. Otherwise choose any z∈IΩ(x,y)\{x,y}. Both endpoint subpairs have strictly smaller supports, and their difference supports are disjoint with union equal to the x,y support. By induction route x→z and then z→y. The supports of the two routes are disjoint, so concatenation does not reacquire a released alternative. Each event meets all coordinates of the same B. Induction terminates because the differing-coordinate count decreases. ∎

The proof works for arbitrary relational acceptance predicates once their full-record domain is authenticated. It does **not** bound the indecomposable sizes independently of inventory size. For example, the accepted family containing only0^n and1^n forces one n-coordinate event. Thus it does not solve arbitrary-relation all-size bounded locality or create missing all-order native charged supply.

As a finite proof sanity check separate from the native census, the audit enumerates all65,535 nonempty accepted subsets of the four-bit cube with nonuniform clock weights(1,3,2,5), and tests six simultaneous slot/clock budgets against full monotone reachability:393,210 domain/budget checks. This supports but does not replace the displayed all-size proof.

## Scope exclusions and reproduction

Strict means no reacquisition of a released normalized alternative; it is stronger here than merely avoiding repeated complete assemblies. Atomic simultaneous replacements are the transition contract. A different microscopic transition contract requires its own proof. The source-frame host is not a new single global Hamiltonian observer. The claims do not free frozen exterior records, preserve ordinary individual marked-port phases, establish odd-primary or exact integer-X currents, or produce an all-size electrical bank.

Run:

    python -B reproduce.py
    python -B verify_manifest.py --relocate

The first command runs the independent threshold audit, unchanged full native authentication and unchanged full filtration replay in normal and optimized Python and requires identical outputs. The second copies **only manifest-listed files** to a fresh temporary directory and runs that complete suite there. No Git, network, producer imports, absolute paths or original checkout are required. SOURCEPINS.json pins every inherited input and historical source snapshot; historical documents and producer snapshots are stored as text archives. MANIFEST.json pins the complete deliverable except itself.
