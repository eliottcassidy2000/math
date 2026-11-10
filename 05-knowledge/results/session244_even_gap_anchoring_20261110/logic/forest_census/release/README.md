# Exhaustive chronological-forest census and a normal-form obstruction

Status: **FINITE-EXACT through order 12, with all stored witnesses independently
recounted.** The fixed-order depth-two scalar-image conjecture is **REFUTED**:
H=97309 requires depth three at order 12. It nevertheless has an explicit
depth-one witness at order 20. The graph-preserving adjacent-swap normal form
below is **PROVED**. Arbitrary allowed-odd scalar coverage remains **OPEN**.

## Universe and inherited mechanism

Vertex 0 is the initial singleton. At birth v>=1, choose a sink (choice 0),
or a unique old vertex r beaten by v (choice r+1). Every other old vertex
beats v. There are product_(v=1)^(n-1)(v+1)=n! chronological histories at
order n. The census enumerates all of them, without a reachability, fixed-root,
strong-connectivity, or connected-forest restriction.

After arbitrary relabeling this is exactly the entire class of
outdegree-1-degenerate tournaments: every nonempty induced subtournament has
a vertex of outdegree at most one. Reverse a low-outdegree elimination order
for the converse. See the proved guarded-grammar-plus-ordinal-joins repair
in `../../ROUND3.md` and THM-5116,
`01-canon/theorems/THM-5116-sharp-outdegree-core-bound-and-six-core-forest-construction-class.md`.
The chronological reversal-forest observation is inherited from Session199,
`session199/hidden-structure/PROOF.md`, Section 3; that older fixed-root class
has chronological caterpillar forests and is narrower.

A code is accumulated as code=code*(v+1)+choice. Decode from the right using
remainders modulo n,n-1,...,2. A sink birth has forest depth zero; a birth
beating r has depth 1+depth(r). The maximum vertex depth is the history depth.
This is a presentation statistic, not a graph invariant unless minimized
over all valid orderings of the same graph.

Each CSV row has four columns:

- H: the number of directed Hamilton paths;
- chronological_code: the first history attaining H in code order;
- min_forest_depth: minimum history depth attaining that scalar at that order;
- min_depth_code: the first history attaining that scalar at that minimum depth.

The minimization is over the entire enumerated history universe, not just
stored representatives. It is not a claim about minimum graph order, the
minimum depth of the first graph, or arbitrary larger orders.

## Algorithm and integer safety

For every induced subset the main DP retains path counts starting at each
vertex, ending at each vertex, and their total H. The empty-set total is one.
A birth changes only subsets containing the new highest bit; DFS can reuse
all smaller subset tables and restore adjacency when backtracking.

At the last layer, let D=H(V), let E_r count paths ending at r, and compute

    U_r = sum_(S contains r) starts(S,r) H(V\S).

The sink child has count D; the child beating only r has count D+U_r-E_r.
This inherited singleton graft identity holds for every root, reachable or
not. It avoids recomputing a full child DP at the final layer. The completed
run explicitly checks that each order's processed total is exactly n!.

The supported range is 2<=max_order<=12. Codes and all path counts fit in
uint32 because 12!=479001600<2^32. For an n-vertex parent, U_r<=n!: for each
subset size k there are at most C(n-1,k-1)(k-1)!(n-k)!=(n-1)! pairs. At the
last layer n<=11, so even D+U_r<=2*11! fits before subtracting E_r. All DP
summands are nonnegative and bounded by the corresponding path/permutation
counts. History totals use uint64. No floating-point quantity influences
any enumerated result; wall times are descriptive only.

## Finite results

| Order | Histories | Distinct H | Maximum H | First allowed miss |
| --- | --- | --- | --- | --- |
| 1 | 1 | 1 | 1 | 3 |
| 2 | 2 | 1 | 1 | 3 |
| 3 | 6 | 2 | 3 | 5 |
| 4 | 24 | 3 | 5 | 9 |
| 5 | 120 | 7 | 15 | 17 |
| 6 | 720 | 16 | 37 | 35 |
| 7 | 5040 | 55 | 123 | 63 |
| 8 | 40320 | 192 | 413 | 311 |
| 9 | 362880 | 706 | 1509 | 991 |
| 10 | 3628800 | 2817 | 6083 | 4585 |
| 11 | 39916800 | 11895 | 25827 | 20407 |
| 12 | 479001600 | 53549 | 115781 | 91081 |

At every order through 11 the full scalar image equals the depth-at-most-two
image. At order 12 there are 53549 scalar values: one has minimum depth zero,
46092 have minimum depth one, 7455 have minimum depth two, and exactly one
has minimum depth three. That exception is **H=97309**, with first and
minimum-depth history code **25112126**. The first order-12 scalar requiring
depth two is 53443. Depth three suffices for every order-12 scalar.

A separately implemented shallow-only census, using the proved adjacent-swap
normalization, exhausts 23670073 normalized shallow histories and obtains
exactly 53548 values. Its complete image equals the full census's depth-at-
most-two image and excludes only 97309 from the complete forest spectrum.
Its independent roots/depth-one/last-depth recurrence checks the history
total; see NORMALIZED_SHALLOW_AUDIT.json. This confirms scalar absence at the
stated order, not absence at larger orders.

All **138488** literal stored witnesses (both codes of all 69244 CSV rows)
pass the independent C++ endpoint-DP audit. The order-one-through-eleven CSVs
also reproduce the prior completed depth-aware checkpoint byte for byte;
in the order-12 run those order-eleven values were computed by full subset
DP rather than the prior run's final-layer graft shortcut. The full
479001600-history run took 1805.13 seconds on this environment.


These are class-specific fixed-order results, not a new global coverage
cutoff. Sink padding makes scalar spectra nested by order but does not show
that every future allowed target eventually appears. “First allowed miss”
means the least positive odd integer absent from that order's spectrum after
excluding 7 and 21. The packet makes no new global nonexistence claim.

## A proved graph-preserving reduction

If consecutive vertices x,y in a chronological history satisfy y->x, then
x is y's unique earlier outneighbor. Swap their positions while keeping the
tournament itself fixed. Only their mutual relative order changes: that
single reversal edge disappears, y becomes a forest root, and x keeps its
previous parent. Every other parent incidence is unchanged. Forest depths
therefore cannot increase.

Repeat while such an adjacent reversal exists. Each swap deletes one forest
edge, so at most n-1 swaps suffice. The resulting history never chooses the
immediate predecessor as parent. Conversely such histories are ordinary
forest histories. Thus this normalization preserves the entire graph class
and every attainable minimum-depth bound, without losing any root response.
There are exactly product_(v=1)^(n-1) v=(n-1)! normalized histories.

The shallow replay adds the depth-at-most-two restriction during generation.
Its counting recurrence tracks r depth-zero vertices, s depth-one vertices,
and the last vertex's depth l. A birth has one sink choice, r-[l=0] depth-one
choices, and s-[l=1] depth-two choices. This independently counts the
normalized search universe; the unnormalized variant instead uses r and s.
This is a partial flattening theorem, not a depth-two flattening theorem.

## A rigorous limit on interpreting shallow scalar witnesses

`NORMAL_FORM.md` proves that the graph with order-nine code 26485, H=217,
has minimum chronological presentation depth exactly three. Sink padding
gives such graphs at every order n>=9, still with H=217. Yet H=217 itself
already has a depth-one order-nine representative.

There is also an independently computed fixed-order transition obstruction.
The same deep parent has a child of count 373, but none of the 458 depth-at-
most-two order-nine histories with H=217 has such a child. All 176850 shallow
histories are checked, using literal DP rather than the graft identity.
Thus scalar equality cannot be used to discard root/attachment information
and reproduce every continuation. These order-nine obstructions alone do not
refute bare scalar-image equality; the separate order-twelve value 97309 does.
Yet 97309=31*43*73: literal depth-one factors of orders 6,7,7 give a depth-one
order-twenty ordinal-join witness. See NORMAL_FORM.md and
DEPTH_EXCEPTION_CHECK.json. The fixed-order refutation must not be presented
as a variable-order scalar obstruction.

## Independent audits and hostile controls

`audit_all_witnesses.cpp` decodes both codes of every published CSV row,
constructs the literal tournament, and recounts it by a direct endpoint subset
DP. It does not use U, the graft identity, or parent-table reuse. It checks
the actual minimum-depth witness depth, strict CSV structure, code ranges,
sorted distinct values, and the stored depth's compatibility with the first
witness. Positive recounting alone does not establish census completeness or
depth minimality; those depend on the exhaustive census.

`audit.py` independently counts deterministic samples in Python, checks every
summary field against the CSV data, and exhausts all labelled tournaments
through order six using low-outdegree elimination instead of history
construction. The full labelled degeneracy-class sizes through six are
1,2,8,64,1000,26720. This different-universe check reproduces the corresponding
complete scalar spectra. It also checks sink-padding nesting. The script is
read-only and uses explicit exceptions, not assertions.

`shallow_census.cpp` is a separate depth-restricted implementation. It never
visits depth-three births and stores no minimum-depth array. Its normalized
run compares the full computed scalar image against every order-twelve CSV
row, with an independent combinatorial history-total check. It shares the
proved singleton graft identity, but not the main census's unrestricted
history universe or depth-minimum accumulation.

`verify_depth_exception.py` independently checks the complete order-twelve
witness profiles, forced reverse deletions, the three depth-one factors,
and every cross edge of their order-twenty join. Counting that join uses the
proved exact product law, not a claimed direct twenty-vertex DP.

`verify_normal_form.cpp` independently recounts every shallow order-nine
history, all literal children in the H=217 fiber, and all 9! vertex orderings
of the deep example. `NORMAL_FORM_CHECK.json` retains all 30 shallow
child-count signatures and the 16 valid ordering depths.

`verify_controls.py` and `CONTROLS.json` record strict warning-free builds,
normal/-DNDEBUG byte equality through order eight, a full small-run witness
recount, Python normal/-O equality, the normal-form independent replay,
rejection of order 13 and malformed/wrapping inputs, and rejection of a
corrupted private CSV by both auditors. No published CSV is corrupted or
modified by the control script. Checks are active under -DNDEBUG and -O.

## Exact source provenance

`RUN_PROVENANCE.json` records the original running-source and binary hashes.
Only command-line validation was hardened after the order-12 run had started:
the bound is now checked before narrowing the parsed integer, and nondigit
inputs are rejected. The complete diff is `ARGUMENT_VALIDATION.patch`.
All enumeration, DP, depth tracking, and output code is unchanged; the valid
argument 12 takes the identical computation path. Reverse the patch to
reconstruct the exact running source. The controls verify its recorded hash,
compile it, and compare its order-one-through-eight outputs byte for byte
against the published wrapper. This source distinction is explicit, not a
claim that the published file was byte-identical to the running source.

## Reproduction

Use a fresh output directory; the census writes one CSV per order.

    g++ -std=c++17 -O3 -Wall -Wextra -Wconversion -Werror -pedantic census.cpp -o /tmp/forest-census
    g++ -std=c++17 -O3 -Wall -Wextra -Wconversion -Werror -pedantic audit_all_witnesses.cpp -o /tmp/forest-audit
    g++ -std=c++17 -O3 -Wall -Wextra -Wconversion -Werror -pedantic verify_normal_form.cpp -o /tmp/forest-normal-form
    mkdir OUTPUT_DIRECTORY
    /tmp/forest-census 12 OUTPUT_DIRECTORY > OUTPUT_DIRECTORY/SUMMARY.jsonl 2> OUTPUT_DIRECTORY/CENSUS_PROGRESS.log
    /tmp/forest-audit 12 OUTPUT_DIRECTORY > OUTPUT_DIRECTORY/ALL_WITNESSES.jsonl 2> OUTPUT_DIRECTORY/ALL_WITNESSES.log
    python3 -I -B audit.py OUTPUT_DIRECTORY > OUTPUT_DIRECTORY/AUDIT.json
    python3 -I -B -O audit.py OUTPUT_DIRECTORY > OUTPUT_DIRECTORY/AUDIT_OPT.json
    cmp OUTPUT_DIRECTORY/AUDIT.json OUTPUT_DIRECTORY/AUDIT_OPT.json
    /tmp/forest-normal-form > NORMAL_FORM_REPLAY.json
    cmp NORMAL_FORM_CHECK.json NORMAL_FORM_REPLAY.json
    g++ -std=c++17 -O3 -Wall -Wextra -Wconversion -Werror -pedantic shallow_census.cpp -o /tmp/forest-shallow
    /tmp/forest-shallow 12 OUTPUT_DIRECTORY 1 > NORMALIZED_SHALLOW_REPLAY.json
    cmp NORMALIZED_SHALLOW_AUDIT.json NORMALIZED_SHALLOW_REPLAY.json
    python3 -I -B verify_depth_exception.py > DEPTH_EXCEPTION_REPLAY.json
    cmp DEPTH_EXCEPTION_CHECK.json DEPTH_EXCEPTION_REPLAY.json
    python3 -I -B verify_controls.py > CONTROLS_REPLAY.json
    cmp CONTROLS.json CONTROLS_REPLAY.json
    sha256sum -c SHA256SUMS

Run the manifest check from this release directory. Replaying into a fresh
directory keeps the published packet read-only. Temporary binaries, abandoned
runs, discovery probes, and duplicate checkpoints are deliberately excluded
from the release packet. All such prior material was preserved outside it.
