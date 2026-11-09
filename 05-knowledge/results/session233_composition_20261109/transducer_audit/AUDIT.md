# Independent audit: bridge-linked block transfer

**Verdict: PASS, with explicit scope and inherited-result qualifications below.** This does not complete Chan–Pak. No mathematical correction to the submitted transfer or its 1,374-value census was necessary.

## Independent executable checks

`python3 -O audit.py` passes. No producer module is imported. Its polynomial arithmetic uses integer coefficient arrays. Its local matrices are built by a parity dynamic program on vertex assignments, rather than the producer's character-sum formulas. A separate directed-simple-path traversal constructs every odd cycle and every open directed path; a disjoint-mask recursion counts packings of the remaining vertices.

- All three polynomial coordinates, not merely their evaluations, agree on all 255 compositions of orders 1–8. This directly checks the claimed meanings of even/odd states, including arbitrary path heads and disjoint completed cycles.
- Independent Hamilton-path subset DP agrees at x=2 on all 511 compositions of orders 1–9.
- All 262,143 compositions of orders 1–18 are enumerated without pruning or merging states. Exactly 1,374 distinct counts are at most 10,000. Their sorted list is `values.json`; `receipt.json` includes its canonical comma-separated SHA-256. Every cumulative count in the producer receipt through order 18 agrees.
- All 176 eligible binary run words through length 10 have their canonical complementary endpoint adjacency exactly equal to the claimed bridge-block graph.
- `python3 -O audit.py --failure-sentinel` exits 1 with a deliberate RuntimeError, recorded in `sentinel-O.txt`.

`run-O.txt` contains the execution transcript. Its first printed receipt precedes the final Hamilton-path/run-window checks; `receipt.json` is the final receipt including those checks.

## Exhaustive transition proof

At any boundary between blocks there is exactly one forward crossing arc. Every cycle crossing the boundary has equally many forward and backward crossings, hence exactly one each. Consequently at most one cycle in a disjoint packing crosses it, and its restriction to the processed prefix is a single directed path ending at the last prefix vertex. For the intermediate states one may equivalently count *all* directed paths ending there, weighted by disjoint odd-cycle packings on the other vertices. No assumption on the initial vertex is needed.

For m>=2, partition configurations by the component using the unique incoming bridge:

1. No component uses it. The old prefix is closed, and either no new path is present or a new increasing path ends at the new last vertex. A transitive block has no internal cycle. Mandatory last vertex gives 2^(m-2) subsets of each parity.
2. The component is a newly closed odd cycle. The old path is extended through an increasing subset containing the first vertex, and returns backward to its old head. Its new subset parity is 1-p. A distinct new terminal path may be absent or present. If present, it contains the last vertex and is disjoint from the first subset. The middle vertices have labels unused/old-cycle/new-path. Their two parity counts give exactly N(1-p,r). The mandatory first/last vertices explain the parity shifts in the character sum.
3. The component remains the terminal path. It crosses forward exactly once. Its tail extension contains both first and last new vertices. Its initial vertex can be old, or can be a new middle vertex with an increasing head segment returning backward to the old head. Every middle vertex is unused, tail, or head. There can be at most one backward crossing, since the path ends on the new side and has only one available forward crossing. Thus this decomposition is exhaustive and uniquely recoverable. Counts by total middle parity are (3^(m-2) +/- (-1)^(m-2))/2.

All backward arcs used in cases 2 and 3 exist: only the incoming first-vertex/old-last-vertex pair is exceptional. That exceptional closure would be a two-cycle, and odd parity excludes it. A head extension uses middle vertices, so cannot use the exceptional pair. Nonempty head extension is essential: (1,3,1) has P=1+6x and H=13.

For m=1 the first and last vertices coincide. There is no disjoint close-and-restart transition or head extension. The only possibilities are leave it unused, start a length-one path, extend an old path by one vertex, or close an old even path with the singleton. These give exactly A_1. This separate treatment must remain in the theorem. Integrality and nonnegativity of the m>=2 formulas follow from the three-label counting, not divisibility assertions alone. Substitution x=2 gives the displayed matrices; the audit checks their first five against the receipt.

## All-order completeness and scope

For two or more blocks and n>=3, every consecutive forward edge and the final-to-initial backward edge exist, making a spanning directed cycle. This proves strongness. One block is transitive at every order, contributing the single value 1; the two-vertex exception is transitive as stated. These exceptions must not be subjected to Busch's strong bound.

Busch's strong-order lower bound is **CITED**, inherited from THM-4744, not reproved here. It is strictly increasing from n=3 onward; g(6)=15, g(7)=25, and g(19)=15,625. Therefore the exact order<=5 inventory excludes 11 for the whole bridge family; order<=6 excludes 19; and orders<=18 cover every nontrivial value <=10,000. The independent unpruned enumeration avoids relying on pruning soundness. The producer's pruning is also sound, either by induced-subtournament monotonicity at x=2 or directly by h' >= h in the nonnegative transfer. Large one-block graphs add no values beyond 1.

For run words with **every run length >=2**, distinct 1-run staircase windows are disjoint. The original endpoint is their ordinal sum with leftover singletons. Each nontrivial window has size>=4, and inherited staircase counts begin 5,9,17, increasing thereafter. Therefore no factor equals 3 or 11, and neither value occurs at this endpoint. The complement is the indicated bridge-block tournament. Original windows have size>=4 and uncovered vertices are singleton blocks; arbitrary compositions with parts 1,2,3,... enlarge that family. The census of 1,374 is for this enlarged bridge family, not an assertion that the restricted run-window endpoints attain all these values.

Since 11 is prime and all component counts are positive integers, ordinal products of either canonical endpoint family cannot equal 11. This argument does not exclude overlapping windows arising from other run hypotheses, noncanonical realizers of a tight word, general ribbons, or broader substitutions. It is a proof of incompleteness of this specific construction strategy, not of Chan–Pak or of its failure.

## Inheritance and novelty

THM-4951 already supplies the canonical staircase/run-window construction and tight packing-number splits. THM-4952 already supplies the odd-cycle packing viewpoint, ribbon finite-state transfer, and staircase polynomial/tribonacci specialization. The x=2 identification with Hamilton paths is the inherited odd-cycle formula, independently checked here but not newly proved. Busch's lower bound and ordinal multiplicativity are inherited.

The concrete new claim supported by this audit is the exact three-coordinate bridge-block compression (including the head-extension transition), plus its family-specific exact bounded census and endpoint obstruction. This audit checks the cited ancestors, not the entire research corpus or external literature, so it does not certify global novelty.

## Producer hardening and provenance

The producer originally used Python assert for its verification loops. The parent subsequently replaced these with explicit RuntimeError checks and added --fail sentinels to both scripts; no mathematical formulas changed. Inspection confirms these changes. The normal/optimized receipts are byte-identical. This auditor independently executed all four producer sentinel combinations (two scripts, normal and -O), each exiting 1 with the expected RuntimeError. See `producer_sentinel_checks.json`. Final source and receipt hashes are in `producer_hashes.json`.

No repository or producer file was modified by this auditor; no commit was made.
