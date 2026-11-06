# Exact H=3 metagraph incidence at every tournament size

Every tournament with exactly three Hamilton paths is Tr_a -> C3 -> Tr_b. This packet explicitly computes the unique ordinary-class line neighbour for every a,b, including its Hamilton count and actual ordinary tiling-fibre degree. Its two-variable generating function has denominator (1-x-x²-x³)(1-y-y²-y³) and an explicit numerator; the counts obey tribonacci recurrences independently in a and b.

The proof classifies the directed cycles created or removed by the relevant chord flip and uses the inherited odd-cycle collection formula. It is not a fitted recurrence or an extrapolation from a census.

Consequences include:

- Complete possible count pairs (3,k), and their exact labelled edge multiplicities, at every order.
- Rigidity of every target, so its ordinary count equals its ordinary fibre degree.
- An exact size-indexed generating function for the sum of target degrees per source.
- Concrete failure of strong lumpability by Hamilton count: at order5 two H=3 sources have all neighbours at H=9 and H=13 respectively.

PROOF.md contains the full statements, proofs, source conventions and provenance. REVIEW.md records independent mathematical review. No global count-spectrum completeness, general connectivity, target converse classification or collapsed-simple-graph degree claim is made.

Reproduce from this directory:

    python3 check.py > replay.json
    cmp replay.json RESULT.json
    python3 -O check.py > replay-opt.json
    cmp replay-opt.json RESULT.json
    python3 independent_audit.py > independent-replay.out
    cmp independent-replay.out independent_audit.out
    sha256sum -c SHA256SUMS

The main checker uses explicit checks that survive optimization:66 direct Hamilton-count cases through order13,21 exhaustive directed-cycle inventories through order8, all tied score cells checked, and rational-coefficient controls on 0<=a,b<=20. The independent checker uses forward path DP through order10 and DFS cycle inventories through order9 (28 cases). It uses assertions and must run without -O. Both scripts need only Python3.
