# Independent adversarial audit: five-run Fermat-factor separation

**PASS.** The exact minimum is five nonempty O/Z runs, within the specified source grammar and excluding the final C1-closing J. This is not minimal graph order, an obstruction to arbitrary tournaments, or a Chan–Pak completeness theorem.

## Independent finite certificate

The auditor's `verify_runs_independent.cpp` does not use the producer's block-power formula. It visits each positive run length 1,...,64 by repeated single-letter source transitions, working modulo M=641·6700417=4294967297. Both starting letters and all alternating patterns with one to four runs are covered. Its counts are:

| Runs | Recipes checked | q-divisible | p-divisible | q-only |
|---|---:|---:|---:|---:|
| 1 | 128 | 2 | 2 | 0 |
| 2 | 8192 | 129 | 135 | 0 |
| 3 | 524288 | 644 | 1360 | 0 |
| 4 | 33554432 | 888 | 53200 | 0 |

Total:34,087,040 phase recipes, zero escapes. The empty word closes to H=3 and is not a counterexample. No assumption of injectivity from recipes to profiles or graph classes enters this proof.

## Why finite phases cover all lengths

For any source s=(D,E,U,d) and k>=1, direct induction gives

    O^k(s)=(D,0,2^k U,d),
    Z^k(s)=(D−E+(2^k−1)U−(2^k−2)d,d,d+2^k(U−d),d).

The producer's d=1 formulas are exactly these specializations. Exact checks give 2^32=−1 modulo each prime, hence order64. Thus k may be reduced to its POSITIVE representative in1,...,64. The representative of residue0 is64, never an empty run. In particular O^64 and Z^64 need not act identically to a zero-length block: they overwrite E. The producer's proof and program preserve this important boundary.

For a word with r maximal nonempty runs there are exactly two alternating letter patterns. Replacing every positive length by its positive phase representative retains that run pattern and all final residues. The finite zero-escape certificate therefore proves the implication q|H => p|H at every word length for r<=4.

## Witnesses and an independent infinite extension

Exact single-letter big-integer arithmetic confirms

    O^31 Z^15 O^9 Z^5 O^9

has closed count1182825546590390812645=176530139331685·6700417, residue275 modulo641. Its source length69, transitive-edge seed of order2, and closing birth give graph order72. Thus five runs suffice and the lower bound is sharp.

The producer extends the final run by64t. Independently extending the FIRST run instead gives another certified five-run family

    H_t=275398033407·2^(32+64t)−275389644827,
    t>=0, graph order72+64t.

Its residues remain q:0 and p:275; the positive leading coefficient proves the values strictly increase. This extension is an audited arithmetic consequence, not needed for the minimum-five proof.

The shorter source prefixOOZJOZJJOJJOZ has state(4909,859,6013,859). Appending20+64t O letters and the closure gives

    H_t=4909+6013·2^(20+64t), order36+64t,

with residue380 modulo641 and zero modulo6700417. Its seed is also the transitive edge. The earlier draft's order35 was wrong and is repaired in the final producer files. Root changes are NOT necessary: the five-run example avoids them in its source. Run complexity, source root-change count, and final graph order remain separate.

## Literal graph and implementation audit

The auditor's graph reconstruction starts from the complete forward chronological tournament, then reverses precisely the root-newborn edge at each Z/J birth and updates the root only at J. This construction is independent of the producer's per-birth edge insertion code. Every edge in BOTH36- and72-vertex published adjacency witnesses matches, and direct reachability checks verify strong connectivity.

A separate Python subset-DP counts literal Hamiltonian paths for every source word of lengths0,...,5 and for its closure:728 checks, maximum graph order8. These counts match source D and closure D−E+U. The source recurrence remains an inherited theorem for the large graphs; no exhaustive36/72-vertex Hamiltonian-path count is claimed.

The producer's18 independent-DP graph inputs through order18 also replay successfully with undefined-behavior sanitization. The complete producer four-run C++ program passes the same sanitizer. Its modular products use signed128-bit arithmetic, necessary because M² exceeds signed64-bit range; other reduced sums fit signed64-bit. The auditor's alternative code uses only sums of residues and doubling, safely below2^64. No assertion removal under optimization can suppress its explicit failure checks.

The producer's full run proof, empty/positive distinction, physical order accounting, exact witness values, arithmetic overflow boundaries, and scope claims all pass. No new error remains open.

## Reproduction

From this directory, with a C++17 compiler:

    c++ -O2 -std=c++17 verify_runs_independent.cpp -o /tmp/session219_runs_audit
    /tmp/session219_runs_audit
    python3 verify_witness_independent.py
    python3 verify_graphs_independent.py PATH_TO_graph_witnesses.json

The last argument is the producer's delivered `graph_witnesses.json`. JSON receipts record the exact universes. Producer implementation sanitizer transcripts are also retained. This audit did not modify the source checkout or publish files.

## Late addendum: independent full36-vertex weighted-DP count

**PASS.** Reviewed `FULL_GRAPH_COUNT.md` and `check_full_36_dp.cpp`, and replayed the program with `-fsanitize=undefined -fno-sanitize-recover=all`. It returns exactly

    B0=10063, B1=16935, B2=6013, full36=6305092397.

This strengthens the initial graph audit: the36-vertex graph now has an independent complete path-count certificate, using a weighted subset DP on its16-vertex outside skeleton. The72-vertex graph still uses the inherited source recurrence for its full count.

The auditor checked the actual published36-vertex adjacency: vertices15,...,34 are the transitive20-vertex O block; all first15 vertices point into it; every block vertex points to closing vertex35. Deleting this block gives exactly the16-vertex skeleton built by the new program. Inserting O vertices does not change the old root, so the closing edges in this skeleton agree with the source-prefix-plus-J construction.

Every maximal block segment is increasing. Any nonterminal segment must exit to the unique closing vertex, because no block vertex points to an old A vertex. Consequently there are at most two segments. Removing the block leaves an ordering of the outside vertices, with each nonempty block segment represented by exactly one marked slot: an initial slot only when the first outside vertex is c; an internal slot only at A-to-c; or a final slot only after A. These are precisely the extra z choices implemented by initialization, transitions, and final accumulation. Native and marked A-to-c transitions are distinct valid choices even when the native edge exists. No zero-length segment is marked. Initial and internal slots cannot both occur because c is visited only once. Thus no degree-three coefficient is omitted by truncation to degree2.

For one marked slot the entire block has one increasing assignment. For two slots, any nonempty proper subset of the m labelled block vertices fills the first slot in increasing order and its complement fills the second, giving2^m−2 assignments. The separating outside vertex prevents segments from merging. This is an exact bijection, so H_m=B1+(2^m−2)B2 for every m>=1. The unweighted B0 correctly handles the separate m=0 graph; it does not contribute when the block is required to be nonempty. The formula establishes the full36 count at m20 and independently recovers the infinite20+64t family.

There is no arithmetic overflow: each outside ordering has at most two binary slot opportunities, bounding all skeleton counts by4·16!<2^64. Every intermediate DP count is no larger than this conservative bound. The only evaluated large-block multiplier is2^20−2, and its product with6013 and addition remain below2^64. The all-m extension is a mathematical identity, not a claim that the fixed-width program evaluates arbitrary m. All shift widths are within the16-vertex skeleton's unsigned32-bit mask range.

Reproduction of the audited late certificate:

    c++ -O2 -fsanitize=undefined -fno-sanitize-recover=all PATH_TO_check_full_36_dp.cpp -o /tmp/session219_full36_audit
    /tmp/session219_full36_audit

Receipt: `full_36_ubsan.out`. No remaining blocker was found.
