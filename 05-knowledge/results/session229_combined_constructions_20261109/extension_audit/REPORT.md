# Independent finite-extension audit

**PASS.** The scoped assertions in the last subsection of `PROOF.md` and in `EXTENSION_PROOF.md` are confirmed. No all-size obstruction or unrestricted coverage result follows.

## Independent method

`check.py` imports no producer module. It constructs every extension arc from the meaning “orientation bit i means appended vertex → i,” verifies the tournament axioms, and enumerates every nonempty directed path sequence recursively. Path counts indexed by vertex subset and first/last vertex give the defining source coordinates of ordered covering path pairs: H, Q, A, G, S, E, T. Empty paths contribute H(empty)=1; the T sum requires both paths nonempty. This is literal path enumeration, not the producer subset-path dynamic program.

Every cut mask is visited for every eligible extension. One global dictionary per experiment groups by the exact ordered tuple (H,Q,A,G,S,E), across both input cores and (in the second experiment) both resulting orders. Each T entry retains **all** physical witnesses, avoiding the producer's representative-only `setdefault` ambiguity. Receipt witnesses are independently checked against their exact core, orientation, adjacency, cut, signature and T. All 40 reported length-two runs are present and unique. The earlier receipt deliberately shows only ten of its 1,772 singleton runs; the independent census covers all of them. No assertions are used, so checks remain active under `-O`.

## Results

- Two H45 cores: all 128 extensions visited. Exactly 48 pass gcd(H,Q)=1, all from the first core. The 6,144 retained cut records give 1,772 interfaces, each with exactly one T value. Thus the claimed absence of a gap-two pair holds, with the stronger finite fact that there is no unequal-T pair at all.
- B95/B389: all 384 extensions visited. Exactly 86 and 154 respectively pass, totaling 240. The 100,864 retained cut records give 97,208 interfaces: 97,056 have one T value and 152 have two. Maximum step-two run length is two, achieved by exactly 40 runs.
- Every one of the 80 physical witnesses in these 40 runs has orientation (core 0, 0), (core 0, 127), (core 1, 0), or (core 1, 255). These mean appended sink or source. There are no hidden nonordinal witnesses behind chosen dictionary representatives.
- The nonordinal no-pair statement is correct even against all physical examples, not merely the displayed ones. No multi-T interface in this census contains a nonordinal witness at all. The text's “may still have unequal T values with a larger gap” is a cautious non-claim; this census in fact finds none for nonordinal orientations. The 112 remaining two-T fibers have larger gaps and only ordinal physical witnesses.

## Scope and reproduction

The checked universes retain all original core arcs, append exactly one vertex, and filter by gcd(H,Q)=1 before grouping all cuts. They say nothing about other cores, changed original arcs, multiple appended vertices, or a weaker interface. The text preserves these boundaries correctly.

Run `python check.py > CHECK.json`, `python -O check.py > CHECK.optimized.json`, then `cmp CHECK.json CHECK.optimized.json`. Both completed successfully with byte-identical receipts. The receipts include hashes of every audited input and canonical complete-fiber census hashes. `SHA256SUMS` covers the checker, report, and both receipts. No repository or git writes were performed.

## Portable release binding

The released checker reads only its sibling `../rich` packet, not an absolute research-workspace path. It imports no producer implementation. Both released receipts are regenerated against these exact packaged inputs; the proof received editorial wording cleanup only.
