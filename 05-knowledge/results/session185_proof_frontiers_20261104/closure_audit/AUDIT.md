# Independent audit of the larger-SCC escape packet

**Verdict: PASS within the stated scope.** Audited the frozen `closure/PROOF.md`, `verify.py`, and `witnesses.json` on2026-11-04, including the strengthened coprime-constant valuation family. No recursive closure, minimum suffix order, native legality, or global spectrum theorem is established or endorsed.

## Independent finite checks

`audit.py` imports no code from the author. It reconstructs all seven signature coordinates for the nine source/new states by enumerating every directed partial permutation, rather than the author's subset recurrence. The six new signatures agree exactly with `witnesses.json`. Literal masks are checked as tournaments. Floyd–Warshall reachability independently verifies all six new suffixes are strong.

A separately written attachment builder constructs the full graphs. Direct recursive enumeration of directed Hamilton paths, with no subset-count memoization, checks the full counts at clocks1,2 and the common-root-clock4 specialization. Full-graph strongness is checked independently. The common counts are exactly

    107→323→973→{2921,2923},
    107→325→979→{2939,2941}.

Exact rational scale comparison independently recovers every stated local shift:

    2→24:−2; 4→44:−3; 24→242:0;
    24→244:−1; 44→442:+1; 44→444:+2.

The constants satisfy c_child=3c_parent+d on all eight edges, including entry edges. The minimum original clock4 is legal for all displayed paths. The full seven-coordinate source fields were recomputed rather than inferred from counts.

## All-clock proof

The inherited attachment formula H_t=c+2^t F reduces each claim to the checked constant and slope equalities. These are genuine identities for every legal source/target clock. A finite graph witness with these exact coordinates proves its all-clock edge; no extrapolation from the two test clocks is involved. Port sets are nonempty on the strong suffixes, establishing the terminal-port hypothesis. The proof correctly distinguishes numerical replacement by a fully specified new graph from preservation of old names or native-operation legality.

## Unbounded normalized valuation with c coprime to3

The final seed rows(24,1,3,6,14), ports{1,2,3}, have signature(11,100,64,49,28,9,6), independently counted by partial permutations. It is strong.

For the20-vertex transitive block, I independently count the joint coefficient without using singleton multiplication. If i is the maximum of the first path and j is the minimum of the second, both ports, the partition count is2^(i−j−1) when i>j; for i<j it is1 only when j=i+1. Summing this over the ports of mask671088 gives J=38626. The ordinary boundary coefficients are F=sum2^i=671088 and G=sum2^(19−i)=60357; the resulting full signature is exactly(1,1048576,671088,60357,38626,1,0).

The integer identity25·671088=16(2^20−1) proves the telescoping geometric sum F(B_k)=2^(6+20k). The block's residue modulo3 is(1,1,0,0,1,1,0); since h remains11, each append adds1 to c modulo3. The seed has c=−17≡1. Therefore c(B_k)≡1+k and infinitely many k have c coprime to3. The first component is marked, the final singleton is marked, and all SCCs have order at most5. No removable initial zero-port prefix is involved. All statements follow for arbitrary k, not merely the100 independently replayed instances.

This refutes a uniform normalized-v2 bound for the enlarged grammar, not a bound depending on the exact target c. The proof retains that essential distinction.

## Reproduction and source boundary

Both author and independent audit scripts run successfully normally and under `python -O`. Each pair's outputs agree. Audit artifacts are `audit.py`, `receipt.json`, `normal.out`, `optimized.out`.

The S182/S184 seven-coordinate law and SCC≤4 obstruction are credited antecedents. Earlier S180 uniform paired constructions are acknowledged, and the packet does not claim to be the first such construction. Targeted current-main exact-tuple/signature searches are useful negative evidence but do not establish historical novelty. The discovery bank sizes and any exhaustive-minimality assertion were not audited; the proof does not depend on them or claim minimality.
