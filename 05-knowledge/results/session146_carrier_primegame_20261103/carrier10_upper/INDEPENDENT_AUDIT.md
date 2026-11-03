# Independent order-ten kernel upper audit

For E={0,...,7}, the 735 pooled-parity signatures on ten old vertices have no isolates and endpoint set E. The file `kernel.txt` contains 56 integer columns indexed by the independent exact-cover representative order in `reps4.txt`. The verifier reconstructs that order from endpoint-pair paths and ordered interior subsets, validates the entire supplied list, and finds rank 56 over both GF(101) and GF(103). Hence the integer columns are rationally independent.

The independently generated canonical diagrams consist of all partial matchings on E; a matched pair is connected by a path of length 2 or 3, and an unmatched vertex has a stub of length 1 or 2. These parity choices give 32,400 diagrams. An independent disjoint-set/forest evaluator computes all 735 row entries for each diagram. In each nonzero row, it verifies exactly that the component count c and number of nontrivial components r are constant across nonzero entries, so their coefficient magnitudes c! 2^r / 2^(c−1) agree over the rationals. It then sums the signed integer coefficients against every one of the 56 candidate columns using integer arithmetic. All sums are zero; 17,042 diagrams are entirely zero. This exact integer test is not a modular vanishing test.

## Reduction to the canonical diagrams

Work in the no-isolate block with endpoint set E. If a selected exterior edge touches one of the two old interior vertices, or if two selected exterior edges meet an endpoint, every block coefficient vanishes because selected degree exceeds two. If an endpoint is untouched, attach a path of two fresh exterior edges to it. For every old forest in the block, this preserves the even-path condition, canonical path sign, and the counts c,r, because the endpoint's old component is already nontrivial and gains two vertices and two edges. Thus the restricted coefficient row is unchanged. The verifier checks this tail rule across all 735 block representatives on a mixed zero/nonzero sample.

After saturation, every endpoint has exterior degree one and the exterior selected graph must have maximum degree two and no cycle. Each component touching E is a path meeting either two endpoints (a bridge) or one endpoint (a stub). The rest are disjoint external paths and isolated vertices. An odd external-only path kills the whole row; every even external-only path or isolate affects only the common coefficient scalar. A bridge's length can be replaced by 2 or 3 and a stub's length by 1 or 2 according to its parity. For either traversal direction, equal-parity replacements differ by a fixed orientation sign; signs caused by external-label permutations are also a single row gauge. All nonzero entries of an individual row have the same c,r because the number of vertices and selected edges is fixed, and the isolated external vertices are fixed. Therefore the arbitrary row is a scalar multiple of one of the tested canonical rows. Loops, branching and cycles make it zero.

Consequently the block has at least 56 rational kernel directions. Relabelling yields 56 independent directions in each of the C(10,8)=45 disjoint endpoint-set coordinate blocks, for 2,520 total. The previously verified saturated-block lower certificate proves dimension at least 345,106, while the all-order carrier dimension is 347,626. The kernel upper bound matches: the full fixed-source-oblivious-attachment dimension at order ten is **345,106**, or 345,105 nonconstant coordinates for normalized laws. This is not an S/B-only or nonlinear claim.

For n≥10, lift each of these 2,520 directions to each ten-set A, with all other old vertices absent from its Fourier support. Any n-order exterior row restricted to that support is an order-ten exterior row, treating the other old vertices as exterior and leaving their old-old edges unselected. Relabelling gives only nonzero sign gauges. Different sets A give disjoint exact nonisolated supports, proving the all-order upper bound `dimension(n) ≤ c_n − 2520 C(n,10)`.

Hostiles cover exterior cycles, disconnected even and odd paths, exterior relabel signs, two-edge-tail saturation, zero and nonzero columns, and corruption of a kernel entry. These checks and the full integer enumeration are built into the verifier. The result inherits the finite coefficient formula and carrier upper bound from Session145 and the independent order-ten lower certificate.

Reproduce:

```sh
python3 verify.py
python3 -O verify.py
sha256sum -c SHA256SUMS
```

If `kernel.txt` is split into numbered parts for distribution, concatenate them in the prescribed order and check its joined SHA-256 before replay.
