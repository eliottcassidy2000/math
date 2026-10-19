# Compact certificates for a consecutive interval of square-sum nice pairs

**Status:** Complete finite certificate. All 2,568 consecutive nice pairs passed the independent decoder in normal and optimized Python execution, with 17,173,589 exact checks. The separate recursion and protected transition tables remain required for an all-size result.

A square-sum nice pair at index \(n\) is a pair of permutations \(S_0,S_1\) of \([1,n]\) and \([1,n+1]\), each with square sums on adjacent entries, such that every shared integer occurs at positions of the same parity. Positions are one-based.

The state needed for the prescribed-seed lift is:

- Both paths begin \(1,3\).
- The even-length member ends at 8.
- The odd-length member ends at 10.

The target base interval is every integer

\[
53\le n\le2620.
\]

Its upper endpoint is \(49\cdot53+23\). Once the separate modulus-49 self-glue table is verified, strong induction gives this state for every \(n\ge53\): for any \(n\ge2621\), write

\[
n=49m+d,\qquad24\le d\le72,\qquad
m=\left\lfloor\frac{n-24}{49}\right\rfloor.
\]

Then \(53\le m<n\). The present base certificate alone does not establish that induction step or the final prescribed-seed theorem.

## Binary format

Each gzip file covers an inclusive consecutive index range \([a,b]\), recorded in `B_base_chunks.json`. Its decompressed bytes consist of:

1. The eight-byte ASCII/version signature `SSBAS1` followed by two zero bytes.
2. Two unsigned 32-bit little-endian integers: the first and last indices.
3. For each index \(n=a,a+1,\ldots,b\), first the \(n-1\) successive edge square roots of \(S_0\), then the \(n\) successive edge square roots of \(S_1\). Each root occupies exactly one byte.

Both paths start at \(x_1=1\). A root byte \(r_k\) reconstructs the next vertex by

\[
x_{k+1}=r_k^2-x_k.
\]

The lengths are fixed by the index, so no separators or extra length fields are needed. Every root is at most 72 because the largest path has 2621 vertices and \(2\cdot2621<73^2\). The encoder compared every decoded path with the original saved permutation before writing its chunk. The independent verifier does not read discovery output.

Chunks contain at most 256 indices. The manifest records both compressed and decompressed sizes and SHA-256 digests. Gzip timestamps are fixed to zero.

## Independent verification

Run from this directory:

```sh
python3 verify_B_base_interval_20261019.py > replay.json
python3 -O verify_B_base_interval_20261019.py > replay.optimized.json
cmp replay.json replay.optimized.json
cmp replay.json B_base_interval_controls.json
sha256sum -c SHA256SUMS
```

The standard-library verifier decodes every path, checks every vertex lies in the required interval and occurs only once, checks the initial edge and required endpoint, and checks equal positional parity for every common vertex. The recurrence guarantees every edge sum equals the square of its encoded integer root. It also checks all chunk ranges are consecutive and cover exactly the target interval. No assertion can be removed by optimized Python execution.

The full interval contains 2568 nice pairs, or 5136 paths, and encodes 6,861,696 edge roots. This is a finite certificate for the specified states, not merely a record of heuristic successes.

The affine recursion into which these bases fit follows the architecture of Robert Gerbicz's 2018 solution of ordinary interval square-sum Hamiltonicity. See [his author comment on OEIS A090461](https://oeis.org/A090461) and [the original discussion](https://www.mersenneforum.org/node/17331/page2). This packet concerns the particular prefix, endpoints, and compatible parity states required by the prescribed-seed extension.
