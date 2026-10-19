# A complete conditional split-lift for the prescribed square-sum seed

## Result and scope

A square-sum path has distinct positive integer labels and a square sum on every consecutive pair. A **nice pair** consists of square-sum paths \(S_0,S_1\) that are permutations of \([1,n]\) and \([1,n+1]\), respectively, with each common vertex occurring at positions of the same parity. Positions are one-based.

**Theorem (complete final transformation).** Let \(n\ge127\), \(24\le d\le72\), and \(N=49n+d\). Put \(j=0\) when \(d\le56\), and \(j=1\) when \(d\ge57\). Suppose a nice pair satisfies:

1. Both paths begin \(69,127\).
2. The path \(S_j\) ends at 10, and the other ends at 8.
3. Both contain the edge \(2\!-!7\), with 2 at an odd position and 7 at an even position. Its traversal direction may differ between the two paths.

Then there is a square-sum Hamiltonian cycle on \([1,N]\) containing the exact unscaled subpath

\[
482,\quad3362,\quad6242.
\]

The construction cuts one edge in each of 49 affine lifts of the input paths and reconnects the resulting segments through the singleton reservoir \([1,24]\). A complete finite certificate covers all 392 combinations of input parity, residue, and the two cut-edge orientations. The independent verifier uses only integer arithmetic and the Python standard library.

The hypotheses on the input nice pair remain substantive. This theorem certifies every final transformation; it does not assert the required input pairs exist for every index. In particular, it is not yet an all-sufficiently-large-\(N\) prescribed-seed theorem.

Ordinary interval square-sum Hamiltonicity is already due to Robert Gerbicz: paths exist for every \(N\ge25\), and cycles for every \(N\ge32\). The affine modulus-49 architecture is attributed to his 2018 construction. The present certificate concerns the additional prescribed seed and the controlled splitting needed to preserve it. See [Gerbicz's author comment on OEIS A090461](https://oeis.org/A090461) and [the original discussion](https://www.mersenneforum.org/node/17331/page2). No claim of novelty for ordinary interval coverage is intended.

## 1. Exact interval coverage

For \(-24\le c\le24\), define the alternating affine lift

\[
T(c,i)[k]=49S_i[k]+(-1)^k c.
\]

Every internal edge remains square-sum: its two offsets cancel, and its sum is 49 times the original square.

For every common vertex \(x\le n\), parity compatibility implies that its two candidate lifts agree for each fixed \(c\). As \(c\) varies, these representatives are exactly

\[
[49x-24,49x+24].
\]

These disjoint intervals cover \([25,49n+24]\). The extra vertex \(n+1\) in \(S_1\) has positional parity \(n+1\), because the counts of odd and even positions in the consecutive-length permutations differ only at that parity. Thus its representative in block \(c\) is

\[
49(n+1)+(-1)^{n+1}c.
\]

Select

\[
i(c)=1\quad\Longleftrightarrow\quad49+(-1)^{n+1}c\le d.
\]

The selected blocks contain each integer in \([25,N]\) exactly once. Adding the singleton vertices 1 through 24 completes the interval. Reversing a block or cutting it into two nonempty segments does not alter this partition.

## 2. The controlled split and its ports

Write \(p=n\bmod2\). The beginning of every selected block is

\[
3381-c,\quad6223+c,
\]

since both inputs begin \(69,127\). Its terminal vertex is

\[
49e_i+(-1)^{p+i}c,
\qquad e_j=10,\quad e_{1-j}=8.
\]

Because 2 has odd position and 7 has even position, the lifted edge being cut has endpoints

\[
98-c,\qquad343+c.
\]

If the input traverses \(7\to2\), the two remaining segments have endpoint pairs

\[
(3381-c,343+c),\qquad
(98-c,49e_i+(-1)^{p+i}c).
\]

If it traverses \(2\to7\), interchange the two cut ports. These formulas involve only \(p,d\), and the two traversal directions, not the length or internal order of the paths. All 196 segment endpoints and the 24 singleton vertices are distinct; the certificate verifier checks this in every state.

The prefix \(69,127\) cannot be severed by cutting \(2\!-!7\). At \(c=19\), that retained prefix is

\[
3362,6242,
\qquad3362+6242=98^2.
\]

The port at 482 occurs at

\[
c_{10}=-8(-1)^{p+j},
\qquad49\cdot10+(-1)^{p+j}c_{10}=482.
\]

The choice rule selects input \(j\) at that offset. For \(j=0\), the extra-vertex threshold there is 57, so every \(d\le56\) selects \(S_0\). For \(j=1\), the threshold is 41, so every \(d\ge57\) selects \(S_1\). Since \(c_{10}=\pm8\ne19\), these are different blocks. Joining the 482 port to the 3362 port is valid because

\[
482+3362=62^2.
\]

## 3. The finite certificate

A certificate word contains:

- Each of the two segments of each offset \(c\) exactly once, in either orientation.
- Each singleton 1 through 24 exactly once.
- A square-sum connection between every consecutive pair of atoms.
- Initial singleton 1 and terminal singleton 8.
- The connection 482–3362.

Closing the word by the edge \(8\!-!1\), whose sum is \(3^2\), gives a cycle. Section 1 proves it spans precisely \([1,N]\). Section 2 proves it preserves the seed subpath, in one of its two orientations. This proves the theorem once every state has a valid word.

The data file contains all

\[
2\times49\times2\times2=392
\]

states. The independent verifier reconstructs the port formulas and input selections from scratch, checks the complete state set without repetitions, verifies segment and singleton coverage, and checks every connecting square with integer square root. It does not call an optimizer, use floating-point arithmetic, or trust solver status. The finite words were discovered with our own generated MILP models using installed SciPy/HiGHS; discovery software is unnecessary for replay.

## 4. Reproduction and receipt

Place the files in this directory together and run:

```sh
python3 verify_final_split_table_20261019.py > replay.json
python3 -O verify_final_split_table_20261019.py > replay.optimized.json
cmp replay.json replay.optimized.json
cmp replay.json final_split_table_controls.json
sha256sum -c SHA256SUMS
```

The two modes return identical receipts, with 234,025 explicit checks over 392 states. No check relies on a removable Python assertion.

The certificate file has 986,763 bytes and SHA-256

`3a035b14c8caa9f1293f66bbb0c9e77e9f1e933b7e5a07e523e2a5e7f59dbffe`.

The compact independent theorem on a single optimal-endpoint cycle \([1,6242]\) is a separate packet. It uses a particular 127-vertex path, a 6247-cycle from this split mechanism, and a checked local repair. That finite example does not supply the missing all-index input-pair hypotheses here.
