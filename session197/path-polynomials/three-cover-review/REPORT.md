# Independent review: three-cover cyclic gate and restricted inverse

**PASS, no correction required.** The positive bilinear cover kernel is exact for actual tournaments, the paired and unpaired inverse algorithms are complete within their named construction languages, and the explicit strong self-converse 16-vertex 7p witness checks independently. Neither inverse is a theorem about unrestricted spectrum coverage.

## Frozen source

Directory: `session197-odd-cycle-path-polynomials/`.

- THREE_COVER_CYCLIC_GATE.md: `0630c4ff7f9a52e69cbce4f18cdac8678609d4565a3e42f0afec8d588afe37cf`.
- verify_three_cover_cyclic_gate.py: `01274c4cc02fbf8afe11917c7fa224dd16710898ea7e3a2a3b526a81a9760475`.
- Both source receipts and both copied-source reruns: `106c9af91fe7ef7701c649cbcabc9c80a686a0dbfecf18f19e9e15bd8e54dca4`.
- three_cover_cyclic_witness.json: `6339e32f9e9d64222e13d55f93857627512ed16f51644158e843885be7d1aba1`.

The source helper verify.py (SHA `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`) and companion PROOF.md (SHA `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`) were copied with the packet. The checker was run in the copy because it writes the witness. Normal and `-O` output and the regenerated witness bytes match the frozen originals. Original sources were not changed.

The named S195 path-cover proof was read at its exact published commit ee07a5aac9cdfd8d3834fa834df4fb357fcaf769. Its Section4 indeed uses C3[A,Tr_m,Tr_m], defines ordered covers, and gives the m=2 row 11,17,4. The new kernel specializes to that row after the stated factorial conversion. This is correctly presented as a specialization, not a new general composition theorem.

## Counting and algebra audit

1. Cutting a Hamilton path into maximal module runs yields a linear word of colors following the directed three-cycle. All three modules are nonempty, and the two-vertex third module has either one or two nonempty runs. Therefore each of the other modules has at most three runs. For each count triple the possible color words are exactly the legal starts and lengths; there is no cyclic rotation quotient because this is a path.
2. The fragments of an unordered local cover are distinct actual paths on disjoint vertex sets. Assigning them to named run positions gives i!j!k! choices. Both one-path and two-path unordered covers of Tr2 have count1. The displayed two color-count matrices and the kernel [[5,6,0],[6,28,24],[0,24,72]] follow exactly. The inverse reconstruction uses intermodule edges at every join, so adjacent same-color fragments cannot be mistakenly merged.
3. The relations B(2)=2h+2c2 and B(3)=3h+6c2+6c3 yield the stated matrix in bin coordinates. A separate integer matrix calculation checks the complete transform, not only evaluations on examples. The binary coefficient is odd because B_U(2) is even and H(U) odd, while all other displayed terms are even. Thus it cannot vanish for a nonempty U. This does not rule out binary rays, paired inputs that vary together, or other cyclic quotients.
4. Taking U=T^op gives the displayed positive quadratic form. The explicit block swap and reversal of Tr2 is an anti-involution on the actual graph. The cyclic quotient is strong and every block nonempty, so the whole graph is strong without any strongness assumption on T.
5. Modulo8 the mixed terms vanish because u(h+u) is even for odd h. Modulo3, integer bin-polynomial coefficients and parity give B_T(2)=(-1)^n h, so u is0 for odd n and h for even n. This proves the stated paired residues. The possible prime or 7p exception involving the prime3 is correctly excluded by the n=1 value5 and the n>=2 lower bound45. The necessary residues 5 mod24 for a prime output and 11 mod24 for its 7p quotient are valid.

## Complete restricted inverse

Every partition of the vertex set into j nonempty parts supplies at least one unordered path cover, by Hamilton-path existence on each induced part. Different partitions cannot produce the same cover. Thus c_j>=S(n,j), with simultaneous equality in the transitive case.

The paired form has nonnegative coefficients and a strictly positive h-square term. It gives exactly the stated Phi_n, and any nontransitive input has h>=3, making the inequality strict. Stirling monotonicity and the positive u dependence make Phi_n strictly increasing for n>=1. The first eight values and the closed exponential expression are correct. Therefore the order cutoff is complete, and enumerating every labelled input at those finitely many orders gives a genuine all-order membership decision in the paired language. It may be extremely expensive; the source makes no stronger efficiency claim.

For the unpaired gate, the first kernel column gives H>=5h+6c2>=3·2^n−1 independently for each input. Both orders are therefore bounded, so the proposed two-input inverse is also finite and complete within that gate. There is no silent strongness, fixed-order or symmetry restriction in this bound.

Through target2000, Phi_4=6069 excludes all orders above three. Direct physical graph counts for every smaller labelled input give exactly5,45,509,621, matching the complete inverse output. Negative or too-small targets yield an empty search rather than an undefined order bound.

## Independent implementation

The standalone `check.cpp` imports no source code.

- It enumerates all 1,622 multiset color words in the relevant count ranges, instead of assuming the three-start description, and recovers the full factorial-weighted kernel.
- It counts path covers by literal vertex permutations and every admissible cut set, then divides by the exact fragment-order factorial. All 1,099 labelled tournaments through order five satisfy the order bounds, equality criterion, paired congruences and binary parity.
- It counts all 121 small physical two-input graphs using inclusion-exclusion over unrestricted walks, separate from the source endpoint subset DP, and checks the exact kernel.
- It checks the bin-coordinate transform by an integer matrix identity, the Stirling formula for Phi through order12, and the complete paired membership set through2000 from actual graph counts.
- It reconstructs the 16-vertex witness directly, checks every arc of the anti-involution and strong connectivity, and counts its Hamilton paths by walk inclusion-exclusion with signed 128-bit intermediates. The count is exactly37,050,461. Its literal input cover vector is (23,273,611), and its adjacency rows agree exactly with the frozen witness.
- A separate sieve generates all342 primes through2300 and tests the quotient5,292,923 against each. Together with the exact square-root inequality this independently certifies primality. The source's stronger all-integer trial-divisor remainder certificate is also reproduced unchanged.

Independent hashes:

- check.cpp: `a1a20b9d6c86a2a9843e62537e0a26fb7bfeb93afa96a3930912354955f236b5`.
- independent-normal.json and independent-debug.json: `574180b541b8c2b3c74bf2229cec7a6867abc073f935995e1c3cba86205ab872`.

Both optimized and assertion-enabled unoptimized independent builds agree. Compile as C++17 and run with no arguments. No external bank or file is required by the independent checker.

Conclusion: the exact physical kernel, sharp paired bound, complete restricted inverses, parity limitations and finite 7p witness all pass at their stated scope. No minimality or unbounded-primality claim is inferred.
