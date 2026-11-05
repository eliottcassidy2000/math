# Independent audit: modular gadget universality

**PASS, no correction required.** The joint residue construction is complete for every even modulus. The dyadic construction has the claimed quadratic word-length bound, and the strong singleton-prefix graphs give infinitely many counts in each parity-compatible residue class. These are congruence results, not exact-target, cofinite-spectrum or prime-generation theorems.

## Frozen source

Source directory: `session197-prime-atoms/modular-gadget-universality/`.

- PROOF.md: `74df1989d75ff743177937c108ac36dce133412e48701e4d412297d8c47b7d63`.
- verify.py: `2b465ffd43d8089cc6de7b7ad41580de1451eadff4d444f47f8d9f0a2dbb8f1a`.
- Both source receipts and copied normal/optimized replays: `30b755a0eb7e14eb39454af32b07a55633542401bd66d4d29f3ae16b045349ab`.

The source was copied before replay, although its verifier does not write files. Originals were not modified. The physical integer gadget and singleton-prefix seam are inherited from the already reviewed ROOTED and ARITHMETIC_INVERSE packets, not new assumptions about realizability.

## General modular loop

For M=2^k q, k>=1 and q odd, choose a>=2 with2^a=1 modulo q. From any state n=0 modulo M, appending1^a leaves f unchanged and makes the current n congruent to2^a-1. Subsequent zero bits add the geometric prefix sum, with one parity correction at the first zero. Thus the exact modular increment is

    Delta=(2^a-1)(2^k-1)+1,

and the new n is again zero modulo M. The proof uses only the state modulo even M, which retains the needed parity information.

Modulo q, Delta is1. If k>=2, Delta is2 modulo4, so its2-adic valuation is exactly one. If k=1, Delta is even and still coprime to q. Therefore gcd(Delta,M)=2 in every case, including M=2 and q=1. This is sufficient to reach every odd f residue from1.

A fixed-width suffix for r=u-1 gives exactly the same additive f shift as starting that suffix from the state(0,1), namely f(r)-1. The dependence on the previous n vanishes modulo M. Since f(r) and z are odd, their difference is divisible by2, and the reduced congruence has a solution modulo M/2. At M=2 that modulus is1 and imposes no condition; the construction handles it explicitly.

Adding multiples of M/2 loop repetitions preserves both residues. The resulting words start with1 and grow in length, so they represent arbitrarily large distinct integers, not merely common-zero padding of one integer. The stated O(M²) word-length bound for the first effective witness follows from a<=M, j<=M-1 and the fixed suffix width. It does not require an efficient factorization algorithm.

## Dyadic short construction

The neutral word10^k returns n to0 modulo2^k and changes f by2^k, so it preserves the modular state while ensuring a nonempty leading1. This remains valid from any loop-return state, not only the initial integer zero.

The a-loop increment is2-2^a=-2(2^(a-1)-1) modulo2^k. Writing the desired correction j as a sum of powers of two, bit i=0 uses one a=2 loop; every higher bit uses a=i+1 and a=2, whose normalized corrections add to2^i. For k>=2 there are at most2(k-1) loops, each at most2k bits. The initial neutral word and final k-bit suffix fit inside the displayed generous bound4k²+2k+1. For k=1, j=0 and no adjustment loop is needed. Extra neutral prefixes give arbitrary largeness separately from the short-witness bound.

## Actual tournament conclusions

The gadget G(L) itself can be weak. The strong conclusion applies to T(L), obtained by adding the singleton-prefix vertex a and reversing only b→a. Deleting b leaves a transitive tournament. In every Hamilton cycle, a's unique incoming arc is b→a; afterward all the remaining vertices must occur in their unique transitive order. The last base vertex points to b because its highest incidence bit is1. Thus for L>1 the cycle exists and is unique, proving strongness and feedback vertex number exactly1.

For u=1 modulo even M, L is odd and H(T(L))=g(L)+L-1 has the prescribed odd residue. The constructed L can be arbitrarily large, and H>=L makes the set of attained counts unbounded, hence infinite. For an arbitrary positive modulus m, use M=2m and the appropriate odd lift; this also includes m=1.

For a fixed physical prefix, h is odd and R=e modulo2. Hence h is invertible modulo2^k, and the required gadget residue computed from an odd target residue is itself odd. Joint modular universality then supplies the physical gadget for any chosen L residue. The statement correctly restricts this fixed-prefix conclusion to dyadic moduli; odd content can obstruct other moduli.

Meeting every allowed congruence class does not prove exact-value coverage, cofinite coverage, or infinitely many prime values. The proof consistently keeps those conclusions separate.

## Independent controls

The standalone `check.py` imports no source code. Its general construction chooses the smallest eligible a>=2 by modular exponentiation, rather than the source's2phi(q), and finds the correction by exhaustive residue search rather than a modular inverse. Gadget counts are independently evaluated by summing the contributions of zero positions in the actual binary word. The physical controls count literal directed Hamilton words and anchored cycles, not subset endpoint DP.

Independent normal and optimized runs agree on:

- 22,880 general joint residue pairs for every even M through64;
- 2,640 arbitrary loop-state translations and3,168 fixed-suffix translations;
- 10,938 dyadic joint pairs, including large modulus bit lengths16,32,64,128;
- 10,922 neutral-prefix arbitrary-largeness controls;
- 43,688 fixed-prefix residue pairs, with all chosen index residues in the small ranges;
- 131 actual gadget/singleton-prefix graph pairs, with direct H and unique-cycle counts and transitive root-deletion checks.

These finite controls support the implementation. The all-order statements rest on the proved loop gcd, suffix translation, correction decomposition and physical cycle argument.

Hashes:

- check.py: `8f5b8a38b3a03cb47214403df2deb5d63e4b80b802a3cac2b4236e236f40efd1`.
- independent-normal.json and independent-optimized.json: `cf3b4407e1dd3bcf208498fc15a4596730ed0e6d6ec18ca5c0c9bf4dbcfc9cb5`.

Run `python3 check.py` and `python3 -O check.py`. Neither writes files or needs external data. The exact source copies and both source receipts remain alongside the independent checker.
