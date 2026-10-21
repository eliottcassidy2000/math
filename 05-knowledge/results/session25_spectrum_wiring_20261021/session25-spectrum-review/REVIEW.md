# Independent Session25 modular Hamiltonian-spectrum review

Status: APPROVED. Frozen producer proof and verifier independently audited; no mathematical revision required.
Reviewer imports no producer code. Repository and Git are read-only in this task; all artifacts are in this review directory.

## Independently recovered proof

Let T(w) have old vertices 1,...,n with i→j for i<j, plus x with x→i exactly when b_i=1. Every Hamiltonian path, split at x, consists of an increasing sequence on the left and an increasing sequence on the right. Conversely each assignment of old vertices to these two sequences is one candidate path, legal precisely when the last left bit is 0 and the first right bit is 1 (a condition omitted for an empty sequence).

If x is first there is one path when b_1=1; if last there is one when b_n=0. For nonempty sides let j be last on the left and i first on the right. If j<i, all vertices can be covered only when i=j+1, giving one path for each adjacent 01. If i<j, the endpoints are fixed, all vertices below i are left and all above j right, and the j-i-1 interior vertices choose either side freely. Thus their contribution is b_i(1-b_j)2^(j-i-1). These cases partition all paths, not just insertions into the unique old path. In particular the weighted inversion sum retains all orphan paths.

Using 1+#10 for the non-orphan contribution, initialize H=1,A=0,p=0 and append bits via

- 1: (H,A,p)→(H,2A+1,1);
- 0: (H,A,p)→(H+A+p,2A,0).

Here A is the binary value of the prefix. The p term accounts for the new adjacent10; A is exactly the sum of weighted inversions ending at the new0. This derives the recurrence for arbitrary finite words, including the empty word.

For L≥1, q=2^L, and a completed block1^L0^L, the block map is

A' = q² A + q(q-1),
H' = H + q(q-1) A + (q-1)² + 1.

This follows by geometric sums over the L ones and L zeros; p=0 at every block boundary. Set M=2^a m with m odd. A positive L≥a with 2^L≡1 mod m exists because2 is a unit modulo m (for m=1 the constraint is vacuous). Thus q≡0 mod2^a and q≡1 mod m. Consequently q(q-1)≡0 modM. Starting A=0, all block-boundary A vanish modM, and H_t≡1+tD modM, D=(q-1)²+1.

Modulo m, D≡1. Modulo2^a, D≡2, so gcd(D,M)=gcd(2,M)=g. Therefore t=1,...,M/g gives exactly all r with r≡1 modg. Positive t is essential for feedback vertex number exactly1; the period endpoint replaces t=0. In particular M=1 and M=2 work with L=t=1 and T(10)=C3.

For positive L,t, the word starts1 and ends0, and n≥2. There are directed routes x→1→i and i→n→x, with repeated endpoints omitted. Every old vertex communicates with x, hence the tournament is strong. Deleting x leaves a transitive tournament, while x→1→n→x is a triangle. Feedback vertex number is exactly1.

The integer closed form is

H_t = 1 + [q²(q^(2t)-1)+2t(q+1)]/(q+1)².

This is an integer equality, not a license to invert q+1 modulo M. It follows by summing the exact block map, and our controls include t=0. The closed form is not needed for modular coverage.

## Scope and hostile boundaries

- L=0 or t=0 collapses the word to empty, leaving a singleton with feedback vertex number0. Strong singleton conventions do not repair that failure.
- Untuned dyadic example: M=8,L=1,t=2 gives H(1010)=9≡1, whereas1+tD≡5.
- Untuned odd example: M=3,L=1,t=2 gives H=9≡0, whereas1+tD≡2. In fact L1 misses residue2 mod3.
- H≡r modM is not H=r. For instance L3,t3 gives H=207127≡7 mod8, not the permanently omitted exact count7.
- Rédei parity gives necessity of the even-modulus restriction; odd moduli impose no parity restriction on representatives because adding the odd modulus changes parity.
- Fixed L gives exponential exact-count growth, not an exact arithmetic progression. The affine progression is only modulo the tuned M. Varying t by the modular period gives infinitely many distinct counts in each covered residue, but proves no density or bounded-gap statement in the integer spectrum.
- The omitted exact counts7 and21, and the global completeness question, are unchanged.

## Canon overlap / novelty boundary

Targeted read-only searches covered the theorem directory by statement synonyms, insertion words, feedback vertex number, modular Hamiltonian spectrum, binary words, and weighted inversions. No prior universal modular-coverage statement was found. This is only a repository-overlap assessment, not a literature-priority claim.

- THM-012b-insertion-decomposition.md and THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md already give the full insertion/orphan decomposition. The word formula specializes that mechanism; it should not be advertised as a new general insertion identity.
- THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md gives general exact boundary response and a finite order-nine interval. Its endpoint/cut data recover the word identity on a transitive base, not the tuned block arithmetic automatically.
- THM-4102-selected-order-ten-strong-ear-solid-interval.md and THM-4104-selected-order-eleven-strong-ear-solid-interval.md are selected finite exact-count interval results. They neither imply universal modular coverage nor are strengthened in their order/prefix claims by the present unbounded-order congruence result.
- THM-4094's two infinite prime lanes concern exact values. Universal modular coverage does not supply a single missing exact prime atom or advance the80407/80689 frontier from THM-4104.
- THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary.md reinforces retention of compositional sidecars. Here (A,p), not H alone, is the necessary sufficient state for this word recurrence.
- THM-1900-the-insertion-response-calculus.md handles general insertion and neutral downsets. THM-1920-the-spectral-insertion-response.md concerns skew characteristic polynomials, a different spectrum.
- THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md concerns numbers of labelled strong tournaments by order, not Hamiltonian-path counts of individual tournaments.
- MISTAKE-499/502 and the October20 completeness correction are respected: orphan paths are included, empty-word conventions explicit, and exact completeness remains open.

## Independent computational audit

Reproduction:

python3 -B session25-spectrum-review/independent_check.py
python3 -B -O session25-spectrum-review/independent_check.py

The normal and optimized outputs are byte-identical. No assert statements gate correctness.

-511 words of lengths0 through8: independently built adjacency and subset-path DP, split-at-apex assignment enumeration, pair formula and recurrence all agree.
-63 words of lengths0 through5: literal permutation enumeration agrees.
-Strongness is checked in both directions for every nonempty word through length8, using graph searches rather than the producer's transitive-closure implementation.
-112 exact closed-form controls for L1..7,t0..15.
-All1024 moduli1..1024:393472 block-state and residue controls;113 direct modular word checks.
-Adversarial boundary controls described above.

These finite checks test implementations and exhibit failure modes. The all-moduli theorem is proved by the split-path bijection, geometric block recurrence, CRT reset, and gcd argument, not inferred from the controls.


## Frozen final-packet decision

APPROVE the main theorem and its explicitly scoped corollaries. The final packet was read in full. The producer confirmed this freeze; both verifier runs reproduce receipt.json byte for byte. No repository files or Git state were modified by this audit.

The additional sparse-family corollary is sound: H_(L,t) >= 4^(Lt)/3 implies Lt <= floor(log_4(3X)), and the harmonic divisor-pair bound counts at most O(log X log log X) distinct family values. This is an upper bound on the two-parameter subfamily, not on all strong tournament counts. The X >= 3 condition avoids the empty-index logarithm boundary. Independently tested 300 lower-bound rows and 30 rows of the exact L1 congruence H_(1,t) = 1+2t² mod3. Infinitely many distinct counts follows from positive block increments; simultaneous prescribed finite valuations at finitely many odd primes follows by CRT using residues p^e mod p^(e+1). None of these claims requires a new external theorem beyond standard CRT and Redei parity.

Required scope safeguards are already present in the frozen packet: L >= max(1,v2(M)), positive t, explicit empty word, integer rather than modular division in the closed form, actual Hamiltonian paths rather than orbit-normalized counts, no exact prime attainment, and no literature-priority claim.

Frozen producer SHA256:

- PROOF.md: 4c32eb0fb63f4a89511eef58b4ccd659ae28db895ab7d1514f66b20dc0ca48a8
- verify.py: dd6c70a46a4ebd671ecd0e85285196571e4e1fd10da35d4903c56c76ddd30443
- receipt.json: 53b649fec4e19bc47324510595fd73424e423acba6d0e92f57c16a8030be9f31

Frozen reviewer SHA256:

- independent_check.py: e6b36bd69a9bf92f0f19dd7318d0f50c869ab59030cef0d908b54c4d71333530
- independent_check.out: 5ff9578cf2118a776bb09133c7c6b5d51d0555308f757338eca997bb3c9ee099

The reviewer script imports no producer module. The producer normal/-O replay outputs both hash to the producer receipt digest above; reviewer normal/-O outputs both hash to the reviewer output digest above. SHA256SUMS records this review file's own final digest separately to avoid a self-reference.
