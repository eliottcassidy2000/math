# Session 174 synthesis: reversal, quadratic norms, and the missing counts

**Status:** linked packets contain proved statements and separately labelled finite computations. This note adds an independently audited quadratic-norm obstruction and closes an entire proposed constructor family. No global novelty claim, giant-integer primality proof, or resolution of missing Hamilton counts 7 and 21 is asserted.

## 1. The supplied numbers become different mathematical objects

- **240** is the ambient E8 root count. Balanced eight-vertex tournament scores give 70 of those roots; reflection completion gives **126 roots of E7**, not all E8. In even dimensions, this score-shell reflection group is finite exactly at 2,4,8 and has dense orthogonal closure in every other even dimension at least 6. See [lattices](../lattices/theorem.md) and [density proof](../lattices/density-extension.md). These are reflections of score vectors, not tournament automorphisms. Independently, decimal powers stabilize modulo 240; the two roles are not identified.
- **7417** is prime but not a palindrome: its reversal is 7147. Its decimal multiplicative order is 2472, and reciprocal trace coordinates have 1237 distinct values. It is also a genuine sieve prime for the giant-palindrome family: divisibility occurs exactly for m=1312 or 2395 modulo 2472. See [clock/count certificate](../counts7417/README.md) and [quadratic sieve](../palindromes/REPORT.md).
- **1001 and 2222** expose cyclotomic factors, while their decimal modular clocks have different periods and transient behavior. Fixed base and length matter; reversal of 240 is the padded word 042, not a three-digit integer. See [reflection algebra](../reflections/README.md).
- **Belphegor's number and the 1,888,529-digit expression** share digit reversal but sharply different complex-root geometry. The latter's digit polynomial has all roots simple on the unit circle; the Belphegor digit polynomial has exactly 2 roots there and 28 elsewhere. Neither claim is an integer primality test.

A discreet source correction: the supplied 2021 giant expression remains a celebrated example, but the checked [PrimePages leaderboard](https://t5k.org/top20/page.php?id=53) lists a [2,718,281-digit leader](https://t5k.org/primes/page.php?id=138383), rather than the [2021 example](https://t5k.org/primes/page.php?id=132851). These are external record claims, not primality certificates proved by this session. Banks–Hart–Sakata's theorem concerns asymptotic scarcity of primes among all palindromes; it does not settle prime infinitude in these thin families.

## 2. Three precise bridges, with their lost information

**Reversal → reciprocal quotient.** A reciprocal polynomial of even degree 2m is X^m Q(X+X^−1). Over the complex numbers this folds inverse root pairs; modulo 7417 it folds the cyclic subgroup <10> by inversion. The same identity gives both coordinate systems, but their fields and predicates differ. Retain the polynomial, base and degree: an evaluated integer alone does not retain digit reversal.

**Alternating evaluation → affine reflection holonomy.** For r_s(x)=s−x, a word of reflections composes to (−1)^length x plus its alternating label sum. Evaluating the label polynomial at −1 is exactly that translation coordinate. Consequently an even palindromic square-sum edge-label word cannot be a nontrivial simple path: it returns to its initial vertex. This preserves endpoint equations, not all internal vertices. The inherited Q_15 path has translation +1 and is therefore not palindromic in its edge labels. See [reflection packet](../reflections/README.md).

**A bit word → two separately observed objects.** The word b encodes both a decimal palindrome with digits 1 or 3 around a central 7 and a self-dual cyclic tournament built from two equal-profile seeds and their mirrors. The tournament's complete profile and chosen reversal-fixed Hamilton count forget b, although nonisomorphic tournaments retain it. Decimal residue automata can select a subfamily without changing those common tournament observables. This is an injective common encoding, not an equality of primality and Hamilton counts. See [spectrum construction](../spectrum/PROOF.md).

## 3. New synthesis theorem: a quadratic-norm obstruction to 7 and 21

Let A be any nonempty tournament, and construct T_A(1) from disjoint A, A* and one feedback vertex z with uniform cyclic orientation
A → A* → z → A.
Let F_A(u)=sum_k p_A(k)(u)_k, where p_A(k) counts unordered covers by k nonempty directed paths. Write a=F_A(1), b=F_A(2).

The arbitrary-A one-feedback skeleton formula in the [spectrum packet, §2–3](../spectrum/PROOF.md) gives

H(T_A(1)) = 3a²−2ab+b² = (b−a)²+2a².                 (1)

**Theorem.** For every odd prime q with −2 nonsquare modulo q, the q-adic valuation of H(T_A(1)) is even. In particular its 7-adic valuation is even, so **no one-feedback mirrored cyclic wrapper, regardless of seed, realizes 7 or 21**.

Proof: put x=b−a. If q divides x²+2a² and q does not divide a, then (x/a)²=−2 modulo q, a contradiction. Hence q divides a and then x. Dividing the norm by q² leaves another integral norm of the same form. Repetition shows that the valuation is twice min(v_q(x),v_q(a)); positivity of H ensures termination. Since −2=5 modulo 7 is not one of the squares 0,1,2,4, q=7 applies. Both targets have 7-adic valuation one. QED.

Equivalently these are splitting restrictions for the imaginary quadratic norm x²+2a², just as the palindrome identity
(2·10^(m+1)−1)²−41=40P_m
restricts prime factors through a real quadratic discriminant. **The shared bridge is the norm identity and its modular obstruction, not a map equating the two families.** The one-feedback norm is stronger than the selected seeds' divisibility by 9: changing those seeds cannot remove this 7-obstruction.

There is also a general residue check. The inherited tournament Hamilton-path oddness theorem gives a odd. From the falling-factorial definition, b=2p_A(1)+2p_A(2) is even. Thus x=b−a is odd, and (1) gives H(T_A(1))=3 modulo 8. This extends the packet's selected-family residue to every A. It does not apply to arbitrary tournaments outside this wrapper.

## 4. What was inherited, and what the session adds

Inherited mechanisms include tournament Hamilton-path oddness, THM-582/583 reversal-fixed/half-data constructions, THM-868's original balanced E8 shell, the previous session's equal-profile five-vertex seeds, standard reciprocal-polynomial identities and modular dynamic programming. Banks–Hart–Sakata is an external cited theorem, not a new result.

Session-specific proved developments include the finite/dense classification of reflected score shells; the dense-palindrome unit-circle theorem and contrasted Belphegor certificate; exact discriminant/orbit sieves and the explicit 7417 connection; the mirrored cyclic storage construction, bounded geometric recurrences and fixed-path specialization; and the all-seed quadratic-norm obstruction (1) above. “Session-specific” records derivation and packaging, not historical priority.

Keep three meanings of “spectrum” separate: complex roots of a digit polynomial, residues/periods of a modular sequence, and attainable integer Hamilton-path counts. A proof connecting two must state its actual map and preserved predicate.

## 5. Next proof targets for {7,21}

1. **Leave the mirrored cyclic wrapper class deliberately.** The corollary below closes arbitrary nonempty feedback tournaments, not only transitive feedback or the tested seeds. Changing the number or internal orientation of feedback vertices cannot supply 7 or 21. Identify the first changed adjacency rule, then derive a literal path-count formula or subset-DP certificate.
2. **Retain the missing profile information.** Before exploring a proposed nonuniform or mirror-breaking gadget, express its Hamilton count in realizable path-cover data. A formal numerical profile is not automatically the profile of a tournament. Use positive skeleton coefficients when proving lower bounds; signed evaluation kernels can hide cancellation.
3. **Audit the inverse problem at its actual target.** A proposed constructor needs a bijection or verified Hamilton count equal to 7 or 21, not equality of a quotient, spectral phase, palindrome residue, or reversal-fixed count. The latter is a different observable even when parity agrees.

### Corollary: every uniform cyclic wrapper A → A* → Z → A excludes 7 and 21

The universal two-feedback positive skeleton kernel is
H(T_A(2))=5q_1²+6q_1q_2+7q_2²+4q_2q_3+2q_3²,
where q_i=i!p_A(i). If |A|>=2, a Hamilton path exists and splitting it once gives q_1>=1 and q_2>=2; hence H(T_A(2))>=45.

Hamilton-path counts cannot decrease on adjoining vertices. Given a directed Hamilton path, any additional vertex can be inserted at an end or between consecutive vertices; choose a deterministic available insertion. Deleting the added vertex recovers the original path, making this map injective. Therefore for |A|>=2 and r>=2 the induced two-feedback wrapper gives H(T_A(r))>=45.

More generally take any nonempty feedback tournament Z, not necessarily transitive. If |A|>=2 and |Z|>=2, any two feedback vertices give the induced two-feedback wrapper, so the same >=45 bound applies. If |Z|=1, the norm theorem applies.

If |A|=1, name the two old vertices a,b. Paths with consecutive a,b are counted by two labelled path slots in Z; paths where a is last must start with b and have a Hamilton path of Z between them. Therefore
H=F_Z(2)+H(Z)=3H(Z)+2p_Z(2).
For every subset S of Z, both induced subtournaments on S and its complement have a Hamilton path (an empty slot counts as one). Thus F_Z(2)>=2^m for m=|Z|, and H>=2^m+1. For m>=5 this is at least 33.

The remaining small cases admit an elementary proof, not an assumed census. For m=1,2 the counts are 3,5. For m=3, the transitive and cyclic cases give 9,15. For m=4, transitive Z gives 17. Otherwise Z contains a cyclic triangle, so insertion monotonicity gives H(Z)>=3. If c>=1 is the number of cyclic triangles, two-path covers with block sizes 1+3 contribute 4+2c, while the three partitions of type 2+2 contribute 3. Hence p_Z(2)=7+2c>=9 and H>=27. None of the cases gives 7 or 21. QED.

For transitive Z_r the singleton formula specializes to 2^r+1: choose the subset appearing before consecutive a,b, plus the one path b,Z_r,a. The stronger arbitrary-Z proof shows that even changing the feedback interior cannot escape. A next constructor must break the mirror relation A*, uniform cyclic intermodule orientation, or the whole architecture.

The finite synthesis test checks the norm identity, inert-prime valuation behavior on bounded integer pairs, all 7-adic local solutions modulo 49, and the odd/even modulo-8 consequence. These validate hostiles and arithmetic implementation; the universal theorem is proved above. Independent literal-wrapper audit of the norm theorem: [audit script](../spectrum/audit_norm7.py), [certificate](../spectrum/norm7_certificate.json). The full synthesis and wrapper corollary have an [independent proof audit](AUDIT.md), [wrapper certificate](wrappers_certificate.json), and a second [literal-feedback audit](feedback.json) using direct permutations with no imported author code. Normal and optimized runs agree.
