# Independent audit: sharp symmetry-slack pole orders

Date: 2026-10-21. Verdict: **PASS in the stated scope**, with separate finite-exact controls. This is a mathematical proof audit, not formal verification, a complete compiler run, or a historical-priority assessment. No Git operation was performed.

## Reviewed statement and immutable inputs

The reviewed candidate is `session23-sharp-slack/SHARP_POLE_ORDER.md`, SHA-256 `1e876eceba203e77031644a6e11c4f31cfebfa6bca649bc6e7eb7a34f3c4326b`.

Its canonical-word compiler dependency is `session23-connections/SYMMETRY_SLACK_RATIONALITY.md`, SHA-256 `22cb6aca899d9a485b2432d37a90d13040dfc6353ee41529a243d8489aa974a1`. That packet was independently reviewed before this follow-on; this audit read its required arguments rather than silently substituting a marked or noncanonical grammar. The actual-orbit charge and exceptional-set proof were checked directly in `math/05-knowledge/results/continuation_width_symmetry_20261020/README.md`, Theorems 3–5.

For every nonempty exact symmetry-slack layer, write uniquely

    d = 3a + 4b,  a >= 0,  b in {0,1,2},
    m = 2a + b + 1,  k = ceil(m/2),
    p(z) = 1 - z - z^3,  q(z) = p(z^2).

The exact reduced ordinary, self-converse, and converse-merged denominators are respectively

    p^m,  q^k,  p^m q^k.

The only empty layers are 1, 2, and 5. The associated minimal eventual homogeneous recurrence polynomials, their exact orders, the ordinary one-dominant-pole asymptotic, and the explicit sufficient recurrence onsets in section 9 follow as stated in the candidate. No mathematical repair to the completed candidate was required.

## 1. Actual-orbit charge accounting

An actual automorphism orbit induces a vertex-transitive tournament and is therefore regular of odd order. In any restricted linear order, the usual final-vertices argument gives internal feedback at least (ell^2-1)/8. This bound remains valid whether the orbit is homogeneous or participates in synchronized actions.

Thus D = sum[3(ell^2-1)/8 - ell] over nontrivial actual orbits is a nonnegative integer, D <= d, and d-D is divisible by 3. The external backward-arc count e in any optimal order satisfies e <= (d-D)/3. The charge is zero for triangle orbits and at least four for every larger orbit. This excludes 1, 2, and 5 exactly as claimed.

The exceptional-set construction includes all R large orbits and every actual orbit touched by an external backward arc. The count of exceptional ORBITS is at most R+2e, including singleton orbits. Replacing this count with a vertex count would lose the sharp result. Replacing actual orbits with score classes would be invalid. Neither replacement occurs in the candidate.

## 2. Rebuilding the exceptional set in the right order

The canonical compiler supplies a weighted quotient with cost C <= d. Its representative-sampling straightening proof gives a genuine optimal order obtained by placing its W blocks consecutively and each block in natural atom order. This choice matters: the proof does not assume that an arbitrary inherited exceptional set is disjoint from the long canonical blocks.

A block of size greater than d cannot touch a backward quotient edge. The exact decorated quotient stabilizer fixes every such block, since the induced long-slot quotient is transitive. Its actual internal orbits are singleton letters and independent triangle atoms. Hence rebuilding B in this particular order leaves every long block disjoint from B.

For an exceptional actual orbit O and a long block M, homogeneity makes the direction from each vertex of O to M well-defined. Transitivity on O and setwise preservation of M make that direction constant on O. Because every cross edge incident with M is forward in the chosen order, O lies entirely to one side of M. This proves the no-straddling assertion without assuming O itself is a module.

The inherited exceptional grammar says each gap between consecutive B vertices is a homogeneous W module. Two distinct long canonical blocks with no B vertex between them would therefore be contained in, or overlap, one common W module and merge by the maximal-W union lemma. Consequently successive long blocks require B separators. Separators in distinct gaps must belong to distinct actual B orbits, because otherwise that orbit would straddle an intervening long block. Therefore

    L - 1 <= number of B orbits <= R + 2e,
    L <= R + 2(d-D)/3 + 1.

This is the genuinely new sharp step. It applies to every actual filling of every valid configuration in the finite compiler, so it bounds each rational summand, not merely a selected displayed family.

## 3. Optimizing the charge and proving attainment

Writing P=2D-3R, each large orbit contributes at least five. If d is 0 modulo 3, P>=0. If d is 1 modulo 3, D cannot vanish and P>=5. If d is 2 modulo 3, either R>=2 and P>=10, or a single orbit must have charge 2 modulo 3; it cannot have order five and hence contributes at least 19 to P. These bounds give exactly the stated residue formula for m.

The ordinary lower family is injective. The slack-three strong component S(W) has two intrinsic fixed endpoints for |W|>=2, and so is never isomorphic to R5, whose action is transitive. Ordered SCC decomposition recovers the sequence of S components, R5 components, and outer singleton/triangle words. Repeated isomorphic components do not introduce a permutation quotient: their SCC order is intrinsic. Feedback and actual moved support add across this order. Its generating function has a positive pole of order m at the positive root of p, establishing sharpness and existence in all remaining layers.

The nested wrapper construction gives a second attainment proof, needed for the self-converse channel. For every tournament U of order N>=2, x in Wrap(U) is the source of its maximum-score class, which consists of x and at most the unique source of U. Likewise y is the sink of the minimum-score class. These two classes are distinct, including N=2. Thus both endpoints are fixed by every automorphism and are intrinsically recoverable. Every automorphism of U extends, so the FULL actual group, not just a visible subgroup, is retained. The weighted triangle has minimum cross cost min(N,N,1)=1. Hence wrapping adds one feedback and three slack while preserving moved support.

Wrap(U) is strong, and converse exchanges the endpoint roles. At each recursion depth, its unique strong component outside W identifies the wrapper, and then its endpoints recover the preceding level. At the base, the ordered R5 components recover every W gap. This proves literal injectivity of all 2a+b+1 words; there is no dependence on a presumed maximal-W quotient of the construction.

Self-converseness is exactly reversal pairing of the outside words at each level and of the base gap list, with a palindrome in the central base gap when present. Since R5 is self-converse, these conditions are both necessary and sufficient. The exact lower-family series in equation (13) therefore has a positive pole at sqrt(rho) of order k. This also works at d=0 after removing only the finite excluded head. In particular the extension does not make the false inference that an arbitrary ordered list of distinct exceptional strong components must be self-converse.

## 4. Exact denominators, not only upper bounds

The scalar compiler puts F_d over p^m. Its positive subfamily forces a pole of full order m at rho. The cubic p has no rational root and is irreducible over Q; this forces the entire reduced denominator to be p^m.

The self-converse compiler puts C_d over q^k. Its positive subfamily forces a pole of order k at sqrt(rho). The candidate's short irreducibility argument for q is valid: modulo 2, the monic polynomial -q reduces to the square of the irreducible cubic z^3+z+1. Any proper monic integer factor must consequently have degree three. Normalized reflection must exchange the two cubic factors, since a reflection-fixed monic cubic would be odd and have zero constant term. Writing the pair as in the candidate, reduction modulo 2 makes A even and B odd, while the zero z^4 coefficient gives 2B=A^2 and makes B even, a contradiction.

As an independent arithmetic cross-check, the checker also proves irreducibility using the incompatible finite-field factor-degree patterns 2+4 modulo 3 and 3+3 modulo 5. It verifies each small factor irreducible by explicit polynomial division, not by trusting a symbolic factor command.

The stated common-root argument proves gcd(p,q)=1. Thus in (F_d+C_d)/2, the function from the other channel is analytic at every root of a given denominator factor; no pole can cancel. The merged denominator is exactly p^m q^k. Reversal of each denominator gives the recurrence polynomials claimed in the candidate. These are minimal eventual recurrences, with no assertion about earliest onset.

The ordinary asymptotic is justified by positivity and the unique smallest-modulus root of 1-z-z^3. The candidate correctly avoids a single all-parity self-converse asymptotic: q has both positive and negative dominant poles.

## 5. Independent finite-exact controls

`check_sharp_slack.py` is a separate implementation using bitset tournaments, degree-constrained backtracking for full automorphism/anti-isomorphism sets, subset dynamic programming for feedback, explicit homogeneous-subset tests, and SCC decomposition. It does not import the producer's checker or the compiler checker.

Normal and optimized Python runs both **PASS 17,182 explicit checks**, with byte-identical JSON results:

- 1,205 actual tournaments for the rebuilt-order and orbit-charge upper bound: all 1,099 labelled tournaments of orders one through five, 35 deterministically selected tournaments at each order six, seven, and eight, and the synchronized two-triangle hostile. These are not claimed to be 1,205 distinct isomorphism classes or an exhaustive order-eight census.
- 1,098 arbitrary-U wrapper controls: every labelled U of orders two through five. Checks include direct exact feedback, the complete actual automorphism group and its restriction, the complete converse criterion and exchanged endpoints, and strong connectivity. Here U is not restricted to W.
- Sixteen explicit long-word families, with a,b drawn from (1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(0,2),(1,2), each in a reflection-symmetric and a broken-reflection variant. The largest has 72 actual vertices and quotient order 17. Full actual automorphisms and anti-isomorphisms are computed on the expanded tournaments. Exact feedback on these large examples uses independent weighted subset DP on the homogeneous quotient plus internal word feedback, not exponential DP on all 72 vertices.
- On those families, the checker separately verifies canonical maximal-W blocks by testing every quotient module, rebuilds actual B from an optimal block order, checks every B orbit against every long block, finds 54 actual separating gaps, and verifies equality in the sharp orbit-count bound. These checks support but are not needed by the cleaner wrapper injectivity proof.
- The relaxed orbit-charge optimization is exhaustive over all charge multisets of total at most 300. It gives precisely the claimed maximum and the three empty layers. This is a finite arithmetic control, not an all-order proof or a realizability assertion for arbitrary charge multisets.
- Two finite-field factorization/irreducibility certificates independently support the algebraic denominator conclusion.

Run:

    python3 session23-sharp-slack-review/check_sharp_slack.py
    python3 -O session23-sharp-slack-review/check_sharp_slack.py

The saved outputs are `receipt.json` and `receipt-optimized.json`.

## 6. Explicit sufficient recurrence onset

The final section 9 was added after the initial sharp-pole review and was independently checked before this final source pin. Its degree calculation passes exactly as written. Let H=d(4d+1), let A be the total short-block size, and let ell be the long-slot count. The finite compiler gives A <= d(q_core-ell) and q_core <= 4d+1. The ordinary long-tail numerator has degree at most d+3. Clearing p^m therefore gives degree at most A+d ell+3m <= H+3m.

For self-converse configurations ell=2r+j, the paired tail numerator has degree at most 2d+6 and the palindromic tail numerator at most d+6. Clearing q^k gives A+r(2d+6)+j(d+6)+6(k-r-j)=A+d ell+6k <= H+6k. The exponents cleared are nonnegative. Bringing the two channels to their coprime common denominator gives merged numerator degree at most H+3m+6k.

Since a coefficient recurrence is homogeneous at every index beyond the numerator degree, these are valid sufficient onsets, with indexing from coefficient zero:

    ordinary: n >= H+3m+1,
    self-converse: n >= H+6k+1,
    converse-merged: n >= H+3m+6k+1.

The bounds exceed the required lag lengths, so no negative-index convention is hidden. For d=0, the directly stated safe starts 3, 6, and 9 are valid; the merged numerator is indeed 1-z^2-z^4-z^6. These onsets are sufficient rather than minimal. Obtaining the required prefix still requires computation; no prefix coefficients or polynomial-time general enumeration procedure are supplied by this degree argument. The 17,182 finite checks above are not advertised as a substitute for the onset proof.

## 7. Scope and remaining open work

The audit establishes the new actual-orbit separator bound, exact pole orders, exact reduced scalar denominators in all three converse conventions, the resulting minimal eventual recurrence orders, and the explicit sufficient onset bounds. It does NOT compute the compiler's finite configuration list or its numerator in uncompiled layers. The leading ordinary constant remains unevaluated. The construction does not establish minimum vertex order for a slack value. No practical complexity bound, complete slack-six numerator, earliest recurrence onset, Hamiltonian-path identification, black-degree classification, or metagraph connectivity theorem follows here.

The frozen compiler packet remained unchanged throughout this audit. All review files are outside the repository, under `session23-sharp-slack-review/`; no repository, remote, or Git state was modified.
