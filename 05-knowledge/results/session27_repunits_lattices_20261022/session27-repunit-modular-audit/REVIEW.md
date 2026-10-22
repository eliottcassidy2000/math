# Session27 independent audit: repunits and modular bridges

## Verdict

PASS for the mathematical mechanisms and finite checks, subject to the exact scope below. Two prose clarifications requested during review: replace '+1 counterparts' with '2^k+1 counterparts'; replace 'integral coset' with 'half-integral coset'. One formal theorem-domain clarification requested: the degree assertion requires attainable h, while unattainable h has identically zero inventory. Final source hashes are in manifest.json, taken after producer freeze confirmations.

No repository edits or Git writes were performed. Review output is a separate packet. A bounded remote-search addendum being prepared by the catalogue producer is NOT covered unless separately named in the manifest.

## Independent verification versus replay

1. Replayed verify_weighted.py, check.py, and check_camion.py. All output JSONs are byte-identical to the supplied result files. The weighted replay covers 33,867 labelled tournaments, orders 1–6, and fans k=4–12. Modular replay checks through q^60; the Camion replay directly computes residues at 21 and 1729.
2. Independently written independent_checks.py does not import those implementations. It counts Hamiltonian paths by permutations, not subset DP, for all 1,099 labelled tournaments of orders 1–5 and reconstructs weighted counts by strong-component convolution using exact fractions.
3. It reads the inherited 6,880-row order-eight census, isolates its exactly two H=63 rows, and independently recomputes each witness's H=63, trivial automorphism group and strong connectivity by permutation enumeration and reachability. This verifies witnesses and consistency of the inherited inventory, NOT a fresh independent exhaustive classification of all order-eight tournaments. The completeness input remains the cited, previously audited Session24 augmentation census and THM-344.
4. Independently enumerates E8 roots in doubled coordinates: 240 roots, including 70 balanced half-integral roots. Independently checks 24-coloured partition coefficients using the logarithmic-derivative recurrence rather than reciprocal products. Independently checks the 21 Camion hostile using exact strong counts and math.comb.
5. Validates all 7,139 catalogue entries carrying content hashes against present source bytes and every recorded match line/700-character preview; none changed. Checks totals: 113,808 = 72,894 text + 40,488 absent + 424 nontext/read failures + 2 oversized; 7,290 entries, 12,379 matching lines. Does not pretend that absent historical blobs have been recovered.

## SCC theorem and measures

Every tournament's strong condensation is a unique total order. A Hamiltonian path must visit components in that order and is precisely a concatenation of component Hamiltonian paths. Deleting singleton components leaves an ordered word A_1,...,A_r. There are binom(n-M+r,r) placements of n-M singleton components in r+1 gaps, for n>=M. Otherwise the combinatorial count is zero.

An automorphism cannot permute condensation positions: a finite total order has no nonidentity automorphism. Therefore |Aut(T)| is the product of component automorphism orders, even when two components happen to be isomorphic. No r! or permutation-of-equal-components correction belongs in this formula. Orbit-stabilizer gives L_h(n)=n! W_h(n). Multiplying W by a fixed h preserves degree. Uniform-class and uniform-labelled coefficient weights differ, but degree comparisons and their limit direction agree.

The sum is finite without any global spectrum conjecture: Camion's Hamiltonian cycle already yields H(A)>=|A|, hence finitely many strong atoms at bounded H. Busch's sharper bound H(A)>=5^((|A|-1)/3) gives order<=1+3 log_5 h. The report's stronger source-dependent use is justified by Session24 PROOF.md and the publisher abstract linked below. Positivity prevents cancellation of the maximum-degree term. For attained h, degree is the maximum strong-factor length; for h=1 use the empty word and degree zero; unattained h has zero inventory rather than a defined positive leading coefficient.

N and W are fixed-h, growing-order quantities. They are not the all-orders strong-atom count S(h), a large-h density, a converse-merged count, or a fixed-order rarity theorem. The report preserves these distinctions.

## Rarity applications and hostile

31 and 127 are attained primes, so degree one. All proper nontrivial multiplicative decompositions of 63 require 7 or 21, which are impossible; 63 is attained, so degree one. Strong 3,5,11,13,43 are supplied by finite inventories and constructions. Thus 33=3·11,65=5·13,129=3·43 have degree two, and each displayed smaller-repunit/larger-neighbor ratio tends to zero for both measures.

15=3·5 has degree two, while 17 has degree one. Likewise 255=3·5·17 has degree three, while 257 has degree one. The primes occur as strong values, for example via the one-reversed-extreme-arc construction. So the 255/257 ratio diverges, a genuine counterexample to universal repunit rarity. These are 2^k−1 versus 2^k+1, a numerical difference of TWO. No assertion that every other odd integer occurs is needed.

For H=63, Busch excludes strong order>=9 because 5^8>63^3; the inherited complete census below order nine supplies only two order-eight atoms, each asymmetric. Singleton insertion therefore gives N_63=W_63=2(n−7), n>=8, and L_63=2(n−7)n!. The order-eight numerator 80,640 is correct. This is not justified from the two witnesses alone: completeness depends on the inherited census and bound, both explicitly identified.

## Modular obstruction

The proposition is correct for a nonzero rational R(q), holomorphic in the unit disk, assumed to give a positive-integral-weight holomorphic SL2(Z) modular form under q=e^(2πiτ).

If constant term is nonzero, approaching each finite rational cusp along a modular transform of iY forces unbounded growth, hence a rational pole at infinitely many distinct roots of unity. If constant term vanishes, the S-transformation makes the behavior along τ=it exponentially flat as t→0, incompatible with the finite Laurent order of a nonzero rational function at q=1. The nonzero, positive-weight and holomorphic conditions are essential to the scope. This does not rule out transformed series, mock forms or arbitrary proposed modular interpretations.

The fixed-H SCC series is a finite positive weighted sum of q^M/(1−q)^(r+1), and is therefore rational. The labelled exponential generating series is the W-series, not the labelled ordinary generating series. F_3=q^3/(1−q)^2 and F_63=2q^8/(1−q)^2 use ordinary class counts. Weighted F_3 instead has coefficient 1/3 because Aut(C3)=3. The report does not conflate them.

## Normalization and lattice typing

With q=e^(2πiτ), theta_L=sum q^(norm/2); eta=q^(1/24) product(1−q^m); Delta=eta^24. The weight-four lattice theta_E8 is E4, not the Eisenstein series E8. The coefficient 1728 in E4^3−E6^2 follows directly from 720−(−1008), independent of 1729. The character theta_E8/eta^8 has shifted coefficients 1,248,4124; 248=240+8, rather than an alleged equality of lattice and Lie dimensions.

For rank 24, theta_N=E4^3+(R_N−720)Delta, so the Leech character is j−720, not j−744. Ordinary binary Golay Construction A necessarily contains 48 coordinate roots and gives E4^3−672Delta; it is the A1^24 Niemeier lattice, not Leech. The cited THM-489 uses a q^2 convention, while the packet consistently normalizes its displayed identities to q=e^(2πiτ).

For n divisible by eight, the tournament score deviation is half-integral with coordinate sum zero. Subtracting (1/2,...,1/2) leaves an integer vector with even sum −n/2, proving membership in D_n^+. Score norm and coset membership are retained, not orientations or H. The E8 score floor is the 70 balanced half-roots; the rank-24 half-coset norm floor is 6, not the Leech minimal norm 4. This distinguishes the particular map from hypothetical unrelated embeddings.

II25,1 has signature (25,1) and rank26, not a positive packing dimension. The primitive-null quotient is rank24 and positive; primitive-null isometry orbits correspond to the 24 Niemeier classes. The 23 deep-hole types correspond to rootful Niemeier lattices. Distinguish 26 critical bosonic spacetime coordinates from 24 transverse oscillators. The stated packing uniqueness is correctly limited to periodic optimal packings. Genus-one theta equality is not an isometry theorem and does not settle higher-genus/Jacobian questions.

The nine Heegner radicands correctly INCLUDE 1, corresponding to discriminant −4. Stark's completeness theorem is explicitly imported, not re-proved or freshly full-text-audited in this packet.

## Fresh 1729 bridge and historical repair

Luca's August2026 Theorem1, accessed directly at the publisher, says the only prime cube in the odd-part Carmichael set K is 27. Thus a Carmichael integer of form 2^n p^3+1 with prime p has p=3. The packet does not falsely assert uniqueness of 1729 or all n at p=3. Its subsequent Camion implication is valid: for each p|N, Korselt gives N≡1 mod(p−1), the inherited quasi-period propagates D(1)=0 to D(N)=0 modulo p, and squarefreeness permits CRT. At 1729 the direct Wright recurrence reproduces s=c=1 modulo1729. At 21 it gives s=c=16 but 2^20=4 modulo21, disproving reversibility. This concerns strong labelled tournament counts, not the fixed-H generating series in the no-go theorem.

The proposed historical correction is justified by reading the actual reflection: its p=19 inert factors refute the proposed universal split-prime pattern, but do not refute all modular interpretations, prove an e-limit or prove a universal Paley maximizer. The proposed repair correctly preserves numerical provenance and changes only unsupported inference. The second reflection's vertex-addition operator is not a bijection of finite tournaments, so it cannot itself be a group generator in the asserted ordinary action. Neither critique should be inflated into an absolute future-map impossibility theorem.

## Primary source access and limits

- Busch, publisher abstract, https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3 : bound independently located in publisher search result. DOI opening failed; alternate full PDF returned403. This review verifies the statement/attribution, not a fresh full-proof audit.
- Borcherds, The Leech lattice, https://math.berkeley.edu/~reb/papers/leech/leech.pdf : accessed full PDF text; §3 explicitly gives lattice/null-orbit correspondences; §§6–7 give Leech/deep-hole geometry.
- Borcherds–Conway–Queen, https://math.berkeley.edu/~reb/papers/shallow/shallow.pdf : full text accessed; Theorem1 explicitly separates 23 deep and284 shallow types.
- Milne, https://www.jmilne.org/math/CourseNotes/MF.pdf : full text accessed; modular normalization, elliptic moduli and CM reference. Finite identities independently checked above.
- Cohn–Kumar–Miller–Radchenko–Viazovska, https://annals.math.princeton.edu/2017/185-3/p08 : publisher abstract confirms optimality and periodic uniqueness, not unrestricted uniqueness.
- Luca, https://link.springer.com/article/10.1007/s11139-026-01450-1 : publisher full HTML accessed; exact Theorem1 and definitions verified. Not an independent reproof of its deep number-theoretic estimates.
- Stark's primary source is pinned in the producer report with an explicit access qualification. This audit does not remove it.

## Catalogue and patch boundary

The catalogue is an indexed discovery pass over tracked filenames, physically available UTF-8 content and four named targeted connector fetches, not all historical content. Nontext, oversized and missing files are quantified. Broad matches and truncated previews are distinguished from relevant proofs. The THM-227 patch correctly preserves Newton trace identities while removing invalid forbidden-H inference. The proposed alpha vector in the old gap synthesis evaluates to 2T−5, not2T+1; an independent triple also forces vertices and pairs, defeating alpha=(0,3,1). These refutations do not prove global spectrum completeness.

## Final repair and provenance record

The producer confirmed the report-domain clarification and both wording repairs were applied before final hash binding. THM-227 corrections.patch passed a dry-run against a copy of the current target, outside the repository. Four connector-fetched source copies initially contained one surplus terminal newline; independent Git blob hashing showed that deleting that exact one byte reproduced every reported remote blob SHA. This is a serialization artifact, not conflicting mathematical source text; final copy status is recorded in manifest.json. The modular producer's later PRIOR_1729_AND_E8_HOLES.md and the repunit producer's REMOTE_SEARCH.md are separate addenda outside this review's frozen scope.
