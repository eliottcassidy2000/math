# Session184: four proof results and the remaining major gaps

2026-11-04. **PROVED:** termination from the common TT3 entry of an entire specified tournament return grammar; infinite arithmetic progressions with exact Collatz first-descent times; complete modular certificates for prescribed odd-perfect supports; and sharp atom classifications for infinite bounded-degree hole sets. **OPEN:** unrestricted tournament-spectrum completeness, Collatz convergence, and odd-perfect-number nonexistence. Historical novelty is not established.

The request for another substantial proof session led to four complete, independently checked scoped results. None supplies a proof of one of those major conjectures. The tournament work closes a natural proposed induction language negatively; the other packets give positive theorems whose remaining global quantifiers are stated below. They should not be presented as major-conjecture solutions.

## 1. The natural singleton–triangle induction is completely decided

The [tournament proof](../spectrum/PROOF.md) continues Session183's actual terminal-triangle idea. Its carrier consists of ordered singleton and directed-triangle strong components, arbitrary actual port subsets, a transitive prefix TT_t, and one nonuniform root. Its count has the inherited form

    H_t = c + 2^t F.

The new inverse theorem decodes F into one-bit singleton and two-bit triangle chunks. Together with a normalization removing initial unmarked singletons, it makes all clock-uniform digit successors decidable with no suffix-length cutoff and with arbitrary fixed integer clock shifts. Applying that complete inverse to the common fully marked TT3 entry gives exactly the language

    {empty, 2, 4, 22, 42, 422, 4222}.

The seven normalized states have maximum depth four. At clock t=1, the two maximal chains are 9→29→89 and 9→31→95→287→863. Neither original sibling admits a +4 successor anywhere in this carrier, even allowing arbitrary replacement inside it. The finite certificate enumerates 112 distinct actual candidates and 280 digit comparisons; completeness comes from the inverse and normalization proofs, not a guessed search bound.

The same packet determines the leading 3-adic coefficient: for k≥2 triangle components, H_t/3^(k−2) is congruent modulo 3 to the negative sum of adjacent triangle port-count products. The coarse divisibility bound is explicitly inherited from Session163; its exact coefficient and cancellation condition refine that bound in this ordered carrier.

The [order-four extension](../spectrum_four/PROOF.md) then permits every suffix whose strongly connected components have order at most four, including all actual port subsets of the unique strong four-vertex tournament. A finite normalization and positive-recurrence search bound cover every suffix length and fixed integer clock shift. The complete return tree from the same fully marked TT3 entry is still exactly the seven states above. Consequently a missing sibling +4 return in this same attachment architecture, with the same clock-uniform requirement, needs a suffix SCC of order at least five. This is necessary, not sufficient; it does not exclude a one-clock return or a different architecture. The [extension audit](../spectrum_four_audit/AUDIT.md) independently checks the enlarged-bank normalization and complete successor certificate.

[Independent audit](../spectrum_audit/AUDIT.md): all-length inversion, arbitrary-shift normalization, complete tree and exact coefficient pass. This is not a global count obstruction, a native-operation certificate, or an old-name-preserving construction. The known global exceptions remain 7 and 21. **Next proof target:** a return carrier with a strong suffix component of order at least five, more nonuniform interfaces, or genuinely non-clock-uniform transformations, followed separately by a coverage argument.

## 2. Exact first descent survives a genuine template change

The [Collatz proof](../collatz/RESULT.md) gives an arbitrary-word theorem. For positive odd-map valuations a_0,...,a_(m−1), let P_j=3^j, Q_j=2^(a_0+...+a_(j−1)), B_0=0, and B_(j+1)=3B_j+Q_j. If the first multiplicative crossing below 1 is at m, then the residue

    n = (Q_m−B_m) P_m^(-1) modulo 2Q_m,
    n > B_m/(Q_m−P_m)

selects exactly the stated valuations and has first strict descent at odd step m. It defines an infinite progression of density 1/(2Q_m) among positive integers.

For W_q=(1,2)^q followed by r(q) letters 2, where r(q)=floor(q log(9/8)/log(4/3))+1, this certifies an unbounded family of genuine template switches and exact first-descents using only valuations 1 and 2. The q=1000 control has first descent at step 2410. The final slopes approach 1 along a subsequence, and the language of these exact words is not context-free by the classical Parikh semilinearity theorem. Thus the small-margin words are realized on positive integer trajectories, rather than merely formal branch strings.

[Independent audit](../collatz_audit/AUDIT.md): residue admissibility, strict thresholds, the S179 mismatch identity, arbitrary expanding-template rotations and the formal-language argument pass. This neither proves a common rank across every possible template switch nor rules out a conservative certificate language covering every source. **Next proof target:** prove that every positive odd source above 1 has a certified descendant strictly below that original source, or supply another well-founded global selector. Infinite progressions and arbitrarily long certified words do not establish that coverage claim. Descent from a later, larger entry value alone would not suffice.

## 3. Prescribed odd-perfect supports have complete modular certificates

The [odd-perfect proof](../oddperfect/PROOF.md) starts from S179's finite primitive-divisor exponent domains for a prescribed finite support S and Euler prime q. Exact divisor sums give finite signed valuation row lists R_p. An odd perfect number on precisely that support exists if and only if one row from each list sums to the integer zero vector.

For each coordinate, bound every selection sum in [L_r,U_r] and choose m_r=1+max(|L_r|,|U_r|). Then a sum is zero modulo every m_r if and only if it is literally zero: the bounds prevent wraparound. Exact nonnegative group-algebra coefficients, or the equivalent finite convolution recurrence, therefore count the actual exponent assignments. Smaller moduli give sound rejection filters but not conclusive acceptance.

[Independent audit](../oddperfect_audit/AUDIT.md): primitive-divisor exceptions, exact row closure, integer selection, no-wrap moduli and positive controls pass. The experiments do not demonstrate an arithmetic advantage over the inherited fixed-support procedure: all tested complete small support domains already contain an empty row list. The integrality distinction is demonstrated with a labelled abstract catalogue. General modular dynamic programming and classical fixed-prime-count finiteness are not claimed as new. **Next proof target:** an obstruction that uniformly controls varying supports. Deciding each given support is not a bound on all possible supports.

## 4. A positive infinite-hole classification closes a separate scoped target

The [hole-atom proof](../nearby/NOTE.md) extends the repository's finite-hole decomposition to arbitrary, possibly infinite H with Omega(h)≤K. Here Omega counts prime factors with multiplicity, K is a positive integer, weak atoms lack any retained proper factor pair, and strict atoms lack any retained unequal proper pair while having a nonempty ambient strict fibre.

The sharp maximum degrees are 2K+1 for weak atoms and 2K+2 for strict atoms. At weak equality, every divisor of degree 1 through K must be a hole. At strict equality, the target must be p^(2K+2), with p,...,p^K holes. The packet gives the complete rank-ideal family H={n:1≤Omega(n)≤K} and proves both associated Dirichlet series have leading growth

    (log(1/(s−1)))^(2K+1)/(2K+1)!  as real s decreases to 1.

For an arbitrary prime-hole set, it also classifies the support exactly: the weak atoms are semiprimes touching H and triple prime products entirely in H; the strict version deletes hole-prime squares and adds hole-prime fourth powers. Its conditional weighted-profile transition occurs when Q(s) is comparable to the square root of the full prime sum P(s). This is a precise infinite-hole result, not an asserted density theorem for arbitrary prime sets.

[Independent audit](../nearby_audit/AUDIT.md): sharp bounds, equality cases, strict-diagonal exclusions, support-level counting and conditional asymptotics pass. The prime-zeta symmetric-function machinery is classical and credited. **Next proof target:** classify other explicit infinite hole families or identify additional hypotheses beyond bounded Omega that still permit exact support control. No consequence for the three major conjectures above is asserted.

## Sources and reproduction

Each primary proof links or identifies its exact predecessors and external prerequisites. In particular, the arithmetic packets retain the classical [Parikh theorem](https://doi.org/10.1145/321356.321364), [primitive-divisor theorem source](https://math.colgate.edu/~integers/k27/k27.pdf), and [fixed-prime-count odd-perfect bound](https://math.colgate.edu/~integers/d14/d14.pdf). The hole-series packet credits [Mathar's almost-prime series treatment](https://arxiv.org/abs/0803.0900v3). These references support the stated ingredients, not a novelty claim for this session.

All producer checkers, including the order-four extension, and the independent packet checkers have matching normal/optimized mathematical receipts. Each packet distinguishes exhaustive finite controls from the proof establishing its unbounded claims. No hosted-CI status or global conjecture verification is inferred. The substantive outcome is four bounded-scope proof results and a clearer next obligation in each lane, not a resolved major headline.
