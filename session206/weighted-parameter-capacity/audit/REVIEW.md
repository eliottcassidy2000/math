# Independent weighted residue-capacity review

Verdict: PASS for `session206-weighted-capacity/WEIGHTED_CAPACITY.md`.

The source and modulus are fixed throughout; the theorem counts encoded tail inputs with multiplicity. Physical attainability, positive width R>d>0, and odd scalar parity are inherited assumptions, not asserted for arbitrary integer tuples.

## Parity coefficients

For input parity b, the first radical-inverse bit is b, so truncation residues range over precisely [b/2,(b+1)/2). The normalization is 1/K per low residue, not 2/K. Consequently the two integrals give exactly c0=log(R/(R-d/2))/d and c1=log((R-d/2)/(R-d))/d, with no missing factor 1/2. Their sum is the established total capacity. The same residue-class counting proof gives fixed-source O(sqrt X) errors. Equal weighting of the input parities under the scalar cutoff would indeed be wrong.

## Joint modular distribution and gcd arithmetic

After middle-block mixing, (a mod2p,q modp) is uniform. A fixed L-bit suffix maps a modp bijectively because 2^L is invertible modulo odd p. Conditional on each old a, q undergoes a translation, leaving q uniform and independent of the final first coordinate. Thus after suffix conditioning, the two residues U=m modp and V=Q(m) modp are independent and uniform, while parity is fixed.

The homomorphism (U,V)->RU-dV has image the multiples of g=gcd(R,d,p), including composite p. Its image has p/g elements and each fiber has pg elements; hence each permitted scalar residue has probability g/p. Outside the coset h-eb modg, the scalar count is exactly zero for that input parity. The p=1 case is immediate and agrees with both parity coefficients.

## Cutoff and finite-bank consequence

The audited endpoint-cylinder argument applies with slopes bounded below by fixed R-d>0 and m<=C_s X. Boundary counting on each low-residue progression bounds the total by O_s(X2^-L+log X+2^L), including after parity restriction. Fixing epsilon and L before taking X large leaves finitely many cylinders, so mixing error sums to o(X). The subsequent epsilon and L limits establish the stated fixed-source/fixed-modulus asymptotic without an independence assumption on the cutoff. Fixed finite-bank sums are legitimate. Replacing all d=1 sources by one root-image bound uses the independently established d=1 collapse and residue-preserving involution; no general-source collision factor is silently applied.

## Physical example and exact inequality

Applying the physical operations gives JJ=(5,6,2,3) and OZJ=(11,12,4,5), with the listed intermediate states. For p=3, JJ omits scalar residue 1 and OZJ has g=1, hence contributes log(12/7)/15. Root distinct outputs contribute at most log2/6. Strict trapezoid bounds yield log2/2+log(12/7)/5 < 3/8+19/168=41/84. The class has X/6+O(1) elements, so the asserted missed-value bound X/252-o(X) is valid. The sharper logarithmic coefficient is also correct. The unweighted bound fails because its root-plus-JJ part already exceeds 1/2. Since gcd(12,5)=1, OZJ has all parity-compatible modular projections by the separately proved exact modular-image theorem. Thus this genuinely strengthens a mere omitted-class criterion, without proving a universal finite-bank obstruction.

## Reproduction

Optimized replay of the author's explicit-check verifier matches RESULT.json byte-for-byte. It checks 2,688 finite homomorphism fibers over 21 generated physical states and records exploratory cutoff counts. These regressions support implementation consistency but do not establish universal quantifiers or effective convergence rates.
