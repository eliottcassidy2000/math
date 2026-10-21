# Independent fixed-degree rationality review

Verdict: PASS for the final proof. Both minor wording cleanups were applied and verified. The original ambiguous odd-degree realizability claim was corrected during review. This audit treats the named finite-atom theorem and center law as inherited inputs, not as newly reproved facts.

## Algebra and quantifiers

- A strong atom of fixed odd h has deficit at most h. The inherited effective strong-order bound therefore makes every P_h finite, including h=1 after separately allowing the singleton. Strong positive-deficit center families are finite by the same theorem. Empty centers are an explicitly distinguished extra object, not strong tournaments.
- Ordered SCC decomposition is essential: deleting rigid components leaves an ordered multiplicative factorization, with r+1 independent rigid gaps. Each factor exceeds one, hence r<=Omega(h). Repeated types and distinct atoms of identical size/signature retain their correct multiplicities.
- A self-converse SCC word is uniquely a half-word, an optional single fixed center, and its opposite reverse. The half is an ordinary type, not a converse orbit. This validates J_h without an extra factor two and validates the inherited S_d center-divisor formula.
- The bound x(x-1)a<=d follows from b<=a. For x>1 it bounds both x and a; for x=1 the strong deficit-d theorem supplies finiteness instead. This separation avoids accidentally using a zero denominator at x=1.
- Fundamental black-degree classes split exactly into S_(2k) and the converse pairs of non-self-converse h=k classes. The latter exist only for odd k and have series (F_k-J_k)/2. The degree-zero formula counts the empty object at n=0, consistently with the stated inherited convention.
- Every x with x²|k is an odd divisor of 2k and satisfies x(x-1)<=2k, proving the claimed joint denominator exponent bound. Denominators have integer coefficients and constant term one, so they yield eventual integer-coefficient recurrences even though the convenient numerator expression contains factors 1/2.
- The positive root rho is the unique minimal-modulus p-root. Its simple-pole normalization is A(z)~[rho(1+3rho²)]^(-1)(1-z/rho)^(-1). For k=3^r the unique maximum-length factor sequence consists of r factors 3, yielding precisely the displayed leading constant. Other factors can affect lower-order terms but not this constant.
- More general rho-scale growth requires F_k nonzero, not merely a self-converse node of degree k. The proof now states this correctly. The checker includes a formal alphabet with degree 7 present but F_7=0 as a diagnostic counterexample to the unqualified implication; it makes no tournament realizability claim.
- The bounded-degree density corollary is valid. At fixed K the numerator is exponential-polynomial in n, whereas the merged-node count is at least 2^(n(n-1)/2)/(2n!). Their ratio tends to zero. The specified uniform merged-node sampling, fixed K, and no-connectivity qualifications are necessary and present.

Final wording verification: the pole statement now uses modulus at least sqrt(rho)>rho, and the incompatible 'or no such factorization exists' clause has been removed under F_k nonzero. Both minor cleanups are complete; neither changes the mathematics.

## Inheritance and novelty boundary

Closest sources checked directly:

1. math/05-knowledge/results/deficiency_three_20261020/centers.md already proves the full arbitrary-positive-deficit divisor decomposition and nonrigid-half bound. Those formulas are inherited here, not new.
2. math/05-knowledge/results/effective_black_defect_order_bound_20261020/PROOF.txt supplies effective finite strong alphabets and centers.
3. math/05-knowledge/results/metagraph_classification_20261021a/README.md already proves the h=3 atom classification, F_3=z^5 A², J_3=z^5 A(z²), rational B_3 with its unknown R_6 polynomial, and the failure of the quasi-polynomial expectation.
4. math/01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md supplies the current all-order B_0 series and the center law.
5. math/05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md records the leaf numerator (2z^4-z^5-z^7)/(p(z)p(z²)).
6. math/05-knowledge/results/session21_recovery_20261021/RECONSTRUCTED_MATHEMATICS.md and RECONSTRUCTED_METAGRAPH_PROOFS.md already carry stronger exact trivalent recurrence/asymptotic details. The present theorem should not claim novelty for the r=1 case.

The defensible added synthesis is rationality for every fixed fundamental degree, its general denominator bounds, the explicit 3^r leading constants beyond inherited special cases, and the immediate fixed-bounded-degree density consequence. Repository phrase searches did not locate an earlier general fixed-degree rationality theorem; that is a bounded search finding, not a claim of global literature priority.

## Independent executable controls

Run from workspace root:

    python session22-fixed-degree-review/check.py
    python -O session22-fixed-degree-review/check.py

Both passed with identical output. The checker does not use assert, so optimization does not erase verification.

- 24,419 literal words through size 10 over a nine-letter formal involutive alphabet, without a presentation-weight pruning cutoff.
- Two distinct letters share the same size/signature; another converse pair is genuinely non-fixed. This tests multiplicities and converse orbit counting independently of tournament realizability.
- F_h and J_h agree with direct word enumeration for odd h<=499; parity of F_h-J_h is separately checked.
- S_d agrees for every even d<=200, and B_k agrees for k<=100.
- SymPy exact denominator divisibility checked for k<=40.
- Known leaf simplified numerator and the inherited trivalent formula with arbitrary R_6 checked symbolically.
- Pole-normalization constant checked symbolically for r=0,...,9.
- Total: 1,004 explicit checks. These are formal algebra controls, not large-order tournament censuses.

The root's separate checker was also read. Its h<=81 pruning is safe for its degree<=10 comparisons: nonrigid halves satisfy deficit>=2h/3, forcing h<=30, and its rigid-half strong centers have h<=9. No source math files or Git/remote state were modified by this review.

Final frozen proof SHA-256: e428906db3b67dbf5c85eef86460a012be8a9574a06bcd62759e73778f23446d
