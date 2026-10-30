# Independent audit of marked periodic content rates

**Verdict: PASS for the stated scoped theorems.** Independent algebraic proof review and separately implemented exact tests. No claim of Collatz convergence, universal finite-word bounds, or formal proof-assistant certification.

## Algebraic review

The affine order formula and cyclic-membership cutoff were checked separately in all three cases.

- When the multiplier has nontrivial order modulo5, its difference from1 is a unit. Integral fixed-point conjugation leaves the multiplier subgroup and the commutator divisibility condition. Its principal-unit depth gives exactly the stated cutoff.
- For multiplier1 modulo5 and translation at least as deep as multiplier-minus1, division by the latter is integral. Requiring the commutator valuation to be at least e+s is necessary, not merely e. The draft retains this factor correctly.
- When translation is shallower, the auxiliary multiplier 1+(m−1)c/b has the required power precisely when v5(c)≥min(h,e). The independent multiplier match is equivalent to commutator depth at least e+h. The draft correctly uses the translation pivot and subtracts r=h from its cutoff.

At each admissible layer, the successful controller exponent is one residue class modulo the exact affine order. The finite/infinite geometric sums consequently give the displayed rational formula.

## Ordinary all-layer mean

The proof does not exchange two unrelated limits. Finite commutator-depth phases are bounded. At an infinite-depth phase, the actual positive valuation words commute, and the inherited faithful affine-word representation plus the elementary commuting-word lemma reduces them to powers of one primitive word. Faithfulness can also be seen directly: length and clock are fixed by the affine slope, while the first exponent is recovered from the valuation of B−3^(d−1), then stripped recursively.

Conjugation by the phase map preserves every identity-modulo-5^e predicate. Its content is therefore that of a fixed affine root to an exponent ak+b. The multiplier-order and geometric-sum LTE formulas bound each observed height by a constant plus log5(ak+b). Counting each successful residue class introduces at most one discrepancy per layer; only O(logN) layers can occur among the first N blocks, and the geometric expected tail beyond them contributes O(1). Thus the ordinary untruncated Cesaro limit and O(logN) discrepancy are justified. The all-finite-depth case has bounded discrepancy and a genuinely periodic charge sequence.

The added fixed-exit corollary is also valid. For each nonempty prefix V of a fixed exit, the same commutator with F_(Vu) either bounds the terminal height uniformly or gives a common-root exponent aN+b and logarithmic height. There are only finitely many exit prefixes. Their sum does not change the asymptotic mean. Entry and exit must remain fixed as N grows.

## The asymptotic bound

The condition ord25(M)=5 ord5(M) is equivalent to r≤1, including the r=0 translation case. It makes each higher-layer return density at most 5^(1−e) times the first-layer density. The inherited all-word first-layer bound applies to the literal entry followed by N controller copies. Removing the fixed entry charge before taking N to infinity yields the claimed asymptotic5/8 bound. The all-one controller gives equality. The primitive-word-root extension is valid because grouping a fixed number of primitive blocks changes neither the mean per clock nor the retained entry. This does not cover every stalled controller or any arbitrary changing-word sequence. The stated stall hostile correctly refutes only the intermediate universal5/4 tail multiplier.

## Independent controls

`check_periodic_independent.py` checks every affine multiplier/translation residue pair modulo25 against literal cyclic subgroup generation at depths1 and2, using positive integer lifts. Additional deliberately stalled unit- and translation-pivot examples run through depths3 and4. It passes **548,000 membership tests and1,016 order tests**.

`check_literal_means.py` independently constructs affine coefficients by the explicit finite sum, then evaluates actual unbounded-integer prefix contents by gcd and5-adic valuation. Ten controllers and four entries, at four block horizons through200, give **160 finite-prefix discrepancy controls**, including the order-stall hostile and repeated controllers. These tests are controls for the proof, not a finite empirical substitute for its limit argument.

Both scripts use explicit exception gates. Normal and optimized outputs match byte for byte. Receipts are `independent_results.json` and `literal_results.json`.

Run both scripts from this directory. Python3 standard library only. The literal checker loads definitions from the independent checker without executing its enumeration; it imports no producer implementation.
