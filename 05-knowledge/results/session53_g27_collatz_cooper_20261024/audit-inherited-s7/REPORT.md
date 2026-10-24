# Independent audit of inherited s7 square-root prime dichotomy

2026-10-24. Audited `05-knowledge/results/cooper_s7_weight_one_root_frobenius_20261018.md`, especially Theorem D, without using its PROVED label as evidence. No repository file was changed.

## Verdict

No mathematical counterexample or invalid implication found in Theorem D. The proof works with the classical modular parametrization and standard modular-form facts stated in its setup. I independently reconstructed the formal-algebra deductions and checked the exceptional-prime Sturm certificates. This is an audit of that proof, not an independent claim of novelty or a replacement for literature review.

For every prime p, the all-n congruence c(pn)=c(n) modulo p^2 holds exactly when p=7 or p mod7 belongs to {1,2,4}. For every other prime there exists n for which even the mod-p congruence fails. The written proof does NOT establish that n=1 always witnesses failure; it establishes existence. n=1 works for every inert prime in my finite range.

## Normalizations and integral rings

The uniquely normalized characteristic-zero square root has G(0)=1. From the cited identity F(t(q))=theta(q)^2, formal uniqueness gives G(t(q))=theta(q), not its negative. Here h(q) and t(q)=h/(1+13h+49h^2) lie in q Z[[q]] with leading coefficient 1. Consequently the compositional inverse q(t) belongs to t Z[[t]], proving G=theta(q(t)) is integral at EVERY prime, including 2. No unjustified division by 2 in characteristic 2 occurs. Similarly s(q)=(1-49h^2)/(1+13h+49h^2), so s(t) is an integral unit with constant 1; its square-root notation is interpreted in characteristic zero before reduction.

The excellent lift t^sigma=t(q(t)^p) has t^sigma/t^p in 1+t Z[[t]], congruent to 1 modp. Thus epsilon_p belongs to t Z[[t]]. Coordinate substitution q<->t is an integral formal-ring automorphism, so coefficientwise congruences transfer both ways.

The p>=5 logarithm argument has p-integral coefficients: log Phi_p = p sum_M sum_{d|M,p∤d} d^{-1}q^M. Denominators are units at p. Dividing the mod-p^2 equality by p in section5.5 is legitimate only after both sides have been shown divisible by p, which the proof does. The factors 2,3,7 are never inverted at a problematic prime: odd split primes are >=11; p=2 and p=7 have separate arguments. In the p=7 Eisenstein identity the displayed denominator is 8, a 7-adic unit.

## The crucial new modular adaptation

Beukers, https://arxiv.org/html/2403.03301v3, equation(5), Lemma3.6 and Proposition4.2 give the product construction, modularity, and coefficientwise divisibility mechanism. Its original single-pole Hauptmodul theorem does not directly prove the inherited claim. The inherited two-pole/Fricke argument is necessary and appears valid:

- Theta has exactly the two order-3 elliptic zeros, with orbifold order1/3 at each, and no cusp zeros. This uses the classical theta transformation, elliptic stabilizers and valence formula.
- The coefficients T_(p,m) have poles of order at most floor(m/3) at each point.
- Conjugation by W7 exchanges the determinant-p cosets and changes their character weights by chi(p). Hence coefficient m has Fricke eigenvalue chi(p)^m.
- The invariant pole-bounded subspace is polynomials in t of degree <=m; the anti-invariant subspace is s times polynomials of degree <=m-1. Integral q-expansions then give integral coefficients within the relevant eigenspace, via the integral inverse coordinate and division by the integral unit s.
- This yields precisely the degree bounds used in Theorem B. The norm-product and trace arguments give the modp and modp^2 coefficient assertions even for p=2. No division by p is needed in Theorem B.

The direct-sum sentence about integral coefficients should be read separately on each eigenspace. Integral q-expansion of an arbitrary mixed invariant/anti-invariant function alone need not imply both components have integral polynomial coefficients. The actual proof only applies the assertion after an eigenvalue has been established, so this wording issue does not invalidate it.

## All-n implication at split primes

Cartier applied to G=G_p G(t^p)+p epsilon_p G_p (vartheta G)(t^p) gives C(G)=G+p C(epsilon_p G_p) vartheta G modp^2. The Eisenstein/Fricke argument shows epsilon_p G_p modp is a polynomial of degree <p with constant0. Cartier therefore kills it, proving equality of complete infinite formal series, not just a truncation.

For clarity, the constant term of the invariant projection A_p is (2+lambda+lambda^{-1})/4 = 1+(lambda-1)^2/(4lambda). This proves A_p(0)=1 modp^2, including the exact strength needed after division by p. The anti-invariant projection similarly has constant (2-lambda-lambda^{-1})/4 = 1-(lambda+1)^2/(4lambda).

## Exceptional primes

At p=2, the modular form of weight12 and level14 in section5.6 is holomorphic: theta^9 clears the largest t-power denominator and theta(2tau)^3 clears the t(2tau) denominator. Its character is trivial. Index24 gives Sturm bound24. The checked congruence implies epsilon_2=t mod2, and epsilon_2 G_2=t mod2; Cartier kills it.

At p=7, f=theta^3 s is holomorphic of weight3 and character chi. Its Eisenstein decomposition gives U7 f-f=-294 A2, divisible by49. The weight30 level49 form clearing the lift identity is holomorphic (theta^27 clears t^9, theta(7tau)^3 clears t(7tau)), with trivial character. Index56 gives Sturm bound140. This yields epsilon_7=2t+4t^2 mod7. The Lucas identity G=(1+t)^2G(t^7) mod7 makes epsilon_7G supported in residue classes1,2,3,4 modulo7. Formal residue invariance justifies c(n)=CT_q(f t^{-n}) and its use with Laurent series. Every constant term involves finitely many coefficients. Therefore the constant-term argument proves the congruence for all n.

## Necessity

For inert p (necessarily odd and not7), assume C(G)=G modp. Theorem B and the Cartier product rule imply C(s H_p)=1. With Q=s^2, Frobenius gives s=Q^{(p+1)/2}s(t^p)/Q(t^p), whence C(s H_p)=s P/Q for a polynomial P. Thus s=Q/P is rational. But Q=(1+t)(1-27t) has distinct simple roots for inert p>=5; at p=3 it is1+t. It is not a rational square in any case. This contradiction proves existence of a failed n for each inert p. It does not identify that n universally.

## Independent reproducible checks

`audit.py` is standard-library-only and does not import the original scripts. Its output is `audit.out` in this directory.

1. Build theta from divisor sums and h,t,s from Euler products; check S4,S5,S6 through q^180 exactly.
2. Check the C7 modular congruence through q^180, well past bound4; its first coefficients are 1,14,91,406,1463.
3. Check the weight3 Eisenstein identity through q^180 exactly, past bound2.
4. Check the p=2 and p=7 modular-form congruences through q^180 modulo4/49, past bounds24/140.
5. Compute F(n) directly from its binomial sum and recursively extract the normalized integer square root through n=180; independently verify the printed second-order recurrence and G(t(q))=theta(q) through degree180.
6. Extend the recurrence with exact integer divisions through n=20000, then verify Theorem D for every prime <=500 and every pn<=20000. All pass. Every tested inert prime fails already at n=1 modulo p. Explicit examples: c(3)=336 versus c(1)=2, difference1 mod3; c(5)=117348, difference1 mod5; c(13)=8123136091556640, difference9 mod13.

The finite range is evidence, not the source of the all-n theorem. The modularity/holomorphy and Sturm argument turn the exceptional-prime bounded checks into universal congruences; the arbitrary-split-prime proof remains symbolic.

## Provenance cautions

MISTAKE-629 concerns adjacent session claims (misidentified sequences, novelty, provenance and mod-p^3 overclaims), not a disproof of this Theorem D. I avoided those failure modes by independently checking the exact s7 definition and normalized root, and I make no mod-p^3 or novelty claim.

The publisher positively verified that the original GP script, proof note, and recorded output are tracked at commit 7818149bf0a1eca778549fdfd1dde34f479421fd and materialized them under session53-publication/inherited-s7-sources. Exact copies and the publisher provenance are frozen under sources/inherited here. The script was merely omitted from the original sparse materialization. Its SHA256 is 9034554ca861d527726ae4436a35be42b977cfae1f77f7d8eb0a527c401aca64; the recorded output SHA256 is f42e870e36b692fc96b3be90d6ecca970133c764d944e2c5296adf71e2bef506. I did not rerun the GP script or certify all its output sections; audit.py independently reconstructs the checks in this audit scope.

## Additional p=3 consequence requested during audit

The root and producer supplied the following consequence for independent review. It is correct in F3[[t]]. From known s7 Lucas congruences, F=(1+t)F(t^3); since F=G^2 and G(t^3)=G^3, cancellation of the unit G^2 gives G^4=(1+t)^(-1). Set P=product over j>=0 of (1+t^(9^j))^2. For each finite product, P_J^4=product(1+t^(9^j))^8=(1+t^(9^(J+1)))/(1+t), by the characteristic-3 identity (1+x)^9=1+x^9. Taking the t-adic limit gives P^4=(1+t)^(-1). Uniqueness of the fourth root with constant1 (4 is a unit modulo3) gives P=G.

Therefore c(n) mod3 is zero unless every base9 digit is0,1,2, and otherwise equals (-1) raised to the number of1 digits. The identity G=(1+t)^2G(t^9) gives C3G=G(t^3). Also G=sG(t^3), by normalized square-root uniqueness applied to F=(1+t)F(t^3), so H=G/s=G(t^3). Hence C3G=H and C3H=G: Cartier swaps the two series. In particular c(9n)=c(n) and c(27n)=c(3n)=d(n) mod3. It is not a scalar-sign relation between C3G and G.

The independent exact integer implementation verifies the digit criterion for all n<=20000. I have not established prior-art or novelty for this compact consequence.

## Acceptance boundary

Accepted: Theorem D with its stated classical modular inputs, its necessary Lemma T/Theorem B machinery, the formal integral normalizations, exceptional-prime Sturm certificates, and the p=3 digit/Cartier consequence. Finite coefficient verification is explicitly separate. Not certified here: other CM tables/Corollary F, Theorem E for every inert prime, higher prime-power conjectures, novelty/publication priority, the complete GP transcript, or every other theorem in the inherited note. Theorem A was not independently audited as a standalone result. No script success is represented as proving the unverified portions.

The manifest validator rejects missing, modified, extra, duplicate, and unsafe-path payloads with checks active under Python optimization. Normal and optimized outputs are independently generated. The portability harness copies only declared manifest payloads into a fresh temporary directory, validates that relocated copy, reruns both modes, and tests rejection controls. MANIFEST.json is the package entry point and hashes every other frozen file; its own hash is reported separately.
