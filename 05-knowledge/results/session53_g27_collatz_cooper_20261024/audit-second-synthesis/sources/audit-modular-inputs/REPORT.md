# Independent closing audit of the inert modular inputs

2026-10-24. Scope: the four remaining inputs to the same-precision two-step Cartier argument. No frozen packet or repository file was changed. This audit does not certify adjacent CM assertions or any priority claim.

## Verdict

**ACCEPTED, with classical modular-form inputs explicitly cited.** The inert Theorem B, the degree/constant-term assertion for epsilon_p H_p, the exceptional p=3 lift, and the companion identity all follow. Together with the separately audited Theorem D and the independently audited two-step formal argument, the statement c(p²n) ≡ c(n) modulo p², for every prime p and every n≥0, can honestly be described as a proved deduction with the cited classical inputs. There is no remaining conjectural hypothesis among these four inputs. This is not a claim that the theorem is previously unknown.

Two wording repairs matter: (i) the anti-invariant projection has constant term congruent to 1 modulo p², not identically 1; its exact value is given below; (ii) a p-integral q-expansion implies integral polynomial coordinates after an eigencomponent is established, not for an arbitrary mixture of the two eigencomponents.

## Sources and what is actually used

The inherited note is preserved in sources/inherited-proof.md, and the specific two-step proof reviewed is preserved in sources/two-step-proof.md. Their hashes are in MANIFEST.json. The earlier independent audit at /workspace/shared/session53-audit-inherited-s7 accepts Theorem D and the general modular machinery; the present proof rechecks the inert extension rather than inheriting that audit's deliberately excluded Theorem E.

Primary source inspected: F. Beukers, *Supercongruences using modular forms*, https://arxiv.org/html/2403.03301v3, equation (5), Lemma 3.6, Proposition 4.2 and its proof, and Appendix C §9.11. These supply the determinant-p coset product, permutation mechanism, trace/norm divisibilities and the s7 parametrization. The single-pole theorem in that paper is **not** being applied directly to the present two-pole situation. The adaptation below is proved here. Other classical inputs are the binary-theta transformation law (or Poisson summation), the modularity and cusp holomorphy of character Eisenstein series and eta products, the valence formula, von Staudt–Clausen, and Sturm's coefficient bound. Numerical matching is not used to assert modularity.

## 1. Geometry, cusps and characters

Use h=(eta(7τ)/eta(τ))^4, D=1+13h+49h², t=h/D, s=(1−49h²)/D, and theta=sum q^(x²+xy+2y²). The theta series lies in M_1(Gamma_0(7),chi), chi(d)=(d/7)=(-7/d). Poisson summation for Gram matrix [[2,1],[1,4]], whose adjugate is integrally equivalent to it, gives theta|_1 W_7=−i theta under the determinant-normalized slash operator. In particular theta has nonzero constants at both cusps infinity and 0. There are no other cusps of X_0(7).

The standard Hauptmodul h identifies X_0(7) with P^1. Fricke sends h to 1/(49h), t to t and s to −s. The two roots of D are the order-3 elliptic points. The identity t theta³=E, where E=(eta(τ)eta(7τ))³, is the s7 eta identity. E has no zeros in the upper half-plane and has order 1 at each of the two cusps. Its weight is 3 and its character is chi: the eta-product character is (−7³/d)=chi(d), and both eta congruence sums are 24. Thus theta has local analytic order 1 at each elliptic point, equivalently orbifold order 1/3, and no other zeros. This also exhausts the weight-one valence total 8/12=2/3.

A modular function with poles of coarse order at most a at these two points and no other poles has invariant part a polynomial in t of degree at most a, and anti-invariant part s times a polynomial of degree at most a−1. Indeed the quotient by Fricke has coordinate t, with one pole downstairs; the two points are exchanged, not fixed. Dividing an anti-invariant function by s is invariant. At either fixed point of Fricke, anti-invariance forces a zero, canceling the simple zero of s, so no new poles arise. At the two poles of t, s has simple poles and division lowers the bound by one. This proves the stated eigenspace description without an unjustified mixed integral decomposition.

All formal coordinates used here are integral: h∈q Z[[q]], D is a unit, t=q+O(q²), its inverse q(t) has integer coefficients, and s∈1+t Z[[t]] is a unit. Therefore an integral q-expansion in an established eigencomponent gives integral polynomial coefficients, using q(t), and dividing by s in the anti-invariant case.

## 2. Inert Theorem B modulo p²

Let p be inert, so p≥3, p≠7, chi(p)=−1. Form T_p(x) as in equation (5) of the primary source, using theta. The determinant-p cosets are represented by diag(p,1) and [[1,b],[0,p]], 0≤b<p. Their character-normalized weight-one quotients are permuted by Gamma_0(7); consequently the coefficients of T_p are modular functions. Each numerator is holomorphic after passage to a common finite-index subgroup, at all cusps as well as in the upper half-plane. The denominator theta has no cusp zeros. Thus there are no cusp poles. At either elliptic point a product of m quotients has local pole order at most m. Invariance forces its Laurent pole exponent to be a multiple of 3, giving coarse pole order at most floor(m/3), not m.

Conjugating [[a,b],[7c,d]] by W_7 gives [[d,−c],[−7b,a]]. Since ad≡p modulo 7, chi(a)=chi(p)chi(d). The scalar −i from theta's Fricke transformation cancels in each quotient. Thus every character-normalized quotient is transformed to chi(p) times another quotient, and T_(p,m)|W_7=chi(p)^m T_(p,m). The geometry above now gives

- T_(p,p)=s alpha, with alpha∈Z[t] of degree at most floor(p/3)−1;
- T_(p,p+1)∈Z[t] of degree at most floor((p+1)/3).

Integrality of the q-expansions follows directly from the cyclotomic product: coefficients are algebraic integers invariant under the roots-of-unity permutations and have integral powers of q, so lie in Z[[q]]. No p-denominator is introduced by passing to the polynomial coordinates.

For completeness the coefficientwise divisibility argument does not require a one-pole Hauptmodul. Factoring the second product in T_p into its leading term gives its roots as chi(p)p theta(q)/theta(zeta_p^b q^(1/p)). All elementary products involving at least two of these roots are divisible by p². The sum of these roots is divisible by p² too, because the trace of 1/theta(zeta_p^b q^(1/p)) is exactly p times an integral U_p series. Hence T_(p,m)≡0 mod p² for m<p. The norm identity product_b theta(zeta_p^b q^(1/p))≡theta(q) mod p, together with theta(q^p)≡theta(q)^p, gives T_(p,p+1)≡chi(p)=−1 mod p. These follow also from the trace/norm proof in Proposition 4.2; no geometric single-pole assumption enters this step.

Write T_(p,p+1)=−1+pB, deg B≤floor((p+1)/3). Substituting the root x=theta(q)/theta(q^p) and canceling its unit pth power gives

G/G^sigma ≡ (1+pB)s alpha mod p².

After division by the integral unit s, H/G^sigma is congruent to a polynomial of degree at most d_p−1, with d_p=floor(p/3)+floor((p+1)/3). Since G^sigma=1+O(t^p) and d_p−1<p, this polynomial equals the truncation H_p modulo p². This proves the inert Theorem B, including p=3.

## 3. The inert epsilon polynomial for every p≥5

Set k=p−1 and lambda=(−7)^(k/2). Euler's criterion gives lambda≡chi(p)=−1 mod p. Let E_k denote the normalized level-one Eisenstein series and

Psi=(E_k(τ)+E_k(7τ))/(2 theta^k).

It is a weight-zero modular function with trivial character, holomorphic at both cusps, with coarse pole order at most a=floor(k/3) at the two elliptic points. Both E_k pullbacks are holomorphic everywhere. Its q-expansion and its Fricke transform are p-integral, because 2 and 7 are p-units and E_k is p-integral by von Staudt–Clausen.

The normalized transformations are E_k|W_7=7^(k/2)E_k(7τ), E_k(7τ)|W_7=7^(−k/2)E_k(τ), and theta^k|W_7=(−1)^(k/2)theta^k. Therefore

Psi+Psi|W_7=(1+lambda)(E_k(7τ)+lambda^(−1)E_k(τ))/(2 theta^k).

Each of the two numerator factors is divisible by p: the first since lambda≡−1, the second since E_k≡1 mod p and lambda^(−1)≡−1. Thus Psi≡(Psi−Psi|W_7)/2=s B_p mod p², where B_p∈Z_(p)[t] and deg B_p≤a−1. Its exact constant is

B_p(0)=(2−lambda−lambda^(−1))/4 = 1−(lambda+1)²/(4lambda) ≡1 mod p².

This verifies the sign and the p² constant precision explicitly; merely knowing Psi(0)=1 would not by itself establish that precision after projection.

To control the lift, put Phi_p(q)=product_n(1−q^(pn))/product_n(1−q^n)^p. Its logarithm is p sum_M(sum_{d|M,p∤d}1/d)q^M. All inner denominators are p-units. Von Staudt–Clausen gives E_(p−1)≡1−2p sum_M sigma_(p−2)(M)q^M mod p²; Fermat identifies the divisor sum modulo p. Hence Phi_p² E_(p−1)≡1 mod p². Logarithm/exponential manipulations are valid for p≥5 in 1+p q Z_p[[q]].

The eta identity gives exactly

t(q^p)/t(q)^p = Phi_p(q)^3 Phi_p(q^7)^3 (theta(q)^p/theta(q^p))^3.

By inert B the last quotient equals theta^(p−1)s H_p modulo p². This is 1 modulo p, as is each Eisenstein series. Expanding products modulo p² therefore gives

p epsilon_p ≡ 3 theta^(p−1)s(H_p−B_p) mod p².

Reduction modulo p first proves H_p−B_p divisible by p coefficientwise, because theta^(p−1)s is an integral unit and 3 is a p-unit. Only now divide by p and multiply by H_p. The relation theta^(p−1)s H_p≡1 mod p yields

epsilon_p H_p ≡ 3(H_p−B_p)/p mod p.

Coefficients of H_p above degree d_p−1 already vanish modulo p² by B; those of B_p vanish above a−1≤d_p−1. Thus the right side modulo p is a polynomial of degree at most d_p−1. Its constant term is zero, since H_p(0)=1 and B_p(0)≡1 mod p². All divisions and degree truncations have the required precision. No assertion about CM values is used.

## 4. Exceptional p=3: a full Sturm certificate, not a numerical fit

Define

M=theta(τ)^9 theta(3τ)^3 [t(3τ)−t(τ)^3]
 =theta(τ)^9 E(3τ)−E(τ)^3 theta(3τ)^3.

The second expression proves holomorphy without treating a rational t-expression as automatically holomorphic. Both terms have weight 12, level 21 and trivial character: theta and E have character chi, their 3-degeneracies have the same character on Gamma_0(21), and each term has character chi^10=1 or chi^6=1. Degeneracy pullbacks of holomorphic forms are holomorphic at every cusp and in the upper half-plane. More explicitly the cusp denominators are 1,3,7,21; E has respective orders 3,1,3,1 and E(3τ) has orders 1,3,1,3. Theta and theta(3τ) are cusp-holomorphic, so neither term has a cusp pole. Interior poles have disappeared in the polynomial expression.

The index [SL_2(Z):Gamma_0(21)] is 21(1+1/3)(1+1/7)=32. The weight-12 Sturm bound is 32. certificate.py constructs theta by direct lattice enumeration and E by its Euler product, independently of the inherited checker. It computes all exact coefficients through q^96 and records the 33 integers [q^0]M through [q^32]M and their exact quotients by 9 in normal.json. Every coefficient in that range, and all the way through 96, is divisible by 9. The standard coefficientwise Sturm congruence theorem works also when the tested prime divides the level; applying it modulo 9 establishes M≡0 mod9 as a complete series. Equivalently one can apply it modulo3, divide the form by3, and apply it again modulo3.

The prefactor theta^9 theta(3τ)^3 has integral unit constant 1, so t(q^3)−t(q)^3≡0 mod9. Dividing out the leading monomial q³ and the unit (t/q)^3 preserves integrality and the modulus. Consequently (t(q^3)/t(q)^3−1)/3≡0 mod3, that is epsilon_3≡0 mod3 after the integral change q↔t.

## 5. The weight-three companion and its character

Let A1=sum_{n≥1}(sum_{d|n}chi(d)d²)q^n and A2=sum_{n≥1}(sum_{d|n}chi(n/d)d²)q^n. The standard character Eisenstein series are E_{3,1,chi}=−8/7+A1 and E_{3,chi,1}=A2, both holomorphic at all cusps in M_3(Gamma_0(7),chi). The first constant is L(−2,chi)/2=−B_(3,chi)/6=−8/7, independently computed from B_3(x) in certificate.py. Thus

f*=1−7A1/8+49A2/8

is holomorphic with character chi; one must not call A1 alone a modular form. The proposed counterpart theta³−3E/4 is likewise in M_3(Gamma_0(7),chi). Eta cusp orders and theta holomorphy were verified above. Their coefficients through q² are (1,21/4,105/4). The index is 8 and Sturm bound floor(3·8/12)=2, so equality is exact. Equivalently compare the integral forms 8−7A1+49A2 and 8theta³−6E. The independent checker obtains (8,42,210) and equality through q^96. Therefore f*=theta³(1−3t/4), with precisely the required character and no overlooked cusp pole. The denominators 4 and 8 are units for every inert prime.

## Reproduction and scope boundaries

Run `python certificate.py` and `python -O certificate.py`; both outputs match normal.json and optimized.json byte-for-byte. Run `python verify.py` or its optimized version to validate the entire declared package and reproduce both modes. The portability harness relocates only declared payloads and checks rejection of missing, modified, extra, duplicate, and unsafe-path payloads, in normal and optimized modes. No Python assert statement controls acceptance.

The exceptional certificate and companion coefficients are finite exact calculations made conclusive by the separate modularity and Sturm arguments. The arbitrary-p inputs in §§2–3 are symbolic proofs, not extrapolations. This closes only the dependencies named in the assignment. It does not certify the inherited CM table, all adjacent theorems, any stronger p³ or p^(2r) congruence, or external novelty. The separately audited two-step constant-term argument is not re-audited here.
