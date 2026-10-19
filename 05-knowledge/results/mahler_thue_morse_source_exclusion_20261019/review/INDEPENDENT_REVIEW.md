# Independent review of Thue–Morse Mahler carry exclusion

2026-10-19. **PASS** for the elementary theorem in [the elementary proof](../THUE_MORSE_FINITE_SOURCE_EXCLUSION.md). **PASS** for the saved [the scalar extension](../THUE_MORSE_SCALAR_MAHLER_EXTENSION.md), including the quantitative exponent bound; source hypotheses are recorded below.

## Elementary proof

1. The carry sum C(w) and iterated affine map are correctly normalized. A difference of two sources whose parities agree has its 2-adic valuation reduced by exactly one under M. This establishes the cylinder isometry and the exact first-mismatch valuation.
2. The periodic source denominator 3^p−2^p is odd; 0≤C(w)≤3^p−2^p. The source is indeed rational and nonpositive, and is assigned the intended periodic carry word by the finite congruence.
3. The common prefix `01101` versus the periodic word `011` has length ratio 5/3. After at least one Thue–Morse substitution, every letter image has the same zero/one count. Applying arbitrary nonempty U,V therefore yields exactly lengths 5k and 3k, even when their lengths differ.
4. The prefix-source estimate B+1≤(3/2)^a(A+1) is correct. The chosen inequality gives 0<(3^(3k)−2^(3k))B+C(w)<2^(5k), with strictness supplied by 3^(3k)−2^(3k)<3^(3k). Thus the first mismatch is genuinely inside the common prefix.
5. The bound is O(log(A+1)) for fixed blocks and prefix. No horizon uniform over positive sources is claimed.
6. For U=000,V=100, every pair of successive ones is at least three positions apart. Every real weighted suffix is at most (2/3)/(1−(2/3)^3)=18/19<1. This is an actual strict-safe formal template.

## Scalar-Mahler extension source check

Primary source checked: Bugeaud–Yao, *Hankel determinants, Padé approximations, and irrationality exponents for p-adic numbers*, https://irma.math.unistra.fr/~bugeaud/travaux/p-adic18.pdf . Theorem 3.1 is on printed page 6; its explicit rational-argument extension is on printed page 13 (PDF indices 5 and 12).

For τ(z)=Σt_nz^n, use τ(z)=z/(1−z²)+(1−z)τ(z²). The nonperiodicity of t makes τ nonrational, hence Kronecker supplies infinitely many nonzero Hankel determinants. Set p=2, w=L, r=1, s=3^L, so b=r p^w/s=(2/3)^L. The factors 1−b^(2^(m+1)) and 1−b^(2^m) never vanish. The cited theorem and remark therefore yield transcendence of τ(b). The specific Thue–Morse Hankel ratio-one verification is on printed page 15 (PDF index 14). In the rational-argument remark, its bound specializes to 2w/(2w−log_2 max{2^w,3^w})=2/log_2(4/3), since 3^w<4^w. This quantitative strengthening also passes. No equality with two at this rational argument is inferred.

The equal-length recoding h(t)=ν(t), ν(0)=UV and ν(1)=VU, is valid because t=σ(t). Consequently the image generating function is an affine expression in τ(z^L). If UV≠VU, its coefficient difference is a nonzero integer polynomial with coefficients in {−1,0,1} and leading coefficient ±1. The rational-root theorem excludes the noninteger root z=2/3. Hence the specialization preserves transcendence. If UV=VU, the recoded blocks coincide, so the image is periodic without needing a separate commutation theorem. Finite prefixes and removal of finite prefixes give invertible rational affine source changes and preserve transcendence.

## Independent exact controls

`verify_independent.py` imports no producer code. Normal and optimized runs agree on **68,642 checks**:

- Every ordered pair of nonempty binary words U,V of lengths at most three.
- Substitution depths one through four.
- Positive sources one through forty and four prefix choices.
- Independent closed-form carry sums, direct actual trajectories, and exact first mismatches.
- Finite suffix controls for the concrete safe family.

The all-parameter conclusions rest on the proofs, not extrapolation from these controls. Reproduce with `python3 -B verify_independent.py` and `python3 -B -O verify_independent.py`.

The earlier receipt-scope wording issue was resolved by the author: regenerated receipts explicitly distinguish the finite computations from the cited analytic theorem. The author subsequently added the exact THM-4072 path and Section 1 dependency to both notes and specified loops from initial state 0 back to 0. Those final clarifications have been checked and resolve the inherited-dependency wording concern. The concrete 000/100 safety claim remains fully self-contained.
