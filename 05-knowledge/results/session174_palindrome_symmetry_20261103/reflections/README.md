# S174: palindromic spectra and reflection holonomy

Status: PROVED elementary results; FINITE-EXACT explicitly identified computations. Independently derived, without historical novelty claims. No enormous-integer primality or current-record assertion. This packet was reconstructed after an ephemeral workspace reset; proofs and tests are complete below.

## Reciprocal folding

For a fixed-length base-b digit word d_0,...,d_n (least significant first), let P(X)=Σd_jX^j and P*(X)=X^nP(1/X). Reversal sends P to P*. Palindromy is P=P*. P−P* vanishes at 1, so integer reversal preserves residues modulo b−1. When n is odd, reciprocal P vanishes at−1, implying divisibility by b+1. Hence the only even-length decimal palindromic prime is 11.

For degree 2m, reciprocal P uniquely equals X^mQ(X+X^−1). Define S_0=2, S_1=Y, S_(k+1)=YS_k−S_(k−1); induction gives S_k(X+X^−1)=X^k+X^−k. Then Q=d_m+Σ_(k=1)^m d_(m+k)S_k. Uniqueness follows from distinct degrees and monic leading coefficients. This quotient pairs reciprocal roots, retaining the complete polynomial when m is retained. Numerical evaluation alone loses the action without base and length sidecars.

Examples: 1001 corresponds to X³+1=(X+1)(X²−X+1), and 1001=7·11·13. For 2222, P=2(X+1)(X²+1), giving 2222=2·11·101. For 7417, P=7X³+4X²+X+7 and P−P*=3X(X−1), hence 7417−7147=270. Trial division establishes 7417 is prime, not palindromic. For 240, P=2X²+4X, fixed-length reversal 042=42 gives difference 198, illustrating the leading-zero sidecar.

## The dense-family unit-circle theorem

**Theorem.** For every integer m≥1, real a>0 and |δ|<a, the polynomial F(X)=aΣ_(j=0)^(2m)X^j−δX^m has 2m distinct roots on the unit circle, none equal to±1.

**Proof.** The quotient is aD_m(Y)−δ with D_m=1+Σ_(k=1)^mS_k, monic degree m. At Y=2cosθ,
D_m(2cosθ)=sin((m+1/2)θ)/sin(θ/2).
At θ_0=0 the limiting value is 2m+1; at θ_j=jπ/m for 1≤j≤m it is(−1)^j. Because δ/a lies strictly between−1 and 1, D_m(2cosθ)−δ/a changes sign in each of the m disjoint intervals(θ_(j−1),θ_j). The degree-m quotient therefore has exactly m distinct roots in(−2,2), one in each angular interval. Each lifts through X²−YX+1 to two distinct conjugate unit-circle roots. These account for all 2m roots. ∎

Each conjugate pair has a certified angular bracket of widthπ/m. Strictness is necessary: m=1, δ=−a yields a(X+1)².

With m=944264,a=9,δ=1, evaluation at 10 gives
N=10^1888529−10^944264−1.
Its word is exactly m copies of 9, then 8, then m copies of 9, total 1888529 digits. Its degree 1888528 digit polynomial therefore has all 1888528 roots simple on the unit circle. This is polynomial root geometry, not integer primality.

For k≥4,10^k≡160(mod 240), because modulo 16 it is 0 and modulo 15 it is 10. Consequently N≡239(mod 240). Also N≡8(mod 11), since m is even. Exact modular computation gives N≡6663(mod 7417) and ord_7417(10)=2472, so 7417 is not a divisor of N.

## Belphegor gives the decisive contrasting spectrum

The integer B=10^30+666·10^14+1 has digit polynomial X^30+6X^16+6X^15+6X^14+1. Its reciprocal quotient is S_15(Y)+6Y+6, whose ascending coefficients are
[6,−9,0,140,0,−378,0,450,0,−275,0,90,0,−15,0,1].
Exact rational Sturm computation gives variation counts 8 and 7 at−2 and 2, respectively, with no endpoint roots and constant final gcd. Hence exactly one quotient root lies in(−2,2), so this polynomial has exactly 2 unit-circle roots and 28 off-circle roots. This is FINITE-EXACT and reproduced independently by the script. Palindromy alone therefore does not force the dense family's spectral property. B≡8(mod 11), sharing N's residue without sharing its spectrum.

## Operational bridge to square-sum paths

Over any commutative ring let r_s(x)=s−x. For a label word s_1,...,s_l,
H=r_(s_l)∘...∘r_(s_1), H(x)=(−1)^lx+A,
A=Σ_(j=1)^l(−1)^(l−j)s_j.
Induction proves the formula. Writing P_s(X)=Σ_(j=1)^ls_jX^(j−1) gives A=(−1)^(l−1)P_s(−1). Alternating digit evaluation is precisely the translation coordinate of affine reflection holonomy. The quotient forgets the word beyond parity and A; nonpalindromic words can have the same map.

An even palindromic word of involutive generators cancels from the middle out to the identity; an odd one is a conjugate of its central generator. For an actual path v_0,...,v_l with edge sums s_j=v_(j−1)+v_j, the dynamics are exactly v_j=r_(s_j)(v_(j−1)). Thus a nontrivial simple path cannot have an even palindromic edge-sum word, since H=id forces v_l=v_0. Taking s_j to be squares gives an actual square-sum obstruction, not a verbal analogy. Odd palindromic words can support positive simple paths: [1,15,10,6] has labels[16,25,16].

The Q_15 path inherited from THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md is
[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9].
Its 14 square labels have A=1 and holonomy x↦x+1, exactly matching endpoints 8→9. Therefore its label word is not palindromic. The map preserves actual endpoint equations; it loses internal vertices unless the full word/start are retained. No Hamiltonicity or LRC result beyond this obstruction is asserted.

## Reproduction

Run python test_reflections.py and python -O test_reflections.py. Explicit exceptions survive optimization. Exact Fraction Sturm sequences independently verify 96 quotient spectra(m=1..24,δ/a∈{−8/9,0,1/9,8/9}), with polynomial reconstruction checks. Exhaustive digit tests cover bases 2..7,lengths 1..5; reflection tests cover alphabet 0..3,lengths 1..7. Hostiles include leading zeroes,7417's symmetry defect, a repeated-root boundary and positive odd-palindrome square-sum path.
