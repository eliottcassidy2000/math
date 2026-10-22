# Independent audit: self-converse tournament census and asymptotics

Date: 2026-10-22. Verdict: **ACCEPTED**, including the all-partition remainder. Audited producer packet: `session29-selfconverse-asymptotics`; exact source hashes are recorded in `receipt.json`. No repository or Git mutation was performed by this auditor. This is a proof and implementation audit, not a claim of literature novelty.

## 1. Exact formula and anti-permutation criterion

The twisted Burnside count is correct: on a self-converse isomorphism class, the number of labelled tournaments times the number of reversing isomorphisms for each is exactly n!. Other classes contribute zero.

For a vertex permutation π, let an unordered pair return after t steps. Consistency requires π^t to swap its endpoints if and only if t is odd. This criterion independently establishes that allowable vertex cycles are all congruent to 2 modulo 4, with at most one singleton. An odd cycle of length greater than one gives an inconsistent internal pair orbit; a cycle divisible by four gives an inconsistent antipodal pair orbit; two singletons give an inconsistent fixed pair. The remaining types satisfy the criterion on every pair orbit.

Writing each nonsingleton cycle as 2a, with a odd, gives internal pair-orbit count a, intercycle count 2 gcd(a,b), and one extra orbit per nonsingleton cycle if a singleton is present. The denominator is 2^ell z_lambda. This yields the producer's two formulas exactly. The edge-free cases n=0,1 are also valid by direct twisted Burnside, without requiring the converse map to be nontrivial.

## 2. Denominator and partition bounds are sufficient

Write nonunit odd parts as a_1,...,a_l; k=sum a_i and d=k-l. Then d is even and k<=3d/2. With S=sum_(i<j) gcd(a_i,a_j), the exact even exponent loss is

L=2md-[k^2-2kl+2k-l+2S].

Using 2S<=(l-1)k gives L>=d(2m-k-1). Since k<=m, it also gives L>=d(m-1).

There is no lost factorial in the weight estimate. Precisely,

m!/z_lambda = (m)_k / product_(a>=3)(a^(multiplicity(a)) multiplicity(a)!) <= (m)_k <= m^k.

The number of nonunit partitions at a given d is bounded above by the number of all compositions of every k<=floor(3d/2), hence by 2^(floor(3d/2)+1). This is deliberately loose but valid; unit parts are determined by k and m and do not create additional choices.

## 3. Uniform remainder, not merely fixed-defect reasoning

Defects 0,2,4 produce 1,A_m,C_m as claimed. At defect 6 there are only (7), (5,3), and (3,3,3), with k<=9 and L>=12m-60. Thus this layer is O(m^9 2^(-12m)).

For 8<=d<=m/6, k<=m/4, so L>=d(7m/4-1). Combining the previous bounds gives at most 2a_m^d per defect, where a_m=(2m)^(3/2) 2^(-7m/4+1). Eventually a_m<=1/2, and summing all integers d>=8 overbounds the even defects by 4a_m^8. Dividing this by m^9 2^(-12m) gives a constant times m^3 2^(-2m), tending to zero.

For d>m/6, overbound the entire partition family by 2^m partitions, each weight ratio by m!, and each term by 2^(-m(m-1)/6). The base-2 logarithm of the resulting bound divided by m^9 2^(-12m) is at most

-m^2/6 + m log_2(m) + (13+1/6)m - 9 log_2(m),

which tends to minus infinity quadratically. Thus this entire tail is negligible too. Neither split requires d fixed. There is no missing intermediate range.

The odd/even normalized term ratio is exactly 2^(-d), so the coefficient divisors 4 and 16 and the odd remainder follow term by term.

## 4. Additional exact defect-six coefficient and sharp error scale

Independently extracting the next layer gives

D_m = 2^(-12m)[2^48 (m)_7/7 + 2^48 (m)_8/15 + 2^60 (m)_9/162].

The terms correspond respectively to (7), (5,3), and (3,3,3). These coefficients are checked exactly, not fitted. Since every defect d>=8 contributes o(m^9 2^(-12m)), the remainder after 1+A+C is asymptotic to

(2^60/162) m^9 2^(-12m)

in even order and to 1/64 of this expression in odd order. Consequently the producer's big-O scale is sharp up to its constant. This additional conclusion follows from the accepted uniform tail proof, not from the finite computation.

## 5. Average-blue interpretation and reciprocal expansion

The fixed triangular-grid reflection has floor((n-1)^2/4) free binary coordinates. Its fixed tilings are exactly the blue tilings, all on self-converse classes. Therefore their total mass is B_n=2^floor((n-1)^2/4), and the arithmetic mean over self-converse classes is B_n/b_n.

For n=2m and 2m+1 this is exactly m! divided by the corresponding normalized census. It is not a maximum, not an average over labelled tournaments, and not an identification of the user's M_k. The mean/m! is strictly less than one for m>=3 and tends to one from below. The producer's first-correction reciprocal formula is valid because A^2=O(m^6 2^(-8m)). Any finer reciprocal expansion must include A^2 at the same exponential scale as C; omitting it would be incorrect.

The existing blue-mass expression is intended for n>=1. Do not use floor((n-1)^2/4) at n=0 to infer a tiling model merely because the census formula includes b_0=1.

## 6. Independent computational controls

Run `python audit.py` or `python -O audit.py`. No producer module is imported. The independent checker:

- computes rational odd-partition formulas through n=71 and matches every producer census value n=0..61;
- checks all unrestricted vertex-cycle partitions through n=30, constructing actual permutations and following ordered endpoints until their unordered pair returns, without assuming the allowable cycle filter;
- independently sums every vertex permutation for n=0..8;
- checks exponent, denominator, defect-count and odd-transfer inequalities on every odd partition through m=35;
- checks exact defect-2, defect-4 and defect-6 coefficients;
- uses explicit RuntimeError guards, not Python assertions.

The normal and optimized receipts are byte-identical. These finite checks validate implementation and coefficients; the uniform proof above establishes the asymptotic theorem.

## 7. Provenance audit and source scope

The exact partition enumeration is classical, and already present in corrected repository THM-283. THM-643 supplies the project's blue-mass interpretation. These claims must not be presented as newly discovered enumeration formulas.

Eplett's original 1979 article, *Self-Converse Tournaments*, Canadian Mathematical Bulletin 22(1), 23–27, was independently opened at the publisher PDF. Page 24 Corollary 1 specifies an odd-partition enumeration for odd order. Page 25 explicitly attributes the even-order computation to Harary–Palmer, *Graphical Enumeration*, p.156. The same page gives census values n=5..11 matching our output. The PDF text extraction partly mangles the displayed formula, so exact formula equivalence is additionally checked through independent derivation and OEIS's readable formulas. No direct inspection of Harary–Palmer p.156 is claimed in this audit.

Primary PDF: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/ED7BC2A1F6FF2D7F31547D771CE69B80/S000843950005743Xa.pdf/selfconverse_tournaments.pdf

DOI: https://doi.org/10.4153/CMB-1979-004-6

OEIS cross-check: https://oeis.org/A002785

McKay's independent census data page agrees for n=2..13: https://users.cecs.anu.edu.au/~bdm/data/digraphs.html

No novelty claim for the asymptotic expansion is established by this audit. The producer's explicit distinction between recovered classical enumeration and new-to-this-session deduction is appropriate.
