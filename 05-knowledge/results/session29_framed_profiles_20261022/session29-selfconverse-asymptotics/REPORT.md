# Self-converse census: classical formula, rigorous correction, and mean blue mass

Date: 2026-10-22. Status: proved from the displayed partition formula, with an independent derivation and finite checks. No claim of literature novelty. This is a defined **average**, not a maximum and not the user's unidentified M_k.

## 1. Results

Let b_n be the number of isomorphism classes of self-converse tournaments of order n. Put (m)_r=m(m-1)...(m-r+1), and

A_m = (m)_3 2^(8-4m)/3,
C_m = [(m)_5 2^24/5 + (m)_6 2^28/18] 2^(-8m).

As m tends to infinity,

b_(2m) = 2^[m(m-1)]/m! · [1+A_m+C_m+O(m^9 2^(-12m))],
b_(2m+1) = 2^(m²)/m! · [1+A_m/4+C_m/16+O(m^9 2^(-12m))].

Consequently the leading terms are 2^[m(m-1)]/m! and 2^(m²)/m!. The first corrections have exponentially small relative size. Every omitted partition term is nonnegative.

Use the previously proved blue-tile mass B_n=2^floor((n-1)²/4). The arithmetic mean of blue tiles over the b_n self-converse classes is B_n/b_n. Hence

mean_blue_(2m)/m! = 1-A_m+O(m^6 2^(-8m)),
mean_blue_(2m+1)/m! = 1-A_m/4+O(m^6 2^(-8m)).

For m>=3 these ratios are strictly below 1, since the partition (3,1^(m-3)) contributes positively; they tend to 1. This does not determine the maximum blue fibre or the distribution. It does imply that the maximum is at least this average.

## 2. Exact formula and anti-Burnside proof

For an odd partition λ=(d_1,...,d_l) of m let z_λ=∏_d d^(a_d) a_d!, where a_d is the multiplicity of d. Then

b_(2m) = Σ_(λ odd) 2^[m+2Σ_(i<j)gcd(d_i,d_j)-l]/z_λ,
b_(2m+1) = Σ_(λ odd) 2^[m+2Σ_(i<j)gcd(d_i,d_j)]/z_λ.

Burnside twisted by converse counts self-converse orbits: b_n=(1/n!)Σ_π |Fix(converse·π)|. To see this directly, an isomorphism class fixed by converse contributes exactly n! pairs (T,π) with π:T→T^op; other classes contribute zero.

A permutation can reverse a tournament precisely when its vertex cycles have lengths 2d with d odd, together with at most one fixed vertex. An odd cycle of length>1 has odd-length internal pair orbits without endpoint exchange and is inconsistent with reversal. An even cycle of length r has antipodal pair orbit length r/2 with endpoint exchange, which is consistent precisely when r/2 is odd. Two fixed vertices are impossible. Conversely these conditions make each pair orbit consistent: non-antipodal internal and intercycle orbits have even length and no endpoint exchange; antipodal orbits have odd length and exchange endpoints.

Inside a 2d-cycle there are d pair orbits; between 2d- and 2e-cycles there are 2gcd(d,e); a fixed vertex contributes one orbit per nontrivial cycle. Each pair orbit permits two orientations. The cycle-type denominator is 2^l z_λ. Substitution yields both formulas, including m=0 with b_0=b_1=1.

## 3. Rigorous remainder estimate

Separate the non-unit parts of λ. Their sum is k, their number is l, and define the defect d=k-l=m-length(λ). The defect is even and k<=3d/2. Let S=Σ_(i<j)gcd(d_i,d_j), now only over the non-unit parts.

Relative to the identity partition 1^m, the even exponent loss is

L=2md-[k²-2kl+2k-l+2S].

Since gcd(a,b)<=(a+b)/2, 2S<=(l-1)k. Therefore

L>=d(2m-k-1)>=d(2m-3d/2-1),
and also L>=d(m-1), since k<=m.

The ratio of partition weights is m!/z_λ <= m^k <= m^(3d/2). The number of possible non-unit partitions at defect d is at most 2^(floor(3d/2)+1): sum the crude composition bound over k<=3d/2. These bounds are intentionally elementary, not fitted to census data.

Defect 0 gives 1. Defect 2 consists only of (3) and gives A_m. Defect 4 consists only of (5) and (3,3), giving C_m. For defect 6 there are finitely many possibilities, k<=9 and L>=12m-60; their total is O(m^9 2^(-12m)).

For 8<=d<=m/6, the loss bound gives L>=d(7m/4-1). Let a_m=(2m)^(3/2)2^(-7m/4+1). Summing even d is bounded by summing all integers:

Σ_(8<=d<=m/6) relative terms <= 2Σ_(d>=8) a_m^d <=4a_m^8

for sufficiently large m. This is O(m^12 2^(-14m))=o(m^9 2^(-12m)).

For d>m/6 there are at most 2^m partitions in all, every weight ratio is at most m!, and L>d(m-1)>m(m-1)/6. Thus their total is at most

2^m m! 2^[-m(m-1)/6] = o(m^9 2^(-12m)).

This proves the even remainder. For each fixed λ, its odd-order relative term is its even-order relative term multiplied by 2^(-d): the odd exponent gains length(λ), whereas the identity gains m. Thus the defect-2 and defect-4 corrections are divided by 4 and 16, and the same remainder bound applies. Inverting 1+A_m+O(m^6 2^(-8m)) proves the stated average formula, since A_m²=O(m^6 2^(-8m)).

## 4. Provenance, distinguishing recovery from new deduction

Classical primary source inspected:
- W. J. R. Eplett, Self-converse tournaments, Canadian Mathematical Bulletin 22(1) (1979), 23–27, DOI https://doi.org/10.4153/CMB-1979-004-6 . Corollary 1 gives the odd-order odd-partition enumeration. The following page cites Harary–Palmer, Graphical Enumeration (1973), p.156 for even order. Eplett also proves conjugacy of reversing involutions and the unique representation used in the fixed-involution approach.
- F. Harary and E. M. Palmer, Enumeration of self-converse digraphs, Mathematika 13 (1966), 151–157, https://doi.org/10.1112/S0025579300003910 . Bibliographic provenance only here; the present derivation does not rely on an unread formula from it.
- OEIS A002785 https://oeis.org/A002785 identifies self-converse tournaments and displays both formulas and counts. Used for numerical cross-check, not proof.

Earlier repository work already contains the exact formula:
- 01-canon/theorems/THM-283-sc-burnside-formula.md (corrected allowed cycles 2 mod4).
- HYP-2074, historical index 05-knowledge/hypotheses/INDEX-HISTORICAL-THROUGH-2026-07-21.md, June 2 S562: both bisections, including even b_(2m)=Σ4^(c2(λ))/z_λ and odd extra factor2^length(λ). That entry explicitly disclaims OEIS novelty.
- 05-knowledge/results/converse_tournament_count_channels_20261021, Session22 independent formula audit and historical-error correction.
- 01-canon/theorems/THM-643-blue-black-line-parity-structure.md proves total blue mass and small self-converse counts.

Accordingly the contribution of this packet is a self-contained rigorous asymptotic expansion and its mean-blue consequence, plus a portable independent check. No priority claim is made for the expansion; a bounded provenance search did not establish its first publication.

## 5. Reproduction and scope

Run `python check.py`. Python standard library only. Uses odd partitions rather than tournament enumeration; integer permutation weights m!/z avoid rational accumulation. receipt.json records:
- exact census n=0..61;
- agreement with OEIS displayed n=1..23;
- an independent anti-Burnside calculation over all permutations for n=1..7, propagating edge orientation parity without assuming the cycle-type formula;
- checks of the loss inequalities for every enumerated partition m<=30;
- exact defect-2 and defect-4 coefficients m=6..30.

Finite checks support implementation only; Section3, not the data, proves the asymptotic remainder.
