# Actual Hamiltonian-path multiplicities: polynomial growth and the 7-adic degree defect

Session 24, 2026-10-21. Status: proofs below plus exact finite controls; independent review still required. No literature-priority claim is made.

## Scope and inheritance

H(T) is the actual number of directed Hamiltonian paths, A(T)=|Aut(T)|, and h(T)=H(T)/A(T). These are not interchangeable. N_n(H) counts ordinary tournament isomorphism classes of order n with H paths; S_n(H) counts the self-converse ones; M_n(H) counts converse-merged classes. Thus M_n(H)=(N_n(H)+S_n(H))/2.

The unique ordered SCC decomposition, multiplicativity of H and A, and converse as reversal with componentwise converse are inherited. So are the SCC-word and square-divisor center architecture from Session 22/23. The existing fixed-h formulas have units K1 and C3 and denominator 1-z-z^3. They are NOT new results here. Fixing actual H changes the unit alphabet to K1 alone, turning exponential coefficient growth into eventual polynomial growth. This packet derives that specialization, the associated arithmetic degree defect, its typical metagraph-degree consequence, and concrete exact multiplicities.

Use the empty tournament only in the auxiliary order-zero H=1 series. All numerical n>=1 statements refer to actual tournaments. Strong means nontrivial strong unless K1 is explicitly included.

Classical ingredients:

* A strong tournament has a Hamiltonian cycle (Camion); its n cuts give n distinct Hamiltonian paths. Hence H(S)>=|S|.
* Busch proves H(S)>=5^((|S|-1)/3). Primary source: Arthur H. Busch, “A Note on the Number of Hamiltonian Paths in Strong Tournaments,” Electronic Journal of Combinatorics 13 (2006), N3, https://doi.org/10.37236/1141 . The publisher abstract explicitly states this bound. It suffices here; no characterization of equality is used.
* H is odd; H=1 exactly for transitive tournaments. The inherited all-order gaps H=7,21 are used, not rediscovered. The repository's THM-1370 explicitly retracts a novelty claim and cites the prior 2016 discussion.

Primary sources for the classical SCC/strong-order context include Moon, “On Subtournaments of a Tournament” (1966), https://doi.org/10.4153/CMB-1966-038-7 . No fresh theorem that the only missing odd values are7 and21 is claimed or assumed.

## 1. Finite cores and exact rational series

For odd q>1 let

    A_q(z) = sum z^|C|, over strong ordinary classes C with H(C)=q.

Let C_q(z) similarly count strong self-converse classes; put C_1^+(z)=1+z (empty or singleton center) and C_q^+=C_q for q>1. All A_q,C_q are finite polynomials: Camion already gives |C|<=q; Busch improves this to |C|<=1+3 log_5 q.

For H>=1 define F_H(z)=sum_n N_n(H)z^n and J_H(z)=sum_n S_n(H)z^n. Then

    F_H(z) = sum_{q_1 ... q_r=H, q_i>1} A_q1(z)...A_qr(z)/(1-z)^(r+1),       (1)
    J_H(z) = sum_{x^2|H} F_x(z^2) C^+_{H/x^2}(z),                          (2)
    M_H(z) = (F_H(z)+J_H(z))/2.                                           (3)

The empty ordered factor tuple occurs only for H=1, giving F_1=J_1=M_1=1/(1-z).

Proof: delete singleton SCCs to leave a finite ordered nontrivial strong core. Insert arbitrary singleton runs in its r+1 gaps. Its H is the product of the core H values. A self-converse SCC word has unique first half U and either no center or one self-converse strong center C; the opposite half is U^op, so H=H(U)^2H(C). The first half is an ordinary class, not a converse orbit. Equation (3) is the involution orbit count. These are specializations of the inherited architecture, not independent novelty claims.

A still more literal coefficient formula is useful. Let B_{H,r,m} be the number of ordered nontrivial strong core words with product H, r components and total size m. Then

    N_n(H) = sum_{r,m} B_{H,r,m} binom(n-m+r,r),                           (4)

where a term is zero for n<m. Each core is an isomorphism-invariant ordered word, so no factorial for rearrangements is silently divided out.

## 2. Exact growth degree and pole; symmetry is lower order

Assume H is attained. Let r(H)=max{r:B_{H,r,m}>0 for some m}, and b_H=sum_m B_{H,r(H),m}. Then the reduced denominator of F_H is exactly

    (1-z)^(r(H)+1),

and

    N_n(H) ~ b_H/r(H)! * n^r(H).                                         (5)

Indeed all numerator contributions at the highest pole are positive at z=1, so no cancellation is possible. Equation (4) is eventually an exact polynomial of degree r(H), with the displayed leading coefficient. A safe elementary onset is n>=H: each core component size is <=its path count and the sum of factors q_i>=3 is <=their product. A sharper effective onset is

    n >= floor(3 log_5 H).                                               (5a)

Indeed the polynomial binomial in (4) is already zero at every integer M-r<=n<M. Thus each term agrees with its polynomial once n>=M-r. Busch gives M-r<=3 log_5 H. Equivalently, after the common denominator (1-z)^(r(H)+1), the numerator has degree at most r(H)+floor(3 log_5 H). This onset sharpening was supplied by the coordinating researcher and checked here. It is an upper bound, not a claim of minimal onset.

For J_H, let t(H) be the maximum r(x) among realized summands of (2). Every such half can be paired to produce a core for H, so 2t(H)<=r(H). Consequently

    S_n(H) = O(n^floor(r(H)/2)).                                         (6)

When J_H is nonzero, its pole at +1 has exact order t(H)+1. Its denominator divides (1-z^2)^(t(H)+1), but cancellation at -1 is possible and DOES occur. Do not infer an exact period-two oscillation or one parity-independent leading constant without checking the center numerator. For H>1, r(H)>=1, so

    M_n(H) ~ (b_H/[2 r(H)!]) n^r(H),
    S_n(H)/N_n(H) -> 0.                                                  (7)

Merged multiplicities are eventually quasipolynomials of period dividing2, and their +1 pole is exactly r(H)+1. If e_- is the actual pole order of J_H at -1, the exact merged denominator is (1-z)^(r(H)+1)(1+z)^e_- for H>1. It can be polynomial rather than genuinely parity-dependent: H=33,41,49 below are examples. Empty J_H is allowed; H=19,35,39,63 have no self-converse classes at any order.

The parity quasipolynomials for S_n(H), and hence for M_n(H), hold for every integer n>=1+floor(3 log_5 H). For a half-core of size M and length r with center size c, the relevant parity binomial extends correctly once n>=2(M-r)+c. Busch bounds this by1+3 log_5 H, including the empty and singleton centers. Consequently the entire all-order fixed-H multiplicity problem reduces effectively to strong tournaments of order at most1+floor(3 log_5 H). Direct labelled enumeration with canonicalization gives the coarse complexity H^{O(log H)} up to bit-operation factors, without implying a speedup for unrestricted A000568. No stronger exact -1 exponent or general symmetry leading constant is claimed.

## 3. The two missing values create an exact two-degree loss

The maximum r(H) is a factorization length over actual strong path counts. It need not equal the arithmetic prime-factor count Omega(H).

For a>=0,

    r(3^a)=a,    b_(3^a)=1,
    N_n(3^a) ~ n^a/a!.                                                  (8)

Every nontrivial factor costs at least one3, and the top-length word is uniquely a copies of C3. Its exact contribution is binom(n-2a,a).

The census and Busch bound below establish that there are exactly two strong classes with H=63 in ALL orders, both at order8, both automorphism-free (Aut order1), and they are a converse pair. Now take a>=2. In any factorization of7*3^a, exactly one factor contains7. Because7 and21 cannot occur, that factor consumes at least two3s. Therefore at most a-1 strong factors occur. Equality is attained by one H63 atom and a-2 C3 atoms. Hence

    r(7*3^a)=a-1,
    top-pole term = 2(a-1) z^(3a+2)/(1-z)^a,
    N_n(7*3^a) ~ [2/(a-2)!] n^(a-1),
    M_n(7*3^a) ~ [1/(a-2)!] n^(a-1).                                    (9)

This loses exactly TWO powers of n compared with the naive Omega(H)=a+1 rule. The gap at7 costs one prime slot; the additional gap at21 forces the7-containing block to consume a second3. This is an all-order implication of the two gaps, not a fit to small tables. For a=0,1 the entire stratum is empty. At H=63 the only possible square-divisor halves in (2) are x=1,3, requiring centers with H=63,7 respectively; both center inventories are empty, proving J_63=0.

Additional exact degree rays, requiring no full-spectrum conjecture:

* r(7*5^b)=b for b>=1. There are exactly two automorphism-free (Aut order1) strong H35 atoms of order7, so the leading ordinary coefficient is 2/(b-1)!.
* r(7^b)=floor(b/2) for b>=2, since no factor7 is possible while actual strong49 and343 witnesses exist. The b=1 stratum is empty. There are exactly eight strong H49 atoms globally, so for b=2m the leading coefficient is 8^m/m!. For odd b the coefficient additionally depends on the total strong H343 inventory, not computed here.

There is a complementary symmetry refinement. If a=2m+1>=3, equation (2) and the gaps force any half to have H=3^t with t<=m-1. The upper bound is attained by the self-converse Paley order7 center of H=189=7*3^3, paired with a C3 half of H=3^(m-1). Thus the self-converse growth degree along odd exponents is exactly m-1=(a-3)/2 in the quasipolynomial sense (at least one parity attains it). For even a=2m, the putative t=m-1 center would have H=63 and is absent; hence the degree is at most m-2 when m>=2, and S_n(63)=0. Equality on the even ray requires additional self-converse center witnesses and is not assumed. These statements do not assert that arbitrary remaining odd integers are attainable.

## 4. Why this changes the typical metagraph node

Here use exactly the inherited fundamental paired-black converse-merged metagraph Q from Statement3 (the sigma-twin theorem) of the current canon `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`. That source calls this graph B; the present Q is the same graph, obtained by identifying each pair of sigma-twin black lines. Its degree is h=H/Aut at a non-self-converse pair and (h-b)/2 at a self-converse class. Loops contribute two to degree, as in the source's multigraph handshake convention. Thus degree-one nodes are genuine leaves with one non-loop edge. This is not a claim about the arc-flip quotient or the literal merged black network with doubled edge multiplicities.

At fixed H and growing n, top-length cores contribute 1-O(1/n) of ordinary classes, because all shorter cores give lower-degree polynomials. The same is true after converse merging.

For H=3^a, a>=1, the dominant cores are a copies of C3. They have Aut=3^a and h=1. Self-converse words are lower-order by (6). Thus a fraction tending to1 of the merged nodes in this fixed-H stratum are degree-one black leaves.

For H=7*3^a, a>=2, the dominant cores contain one automorphism-free (Aut order1) H63 atom and a-2 triangles. They have Aut=3^(a-2), hence h=63. They are NEVER self-converse: the unique H63 component would have to be central and self-converse, but neither H63 atom is. Their exact ordinary count is

    2(a-1) binom(n-2a-3,a-1),

with zero when n<3a+2. Thus a fraction 1-O(1/n) of the merged fixed-H nodes have black degree63. At a=2 this is exact for every n>=8.

The actual-H ray3^a and its7-multiple therefore contrast degree1 versus63 in the dominant metagraph population, simultaneously with the two-degree multiplicity deficit. This is a frequency/degree result, not a new adjacency, connectivity, or edge-placement theorem.

## 5. Complete atom inventory through H=73 and examples

The frozen n<=7 atlases are reused, not re-enumerated. A bounded augmentation of all456 ordinary n7 representatives by the128 orientations of a new vertex yields all6880 ordinary n8 classes. Score-partition canonicalization is exact: isomorphisms preserve degrees, and every permutation within equal-score cells is checked. Hamilton paths use exact integer subset DP. Converse canonicalization determines self-converseness; automorphism order is the number of score-preserving permutations attaining the canonical code.

The n8 result matches EVERY one of the320 bins of the independently frozen full2^28 labelled census after weighting each ordinary class by8!/Aut. Additional totals:6008 strong,176 self-converse,3528 converse-merged; sum H/Aut=2^21. Independent SCC convolution reconstructs every ordinary and self-converse multiplicity through n8; direct permutation enumeration checks all labelled tournaments through n5.

Completeness in H, not merely in n: 5^8>73^3, so Busch rules out H<=73 for any strong tournament of order>=9. Thus the computed strong/strong-self-converse inventories are complete globally for every odd H<=73. In particular the global strong H63 count is exactly2. All formulas in low_H_formulas.json are consequently all-order formulas, not extrapolations.

Selected ordinary multiplicity polynomials (valid at least for the stated n):

| H | N_n(H) | validity | S_n(H) |
|---|---|---|---|
| 1 | 1 | n>=1 | 1 |
| 3 | n-2 | n>=3 | 1 for odd n,0 for even n |
| 5 | n-3 | n>=4 | 1 for even n,0 for odd n |
| 9 | (n-4)(n-1)/2 | n>=5 | floor((n-4)/2)+2*[n odd] |
| 15 | (n-4)(n-3) | n>=5 | 2*[n odd] |
| 19 | 2(n-5) | n>=6 | 0 |
| 25 | (n-5)(n-2)/2 | n>=6 | floor((n-5)/2) |
| 27 | n(n-5)(n-4)/6 | n>=6 | (n-5)/2 for odd n>=7,0 for even n |
| 33 | n^2-5n-1 | n>=6 | 1 |
| 35 | 2(n-6) | n>=7 | 0 |
| 39 | (n-6)(n-5) | n>=7 | 0 |
| 41 | 6n-35 | n>=6 | 1 |
| 49 | 8n-49 | n>=7 | 1 |
| 63 | 2(n-7) | n>=8 | 0 |

H=7,21 have zero multiplicity at every n. Divide N+S by2 for merged counts. In particular M_n(33)=n(n-5)/2, M_n(41)=3n-17, M_n(49)=4n-24, and M_n(63)=n-7 on the stated ranges.

The H39 row is deliberately quadratic: besides its two strong order7 atoms it has the SCC factorization3*13. Merely padding the strong atoms would incorrectly give2(n-6). This is a useful hostile to incomplete counting.

At n8 the modal actual-H value by ordinary class multiplicity is H345:49 ordinary classes,25 merged,1 self-converse. H281 occurs48 ordinary/24 merged times. H63 occurs2 ordinary/1 merged times, while H49 occurs15 ordinary/8 merged times. These numerical observations are finite-exact; no persistent unimodality is claimed.

## 6. Reproducibility and boundary

Files: extract.py reads the original n<=7 atlases and frozen labelled n8 histogram; augment8.cpp produces n8_classes.csv from n7_representatives.txt; analyze.py produces the complete strong inventories, all count tables, and exact rational functions for H<=73; check.py independently reconstructs all census channels by SCC and center convolution and verifies direct Hamilton-path witnesses. receipt.json and optimized-receipt.json are byte-identical with2063 active checks, including under python -O.

The n8 isomorphism census is a bounded new computation; the all-order H<=73 completeness is theorem-backed. Existing SCC/center word identities are inherited. No raw nonintegral division of labelled counts by n! is used. No full odd-spectrum completeness theorem, general unknown atom inventory, or metagraph wiring theorem is asserted. The exact source fingerprints are in census.json. This directory was written outside the repository and no Git operation was performed.
