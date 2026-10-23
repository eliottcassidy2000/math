# Fixed root-shell genus cannot determine H, even at fixed strong scores

Verdict: PASS. This addendum audits the all-order pigeonhole argument, its stronger Moon-chain construction, its eventual-every-order version, and an explicit small genus-two hostile. It does not concern full unbounded-norm theta coefficients.

## 1. Signature capacity

Fix g≥1. For j ordered norm-two roots, each off-diagonal Gram entry lies in {-2,-1,0,1,2}, and each alternating entry lies in {-4,...,4}. Thus there are at most 45^(j(j−1)/2) possible joint bins. Put

D=Σ_(j=1)^g45^(j(j−1)/2).

The number of roots is n(n−1), so every genus≤g coefficient is at most n^(2g). At fixed n, at most (n^(2g)+1)^D complete root-shell signatures can occur. This loose bound is sufficient; it counts many impossible bins and signatures but never undercounts the capacity.

## 2. Independent H coordinates from Moon chains

The full product proof in THM-4744 was inspected, including both exceptional cases in its path-merge argument. For a composition with all parts≥2, IC(parts) is strong, has

H=∏_(part a)(2^(a−1)+1),

and has the same labelled outdegree vector at fixed total part length N: d_0=N−1, d_i=N−i for1≤i≤N−1, and d_N=1. The product law is inherited Moon mathematics, not a Session36 claim of novelty.

Choose k distinct odd primes p_i>3, and let b_i=(2^(p_i)+1)/3. These b_i are >1, pairwise coprime, and not divisible by3. For completeness: if d divides both 2^p+1 and2^q+1, the order of2 modulo d divides gcd(2p,2q)=2, hence d divides3. Also p>3 prime means p≡1 or5 mod6; direct reduction modulo9 gives 2^p+1≡3 or6, so the 3-adic valuation is exactly1.

A part p_i+1 contributes3b_i. Replacing it by (p_i+1)/2 parts2 preserves total part length but contributes only a power of3. Take B such slots for each i and choose m_i∈{0,...,B} slots to be large parts. All resulting chains have the same order and labelled score vector. Their H-values are distinct because a prime divisor of b_i isolates the exponent m_i. No unique factorization of a complicated unknown number is required for this proof.

Let L=Σ_i(p_i+1). At n=1+BL there are (B+1)^k distinct H-values. Choose k>2gD and B with

B^(k−2gD)>2^D(L+1)^(2gD).

Since n≤B(L+1), the number of signatures is strictly less than B^k<(B+1)^k. Some two strong same-score chains have identical genus≤g root-shell data but different H.

## 3. Every sufficiently large order

For arbitrary sufficiently large n, set B=floor((n−1)/L)−1. The remaining part length r=(n−1)−BL lies in [L,2L−1]. Since L≥6, r has a composition into2s and3s. Append this common padding to every member of the family. It changes H by one common factor and preserves the distinction of the b_i exponents.

For B≥1, n≤3BL. The sufficient inequality is now

B^(k−2gD)>2^D(3L)^(2gD).

It eventually holds because k>2gD. Thus the collision exists at every sufficiently large order, still inside the strong fixed-labelled-score family.

## 4. Complete genus-two root data is exactly (c3,c4)

The conference moment proof generalizes. Write A=n(n−1)(n−2), and let c3,c4 be the numbers of directed simple3- and4-cycles. For ordered norm-two roots, each g=±1 has N(w=±3)=6c3 and N(w=±1)=A−6c3. At g=0:

N(w=±4)=4c4,
N(w=±2)=A(n−3)/6+8c3(n−3)−16c4,
N(w=0)=2A(n−3)/3−16c3(n−3)+24c4.

Each g=±2 has n(n−1) occurrences at w=0. These formulas follow from the same tight-frame moment calculation, using ||M1||²=(n³−n)/3−8c3. Conversely the ±3 and±4 bins recover c3 andc4, so this is an exact characterization at fixed n.

Independent verification checks every tournament through n5 and160 deterministic samples of orders6..10 in `independent_rootpair_formula.py`.

## 5. Explicit strong fixed-score genus-two hostile

IC(2,2,4) and IC(2,3,3) have order9, identical labelled outdegrees, c3=7 andc4=10, and identical complete root-pair distributions. Their Hamiltonian counts are81 and75 respectively. This is an explicit hostile, with no global minimality claim.

`independent_interval_chains.py` directly checks all232 chains of orders3..13: H by independent path recursion, strong connectivity by Tarjan, labelled degrees, and product law. It finds26 same-c4/different-H pair witnesses and explicitly compares their root-pair histograms. It also checks the four order15 prime-replacement examples with H=2187,2673,3483,4257, and pairwise coprimality for every chosen prime below100.

Normal and optimized outputs are identical. The finite tests verify constructions; the unbounded-genus impossibility is proved by the preceding inequalities, not inferred from these small examples.
