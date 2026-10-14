# The full signed Collatz coupling at every ternary level

Date: 2026-10-14. Status: PROVED, with an independent audit of the all-level matrix theorem and exact finite checks. New derivation in this session; external novelty is unclaimed. The result concerns the two accelerated outputs at the same source and their laws conditioned on native path prefixes, not an independent temporal transition model.

## 1. Result and inherited interface

For positive odd n let

    a=v₂(3n+1), b=v₂(3n−1),
    u=U₊(n)=(3n+1)/2ᵃ, v=U₋(n)=(3n−1)/2ᵇ.

The [signed-pairing audit](signed_collatz_pairing_and_exponent_bias_20261014.md) proves

    min(a,b)=1, max(a,b)≥2,  2ᵃu−2ᵇv=2,
    density(a=k,b=1)=density(a=1,b=k)=2^(−k), k≥2.

It also proves that the output pair modulo 3 has exactly three possibilities, (1,1),(2,1),(2,2), each of density 1/3. The present result retains the higher ternary digits jointly, rather than only their separate marginals.

Fix j≥1, q=3ʲ, t=3^(j−1) and d=2t. Let Pⱼ(u,v) be the relative natural density among odd n of the pair (U₊(n),U₋(n))≡(u,v) mod q, with rows and columns indexed by units modulo q. Then:

1. Pⱼ=Mⱼ/[2t(4ᵗ−1)] for the explicit nonnegative integer matrix below.
2. Pⱼ(u,v)>0 exactly when it is not the case that u≡1 and v≡2 mod 3. Thus its support has 3t² of the 4t² possible unit pairs.
3. Every Mⱼ is invertible, and

       |det Mⱼ|=2^(2t)(4ᵗ−1)^(2t−2),
       |det Pⱼ|=1/[t^(2t)(4ᵗ−1)²].

4. The matrices project coherently from q=3^(j+1) to q=3ʲ.
5. Their inverse-limit joint probability measure has full topological support on the three allowed mod-3 rectangles, but is singular with respect to product Haar measure. It is supported by a countable union of affine graphs.
6. Conditioning the source modulo 3ʰ gives exact rank 2·3^(j−h−1), for 0≤h≤j−1. A positive mixture over m such source classes has rank m times this value.
7. Conditioning on an actual fixed-sign valuation prefix of length L gives endpoint paired-output rank 2·3^max(j−L−1,0). A mixture over native words counts the distinct surviving ternary affine-carry classes by the same rank formula.

The finite-level full rank does not imply independence, absolute continuity in the limit, an orbit mixing theorem, or convergence of either fixed-sign map.

## 2. Exact finite compiler

For k≥2 there are two valuation branches:

    Aₖ: (a,b)=(k,1),  v=2^(k−1)u−1,  u≡(−1)ᵏ mod 3;
    Bₖ: (a,b)=(1,k),  u=2^(k−1)v+1,  v≡−(−1)ᵏ mod 3.

Each branch has total relative source density 2^(−k). Conditioned on that branch, the free output coordinate is uniform over its t allowed residue classes modulo q. This follows by fixing the exact dyadic valuation cylinder and using the Chinese remainder theorem for the source modulo 3^(j−1). Each point of its affine graph modulo q therefore has mass 2^(−k)/t.

The element 2 has exact order d=2·3^(j−1) modulo 3ʲ. One elementary proof is v₃(4^(3ʳ)−1)=r+1, proved by induction using the binomial theorem. Hence 4 has order t and generates H={h:h≡1 mod 3}; because 2≡−1 mod 3, it has order 2t in the full unit group.

Since d is even, both the graph and its required parity class depend only on k mod d. Choose representatives k=2,…,d+1. Summing the geometric tail in each class gives

    Σ_{ℓ≥0} 2^(−k−ℓd)/t
      = 2^(d+1−k) / [2t(2ᵈ−1)].

Thus set D=2t(2ᵈ−1)=2t(4ᵗ−1), initially set Mⱼ to zero, and for each k=2,…,d+1 put cₖ=2^(d+1−k). Perform these additions modulo q:

- For each u≡(−1)ᵏ mod 3, add cₖ to Mⱼ(u,2^(k−1)u−1).
- For each v≡−(−1)ᵏ mod 3, add cₖ to Mⱼ(2^(k−1)v+1,v).

There are t entries contributed by each branch and Σ_{k=2}^{d+1}cₖ=2ᵈ−1. Therefore all entries sum to D. A repeated cell receives both contributions; the two branch mechanisms must not be treated as disjoint output cells.

These probabilities are natural densities, not just formal geometric weights. Finite valuation truncations are finite unions of arithmetic progressions, and the omitted sources with max(a,b)>K have total relative density 2^(1−K), tending to zero.

## 3. Exact support and marginals

Every Aₖ point has v≡1 mod 3, since b=1. Every Bₖ point has u≡2 mod 3, since a=1. This excludes the rectangle (u≡1,v≡2).

Conversely, suppose v≡1 mod 3 and u is any unit. Then (v+1)/u is a unit modulo q. Because 2 generates all units, it equals 2^(k−1) for a unique representative k∈{2,…,d+1}. Reducing the relation modulo 3 forces u≡(−1)ᵏ. Thus this pair occurs in Aₖ. Likewise every pair with u≡2 mod 3 and any unit v occurs in one Bₖ. These are exactly the three allowed rectangles.

The row marginal is 1/q for u≡1 mod 3 and 2/q for u≡2 mod 3. The column marginal is 2/q for v≡1 mod 3 and 1/q for v≡2 mod 3. This rederives the preceding audit's all-level marginals, but the joint entries generally are not uniform within a mod-3 rectangle and are not products of the marginals.

For q=9, in the common row/column order (1,2,4,5,7,8), the exact matrix is P₂=M₂/378, where

    M₂ = [32  0  2  0  8  0]
         [ 2  2  8 32 32  8]
         [ 2  0  8  0 32  0]
         [32 32  2  8  8  2]
         [ 8  0 32  0  2  0]
         [ 8  8 32  2  2 32].

For example P₂(1,1)=32/378, whereas the product of the marginals would be (1/9)(2/9)=2/81. The distinction is measurable already at the first digit: the product law incorrectly assigns positive mass to the forbidden (1,2) rectangle.

## 4. Full rank and determinant at every level

Order rows and columns by H={1 mod 3} first, then K={2 mod 3}. Each has size t. The matrix takes block form

    Mⱼ = [ A  0 ]
         [ B  C ].

Only the even high valuations k=2r, 1≤r≤t, contribute to A and C. Their weights are

    c₂ᵣ=2·4^(t−r).

For A, the source branch is a=2r,b=1 and

    v+1=2·4^(r−1)u.

The column change h=(v+1)/2 is a bijection H→H. Its inverse is v=2h−1, also in H. Write both h and u in the cyclic coordinates of H generated by 4. The relation becomes h/u=4^(r−1). After these explicit independent row/column permutations, A is the weighted cyclic-convolution matrix

    W=2Σ_{ℓ=0}^{t−1}4^(t−1−ℓ)R^ℓ,

where R is a cyclic permutation matrix of order t. Direct multiplication gives the especially simple inverse identity

    W(4I−R)=2(4ᵗ−1)I.

Equivalently, in these reordered coordinates W has eigenvalues

    2(4ᵗ−1)/(4−ζ),  ζᵗ=1.

All are nonzero. Since Π_{ζᵗ=1}(4−ζ)=4ᵗ−1,

    |det A|=det W=2ᵗ(4ᵗ−1)^(t−1).

The formula concerns W after coordinate permutations; it need not be the eigenvalue list of A in its original numerical residue ordering.

For C, put h=u−1 and g=v/2. These are bijections from the respective K coordinates to H. The relation u−1=2·4^(r−1)v becomes h/g=4ʳ. Its weights are the same c₂ᵣ. Compared with the previous convolution this merely shifts the cyclic index, and may transpose it depending on the chosen row convention. These permutations preserve absolute determinant. Hence

    |det C|=2ᵗ(4ᵗ−1)^(t−1).

Block triangularity now proves

    |det Mⱼ|=|det A det C|=2^(2t)(4ᵗ−1)^(2t−2).

Dividing by D^(2t) yields the stated normalized determinant. This result is about a coupling matrix, not the transition matrix governing U₊(n)→U₊²(n), and no dynamical invertibility is inferred.

## 5. Coherence and the singular inverse limit

The probability at level j is defined by the same exact input events as at level j+1, with nine finer pairs above each coarse pair. Summing the nine probabilities therefore recovers the coarse probability. The finite compiler also satisfies this identity, because its geometric-series derivation is exact.

A direct construction of the limiting joint measure avoids any ambiguity about the limiting weights. For each k≥2:

- With total probability 2^(−k), choose u with normalized Haar distribution on u≡(−1)ᵏ mod 3 and set v=2^(k−1)u−1 in the 3-adic integers.
- With total probability 2^(−k), choose v with normalized Haar distribution on v≡−(−1)ᵏ mod 3 and set u=2^(k−1)v+1.

The total weight is 2Σ_{k≥2}2^(−k)=1. Its finite projections are exactly Pⱼ.

For any fixed k and orientation, its affine graph has only t points modulo 3ʲ. The normalized product Haar measure of the corresponding cylinder cover is t/(2t)²=1/(4t), tending to zero. Therefore each graph has product Haar measure zero. A countable union of these graphs has product Haar measure zero, but carries the entire coupling measure. The measure is singular.

At the same time, every cylinder inside any of the three allowed mod-3 rectangles has positive probability by the support proof in section 3. Its topological support is exactly their union. Thus all finite matrices can have full rank, all admissible finite cells can have positive probability, and the infinite joint law can still be singular. Retaining the valuation label explains this apparent tension: each label selects a one-dimensional affine graph, while arbitrarily high labels are required to see all finer cells.

No statement about the empirical distribution along one integer orbit is used here. The limit is the coherent residue law obtained by sampling starting integers and then taking their residues.

## 6. Source conditioning: exact rank budgets

The companion [Mullin–Collatz native-gate investigation](mullin_collatz_native_gate_bridge_20261014.md) suggested the following sharper use of the full-rank theorem: a native arithmetic gate can preserve all dyadic valuation options while restricting the ternary source coordinate. This produces an exact rank loss, rather than an analogy based only on residue labels.

Let 0≤h≤j−1, p=3ʰ and r=t/p. Condition the source on n≡s mod p. Write Pⱼ,ₛ for its paired-output probability matrix. In the Aₖ branch, b=1, so

    v=(3n−1)/2 ≡ (3s−1)/2 mod 3^(h+1).

In the Bₖ branch, a=1, so

    u=(3n+1)/2 ≡ (3s+1)/2 mod 3^(h+1).

The inverses of 2 are taken modulo 3^(h+1). Thus Pⱼ,ₛ is supported on r selected columns plus r selected rows, giving rank at most 2r. The source class is independent of the dyadic valuation cylinder, and the two branch contributions must be retained separately where their row and column strips overlap.

Uniform sources modulo p give

    Pⱼ = (1/p)Σₛ Pⱼ,ₛ.

Since rank Pⱼ=2t and Σₛ 2r=2t, subadditivity of rank forces rank Pⱼ,ₛ=2r for every s. The same argument with a subset and its complement proves that any unweighted nonempty m-class sum has rank 2mr.

In fact this holds for arbitrary nonzero weights. Choose rank factorizations Pⱼ,ₛ=UₛVₛ with inner dimension 2r. Concatenating the Uₛ as columns and stacking the Vₛ as rows gives square 2t-dimensional matrices U,V. Their product ΣₛPⱼ,ₛ=pPⱼ is invertible, so both U and V are invertible. For any coefficients αₛ,

    Σₛ αₛPⱼ,ₛ = U diag(αₛI₂ᵣ) V.

Its rank is exactly 2r times the number of nonzero αₛ. In particular, any probability mixture supported on m source classes has rank 2mr, regardless of its positive weights. This is a direct-rank decomposition, not a consequence of probabilistic independence between the two outputs.

## 7. Native path composition: ternary carry controls rank

Fix a sign ε and an ordered, native valuation prefix w=(a₀,…,a_{L−1}). Put A=Σaᵢ and let B be its usual affine carry, so its actual endpoint is

    nL=(3ᴸn₀+εB)/2ᴬ.

The initial sources with precisely that prefix are one odd residue class modulo 2^(A+1). The native-cylinder statement follows inductively: after an existing prefix, sources n₀=r+2^(A+1)z have endpoints c+2·3ᴸz; specifying the next positive valuation fixes exactly the required additional dyadic bits, because 3 is odd. This preserves the actual ordered word, not only the final affine fixed-point equation.

Consequently, as the sources run over that native cylinder, the endpoints run over the progression

    nL=c+2·3ᴸz.

Any initial finite cutoff needed for positivity has no effect on the density law. The endpoint's ternary class is

    nL≡r_w:=εB·2^(−A) mod 3ᴸ,

while its odd dyadic tail and its higher ternary digits remain uniform in the exact CRT sense. This is an assertion about the ensemble of actual sources following the given prefix. It does not randomize the future of an individual source.

At the observation level q=3ʲ, set h=min(L,j−1). Only the endpoint class modulo 3ʰ is relevant to the conditional paired-output law. Applying section 6 gives the exact native-path rank

    rank Pⱼ conditioned on w = 2·3^max(j−L−1,0).

In particular it is rank 2 when j≤L+1. This is a genuine finite-path composition statement: the changing conditional law remembers the ordered affine carry. Iterating the original unconditioned full-rank Pⱼ as if it were an independent transition kernel would miss the entire effect.

For a collection of native words of the same length L and arbitrary positive mixture weights, group the words by r_w mod 3ʰ, where h=min(L,j−1). All words in the same class have the same endpoint paired law at level j. If R distinct classes remain, then

    rank of their mixture = 2·3^(j−h−1)R.

Thus this rank is an exact counter of distinct visible ternary carry classes. At j=L+1 the count is simply half the rank.

There is also a pointwise decoding sidecar valid on every actual positive orbit, without an ensemble assumption. If u=U₊(n)>v=U₋(n), then a=1 and n=(2u−1)/3. If u<v, then b=1 and n=(2v+1)/3. The tie occurs only at n=1. Hence the pair modulo 3^(j+1), together with the real-order comparison bit, recovers n modulo 3ʲ. The mixed residue cell (2,1) is precisely where the two high-valuation orientations alias if the comparison bit is dropped.

**Forward research question, OPEN.** For source-preserving native prefix families selected by a no-descent condition, how many distinct carry classes r_w survive, and are there quantitative collisions beyond what ordinary word counts and dyadic cylinder bounds already force? The rank identity supplies an exact observable for this question. It does not assert a deficit, and a deficit in a finite quotient would still need a source-preserving lift to imply actual descent. Arbitrarily long Mersenne valuation-one runs already rule out a naive uniform-time equidistribution hypothesis.

## 8. Independent exact verification

The checker uses four genuinely distinct checks:

1. Compile Mⱼ from the finite affine-line/geometric formula.
2. Independently enumerate exact source CRT cylinders for each valuation and sign; compute both accelerated outputs as integers and compare their residues with the compiler.
3. Compute determinants using fraction-free Bareiss elimination, checking every exact division.
4. Project the level-j+1 matrix and compare its rational probabilities with level j, by cross multiplication.

All 8,064 checks pass for j=1,2,3 in ordinary Python and with optimization enabled. The two JSON outputs are byte-identical. These include every source class at each h≤j−1, a proper two-class mixture, and fourteen native-prefix examples on the two signs, with actual integer replay of each word and its endpoint progression. The exact small summaries are:

| q | Matrix size | Integer denominator D | Support size | det M |
|---:|---:|---:|---:|---:|
| 3 | 2×2 | 6 | 3 | 4 |
| 9 | 6×6 | 378 | 27 | 1,008,189,504 |
| 27 | 18×18 | 4,718,574 | 243 | 2¹⁸·262143¹⁶ |

The checker saves every matrix entry and the complete determinant at q=27. The all-j proofs are in sections 2–7; a three-level check is not presented as the proof of an infinite statement. A separate independent reviewer read the matrix, density, coherence and singular-limit proofs and found no substantive gap; its independent checker passed 12,039 conditions through j=4 in normal and optimized Python. The later native-path addition is justified explicitly in sections 6–7 and checked by the current 8,064-check script.

Reproduce from the repository root:

    python3 04-computation/experiments/signed_collatz_joint_tower_exact_check_20261014.py > 05-knowledge/results/signed_collatz_joint_tower_exact_check_20261014.json
    python3 -O 04-computation/experiments/signed_collatz_joint_tower_exact_check_20261014.py > 05-knowledge/results/signed_collatz_joint_tower_exact_check_20261014.optimized.json
    cmp 05-knowledge/results/signed_collatz_joint_tower_exact_check_20261014.json 05-knowledge/results/signed_collatz_joint_tower_exact_check_20261014.optimized.json

Canonical output-payload SHA-256, excluding its hash field:

    669bd70246e43674c4b4e8ca085f8b48cb180904f17f84bc0cbf5de5a04b0fde

## 9. Scope and next boundary

This extends the same-source signed pairing to a complete exact family of joint residue laws. The closest inherited mechanism is the valuation-labelled affine identity and CRT source-density calculation. The hostile is the false independent-output product law, which already puts positive mass in the forbidden block. The needed sidecar is the high-valuation branch label, retained in the affine-graph mixture.

The common character calculation with prime powers remains valid, but no prime-race zero distribution or primitive-orbit count has been transferred to Collatz. No classification of fixed-sign cycles, no proof of either fixed-sign convergence statement, and no claim of novelty is made.
