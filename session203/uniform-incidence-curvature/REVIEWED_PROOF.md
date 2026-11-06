# Highest curvature of uniformly copied binary incidences

Status: **PROOF CANDIDATE, pending independent audit.** This is an all-order obstruction for a specified uniform-incidence class, not an odd-spectrum obstruction.

## Context class

Take a fixed tournament core A, with k designated sensitive vertices. Adjoin a transitive spine p1 -> ... -> pm. At each sensitive vertex s choose a polarity: either s -> pi iff bit bi=1 for every i, or s -> pi iff bi=0 for every i. Each remaining core vertex is uniform toward the entire spine, either beating all of it or losing to all of it. All core and spine edges are independent of the bits. Let H_A(b1,...,bm) be the Hamilton-path count of the resulting actual tournament.

Every Boolean function has a unique multilinear polynomial over the rational numbers. For a set I of positions, its full mixed difference is the coefficient of the product of their bits after fixing all other bits. The core size and orientations are unrestricted. Copy/complement choices may differ between sensitive vertices, but not between spine positions at one vertex.

## Theorem

Suppose k>=1 and m>=2k. The Boolean polynomial H_A has degree exactly 2k. More precisely, fix any 2k consecutive spine positions as free and fix all remaining bits arbitrarily. Delete the k sensitive vertices, and absorb the fixed spine vertices into a tournament W. Every vertex of W is uniform toward the free consecutive block: original nonsensitive core vertices are uniform by definition; fixed spine vertices lie wholly before or after the block.

For r>=0 let a_r be the following ordinary transitive-module gap coefficient of W. Sum over every ordering of V(W), choose r of its linear gaps to be occupied, require that every unoccupied interior gap is a forward core edge, and require that every occupied gap has its left endpoint beating a module vertex and its right endpoint losing to a module vertex. A missing end endpoint imposes no restriction. Put a_r=0 for r>|W|+1. In particular a_1=H(W+x), where x has the specified uniform module incidences, so a_1 is odd by the standard oddness theorem for tournament Hamilton paths. When W is empty, a_1=1.

The mixed difference on the free block is exactly

    (-1)^k C_k Q_k,
    C_k = k! (2k)! / 2^(k-1)
        = 2 (k!)^2 (2k-1)!!,
    Q_k = sum_(r=1..k) a_r binom(k-1,r-1) 2^(r-1).     (1)

Thus Q_k is a positive odd integer. In particular this difference is nonzero and its exact 2-adic valuation is

    1 + 2 v2(k!).                                      (2)

This valuation does not depend on the core edges, the uniform attachments, the sensitive polarities, or the fixed exterior bits.

## Proof

Express the count as a sum over all vertex permutations of the product of the indicators of its consecutive directed edges. Only edges incident to sensitive vertices depend on bits. Each sensitive vertex has at most two incident edges in a path, so every summand has at most 2k variable factors. Boolean reduction cannot increase degree. Hence the count has degree at most 2k.

A term can contribute to the coefficient of all 2k free bits only if every sensitive vertex is internal in the path and both its neighbors are distinct free spine vertices, with every free vertex used at exactly one of these variable adjacencies. There are exactly 2k available sensitive incidences and 2k distinct required variables. Any missing incidence or repeated free vertex lowers the multilinear degree. Consequently such a path decomposes uniquely into k blocks

    (free vertex, sensitive vertex, free vertex)

and the remaining vertices W. No sensitive-sensitive edge occurs. For either polarity, the product of the two variable edge indicators in a block has top coefficient -1: it is (1-bi)bj or bi(1-bj). The total sign is therefore (-1)^k, independently of all polarities and core edges touching sensitive vertices.

Fix an ordering of W and its r occupied gaps. Each occupied gap contains a nonempty run of blocks. Their empty/occupied permissions are exactly those counted by a_r. There are binom(k-1,r-1) compositions of the ordered k blocks into r nonempty runs, and k! assignments of the distinct sensitive labels to their centers.

The free labels can be assigned to the 2k surrounding slots in (2k)! ways, subject only to the required spine edge between consecutive blocks inside each run. There are k-r such inequalities, each between the last free vertex of one block and the first free vertex of the next. These inequalities involve disjoint pairs of slots. Independently swapping the labels within each such pair reverses precisely that inequality, so exactly a fraction 2^(-(k-r)) of assignments satisfies all inequalities. Uniform W incidences impose no further label restrictions.

The coefficient is consequently

    (-1)^k k! (2k)! sum_(r=1..k)
        a_r binom(k-1,r-1) / 2^(k-r),

which is (1). Since a_1 is odd and all r>=2 terms in Q_k are even, Q_k is odd. Its a_1 term is positive. The displayed factorization of C_k proves (2). Nonvanishing together with the degree upper bound proves exact degree 2k.

## Even-multiplier obstruction

Let P and C be two contexts of the stated class on the same m-bit spine, with respectively k_P and k_C sensitive vertices. Suppose at least one of these counts is positive and m>=2 max(k_P,k_C). Let a be any nonzero even integer and c any constant independent of the spine bits. Then the identity

    H_C(b) = a H_P(b) + c                          (3)

cannot hold on the full Boolean cube.

If k_P differs from k_C, the exact polynomial degrees differ, so (3) is impossible. If the positive counts agree at k, take the mixed difference on any consecutive 2k-block. Formula (1) makes both coefficients C_k times signed odd integers. Multiplication by even a changes their 2-adic valuation, contradicting equality. This also excludes a negative even multiplier. No assertion is made when both sensitive counts are zero.

For fixed contexts, the condition m>=2 max(k_P,k_C) holds at every sufficiently large length. Contexts may also depend on m: the same argument applies at each length meeting this inequality. Arbitrary additional uniform core vertices are allowed and do not alter the conclusion.

## Scope and surviving routes

This theorem extends the familiar one-sensitive-root coefficient -2 times an odd count to any number of uniformly copied/complemented sensitive vertices. It excludes doubling plus a bit-independent correction in this class, not all physical doublings or all target-dependent constructions.

It does not cover position-dependent incidence tables, nonlinear functions of the full bit word, changing spine edges as bits vary, a restricted sublanguage of bit words, or a length-dependent number of sensitive vertices larger than m/2. Odd multipliers are not excluded by the valuation argument. It also does not forbid a controlled, nonconstant correction or equal-profile substitutions with multiplier 1.

## Literal controls

check.cpp independently counts each literal tournament by subset dynamic programming and compares its alternating cube sum with the ordinary core-gap coefficients in (1). NORMAL.txt records 38 fixed-seed contexts with k=1,2,3,4 and 1,520 actual graph counts, with both polarities and varied uniform cores. Their top coefficient constants are 2,24,1080,120960 and their valuations are 1,3,3,7. These are finite controls of the formula, not the proof of the arbitrary-order statement.

    g++ -O3 -std=c++17 check.cpp -o /tmp/uniform_curvature
    /tmp/uniform_curvature > replay.txt
    cmp replay.txt NORMAL.txt

All checks use explicit failure branches, not assertions. The only general dependency beyond the displayed bijection is the standard oddness of tournament Hamilton-path counts.
