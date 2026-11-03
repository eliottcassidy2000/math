# MRB-weighted Hamilton paths and a quadratic-character 145 seed

## Status and scope

Proved elementary weighted partition identity and exact logarithmic identity; finite-exact seven-vertex tournament counts, checked independently by enumerating every permutation. No literature-priority claim or irrationality consequence. The MRB connection is an encoding, not an identity with an unweighted tournament path count. The role of 290=2·145 is motivation only.

## Strong wrapped-transitive carrier

Let n≥1. Take a transitive tournament A={a₁,…,aₙ}, with aᵢ→aⱼ when i<j, and add b,c with A→b→c→A. This tournament is strongly connected. Its Hamiltonian paths are exactly

- a_S,b,c,a_{A\S}, for all S⊆A, with each of the two A subpaths increasing, and
- the additional path c,a₁,…,aₙ,b.

Indeed b must immediately precede c unless b is last, and c must immediately follow b unless c is first. If b is last, then c is first. Weight a path by the product of wᵢ over all aᵢ before b. Then

    Zₙ(w)=∏ᵢ(1+wᵢ)+∏ᵢwᵢ.

In particular at wᵢ=1 the count is 2ⁿ+1. If wᵢ=t, the generating polynomial by the number before b is (1+t)ⁿ+tⁿ, whose value at −1 is (−1)ⁿ. Thus its even/odd path split differs by exactly one. More generally a root-of-unity filter on this polynomial gives exact counts in each residue class of the number before b.

The generalization uses any tournament H on A. Let P(H) be its number of Hamilton paths and Q₂(H) the number of ordered pairs of vertex-disjoint directed paths covering A, with empty paths allowed and each nonempty path ordered. Then the cyclic wrap has P(H)+Q₂(H) Hamilton paths. If H is the transitive join H₁→H₂, then P(H)=P(H₁)P(H₂) and Q₂(H)=Q₂(H₁)Q₂(H₂); each directed path meets the two blocks in that order. Empty paths are essential to these exact formulas.

## MRB specialization

Take n=N−1, index A by i=2,…,N and set f(i)=i^(1/i)−1 and wᵢ=(−1)ⁱf(i). All 1+wᵢ are positive. Write C_N=Σ₂≤i≤N wᵢ. Taking logs yields the exact finite identity

    C_N = log(Zₙ(w)−∏wᵢ) + Σ₂≤i≤N [wᵢ−log(1+wᵢ)].

Since f(i)=O(log i/i), the correction is absolutely convergent by comparison with Σ(log i)²/i²; since C_N converges by alternating-series convergence, the product has a positive finite limit. Therefore the limit identity represents the MRB alternating sum exactly as a logarithm of a corrected weighted Hamiltonian-path partition function plus an absolutely convergent correction. The exceptional path term must be removed. The signs are weights, not negative path counts.

## F₇ Paley modification with exactly 145 paths

On residues 0,…,6 orient i→j iff j−i mod 7 is in {1,2,4}, and reverse 1→2, 2→3, and 5→6. The resulting out-neighbor masks, with bit j representing i→j, are [22,40,82,53,98,5,43]. Endpoint counts from independent path enumeration are [21,24,20,17,18,28,17], summing to 145. The explicit directed cycle 0→1→3→2→4→6→5→0 proves strong connectivity. For flip indicators x,y,z of those three arcs, all eight exact counts are encoded by

    189−30(x+y+z)+30xy+12xz+12yz−8xyz.

This response polynomial is a finite enumeration certificate, not a general theorem about arbitrary flips. The even value 290 cannot itself be a directed Hamilton-path count of a tournament, by Rédei's oddness theorem; that fact does not make 145 the MRB constant or identify it with a constant derived from it.

Run `python certificate.py` for an independent permutation certificate of the wrap identity and a relabelled 145 seed. The quadratic-character version is verified by the companion `paley_flips.py`.

## A marked-cut derivative on an unweighted strong tournament

For n≥2, consider the same unweighted wrap G_n and fix a vertex a_i. Deleting a_i from any Hamilton path in which a_i appears before b is a bijection to the Hamilton paths of G_{n−1}: for each ordinary subset path the deletion removes i from its pre-b set, while the exceptional path maps to the exceptional path. Consequently the marked cut has exactly

    M_n=2^(n−1)+1=P(G_{n−1})

paths. This equality is also ∂Z_n/∂w_i at all w_j=1. For an arbitrary direction v, the derivative of t↦Z_n(1+tv₁,…,1+tvₙ) at zero is therefore M_n Σv_i. For n=N−1 and v_i=(−1)^i(i^(1/i)−1), the normalized derivative is the finite MRB partial sum C_N. The same sum is a directional observable of the **unweighted** path set; it is not itself a path count.

Uniformly choose a directed Hamilton path in G_n and write X_N=Σ v_i 1{a_i before b}. Then

    E X_N = [(2^(n−1)+1)/(2^n+1)] C_N → C/2.

There is one path for each pre-b subset S and one extra path at S=A. Therefore the distribution of the subset S differs from the uniform distribution on all subsets in total variation by exactly (1−2^(−n))/(2^n+1), which is at most 1/(2^n+1). The induced scalar distribution of X_N differs from the corresponding weighted uniform-subset sum by at most this amount; equality need not hold. The finite variance follows from

    E X_N² = [2^(n−2)(C_N²+Σv_i²)+C_N²]/(2^n+1)

and tends to (1/4)Σ_{i≥2}f(i)². These facts show exactly what source data the marked-path statistic preserves, and they do not imply that the MRB constant is the unweighted path count of any tournament.

## Exact integer marked-count approximation

Let p_i=floor(Q f(i)), Q a positive integer. Compute p_i with integer arithmetic only by the equivalent condition (Q+p_i)^i≤iQ^i<(Q+p_i+1)^i. Put E and O for the numbers of triples (P,i,r) with P a Hamilton path of G_{N−1}, a_i before b, i respectively even or odd in 2,…,N, and 1≤r≤p_i. Then E and O are literal finite unweighted combinatorial counts and

    (E−O)/(2^(N−2)+1) = Σ_{i=2}^N (−1)^i p_i.

Dividing by Q yields a computable approximation to C_N, with absolute rounding error at most (N−1)/Q. The alternating remainder is at most f(N+1), so the bound is fully rigorous using integer root comparisons. Note carefully that E and O count **marked paths with token multiplicity**, not Hamilton paths of two ordinary tournaments.

For a tighter example, f''(x)>0 for x≥5, because writing g(x)=log(x)/x gives f''=e^g[(g')²+(2log x−3)/x³]. If N≥6 is even, the alternating remainder after the midpoint (C_N+C_{N+1})/2 lies between −(f(N+1)−f(N+2))/2 and 0. This follows by applying the alternating-series estimate to the decreasing differences f(k)−f(k+1). With N=200 and Q=10^9, exact integer comparisons give D=201254860, p_201=26735756 and p_202=26626878, hence

    375664886/2000000000 ≤ C ≤ 375774164/2000000000.

The interval is a checkable construction, not evidence of a new analytic value or exact-spectrum result. Run `python derivative_certificate.py` and `python marked_mrb.py` for the finite checks.
