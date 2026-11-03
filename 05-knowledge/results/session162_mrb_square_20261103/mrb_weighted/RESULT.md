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
