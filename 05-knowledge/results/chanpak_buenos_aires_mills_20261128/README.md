# A two-channel Buenos Aires criterion for Chan–Pak, and a Mills obstruction to fixed doubling rails

2026-11-28. Two different uses of nested arithmetic are at issue. The Buenos Aires constant retains **every first failure of divisibility**; a Mills constant retains **one infinite chain of intervals**. For Hamiltonian path counts the first becomes an exact two-channel reformulation of Chan–Pak. The second rules out taking the existing four +2/doubling rails as a closed cubic recursion. Neither proves Chan–Pak.

## 1. The roles of 2 and 3

Use 2 and 3 as integer stand-ins for e and π (their respective floors), not as approximations inside analytic identities. They have distinct exact combinatorial roles. In the [Grinberg–Stanley formula](https://arxiv.org/abs/2307.05569) a collection of j vertex-disjoint directed odd cycles contributes 2^j to the path count. The smallest such cycle is the directed triangle, with count 1+2=3. Ordered joining with that triangle multiplies the host count by 3. The dyadic gate instead doubles a *difference* and leaves the independent +2 contact intact. Cubing a count by joining three copies is yet another operation: it produces the left boundary of a Mills interval, not an interior successor.

[OEIS A249270](https://oeis.org/A249270) is the mean of the least prime not dividing an integer. Its [prime decoder](https://arxiv.org/abs/2010.15882) starts at 2.920050977316… and uses x ↦ ⌊x⌋(1+{x}); its first two floors are 2 and 3. Starting instead at the integer 2 freezes the decoder. The remainder, not just the floor, carries the later primes. Generalizing the *divisibility experiment* rather than the decoder is what preserves enough information for the tournament problem.

## 2. Two divisibility experiments

Let S be the set of Hamiltonian-path counts. An ordered join multiplies counts since every path concatenates a full path in each summand, so S is a multiplicative monoid. Write S* for the positive odds other than 7 and 21. Use the known exclusion of 7 and 21, as recorded in [Grinberg's lectures, Remark 1.1.11](https://www.cip.ifi.lmu.de/~grinberg/t/22s/lec12.pdf).

* Channel 0 uses all d in S. Its benchmark uses S*.
* Channel 1 uses d with 7∤d and 7d in S: divide the **exactly-one-7** layer by 7. Its benchmark uses all odd d with 7∤d other than 1 and 3.

Given either set U of admissible indices, let D_U(x) be the least d in U not dividing x. Put L_U(t)=lcm{d∈U : d≤t}, with empty lcm 1. The crucial exact event is D_U(x)>t iff L_U(t) divides x. Thus

    mean D_U = 1 + Σ_{t≥1} 1/L_U(t).

This is a cyclic-quotient version of the [Bou-Rabee–McReynolds residual averaging framework](https://arxiv.org/abs/0909.1343): x is separated from zero by the quotient ℤ/dℤ, and the surviving integers are in the intersection of its kernels up to index t. Taking U to be the primes recovers the Buenos Aires mean, since L_U is then a primorial. For the full tournament spectrum the filtration is by lcm rather than product, so imposing the 7-adic layer *before* taking lcm matters. This complements the repo's [primorial crossing-gap and prime-power inverse](../../../01-canon/theorems/THM-5110-primorial-crossing-gap-repair-and-brown-prime-power-inverse.md); that work likewise isolates information lost by passing to an arithmetic projection.

For both channels and both benchmarks this formula is an ordinary Cesàro limit, not a formal expectation. In channel 0 the counts 3 and 5 and the join law supply their powers. In channel 1, counts 35 and 63 supply indices 5 and 9; joining with counts 5 and 3 supplies all 5^b (b≥1) and 3^a (a≥2) as indices. The replay verifies the starting counts by actual arc masks. For t≥9 these two relatively prime powers give L_U(t)>t²/15. The finite average over 1≤x≤X has tail summands ⌊X/L_U(t)⌋/X ≤ 1/L_U(t), so summable domination proves the limit. One can equivalently use Haar probability on the profinite completion of ℤ.

## 3. Exact detection theorem

**Chan–Pak holds if and only if the two tournament means equal their respective benchmark means.** More precisely:

* Channel 0 equality holds iff every odd prime other than 7 is realized.
* Assuming channel 0 equality, channel 1 equality holds iff every 7p for prime p≥5, p≠7, is realized.

To prove detection, benchmark L is always a multiple of actual L. All mean differences (actual minus benchmark) are therefore nonnegative, and equality holds exactly when the filtrations agree at every integer threshold. A missing prime p in channel 0 changes the lcm at t=p: it divides no earlier index. A missing 7p in the layer changes the channel 1 lcm at t=p by the same argument. The mean increase is **at least (p−1)/L_benchmark(p)** in its corresponding channel, because the p factor occurs once at that threshold in the benchmark and not at all in the actual lcm.

For the first converse, all prime-power jumps away from 7 are supplied by ordered joins of prime witnesses. In the benchmark the first jump at 7 is at 35, then at 49, then at 343 and each higher pure power of 7. Counts 35, 49, and 343 are literal witnesses in this packet. Higher powers follow by joining 49 and 343, since every integer exponent ≥2 is a nonnegative combination of 2 and 3. This proves equality in channel 0. For the second converse, take any admissible benchmark index d. If a prime p≥5 divides d, build 7d by joining a witness for 7p with witnesses for the prime factors of d/p. If no such prime divides d, it is a power 3^a, a≥2, and 63 followed by copies of the triangle suffices. Finally, every target with v₇≥2 is a join of its pure 7-part and its 7-free part. This proves the full equivalence and gives an explicit assembly rule from the two families of seeds.

In particular the unfiltered mean cannot detect a missing 77 or 63 if all other conjectured values are present: 77 is redundant after 35 and 11, and 63 is redundant after 35 and 9. The channel 1 quotient makes both visible. A separate 7-channel is not decorative; it recovers exactly the prime tests lost in the unfiltered lcm. Conversely, agreement to many decimal places does not prove equality. For a concrete demonstration, the replay verifies 91 literal witnesses (moving-root words, sometimes with one arc reversed), the common lcm filtrations through 200, and 83 and 82 common decimal places respectively:

    D₀ mean: 3.982460042366362357753322...
    D₁ mean: 5.849126709035833817467332...

The certification uses rational arithmetic. For N=200, each common L(N) divides L(t) at t>N; additionally L(t) contains 3^floor(log₃t) and 5^floor(log₅t). Since their 3- and 5-parts in L(N) are at most N each, L(t)/L(N)>t²/(15N²) and the remaining sum is at most 15N/L(N). Thus the decimal control does not assume unproved spectral coverage above 200. The point is not precision but the bottleneck: any exact construction that proves all ordinary prime and 7p seeds finishes the problem; higher composites need no independent descent.

## 4. A fixed dyadic bank cannot carry a Mills chain

[Mills' theorem](https://doi.org/10.1090/S0002-9904-1947-08849-2) gives A>1 with all floors of A^(3^j) prime. The necessary step from floor h to floor h' is h³≤h'<(h+1)³. Earlier repo work on [Mills prime shells and remainder decoding](../seam_mills_20260925.md) stresses the stored gap inside that shell. The following obstruction asks whether the existing +2/doubling mechanism can even inhabit the shells.

**Finite-rail theorem.** Fix finitely many integer pairs (a_i,b_i), a_i>0. No infinite sequence of odd floors ⌊A^(c^j)⌋, for A>1 and integer c≥2, can eventually lie in the union of the rails a_i2^k+b_i (k≥0).

Write x_j=A^(c^j), and L_j=a_i2^k for the leading term of a chosen rail at time j. Then x_j−L_j is bounded, so log₂ x_j−log₂ L_j tends to zero. Also log₂ L_{j+1}−c log₂ L_j tends to zero. The latter difference lies in a finite union of translates of ℤ, one for each pair of slopes, and therefore equals zero from some point on. Thus L_{j+1}=L_j^c. The logarithmic error now multiplies by c at each step and tends to zero, forcing it to be identically zero on the tail. The odd part of L_j also raises to the c-th power at each step, but is bounded by the largest slope. Hence it is eventually 1. Therefore x_j is eventually a positive power of 2 and therefore even, contradicting the odd floors. This uses real interval nesting, not just agreement of leading residues.

For the [stored +2/dyadic gate](../chanpak_unit_gap_dyadic_gate_20261128), the conclusion is stronger: **there is not even one cubic-shell edge between its four rails** h=a·2^k+b, a∈{495,567}, b∈{303,305}, k≥1. For a candidate edge h→h', set L=h−b and L'=h'−b'. The shell implies

    1 < L'/L³ < (1+306/L)³.

The lower inequality follows from (L+303)³−305>L³. The ratio is b_s·2^m/a_s³ for a pair of slopes a_s,b_s and an integer m. The least such ratio above 1 is 495·2^18/495³ > 1.0698. For k≥5 the right-hand side is less than 1.060. Direct integer checks at k=1,2,3,4 exclude the remaining cases. A finite enlargement by the [fixed-root restart](../chanpak_dyadic_gate_restart_costs_20261128) can create new single edges, but as long as its released counts stay on finitely many dyadic affine rails the infinite-chain obstruction still applies. The [odd CRT chart](../chanpak_joint_crt_and_dyadic_boundary_20261128) addresses residues and does not affect this real-size constraint.

The constructive takeaway is different for the two models: for complete coverage, target the prime and 7p seed families; for a genuinely Mills-like infinite cubic recursion, give the leading coefficient infinitely many choices rather than trying to tune the bounded +2 offset. Neither obstruction applies to general tournaments.

## Replay

Run python3 replay.py. The 91 arc descriptions are in witnesses.json. The replay uses subset Hamiltonian-path DP, join closure, exact rational intervals and exact cubic-shell inequalities.
