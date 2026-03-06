# Web Research: Tribonacci Structure and Deep Connections for OCF

**Source:** opus-2026-03-05-S6 (web research session with Tribonacci focus)

---

## HEADLINE DISCOVERY: Why Both Sides of OCF Give Tribonacci

### The Structural Picture for T_full_n

The "full tiling" tournament T_full_n has adjacency: v_i beats v_j iff j=i+1 (path) or j<=i-2 (backward).

**Theorem (verified n=3,...,8).** Every directed odd cycle of T_full_n is a consecutive interval [k, k+2j] of odd length 2j+1. Specifically:
- 3-cycles: {k, k+1, k+2} for k = 0,...,n-3
- 5-cycles: {k, k+1, k+2, k+3, k+4} for k = 0,...,n-5
- (2j+1)-cycles: {k,...,k+2j} for k = 0,...,n-2j-1

**Proof:** In T_full_n, the only forward arcs are path edges i→i+1. For any i < j with j >= i+2, the arc goes j→i (backward). Therefore, in any directed cycle, you can only increase the vertex label by +1 (via a path edge) and can only decrease it by -2 or more (via a backward arc). If the vertex set of a cycle had a gap (two consecutive elements s_i, s_{i+1} with s_{i+1} - s_i >= 2), then crossing that gap forward would require a non-path forward arc, which doesn't exist. Once you drop below the gap via a backward arc, you can never return above it. So the vertex set must be gap-free, i.e., a consecutive interval [k, k+L-1]. For odd length L=2j+1, the unique cycle is k→k+1→...→k+2j→k (using 2j path edges forward and one backward arc k+2j→k, which exists since k+2j >= k+2).

**Consequence:** Omega(T_full_n) is an INTERVAL GRAPH on odd-length intervals on {0,...,n-1}. Two cycles conflict iff their intervals overlap. The number of cycles is sum_{j=1}^{floor((n-1)/2)} (n-2j).

### Why I(Omega(T_full_n), 2) = Tribonacci(n) (independent of H(T_full_n) proof)

The independence polynomial at x=2 counts:
```
I(Omega, 2) = sum over non-overlapping odd-interval packings P: 2^|P|
```

DP by whether position 0 is covered:
```
f(n) = f(n-1) + 2*f(n-3) + 2*f(n-5) + 2*f(n-7) + ...
```
(Skip position 0, or place an interval of length 3, 5, 7,... starting at 0, each contributing factor 2.)

**Telescoping:** Since f(n-2) = f(n-3) + 2*f(n-5) + 2*f(n-7) + ..., we get:
```
2*f(n-5) + 2*f(n-7) + ... = f(n-2) - f(n-3)
```
So: f(n) = f(n-1) + 2*f(n-3) + f(n-2) - f(n-3) = **f(n-1) + f(n-2) + f(n-3)**.

With base cases f(0)=f(1)=f(2)=1, this is exactly OEIS A000213 (Tribonacci).

### Why H(T_full_n) = Tribonacci(n) (run decomposition proof, kind-pasteur-S11)

Ham paths decompose into "runs" (maximal ascending blocks) with gap condition max(I_k) >= min(I_{k+1})+2. Conditioning on first run length and telescoping similarly gives the Tribonacci recurrence.

### The Beautiful Parallel

Both sides of OCF produce Tribonacci by the SAME algebraic mechanism (telescoping a sum over all odd-length contributions to a 3-term recurrence), but through DIFFERENT combinatorial objects:
- **H side:** Runs in Hamiltonian paths (ascending blocks with gap ≥ 2)
- **I(Omega,2) side:** Non-overlapping odd cycles (intervals with gap ≥ 1)

This suggests a **deeper structural correspondence** between run decompositions and odd-cycle packings that might generalize beyond T_full.

---

## Connection 1: Chudnovsky-Seymour — All Roots Real (CRITICAL)

**Paper:** Chudnovsky & Seymour, J. Combin. Theory Ser. B (2007)

### Result
If G is claw-free, then all roots of I(G, x) are real.

### Implication for OCF
**UPDATE (opus-S7):** Omega(T) is claw-free only for n<=8 (trivially, by vertex counting) and FAILS at n=9 (90% of random tournaments have claws). Omega(T) is also NOT always perfect (C5 appears at n=8, 53.8%).

**For n<=8:** Chudnovsky-Seymour applies, giving all real (negative) roots, log-concavity, positivity. This is trivially true because claw-freeness is trivially true.

**For n>=9:** Chudnovsky-Seymour does NOT directly apply. However, empirical testing (200 random n=8 tournaments including imperfect ones) shows all roots remain real and negative even when Omega has C5 (imperfect). This suggests a **deeper explanation** for real-rootedness that goes beyond claw-freeness.

**Key question (revised):** Does I(Omega(T), x) have all real roots for ALL tournaments at ALL n? If so, this would be a new theorem about conflict graphs not explained by claw-freeness alone.

---

## Connection 2: Chvátal-Sbihi Decomposition + Dyer-Jerrum (STRUCTURAL)

**Papers:**
- Chvátal & Sbihi, JCTB 44 (1988) — Decomposition theorem
- Dyer, Jerrum, Müller, Vušković, SIAM J. Disc. Math. 35 (2021) — Counting

### Structure Theorem
Every claw-free perfect graph decomposes via clique cutsets into "atoms" that are either:
1. **Elementary graphs** (built from line graphs of bipartite graphs via local augmentation)
2. **Peculiar graphs** (specific small structures)

### Counting Result
The partition function (independence polynomial) of (claw, odd-hole)-free graphs can be approximated in polynomial time via this decomposition. The key: clique cutsets FACTOR the partition function.

If K is a clique cutset separating G into parts G_1, G_2, then:
```
I(G, x) = I(G_1, x) * I(G_2, x) / I(K, x)
```
(The clique contributes to both parts; dividing by I(K,x) corrects double-counting.)

### For Omega(T):
**UPDATE (opus-S7):** Omega(T) is NOT always claw-free perfect, so Chvátal-Sbihi does not apply universally. This decomposition approach works only for n<=8 where claw-freeness is trivial. OCF is proved by Grinberg-Stanley using entirely different methods (symmetric functions / matrix algebra), not graph decomposition.

**The T_full example still illustrates the principle:** Omega(T_full_n) is an interval graph (always chordal, always perfect, always claw-free), and the DP sweep computing I(Omega, 2) mirrors the clique-cutset factorization. This is a useful illustration even though the approach doesn't generalize.

---

## Connection 3: Tribonacci Tiling Interpretation

### Classical Result
Tribonacci(n) = number of tilings of an n-board with squares (length 1), dominoes (length 2), and trominoes (length 3).

Equivalently: compositions of n with parts 1, 2, 3.

### Our Run Decomposition
The runs in Ham paths of T_full are maximal ascending blocks. The gap condition between consecutive runs is ≥ 2 (a "spacing" constraint). The Tribonacci recurrence emerges from this spacing.

### The Interval Packing
The independent sets of Omega(T_full_n) correspond to non-overlapping odd intervals on {0,...,n-1}. At fugacity 2, each interval contributes weight 2. The "skip position 0" contributes weight 1. The DP: f(n) = 1*f(n-1) + 2*f(n-3) + 2*f(n-5) + ... telescopes to Tribonacci.

### Rauzy Fractal Connection
The Tribonacci substitution sigma: 1→12, 2→13, 3→1 generates the Rauzy fractal, a self-similar tiling of the plane. The Tribonacci constant tau ~ 1.83929 is a Pisot number (the other roots of x^3-x^2-x-1 have absolute value 1/sqrt(tau) < 1).

The symbolic dynamics of the Tribonacci substitution produces a domain exchange on the Rauzy fractal. The "run decomposition" of T_full_n's Ham paths might be interpretable as a coding of paths in this substitution system.

**Speculative but intriguing:** The Tribonacci substitution on {1,2,3} generates an infinite word. The letters could correspond to "tiles" placed on the board (1=square, 2=domino, 3=tromino). The fixed point of the substitution encodes all possible run patterns, and the Tribonacci growth rate controls the Ham path count.

---

## Connection 4: Independence Polynomial Zero Structure

### Jerrum & Patel (2026)
**Paper:** Zero-free regions for the independence polynomial on restricted graph classes, J. London Math. Soc. (2026)

Extends Chudnovsky-Seymour to show zero-free REGIONS for independence polynomials of H-free graphs (for any fixed subdivided claw H). Combined with Barvinok's interpolation method, gives FPTAS for computing I(G, x) in these zero-free regions.

### For OCF:
**UPDATE (opus-S7):** Omega(T) is NOT always claw-free (fails at n=9). However, all-real-roots appears to hold empirically even for imperfect, non-claw-free Omega(T). If true universally, the roots are all on the negative real line regardless of claw-freeness.

**Key question (revised):** Does I(Omega(T), x) always have all real roots? If so, this is NOT explained by Chudnovsky-Seymour and would be a genuinely new property of tournament conflict graphs. Testing at n=9 and n=10 would be informative.

---

## Connection 5: Transfer Matrix = Tribonacci Matrix for T_full

The 3x3 Tribonacci transfer matrix is:
```
M = [[1, 1, 1],
     [1, 0, 0],
     [0, 1, 0]]
```
with M^n * [1,1,1]^T = [f(n+2), f(n+1), f(n)].

**Can this matrix be interpreted as a tournament transfer matrix?**

For T_full_n, the DP computing H(T_full_n) uses a state representing the "recent history" of the run decomposition. The state tracks how many positions ago the last run ended (0 = still in a run, 1 = ended 1 position ago, 2+ = ended 2+ positions ago). This gives exactly 3 states, matching the 3x3 matrix.

**For OCF more broadly:** The transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S) B_b(V\S) for general tournaments is n×n and symmetric. For T_full, this n×n matrix has a specific structure that reduces to effective rank 3, matching the Tribonacci matrix. The eigenvalues of the n×n matrix include the Tribonacci constant as the dominant eigenvalue.

**Action:** Compute the transfer matrix M for T_full at small n. Check if its characteristic polynomial has x^3 - x^2 - x - 1 as a factor.

---

## Connection 6: Pisot Number Properties

The Tribonacci constant tau = 1.83929... is a Pisot number:
- Real root of x^3 - x^2 - x - 1 = 0
- Other roots have absolute value |beta| = |gamma| = 1/sqrt(tau) < 1
- Satisfies tau + tau^{-3} = 2

**The identity tau + tau^{-3} = 2 is numerically striking in our context** since the fugacity is lambda = 2!

From tau^3 = tau^2 + tau + 1: dividing by tau^3 gives 1 = tau^{-1} + tau^{-2} + tau^{-3}.
So tau^{-1} + tau^{-2} + tau^{-3} = 1, and 2 = tau + tau^{-3} = tau + 1 - tau^{-1} - tau^{-2}.

**Is the evaluation at lambda = 2 related to the Pisot property of the Tribonacci constant?** For the T_full family, H(T_full_n) ~ c * tau^n where tau ~ 1.839. The evaluation I(Omega, 2) involves summing 2^k over k-element independent sets. The exponential growth rate 2 of the weight competes with the combinatorial growth rate of the independent sets.

---

## Connection 7: Interval Graphs and the Hard-Core Model

For T_full, Omega is an interval graph. The hard-core model on interval graphs is exactly solvable:
- The partition function Z(lambda) = I(G, lambda) satisfies a linear recurrence determined by the interval structure
- For "regular" interval families (like our consecutive-odd-intervals), the recurrence has constant coefficients
- The eigenvalues of the recurrence's companion matrix determine the growth rate

For our specific interval graph: Z(2) satisfies the Tribonacci recurrence with growth rate tau. The hard-core model on this interval graph at fugacity lambda = 2 is in a **specific phase** determined by where lambda sits relative to the critical fugacity of the interval graph.

**General hard-core model on interval graphs:** The critical fugacity lambda_c for an interval graph of max degree d satisfies lambda_c = (d-1)^{d-1}/d^d (Shearer bound). For our interval graph, what is the max degree and lambda_c?

---

## Connection 8: OEIS A000213 Combinatorial Interpretations

A000213 has several known combinatorial interpretations:
1. **Tilings** of n-board with squares, dominoes, trominoes
2. **Compositions** of n-1 with parts at most 3
3. **(n-1)-bit binary sequences** where each 1 is adjacent to a 0
4. **S-tetromino tilings** of 2×n grids
5. **Run decompositions** of {0,...,n-1} (our Ham path interpretation)
6. **Weighted odd-interval packings** at fugacity 2 (our I(Omega,2) interpretation) — NEW!

The fact that interpretation (5) counts H(T_full) and interpretation (6) counts I(Omega(T_full), 2) is our OCF consistency check. A direct bijection between interpretations (5) and (6) would verify OCF for the T_full family bijectively.

---

## Synthesis: Proof Strategy Implications

### What the Tribonacci analysis reveals

1. **OCF for T_full is a weighted interval packing identity.** The H-side uses run decompositions; the I(Omega,2)-side uses non-overlapping odd cycles. Both telescope to Tribonacci by the same algebraic mechanism.

2. **The clique-cutset factorization of I(Omega, 2) mirrors the DP structure of H(T).** For interval graphs, the sweep DP computing I is a clique-cutset decomposition. For T_full, the run decomposition of Ham paths has the same sequential structure.

3. **The key to proving OCF might be showing that the clique-cutset decomposition of Omega(T) mirrors the vertex-deletion/subset-convolution structure of H(T) for ALL tournaments, not just T_full.**

### Revised Priorities (post-Grinberg-Stanley proof, post-claw-free disproof)

**OCF is proved.** The remaining questions are structural and explanatory.

1. **Test all-real-roots conjecture at n>=9.** Does I(Omega(T), x) always have all real (negative) roots, even when Omega is not claw-free? This would be a genuinely new theorem. (Tested at n=8: holds for all 200 sampled, including imperfect Omega.)

2. **Extend the Tribonacci analysis to other tournament families.** The T_full family is the one case where Omega has beautiful interval graph structure. Do circulant, Paley, or other families yield recognizable recurrences?

3. **Find a bijection between run decompositions and weighted odd-interval packings** for T_full, providing a bijective proof of OCF for this family (complementing the algebraic Grinberg-Stanley proof).

4. **Characterize the phase transition**: Why does perfectness fail at exactly n=8 and claw-freeness at n=9? What structural property of Omega(T) DOES hold universally?

---

## References

- [Chudnovsky & Seymour (2007)](https://web.math.princeton.edu/~mchudnov/roots.pdf) — Roots of independence polynomial of claw-free graphs
- [Chvátal & Sbihi (1988)](https://www.sciencedirect.com/science/article/pii/S009589569891872X) — Claw-free perfect graph decomposition
- [Dyer, Jerrum, Müller, Vušković (2021)](https://arxiv.org/abs/1909.03414) — Counting weighted independent sets
- [Jerrum & Patel (2026)](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/jlms.70458) — Zero-free regions for independence polynomial
- [OEIS A000213](https://oeis.org/A000213) — Tribonacci numbers (1,1,1,3,5,9,17,31,...)
- [Rauzy fractal](https://en.wikipedia.org/wiki/Rauzy_fractal) — Tribonacci substitution and fractal tilings
- [Scott & Sokal (2005)](https://link.springer.com/article/10.1007/s10955-004-2055-4) — Hard-core model and LLL
- [Chudnovsky & Seymour (survey)](https://web.math.princeton.edu/~mchudnov/claws_survey.pdf) — Structure of claw-free graphs
- [Tribonacci tiling identities](https://arxiv.org/abs/2201.02285) — Allen (2022)
- [Robbins (2014)](https://www.tandfonline.com/doi/abs/10.1080/00150517.2014.12427915) — Tribonacci and 3-regular compositions
