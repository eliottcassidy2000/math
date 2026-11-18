# A recursive 3+1 tournament wheel at 15, and why it does not close Chan–Pak

The [least-order selector](../chanpak_prime_frontier_selector_20261118/README.md) always terminates. An input gives a witness, or a finite absence certificate (by stripping singleton strong components and using Busch's bound). Showing that the latter branch occurs exactly at nonpositive integers, evens, and 7,21 is **the conjecture**, not a consequence of termination. Zero belongs with the nonpositive inputs.

The [local-reduction packet](../chanpak_weierstrass_reduction_levels_20261118/README.md) proved it suffices to cover odd units modulo 15 other than 7; products by 3 and 5, with 35 and 63 for the exceptional remaining core, recover the rest. The still sharper ordinal atom reduction leaves only strong p and 7p sources for primes p≥11 once the fixed seeds are supplied. In the modularity analogy this is the step where a restricted theorem is identified as sufficient. The missing step remains an absolute lift inside that class.

## 1. Parity and sign are two projections, not two operations on H

The positive odd quadrant is the only one where path counts can occur. There is no unsigned tournament operation mapping H to −H. On the unit chart one can be more precise: the eight units modulo 30 correspond by CRT to (nonzero mod 3, nonzero mod 5), or two phase squares. Multiplying by 7 leaves the mod-3 sign unchanged and turns the mod-5 phase by one quarter (7≡2 mod 5). Negation reverses the mod-3 sign and turns the phase by two quarters. This is a **finite residue plane**, not a geometric placement of absolute targets. The 21 hole is outside the chart and reduces to the 7 hole on stripping its factor 3. The two holes are not entire residue cells.

A second parity lives on labels, not on counts. Every tournament automorphism has odd order as a group: an involution would interchange some pair of vertices and reverse an arc, impossible. Hence there are exactly two even-permutation orbits within an ordinary isomorphism class (order≥2). In a self-converse class, an antiautomorphism can be raised to an odd power to become an anti-involution: its square is an automorphism of odd order. The anti-involution fixes at most one vertex, so its sign is (−1)^floor(n/2). The two labeling sheets are preserved by converse at n≡0,1 mod4 and exchanged at n≡2,3 mod4. This statement separates **orientation**, **permutation sign**, and **oddness of H**.

## 2. Three physical 3+1 probes act on one single polynomial

Use the inherited **slot polynomial** f_T(t), giving the number of ordered covers by t optional directed path slots at nonnegative integer t. Equivalently the [Grinberg–Stanley odd-cycle formula](https://arxiv.org/abs/2307.05569) makes it a monic parity polynomial: a nontrivial directed odd cycle has weight 2 and contracts to one atom. In particular f(0)=0 for nonempty T, f(1)=h and f(2)=2h+q, where q counts ordered nonempty two-covers. Disjoint palette sets give the convolution

    f_T(t+s) = sum_(A⊆V) f_A(s) f_(V−A)(t).

Also f_A(−1)=(−1)^|A| H(A). Thus sums over odd or even size selected subsets are exactly the half-difference or half-sum of the two reads f(t+1), f(t−1), and the even positive sum additionally excludes the empty subset.

In each probe replace a single marked vertex of a four-vertex host by T, with uniform contacts, so three actual vertices are added. Use a cyclic triangle with an ordinal singleton for P; Q and R use the strong H=5 tournament with outgoing masks (2,4,9,3), marks 0 and 1 respectively. The exact transformations, with f±=f(t±1), are

    P f = t³ f + t²(f+−f−),
    Q f = t³ f + t[(t+1)f+ + (1−t)f−],
    R f = t³ f + 2t²(f+−f−) + t(f++f−−2f).

Proof from crossing cycles. For P an external cycle through T must use the two cyclic exterior vertices, a single odd-length T path, and no other exterior vertex, leaving the ordinal singleton. For Q the exterior itself is a directed triangle. Any crossing cycle must leave T at 1 and return through 2 or 3; it can use exterior word 12 or 123. The former leaves one singleton, requires odd core size, and the latter uses all exterior, requires positive even core size. Thus the term from avoiding crossing cycles is (t³+2t)f and the other two terms combine to the displayed Q. For R the exterior is transitive 2→3→0; a crossing cycle starts its outside word at 2, with two two-vertex choices 20,23 and one three-vertex choice 230. All crossing cycles for any given probe share their exit vertex, so no packing can contain two. This makes the formulas exact for arbitrary T, not only transitive padding.

At t=1 they read (3h+q,5h+2q,5h+3q). In particular P,Q are a **unimodular** pair: h=2H(P)−H(Q), q=3H(Q)−5H(P). This compresses both lowest coordinates without losing an integer digit. It does not imply one can freely choose both measurements or physically subtract them. Under further composition the full polynomial is the exact state. After d probes the closed count can depend on f up to d+1; the coefficient of the top read is a positive product of the upward coefficients at 1,…,d. At fixed modulus, higher reads do have a cutoff: the coefficient of f(j+1) is divisible by j! since every evaluation walk from 1 must cross upwards once from each of 1,…,j, and each upward coefficient has that level as a factor. The exact polynomial therefore gives a recursive compression with a genuine prime-power precision budget, not a frozen rank-two state.

## 3. The physical quarter-turn and an indefinitely reusable local wheel

For h a 3-unit, reduction at 3 forces P for an even-order module and Q for an odd-order module if using only P,Q; either safe step negates h modulo 3. The other choice gives zero, which subsequent P,Q,R probes cannot undo, by parity and f(2)≡(−1)^|T| f(1) mod3. R is always safe at 3 and also negates h there. Consequently (-1)^|T| h mod3 remains fixed along every safe word.

For a 5-unit state record (order parity,evaluations u=f(1),v=f(2) mod5): f(3)=±v and f(4)=±u by parity. For the two safe choices at each parity the mod-5 matrices on (u,v) are

    even P  [[1,1],[1,2]]     odd Q [[1,2],[3,2]]
    either R [[4,3],[4,4]].

With only P then Q on an even position, the two-step matrix is diag(3,2). The even-layer first coordinate therefore rotates by ×3 modulo 5, while its mod-3 sign is unchanged. By CRT and oddness, this is exactly ×13 modulo 30, the **inverse** of the ×7 target-wheel rotation. Four probe pairs return the residual state. With R admitted, discard any move landing at u=0; the 40 states with u≠0 have four traps (odd, u, 2u), and the other 36 form a single strongly connected component with diameter 11. This finite calculation is reproduced exactly, including every start state, in replay.py. Starting at TT2, one stays in the live component and can reach **each** coprime residue modulo 30 with its required order parity, in at most 11 steps from any live state in that same signed chart. The eight literal seed routes and integer counts are in replay.out. In particular one can repeatedly request residues without falling into a trap. This is a widening from one forced direction to a reusable local selector, not a claim of absolute coverage.

## 4. A constraint on an attempted completion: growth beats the alphabet

Even with arbitrary words in all three probes and arbitrary ordinal joins with any fixed finite library of satellite tournaments between probes, a finite set of seeds produces only X^{o(1)} distinct counts up to X. This includes arbitrary singleton padding even though it can flip parity at zero immediate cost in H.

Indeed for a core of order n, q≥(n−1)h by cutting every Hamilton path at each of its internal gaps; the resulting ordered covers identify the original path and cut. Each probe therefore multiplies h by at least n+2, and adds three vertices. So a word reaching X uses d=O(log X/log log X) probes. Every probe output is also at least the current f(2). For each ordinal satellite B before that probe, f(2) multiplies by f_B(2)≥2^|B|, since every partition of vertices into two slots has at least one internal directed path in each part. Thus the sum of the sizes of all satellites before probes is at most log₂ X. The order of satellites between any two probes is invisible to f, by multiplicativity of the polynomial. For a library of M fixed types, the multiplicity schedules are bounded by binom(floor(log₂ X)+Md,Md); combine with at most 3^d probes. The logarithm is O(log X·log log log X/log log X)=o(log X). Final ordinal factors contribute at most a fixed power of log X (ignore count-one factors). Summing over the bounded seed set and all d preserves subpower growth. No “widening” within this fixed recursive architecture can by itself cover the unit chart in absolute values.

So the next test is not a bigger fixed alphabet of these same three-vertex returns. It is a target-dependent incoming core carrying enough adjustable height, or a nonuniform contact that preserves the readout without paying the factorial growth. A local cycle, a reachability theorem, and an absolute selector are three distinct things.

## 5. Earlier fifteens with a mechanism, and those without one

- The [Fibonacci–Farey q=15 study](../../../07-reflections/fibonacci-farey-q15-four-state-and-mod6-24-state-bridge-codex-20260814.md) already identifies units mod15 up to sign with a four-cycle on P¹(F3). Its vector action lifts the projective period from four to eight. Our new CRT phase is an actual tournament probe word acting on slot values; neither cycle determines an absolute integer.
- The [Hamming / projective atlas](../relation_gluing_7_15_20261028/ATLAS.md) follows seven nonzero vectors in F₂³ to fifteen Hamming codewords in F₂⁴ minus zero. The latter give the fifteen points of PG(3,2). A chosen symplectic form distinguishes fifteen matching lines from twenty triangle lines; the related [Bolza packet](../session52_saturation_residues_metagraph_20261024/fifteen-kernels/README.md) has fifteen Richelot kernel choices but just four observable classes with weights 6,4,4,1. The transfer is the distinction between a chart of choices and what survives a reader, not a map from those kernels to tournament residue states.
- The [earlier tournament fifteen](../outdeg_converse_descent_fifteen_20261116/FIFTEEN.md) has 15 paths on the cyclic five-tournament, three orbits of five, also the orbit counts of the fifteen A5 double transpositions under a five-cycle. The displayed path words through a fixed midpoint have permutation signs −,−,+, so globally 10 negative and 5 positive, versus all even involutions. There is a cyclic matching but not a sign-preserving one with those choices. Replay computes the signs directly.
- [The earlier ledger](../analogy_reduction_atlas_20261025/LEDGER_15_25_MULTIPERFECT.md) separates Ogg's fifteen primes, the ternary exponent 15, and the 15-state machine for a different ternary question. The square-sum G15 threshold instead has actual endpoint constraints. In each case a genuine connection requires compatible transitions, not just a state count. None supplies the missing co-prime target constructor.

Reproduce the inherited-or-independent literal checks with `python replay.py | diff -u replay.out -`. The requested all-order absence classification remains open.
