# Valuation cycles, divisor heights, and branched-cover completions

Session43, 2026-10-23. Status: elementary proofs plus independently replayable finite checks. The prior prime-only quotient obstruction is retained, not contradicted. The positive divisor-height extension is a different, explicitly defined domain. No Monster-to-tournament equivalence, Collatz implication, or novelty claim for Belyi maps is asserted.

## 1. Recover the operation before interpreting the labels

The repository source is `05-knowledge/results/monster_levels_square_sum_7_21_20261022.md`, section3, and its companion `05-knowledge/results/session32_square_monster_20261023/VALUATION_CYCLE_AND_QUOTIENT.md`. The input is

N = 2^5 3^7 5^3 7^4 11 13^2 = 39,046,339,332,000.

Its prime-to-exponent arrows are

13 → 2 → 5 → 3 → 7 → 4; and 11 → 1.

They arise from the small-prime core of the ratio |Monster|/|Baby Monster|. The cycle is not initially a Collatz operation, prime-index operation, or residue multiplication. Identifying the terminal4 with the initial13 turns the six-symbol chain into a FIVE-state cycle.

The source already proves that this map does not descend modulo9: inputs11 and2 coincide modulo9 but outputs1 and5 do not. Any modulus identifying13 with4 divides9. Modulo3, inputs2 and5 agree but outputs5 and3 disagree. Thus no nontrivial quotient ring realizes the full displayed prime map. Even the selected five arrows cannot be an integer-polynomial map modulo9, since such a polynomial respects congruence modulo3.

Those obstructions concern arithmetic quotients of the prime-only map. They do not rule out the following natural extension.

## 2. A literal integer cycle through the composite-base divisor height

For any integer N>1 and base b≥2, define

h_N(b) = max{k≥0 : b^k divides N}.

This is a divisibility height, not a p-adic valuation when b is composite. If N=∏p^e_p and b=∏p^f_p divides N, then

h_N(b) = min_(f_p>0) floor(e_p/f_p).

Proof: b^k divides N exactly when k f_p≤e_p for every prime appearing in b. In particular h_N(p)=v_p(N), so this extends the source operation without fitting an arbitrary table.

For the displayed N,

h_N(4)=floor(5/2)=2=h_N(13).

Consequently

4 → 2 → 5 → 3 → 7 → 4

is already a literal five-cycle on integers. The prime13 is a predecessor of2. The equality13=4 can now be interpreted as folding two distinct states with identical forward images, rather than as pretending that13 is an integer equal to4 or as an invalid mod9 quotient. The fold is a valid forward dynamical congruence. It does not preserve incoming neighborhoods or arithmetic labels: marking the prime13 and the composite4 remains necessary when those are observables.

No value is defined by the maximum formula at b=1: all powers of1 divide N, so the maximum is infinite. In this note1 is a TERMINAL state, not a silently assigned finite valuation. One may add a bookkeeping self-loop at1, but it is an explicit convention. Likewise bases not dividing N have height0 and are outside the nonunit-divisor domain.

### Exact all-divisor census

For k≥1,

A_N(k) := #{b≥2 : b^k divides N} = ∏_(p|N)(floor(e_p/k)+1) − 1.

Proof: choose each exponent of b independently from0 through floor(e_p/k), and discard b=1. Therefore the number of bases of exact height k is A_N(k)−A_N(k+1). This is a general formula, not an enumeration-dependent claim.

Here the tail counts A_N(1),...,A_N(8) are

5759, 143, 23, 7, 3, 1, 1, 0.

The exact-height generating polynomial is

Φ_N(z) = 5616z + 120z² + 16z³ + 4z⁴ + 2z⁵ + z⁷.

All outputs lie in{1,...,7}. Among2,...,7 the five-cycle is as above, and6→5. Thus every nonunit divisor of height at least2 enters that cycle; every divisor of height1 terminates at1. There are exactly143 cycle-basin bases and5616 terminal-basin bases. The cycle-basin condition is precisely b²|N. This is a concrete role for143, not evidence that every prior appearance of143 has the same origin.

The exact census is checked by two independent calculations for all5759 nonunit divisors: exponent minima and repeated integer divisibility. Normal and optimized Python outputs agree byte-for-byte.

### Why this yields larger families without proving exceptional significance

For an arbitrary list of distinct primes p_1,...,p_r and permutation π, the integer

N_π = ∏ p_i^(p_(π(i)))

has h_(N_π)(p_i)=p_(π(i)). Thus ANY prescribed permutation of prime states occurs inside a divisor-height graph. More general directed prime-state maps are realizable the same way. This family proves that cycle shape alone is not a selective fingerprint of the Monster. The meaningful extra data here are the actual group-order origin, the composite-base extension, the full divisor profile, and the effect of changing the group normalization.

For example replacing |B| by the actual involution centralizer order|2.B| replaces N by N/2, so2→4→2 becomes a literal two-cycle. The original positive repair therefore preserves, rather than evades, the source's normalization warning.

For any N with largest prime exponent E, all nonunit-divisor outputs lie in{1,...,E}. Consequently its recurrent divisibility dynamics can be determined on at most E−1 nonterminal integer bases, regardless of the magnitude of N or its divisor count. If an intermediate output fails to divide N, the next height is0, which must be treated as another terminal state; for this particular N every integer2,...,7 does divide N. This proviso is essential for the general statement.

## 3. A different positive completion: finite monodromy

Return to the selected quotient-state graph, with five-cycle C={ [4=13],2,5,3,7 } and side arrow11→1. The partially prescribed map is injective. On exactly these seven states it has a unique permutation completion: add1→11. The result has cycle type(5,2).

This completion is OPTIONAL additional structure. It is not the divisor-height extension of section2, where1 is terminal. It is useful precisely as a reversible alternative, not as a claim that the two models are equal.

A connected finite unbranched cover of C* has cyclic transitive monodromy, hence one cycle on its fibre. The permutation(5)(2) instead represents a disconnected degree5-plus-degree2 cover. To connect its two monodromy orbits, add a second generator that transposes one state in each orbit. This gives a concrete connected branched-cover model.

### Explicit rational map

For positive integers a,b, n=a+b, define

P_(a,b)(z) = n^n/(a^a b^b) · z^a(1−z)^b.

The derivative is

P'_(a,b)(z) = n^n/(a^a b^b) · z^(a−1)(1−z)^(b−1) · (a−nz).

The extra critical point is c=a/n and P(c)=1. The zero fibre has multiplicities(a,b); the fibre over1 has one double point and n−2 simple points; infinity is a single pole of order n. Hence the branch passports are

(a,b), (2,1^(n−2)), (n).

For(a,b)=(5,2),

P(z)=7^7 z^5(z−1)^2/(4·5^5).

The explicit polynomial proves existence without invoking a general existence theorem. Its dessin is the double-star tree, and its monodromy has two cycles joined by one transposition. A recent primary exposition of this known family is Dinamo, https://arxiv.org/abs/2609.19053. We independently prove the statements needed here; this is a connection and realization, not a claim to invent the family.

### Monodromy group and exact obstruction to full symmetry

Let σ be disjoint cycles of lengths a,b, and τ a transposition joining one element of each. Put c=στ. Cutting and joining the cycles shows that c is an n-cycle. In the cyclic order supplied by c the transposition endpoints have separation a (or b in the reverse convention). The conjugates c^j τ c^(−j) transpose pairs at cyclic distance a. Their graph has d=gcd(a,b)=gcd(a,n) connected components, each of size m=n/d.

Transpositions along the edges of a connected graph generate its full symmetric group: along any simple path, conjugation transports an endpoint transposition to any desired pair. Thus the subgroup generated by these conjugates is (S_m)^d. The n-cycle c rotates the d components transitively, and c^d belongs to that base subgroup. Therefore

<σ,τ> = (S_m)^d ⋊ C_d = S_m wreath C_d,

of order (m!)^d d. In particular it is S_n exactly when gcd(a,b)=1. For(5,2) it is S_7, not the Monster or Baby Monster.

The same obstruction is visible analytically:

P_(a,b) = (P_(a/d,b/d))^d.

The factorization is exact because normalization constants match. Thus gcd(a,b) is a genuine common structure between polynomial decomposition and invariant monodromy blocks. It is an example of the kind of structure-preserving bridge sought in the research program.

### Sphere source is not spherical triangle monodromy

The source of every P_(a,b) is the Riemann sphere. That fact must not be conflated with the sign of the triangle-orbifold characteristic based on branch orders.

For(5,2), the branch permutation orders are(10,2,7), not(5,2,7). The corresponding triangle characteristic is

1/10 + 1/2 + 1/7 − 1 = −9/35 < 0.

So the natural triangle orbifold is hyperbolic despite the genus-zero source. The induced orbifold upstairs has cone orders2 at z=0,5 at z=1, and2 at each of the five simple preimages of1; its characteristic is

2 − (1−1/2) − (1−1/5) − 5(1−1/2) = −9/5 = 7(−9/35).

The polynomial has a spherical underlying surface and hyperbolic marked-orbifold structure simultaneously. These are compatible statements about different retained data.

The genus-zero passport check is also elementary: σ has2 cycles, τ has n−1 cycles, and στ has1 cycle, for a total n+2; the usual permutation Euler count gives Euler characteristic2. The verifier checks the passports for2≤a,b≤7 and enumerates generated groups independently whenever a+b≤8.

## 4. Logarithmic sheets: the exact analogy and its limits

The complex plane itself is single-sheeted. The logarithm becomes single-valued on the universal covering exp:C→C*, whose sheet coordinate changes by2πi under one loop around0. Cutting C* along the negative real axis chooses a branch of log; it does not make the negative axis an intrinsic singular line of the underlying plane. DLMF4.2 explicitly gives this branch convention: https://dlmf.nist.gov/4.2.

The selected five-cycle can be realized as finite monodromy of w^5=z, with five sheet states corresponding in order to[4=13],2,5,3,7. This preserves the cycle and winding action, but not integer multiplication, prime status, or group valuations. The seven-state completion above requires an extra branch generator to make its cover connected. These are precise positive versions of the spiral analogy.

The common mathematical operation is: retain a fibre/state address, prescribe its transport around loops, then ask whether a quotient forgets an address needed by another incidence. This is relevant to shared-port compilation as a design pattern. No map from these particular arithmetic covers to native tournament ports has been constructed, so it is not yet a theorem about tournament wiring.

## 5. Prior modular geometry and exact scope boundaries

The repository's THM-2996 already gives a genuine spherical/Euclidean/hyperbolic trichotomy. For integral matrices S=[[0,1],[1,0]], T=[[0,−1],[1,−1]], affine lifts of C2*C3 on F_p² have defect λ=2b_2−b_1. A nonzero defect produces cusp translation(st)²=λe_1, cusp order2p, and image F_p²⋊S3. The triangle groups(2,3,2p) are spherical at p=2, Euclidean at p=3, hyperbolic for p≥5. The output at p=13 has169 affine points, not4. Its bad-prime3 root/weight-lattice distinction is structural; it is not evidence for an arithmetic13=4 rule.

THM-4076 supplies another genuinely structural modulus: for odd q, the Stern inversion packet S(q) is −φ(q)+4·1_(ω(q)=1) modulo8. The next bit comes from fixed points of modular inversion. Its q=13/q=29 hostile shows why the fixed-point parity does not determine modulo16. Thus the correct connection between the mod8 and mod9 investigations is the requirement that a proposed quotient preserve the actual operation and its fixed or shared incidence data. Merely seeing8 and9 near each other gives no identification.

These sources and the two valuation notes are pinned in PROVENANCE.json. No claim here depends on treating a branch cut as a physical geometric defect, the divisor fold as modular reduction, or the Belyi completion as canonically selected by the Monster.

## 6. Exact realizability of prescribed divisor-height wiring

This extension and its proof were supplied by the session root, then independently tested here. It turns the constructive family into a complete finite constraint solver rather than a pattern-fitting exercise.

Let B be a finite set of bases b≥2, and prescribe integers k_b≥0. Allow N≥1. Write f_(p,b)=v_p(b), and define

L_p = max_(b in B) k_b f_(p,b),

where only primes dividing some prescribed base need be included. Put N_min=∏p^L_p. Then the following are equivalent:

1. Some integer N≥1 satisfies h_N(b)=k_b for every prescribed base.
2. h_(N_min)(b)=k_b for every prescribed base.
3. For each b there exists some prime p|b with L_p < (k_b+1) f_(p,b).

Moreover, whenever realizations exist, N_min is their unique least member in the divisibility order.

Proof. The lower bound h_N(b)≥k_b means e_p≥k_b f_(p,b) for every p|b. Combining all bases gives e_p≥L_p, so every realization is divisible by N_min. These inequalities already guarantee h_(N_min)(b)≥k_b. Exact equality requires one upper witness p|b with e_p<(k_b+1)f_(p,b). If N_min has no such witness, increasing exponents can never create one, so no realization exists. If each base has a witness at N_min, the displayed minimum realizes all targets. This proves all assertions.

The full family is explicit: choose e_p≥L_p subject to one of these upper witnesses for each base, and multiply by arbitrary powers of primes absent from every base. This is an intersection of upward lower bounds and disjunctive upper-witness conditions. It is not closed upward under arbitrary multiplication: multiplying by a prime used by a witness can destroy exactness.

For pairwise coprime bases all prescriptions are feasible, realized by N=∏b^(k_b). Shared prime coordinates create the actual compatibility obstruction. For instance h_N(2)=3 and h_N(4)=2 are individually possible but jointly impossible: the first fixes v_2(N)=3, while the second requires v_2(N)≥4. Replacing4 by6 repairs it: N=2³3² has heights3 at2 and2 at6. A second prime supplies a genuine alternative upper witness.

For the five-cycle of section2, the unique least realization is

N_min = 2^5 3^7 5^3 7^4 = 21,003,948,000.

All realizations have exactly these four prime exponents, and may contain arbitrary other prime factors. Adding the prescribed arrows13→2 and11→1 fixes those two exponents and recovers the full input N as the least realization. This explains precisely which part of the large group-order ratio is needed by the cycle and which part controls additional states.

The verifier compares the criterion against independent exhaustive exponent-box search for every three distinct bases from{2,3,4,6,8,9,12}, with each target in{0,1,2,3}:2240 systems,554 feasible and1686 infeasible. The search independently checks every pair of exponents0≤v_2≤11,0≤v_3≤8. These bounds include every potential least solution because the lower envelopes cannot exceed9 and6. It also verifies that every solution in the box is divisible by the proposed least solution.

All-zero targets have N_min=1. If a separate application insists on N>1, one may add a prime outside all bases, but then the least-N>1 assertion need not hold. This boundary is not suppressed.

### Connection to shared-port research: what is and is not proved

The common-prime coordinates are literal shared resources: separate feasible height assignments need not be jointly feasible. The theorem solves their compatibility by a common lower envelope plus per-base upper-witness disjunctions. This is a useful prototype for compiler design: carry the common coordinate instead of independently picking local witnesses. It does not prove that native tournament-port relations possess this minimum structure or these witnesses. A structural transfer would require a map identifying the actual port constraints with this exact constraint language.
