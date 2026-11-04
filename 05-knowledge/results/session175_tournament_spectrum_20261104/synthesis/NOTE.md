# Session 175: exact enumeration quotients and the actual spectrum frontier

**Status:** proved structural statements, independently audited implementations, and explicitly bounded computations. No polynomial-time enumeration theorem, robust speedup, new globally attained Hamilton-path value, or resolution of spectrum completeness is asserted.

## 1. Correct the target before constructing

The suggestion at the end of [Session 174](https://github.com/eliottcassidy2000/math/blob/a55a46e2a9fd3880980b26149c0c6680f3ec8909/00-navigation/S174.md) to seek ordinary Hamilton-path counts 7 or 21 by changing the wrapper architecture was misguided. **Both values are already globally impossible.** Changing a constructor cannot produce either. S174's proved wrapper identities and norm obstruction remain valid; its proposed global target does not.

The open question is whether the ordinary Hamilton-path spectrum is exactly the positive odd integers other than 7 and 21. This was already carefully distinguished in the repository's [pinned spectrum frontier](https://github.com/eliottcassidy2000/math/blob/a55a46e2a9fd3880980b26149c0c6680f3ec8909/05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md), and is explicitly Conjecture 6.2 in [Chan–Pak, §6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf). [Grinberg's April 27, 2026 notes, Remark 1.1.11](https://darijgrinberg.gitlab.io/t/22s/lec12.pdf) still describe the higher-value question as unknown.

The distinction is supported by a short exact argument. Hamilton counts multiply over ordered strong components. [Busch's exact strong minimum, Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), starts 3,5,9,15,25 at orders 3 through 7. Therefore a strong witness to 7 would have order at most 4, and one to 21 order at most 6. The small spectra exclude both. A component product of 7 or 21 necessarily requires a whole factor 7 or 21, so neither occurs at any order. This is inherited theory, not a result of the new wrapper experiments.

[Royle's primary computation report](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament) covers every allowed odd value **strictly below 80,557**, equivalently through 80,555, by order 11. We do not include 80,557 in that certified prefix and do not infer that it is missing.

The inherited completeness reduction has two genuinely infinite obligations:

1. Every odd prime p other than 7 has a strong witness.
2. Every 7p, for an odd prime p other than 3, has a strong witness.

They are necessary by strong-component factorization and the absence of a factor 7. Together with known strong carry witnesses 49,63,343 they are sufficient: factor an allowed odd number as 7^a m. Even a uses 49; odd a≥3 uses 343 and 49; a=1 either uses 63 when m is a power of 3 with exponent at least two, or uses 7p for a prime divisor p other than 3,7. The remaining factors use the first lane. This does not supply the missing infinite witnesses.

Keep three quantities separate: A000568 counts tournament isomorphism classes; H(T) counts directed Hamilton vertex words; H(T)/|Aut(T)| counts path orbits. A result about one is not automatically a result about another.

## 2. A000568: a stricter exact quotient, not a demonstrated faster algorithm

The [enumeration packet](../a000568/README.md) starts from the [Davis odd-partition Burnside formula recorded by OEIS](https://oeis.org/A000568). The singleton-tail split, gcd-response quotient and odd-factorial integral transfers were already present in the repository. The change here is to exploit the request for **one fixed target coefficient**.

At a descending-cycle-length stage, let a past history have mass d, cycle count t, and future response R_j=sum gcd(a,j) over its past cycles. With exact remaining mass M=n−d and singleton cycles still available, two histories have proportional interaction factors on every individual completion precisely when their response differences lie in the mass direction:

    R'_j−R_j = j c.

Thus S_j=j t−R_j for j≥3 is a complete invariant of this wordwise projective equivalence. Absorb the removed scalar 2^(M t) into the historical weight. The exact forward recurrence then permits negative binary exponents; every corresponding division is proved integral and checked. Omitting the compensation would be wrong.

A genuine merge occurs for histories (27) and (9,9,9), at mass 27 with future lengths 1 and 3: different cycle counts and old responses, but the same projected response. The theorem is minimality for individual completion responses up to a common scalar, not optimality among all possible symbolic algorithms. It depends on the exact remaining mass and the available singleton completion; it must not silently be applied to an all-degree output or another species.

The [independent orbit packet](../orbits/README.md) supplies a backward factorial-scaled implementation with a different projective normalization, plus a signed permutation action on edges that does not assume the odd-cycle/gcd formula. Controls cover six representations through n=30, all 77 cached OEIS values through n=76 in the main suite, and exact comparisons at n=100,150,200. The new quotient preserves aggregate enumeration weights; it discards tournament representatives and Hamilton-path data.

**Measured outcome:** at n=200, the fixed-target peak state count falls from 96,836 to 95,709 and transitions from 2,770,773 to 2,303,050. Those reductions are reproducible. The ratio of ungauged time to projective time is 1.028 in one normal run and 0.954 in the optimized-mode replay, so the apparent advantage changes sign. There is no demonstrated robust runtime win. At n=150, traced Python allocation peaks are 10,061,760 versus 9,719,012 bytes, also not a memory improvement. These are output-matched fixed-target comparisons, not a claim against the inherited whole-prefix API or the external state of the art.

The [orbit packet](../orbits/README.md#3-no-finite-order-polynomial-coefficient-recurrence) also proves a recurrence-class obstruction: A000568 is not P-recursive, and its ordinary generating series is not D-finite. Indeed a(n)≥2^(n(n−1)/2)/n!, whereas a fixed-order polynomial-coefficient recurrence permits at most exp(O(n log n)) growth. This excludes that recurrence class, not exponential-coefficient, nonlinear, growing-order, or state-dependent recurrences. No polynomial bound on the number of projective states has been proved.

## 3. Constructor progress: explicit families and exact limitations

The [construction packet](../construction/README.md) treats every tournament made by adding one vertex to a transitive order. Encode the new vertex's wins by a binary word. A three-register streaming recurrence gives H exactly. Removing initial zero sources and final one sinks preserves H, leaving a strong core beginning in 1 and ending in 0. For a core of length n,

    H ≥ 2^(n−2)+2.

Consequently a finite word search is a complete inverse decision procedure **for this entire constructor class**, even with unlimited source/sink padding. Its exclusions include 13,29,35,41,43,49,55,57, which are not global exclusions. This is a useful restriction on a proposed algorithm rather than another fixed-order failure to find examples.

The same packet derives the asymmetric uniform wrapper A→B→z→A:

    H = (F_A(2)−H(A))(F_B(2)−H(B)) + 2 H(A)H(B).

Choosing a transitive r-vertex A and the literal four-vertex seed B=(8,1,3,6) yields the strong family

    H_r = 27·2^r−17.

Examples include 37,91,199. For r=3k+2, every value is divisible by 7; its 7-adic valuation is exactly one whenever k is not 5 modulo 7. This proves infinitely many escapes from the mirrored norm's even-valuation restriction. It does **not** prove infinitely many prime cofactors or all 7p values.

The [cyclic-module obstruction packet](../obstructions/PROOF.md) supplies a positive ordered-path-cover kernel for arbitrary uniform C3[A,B,C], with no mirror assumption. At fixed positive module sizes, transitive interiors uniquely minimize H. Its complete low spectrum is

    {H(C3[A,B,C]) ≤45} = {3,5,9,11,15,17,23,27,33,37,45}.

This is an all-order constructor classification because the proof bounds possible module sizes before completing a finite case analysis. It shows that breaking the mirror alone cannot make this architecture realize even such permitted small targets as 13,19,25,29,31,35,39,41,43. There is no contradiction with the asymmetric H=91 family: 91 is above the classified range, and the norm restriction is stronger than the low-spectrum restriction in a different direction.

## 4. Nonuniform interfaces: retain endpoints and orphans

The [interface packet](../interface/PROOF.md) changes actual intermodule edges. For endpoint counts E(S,u), B(S,u), define

    C(u,v)=sum_{S subset V\{u,v}} E(S+u,u) B((V\({u,v} union S))+v,v).

Reversing u→v changes H by C(v,u)−C(u,v). Each factor omits one of the changed edge's endpoints, so both contributions can be computed from the original endpoint tables. This is a literal split of paths that use the changed edge.

For vertex insertion, legal positions on old Hamilton paths are not enough. Some new paths delete to a broken shortcut. These orphan paths have x→new→y but old y→x and are counted by the corresponding endpoint subset product. The construction and interface packets derive the same correction with opposite bit conventions, explicitly defined in each. The transitive three-vertex example has two legal-slot paths and three orphans, hence H=5.

Three fully enumerated seven-vertex interface cubes supply 1,536 labelled states, not an order-seven census. Every state has an independent literal 5,040-permutation count. Additional independent review verifies strong connectivity, the exact maximal response-component partitions, and all 30 prime/7p lane values. Literal strong witnesses include 49,77,91,133,157. These values were already globally known; the deliverable is their explicit interface construction and verified local responses.

Restrict cube edges to |ΔH|≤2. Each connected component realizes a solid odd interval by the inherited discrete intermediate-value argument. The substantive findings are the actual component sizes and ranges, including an eleven-state strong TT component covering 153,155,157,159. Fragmentation prevents treating a whole cube as one interval and supplies no unbounded coverage theorem.

## 5. Why cloning these successes is insufficient

The [clone obstruction](../clone_obstruction/PROOF.md) closes a natural scaling attempt. Fix K quotient tournaments of order at most d, allow arbitrary nonempty transitive inflation of their vertices, and then allow reversal and any number of ordered sums. The number of distinct prime H values at most N produced by this whole class is at most

    K d floor(log_2(N−1))^d,  for N≥3.

Each vertex of a nontrivial strong quotient lies in a directed triangle. Inflating that vertex to size r and retaining one vertex in the other two modules gives an induced tournament with 2^r+1 paths. Insertion monotonicity bounds r by floor(log_2(N−1)). Only polylogarithmically many parameter tuples remain. For prime totals, component multiplication allows just one nonunit factor. Reversal adds no values.

An elementary central-binomial-coefficient argument gives at least a constant times N/log N primes up to N. Therefore any fixed finite clone library covers a proportion of primes tending to zero. This applies in particular to the finite interface library above. It does not exclude infinitely many prime outputs, prove any particular prime globally missing, or obstruct a construction with a growing number of interface roles or genuinely richer internal modules.

The positive interface experiments and negative clone theorem fit together: local norm escape is real, but it cannot be promoted to completeness merely by lengthening transitive twin blocks.

## 6. Next targets justified by the results

1. **Enumeration:** keep the exact projective equivalence, but investigate a change that reduces integer-operation cost or state growth. Benchmark one coefficient against one coefficient, with repeated timings and allocation measurements. State/transition reduction alone is not a runtime theorem.
2. **Construction:** grow the number or structure of nonuniform interface roles. Derive a closed endpoint-sensitive recurrence including orphan paths, rather than adding an unproved scalar transition law.
3. **Interval coverage:** exhibit connected realized strong response components across increasing orders and prove that their odd intervals overlap. A response gcd of two or isolated finite intervals is insufficient.
4. **Prime lanes:** supply an explicit strong witness mechanism for every required prime and 7p. Infinite exponential families, finite seed products, and fixed-template transitive inflation do not close these obligations.

## 7. Verification map and limits

- [Enumeration derivation and measurements](../a000568/README.md); [independent orbit audit and recurrence obstruction](../orbits/README.md).
- [Constructors](../construction/README.md) and [independent audit](../construction/AUDIT.md).
- [Uniform cyclic classification](../obstructions/PROOF.md) and [independent audit](../obstructions/AUDIT.md).
- [Endpoint/interface proof](../interface/PROOF.md) and [independent reviewer audit](../interface/REVIEWER_AUDIT.md).
- [Clone obstruction proof](../clone_obstruction/PROOF.md) and [reviewer replay receipt](../clone_obstruction/reviewer_replay.json).

Normal and optimized-mode mathematical receipts agree within the packets; timing fields are not expected to be byte-identical. General statements rest on proofs, finite experiments validate only their declared universes, and external published results retain their attribution. None of the enumeration quotients identifies Hamilton witnesses, and none of the constructor packets settles the unrestricted spectrum.
