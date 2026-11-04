# A fixed-target projective quotient for exact tournament enumeration

**Status:** PROVED identities and quotient criterion; FINITE-EXACT implementation controls. No historical-priority or polynomial-time claim. This packet computes one requested A000568 coefficient, not an entire coefficient vector in one run.

## Inheritance and what changes

The primary [OEIS A000568 entry](https://oeis.org/A000568/internal) defines the number of unlabelled tournaments and records Davis's odd-partition Burnside formula. We use that formula as a cited input. The repository already went substantially further:

- `01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`: the soluble singleton tail.
- `01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md`: the gcd/divisor interaction.
- `05-knowledge/results/tournament_future_boundary_20261020/README.md`: merge histories with equal degree and all affordable future responses; exact n=200 computation.
- `05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py`: degree-local odd-factorial integral transfers, retained verbatim here as `inherited_local_scale.py`.

The fixed-target quotient itself is also inherited. [Session98, Sections3–4](../../session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md), already proves the identical if-and-only-if criterion for individual continuations, with s_j=R_j-jR_1, and absorbs2^(M R_1) into the weight. [Session105, Sections2–3](../../session105_metagraph_recursion_20261029/enumeration/INTEGER_PROJECTIVE.md), already proves the odd-local integer scale, exact negative shifts, and the same narrow packing bound. This packet independently implements and measures those mechanisms; it does not introduce a new projective quotient or a new exact-right-shift theorem.

The [October25 production defect recurrence, Section2](../../lca_ports_reverse_20261025/counting/README.md), retains exactly D(j)=j*t-R_j, omits t from the state, and already merges histories(23) and(9,9,5) before stage3 at target26. Its weight Y=O(d)*2^[M(t-b)]*w and this packet's W=O(d)*2^(Mt)*w satisfy W=2^(Mb)Y, where b=max_(j>=3)ceil(D(j)/(j-1)) is state-determined. That earlier gauge makes all shifts nonnegative; its reverse valuation potential additionally supports normalized low-bit queries. Thus the present state reduction relative to the ungauged comparator reproduces an existing quotient. The formulas, controls and measurements below remain valid. Completing the singleton tail analytically is inherited and is applied identically in both benchmark modes.

## 1. Exact-target continuation equivalence

Fix n. At a stage let mu be the processed odd cycle multiset, d its mass, t its cycle count, and M=n-d. Process cycle lengths in descending order, leaving singleton cycles until the final step. For a remaining affordable odd length j, write

    R_j(mu) = sum over a in mu of gcd(a,j),
    S_j(mu) = j*t - R_j(mu).

Here R_1=t and S_1=0. Only coordinates j>=3 with j<=M and not already processed are stored.

**Theorem 1.** For histories mu,nu at the same stage and same d, their interaction factors with every individual mass-M continuation are proportional by one common factor if and only if S_j(mu)=S_j(nu) for every retained j. The proportionality factor is 2^[M(t_mu-t_nu)].

**Proof.** A continuation with multiplicities k_j satisfies sum j*k_j=M. Its past/future cross exponent is

    sum k_j*R_j = M*t - sum k_j*S_j.

Equal S therefore gives the claimed common factor. Conversely compare the all-singleton completion with the completion containing one j-cycle and M-j singletons. If all continuation factors are proportional, subtracting their exponents gives Delta R_j=j*Delta t, equivalently Delta S_j=0. These completions are valid because singleton cycles remain available and j<=M. This is minimality for individual continuation responses up to a scalar, not minimality among arbitrary symbolic summation algorithms. QED.

A genuine merge appears at d=27 and remaining lengths {1,3}: mu=(27), nu=(9,9,9). Their old profiles (R_1,R_3) are (1,3) and (3,9), but both S_3=0. For n=30 their normalized continuation vectors coincide. This is not merely dropping a useless coordinate: both mass-three completions, 1^3 and 3, remain possible.

**Hostile control.** Histories (11,5) and (9,7) have equal mass and equal cycle count but S_3=4 and 2. Their three-cycle completion responses differ. Dropping S_3 is invalid.

**Boundary.** The argument uses exact target mass and the available singleton completion. It does not directly apply to simultaneous all-degree output, no-fixed-point restrictions, or signed core-substitution excess. In particular, excess contributes degree but is not an actual cycle; silently using this recurrence for pair-free cores would be wrong.

## 2. Integer recurrence and exact right shifts

Let O(d) be the odd part of d!. For a state of degree d and projected response S store

    W[d,S] = O(d) * sum_mu w(mu)*2^[(n-d)*t(mu)],
    w(mu) = 2^e(mu)/z(mu).

The sum is over histories in that state. Set W[0,0]=1. For k cycles of current odd length ell>=3, let q=d+k*ell and M=n-d. The updated signature is

    S'_j = S_j + k*(j-gcd(ell,j)),

followed by deletion of coordinates j>n-q and of the processed coordinate. Define

    T = O(q)/(O(d)*ell^k*O(k)),
    E = k*(ell-1)/2 + ell*binom(k,2) - v2(k!)
        + k*(M-k*ell) - k*S_ell.

Then add T*W*2^E to the successor weight. Unlike the inherited ungauged recurrence, E can be negative. This is an exact division by a power of two, not rounding.

**Proof.** In the ungauged Burnside weight the exponent increment is k*(ell-1)/2+ell*binom(k,2)+k*R_ell. Replacing the old normalization M*t by the new one (M-k*ell)*(t+k) changes this by

    -k*ell*t + k*(M-k*ell).

Use R_ell=ell*t-S_ell and remove the 2-part of k!. The displayed recurrence follows. Integrality of T is inherited from q!/(d!*ell^k*k!) being a permutation-class extension count; taking odd parts gives its stated divisibility.

Each individual O(d)*w(mu) is an integer: z(mu) divides d!, while the 2-part of z(mu) is cancelled by 2^e(mu), because e(mu)>=sum_j binom(m_j,2)>=sum_j v2(m_j!). Thus each normalized successor is integral. The multiplier is common to all histories in a projected state by Theorem 1, so a negative shift divides the aggregate exactly as well. The executable checks every such division.

After processing ell=3, complete the unique singleton tail r=n-d. With no S coordinates left, add

    W * [O(n)/(O(d)*O(r))] * 2^[binom(r,2)-v2(r!)]

to a common numerator. Divide the final numerator by O(n). This equals a(n). All displayed divisions are checked, and no floating-point arithmetic participates in the answer.

## 3. Packing and complexity

The production `packed_count` stores nonnegative S_j coordinates in fixed-width binary digits. Every processed part a is larger than every remaining length j. Therefore j-gcd(a,j)<a and

    0 <= S_j = sum_a (j-gcd(a,j)) <= d <= n.

A digit width bit_length(n), at least one, prevents carries between coordinates. This sharper bound is essential: the naive bound n^2 wastes bits. Coordinates are ordered by increasing odd j>=3, updated by one precomputed packed vector, and masked by the remaining mass. The ungauged control uses the inherited R_j coordinates including R_1 with the same digit width.

The exact negative-shift check uses a low-bit mask, not long integer division. Since every T is odd, a negative shift is applied to W before multiplying by T. This preserves exactness and reduces the size of the integer product. Both optimizations are included in the final executable and receipts.

At every stage and fixed d, the projective signatures are images of the ungauged signatures under R->S. Hence the projective number of states is no larger, and using the same allowed k and tail closure the number of mathematical transitions is no larger. This does **not** prove lower wall-clock time: weights can be larger, and right-shift divisibility checks cost work.

Writing P_ell for encountered projected states before ell, the transition bound is the explicit sum over ell and (d,S) in P_ell of 1+floor((n-d)/ell). Transfer-table size is O(n^2 log(n+2)), as in the inherited ordinary local-scale method. Integer weights use O(n^2+n log n) bits. No polynomial bound on the number of response states is proved.

## 4. Reproduction and validation

Run from this directory:

    python3 check.py receipt.json
    python3 -O check.py receipt-optimized.json

The checker uses explicit exceptions, so optimization does not disable controls. Timings are intentionally not expected to match between runs. Compare answers, state counts and transition counts instead.

Controls include all 77 values of the cached primary OEIS table (0..76); direct independent rational odd-partition summation through 35; tuple versus packed implementations through 76; same-target ungauged versus projective through 76 and at every benchmark; inherited all-degree integer DP at n=100,150,200; genuine merge and hostile witnesses; and individual completion identities. The supplied OEIS cache is inherited from the Oct-20 packet; the current entry was checked again for the definition and formula, not claimed as a newly downloaded table.

Benchmarks alternate projective and ungauged execution order over three repetitions and include setup, exact divisibility checks and singleton-tail closure. They compare **one exact coefficient against one exact coefficient**, using the same normalization code and closure. Comparison with the inherited all-degree API is only an equality control, not a fair output-matched timing claim.

The quotient preserves the aggregate a(n) and every normalized individual completion weight. It loses cycle histories, automorphism identities, tournament representatives and Hamiltonian-path counts. It cannot by itself recover the Hamiltonian-path spectrum or identify which isomorphism classes realize a given path count. Those require separate structure or an enriched enumerator.

## 5. Recorded measurements and independent audit

Both normal and optimized-mode runs passed 537 checks with identical exact answers, state/transition tables and source hashes. Three-sample alternating timing ratios (ungauged fixed-target / projective) were:

- n=80: states 1,227 -> 1,012; transitions 21,833 -> 18,031; timing ratios normal 0.875, optimized-mode 1.050.
- n=100: states 2,953 -> 2,721; transitions 61,947 -> 49,946; timing ratios normal 1.006, optimized-mode 1.115.
- n=150: states 20,052 -> 18,207; transitions 506,659 -> 408,188; timing ratios normal 0.966, optimized-mode 1.087.
- n=200: states 96,836 -> 95,709; transitions 2,770,773 -> 2,303,050; timing ratios normal 1.028, optimized-mode 0.954.

Ratios below one mean slower. The two replay measurements change the sign of the n=200 timing comparison, so there is no demonstrated robust runtime advantage. The transition reduction is exact and reproducible; timing is host-dependent. The final implementation already incorporates narrow packing, low-bit divisibility tests, and shifting before odd multiplication. Earlier wider-packing spot timings are not the final benchmark.

The inherited all-degree engine remains appropriate when a whole prefix is requested. The inherited projective quotient is structurally coarser than its ungauged comparator. This independent implementation and measurement packet establishes neither a fastest algorithm nor a solution to the partition-growth barrier.

Independent audit: `../orbits/independent_results_normal.json` and `../orbits/independent_results_optimized.json` pin the current projective source SHA-256 `1a37aa21e66f820473b14f38c544c021642b8a7436d3e970588cfb8565a0f63a`. Those audits check six implementations through 30, a direct permutation action on oriented edges through 7, 600 projective-equivalence cases, and reverse-recursion controls at 50,100,150,200. No defects were found.

`python3 memory.py` measures Python allocation peaks with tracemalloc (not process RSS). At n=150 the final projective kernel used 10,061,760 traced bytes versus 9,719,012 for its ungauged fixed-target control. Thus fewer states did not yield a memory improvement in this measurement either. Exact answer equality and source hash are in `memory-receipt.json`.
