# Tournament growth: faster preparation, exact weighted recursions, and a compressed deletion atlas

**Session105, 2026-10-29.** This session improves actual generator preparation and the arithmetic realization of an existing unrestricted A000568 recurrence, recovers a stronger older deletion theorem, and connects it to exact refresh spectra and a compact torsion atlas. Independent audits preserve the distinction between inherited results, new deductions, and finite measurements. It does not complete the unrestricted prime-tournament atlas or close the native metagraph on endpoint classes alone.

## 1. Practical enumeration: improve the implementation, keep the strongest baseline

### Full rooted bank: avoid redundant canonical work

The [transport-seeded generator](generator/TRANSPORT_SEEDED_BANK.md) changes preparation in the existing Session88 adaptive block-fill implementation. A complete set of parent-to-canonical transports already determines the canonical parent's full automorphism group, so a second search is unnecessary. Combined with Session87's already-established mask-orbit traversal, this retains every rooted transport and all six old output fields.

At child order7:

- canonicalizer requests: **6,289 → 2,649**;
- actual canonical searches: **576 → 525**;
- adaptive root queries: unchanged;
- observed median runtime: about **1.05× faster**, not a dramatic speedup.

The complete bank, groups, alignments, pointed states and transporter sets agree through order7. The 81/27 hostile is retained. These are unrestricted ordinary tournament classes, not a regular/prime/converse-merged subfamily.

### Projective A000568 recurrence: exact integers instead of Fractions

The [integer recurrence](enumeration/INTEGER_PROJECTIVE.md) realizes Session98's fixed-degree projective gcd-response state with exact integer weights. For fixed N and remaining mass r, scaling the old weight W by N!/r! makes every history integral. The negative powers of two that arise after projection become **provably exact divisions**, checked before shifting. The smaller inherited odd-factorial scale also works. Packed nonnegative response deficits are carry-safe in the actual descending-length order; shared transition templates avoid repeated metadata construction.

This preserves every odd cycle type and the exact projective quotient. It is not a newly discovered Burnside formula or a stricter state quotient.

At N=200, three-run medians were approximately:

| Implementation | Output | Seconds |
|---|---|---:|
| Session98 Fraction projective | a_200 | 19.97 |
| New packed odd-local projective | a_200 | 1.95 |
| Session87 integer degree buckets | all a_0,...,a_200 | 1.38 |

Thus the new implementation is **10.2× faster than its Fraction predecessor but 41% slower than the stronger existing prefix engine**. Keep Session87 for production prefix queries. Fewer projective states alone do not guarantee faster arithmetic. Timings are observations on a shared machine, not universal performance claims.

## 2. Recovering the older deletion structure gives a much larger computational payoff

A targeted history search recovered [October19's Hasse/deletion theorem, §5](https://github.com/eliottcassidy2000/math/blob/c49370cdbf6c331a735a17c690a8e865beb3ac2b/05-knowledge/results/tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md). It already proves a source-free induced-subtournament count basis, complete integral deletion cokernels, and finite-field ranks. These are **inherited**, including the order-five Z/10 example, not new classifications based on today's small tests.

Write a_n=A000568(n), a_0=1, a_(-1)=0, and h_r=a_r−a_(r−1). The number h_r counts source-free classes: universal-source adjunction bijects all order-(r−1) classes with source-containing order-r classes. For the integer class-vector map D_n that sums all vertex-deleted cards, the recovered theorem implies

    coker D_n ≅ direct-sum_(r=0)^(n−1) (Z/(n−r)Z)^(h_r).

The displayed cyclic factors are not necessarily in Smith divisibility order. Our [compressed Smith implementation](updown/compressed_smith.py) sorts prime valuations as weighted runs and merges their breakpoints, without generating tournaments or expanding a diagonal of length a_(n−1). It uses only the exact A000568 prefix.

At order200, the deletion matrix rank has **5,558 decimal digits**, yet its complete invariant-factor list compresses into **60 runs**. The proved bound is

    1 + sum_(primes p≤n) floor(log_p n),

which is61 at n=200. The [full atlas](updown/smith_atlas_200.json) contains exact factors and multiplicities; a [compact summary](updown/atlas_summary.json) provides sizes and hashes. This is a compressed realization of the recovered classification, not a claim to have explicitly built the astronomical matrix.

## 3. A growth operator with a complete spectrum

Let U attach a new vertex in all 2^n ways and D sum all n vertex deletions, keeping multiplicities. With E_n=U_n/2^n,

    D E − E D = I.

Deleting the new vertex contributes the identity after normalization; deleting an old one loses one independent edge bit. The [up/down proof](updown/UPDOWN_HARMONIC_AUDIT.md) gives a rational decomposition by deletion harmonics H_j=ker D_j, with dimension h_j. This Weyl-algebra mechanism has classical graded-graph precedents, and the older source/deletion counterpart is explicitly credited.

For the genuine process “delete a uniformly chosen vertex, then reinsert with fair incident edges,” the class transition matrix is K_n=E D/n. Its full spectrum is

    1−j/n, with multiplicity h_j, 0≤j≤n.

Its stationary law is proportional to **1/|Aut(T)|**, not uniform on isomorphism classes. The spectral gap is exactly **3/n** for n≥3. The first nontrivial mode is the number of cyclic triples minus binom(n,3)/4: a refreshed vertex randomizes each affected triple to cyclic with probability1/4.

This makes binary intransitivity an exact spectral statement. It also gives the counting identity

    a_n=(n+1)trace(K_(n+1))−n trace(K_n).

That identity is a consistency relation and a new interface to counts, not a cheaper closed scalar recurrence unless the required traces are independently accessible.

The integral side must not be erased. For example D_5 has Smith factors1,1,1,10. The extra factor2 is visible in cyclic-triangle card parity, while the full all-degree answer comes from the recovered theorem. Rational full rank does not imply an integral or finite-field isomorphism. Source-free class indicators and signed vectors in ker D are also different carriers, despite their equal counts.

## 4. The normalized generating function acquires an exact doubled pole ladder

The [harmonic generating transform](harmonic/HARMONIC_GF.md) combines the recovered source-free count with Session98's global meromorphic theorem. Put

    beta_n=n! a_n / 2^binom(n,2),
    eta_n=n! h_n / 2^binom(n,2),
    B(z)=sum beta_n z^n,  J(z)=sum eta_n z^n.

Then

    J(z)=B(z)−zB(z/2)−(z²/2)B'(z/2).

Since B has exactly the poles4^k of order3k+1, J has exactly the poles **2^r**, of order **floor(3r/2)+1**. The old poles and differentiated doubled poles are disjoint, so the top orders cannot cancel.

Each paired core term has nonnegative Taylor coefficients. Its factor1−f·2^(1−c−f) has a direct Burnside interpretation: exclude the mutually exclusive universal-source events among f fixed vertices of a permutation with c nontrivial cycles. The zero/empty boundaries are explicit. This is an exact refinement of the symmetry-core enumeration, not a new modular-prime atlas or a Collatz convergence argument.

## 5. Native metagraph recursion: the missing alignment cannot be ignored

The [weighted native operator](metagraph/WEIGHTED_NATIVE_RECURSION.md) retains an actual directed Hamilton path and reverses only its nonpath arcs. Its symmetric path-orbit mass matrix M and fixed-source witness matrix W satisfy

    W_AB=|Aut(A)| M_AB,   M_AB=M_BA.

The H/Aut fibre identity is inherited. The new strip recurrence retains the complete normalized path word and commutes exactly with nonpath complement under bitwise-complemented extensions. An SCC specialization likewise retains the exit/entry bridges.

The proposed smaller state fails early: at order4, words4 and18 have the same source and target isomorphism-class pair, but different eight-choice extension distributions at order5. A second witness survives merging converse classes and forgetting endpoint order. This is a smallest tested/proved class-pair lumpability obstruction in the declared convention, not a proof that every aggregate recurrence is impossible.

Retaining the entire jointly aligned native pair does recover its directed path from the agreement arcs. It therefore retains the whole path word. That is a reconstruction-state lower bound, not a computational lower bound for A000568.

### A compatible exploration process, with its real cost stated

Conditioning on selecting a legal Hamilton path gives stationary mass proportional to H/Aut. The refresh process instead has mass proportional to1/Aut. These clocks cannot be mixed without accounting for the changed measure.

The repair is to propose uniformly from all n! vertex words and stay when the path guard fails. This native proposal process shares the refresh law and is at least half-lazy by path reversal. A mixture using refresh with probability alpha has absolute spectral gap at least **3alpha/n**. It is a genuine source-weight-preserving exploration operator, not proof that the native graph alone is connected.

There is no hidden native-only speedup: under labelled-uniform stationarity, proposal acceptance is exactly2^(−(n−1)). The gap guarantee comes from the refresh fraction. Uniform isomorphism-class sampling is not claimed.

## 6. Reproduce and interpret the evidence

All inherited baselines needed by the generator are vendored byte-for-byte with hashes. No external checkout or live network is needed. Python3.11+ and SymPy1.14.0 are used; the enumeration, native transport and compressed-atlas engines themselves use the standard library.

Run `python3 verify_all.py` for isolated normal/optimized replays and `python3 verify_manifest.py` for file hashes. The aggregate receipt records exact scope; benchmark times are deliberately excluded from structural equality comparisons. The [native/generator audit](audit_native/REPORT.md), [recurrence audit](audit_recurrence/REPORT.md), and independently enumerated up/down matrices provide separate checks.

The strongest improvements are exact preparation savings, a usable integer projective engine, a matrix-free compressed deletion atlas, a complete refresh spectrum, and a sharp obstruction preventing an incorrect native recurrence. None is advertised as a complete structural classification of all tournaments or a uniform replacement of the best existing counter.
