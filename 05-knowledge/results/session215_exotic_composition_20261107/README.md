# Session215: positive tournament constructions and their composition

**Status: PROVED scoped construction laws, with independently audited finite-exact certificates. The Chan–Pak conjecture that every odd integer at least 23 is a tournament Hamiltonian-path count remains OPEN.** No enlarged unrestricted finite-coverage range or literature-priority claim is made.

This session moves beyond Session214's fixed source pairs and restricted two-apex surgery. The useful discoveries are new persistent paired sources, an all-order explanation of homogeneous-block transfers for arbitrarily many apexes, and a common formula connecting those blocks to cyclic substitution with arbitrary tournament seeds.

## 1. A sink repairs two new paired sources permanently

[Construction and proof](carriers/PROOF.md), [literal witnesses](carriers/NEW_CARRIERS.json), [independent audit](audit/CARRIERS_AUDIT.md).

The marked state is (X,e,V,d), where X is the adjusted base, e is the root-ending path count, V=R+b retains two-path-cover data, and d is the count after deleting the root. Two new physical pairs have common states

    (4225,0,8910,195), with starts 25 and 11;
    (4357,0,9836,251), with starts 13 and 14.

Their start labels satisfy the required staggered budget j<=8. Each pair therefore gives both actual counts 3X+2 and 3X+4, initially 12677/12679 and 13073/13075.

The mechanism is constructive. Different sources agree in X,V,d after J-containing prefixes, but disagree in e. Appending a global sink O resets e to zero without changing X or d. The two resulting states now agree and remain synchronized after **every common word over O and Z**. The new pairs are not obtained by forcing a common J through the four old pairs, which Session214 ruled out.

In particular, every k>=1 gives paired actual counts 3X_k+2 and 3X_k+4 at

    X_k=8715*2^k-4295,
    X_k=9585*2^k-4977.

The two initial bases are absent from the four previous O/Z scalar orbits, by the inherited finite-order cutoff and an exact depth-six check. Every new descendant **four-state** stays outside the four previous state orbits because d is invariant and different. This does not claim that all descendant scalar values are new, or that these values were previously missing from the unrestricted spectrum.

The packet also preserves earlier finite asynchronous collisions, including X=8957, and explains exactly why matching X alone did not provide persistent continuation. Widening the source bank, rather than extending those failed equalities indefinitely, produced the positive families.

## 2. Arbitrarily many apexes have an exact positive block compiler

[Proof and verifier](multi-exit/README.md), [independent audit](audit/MULTI_AUDIT.md).

Fix any k-apex tournament A and a transitive chain. Every chain vertex may have any of the 2^k incidence masks against the apexes. Hamilton paths are counted by assigning chain vertices to the k+1 gaps of an apex permutation, retaining each gap's first and last boundary conditions.

For a fixed incidence letter c, the individual gap-assignment operators E_r are commuting idempotents. Consequently its length-m homogeneous block has the **nonnegative exact expansion**

    M_c^m=sum_U |U|! S(m,|U|) product_(r in U) E_r.

This proves, for every k, that a homogeneous block has only exponential bases 0,1,...,k+1. There are no polynomial factors in m. Different incidence letters need not commute; arbitrary blocks compose through their complete boundary states. The transfer method and the two-apex eigenvalues were inherited. The structural all-k idempotent/Stirling explanation is the contribution here.

For any k>=1 and any k-apex tournament, the word consisting of k vertices dominated by all apexes followed by m vertices dominating all apexes has highest exponential term

    (k!)^2 (k+1)^m.

For a directed-triangle apex and k=3 the exact formula, for m>=1, is

    H=87-189*2^m+72*3^m+36*4^m.

Its m=0 value is separately 3. The fourth exponential mode separates this parameterized family from every fixed-context homogeneous family with at most two apexes. It does not establish that any individual value needs three apexes. This symmetric example is always divisible by 3, so it is not itself a completeness route.

## 3. Cyclic composition replaces integer bases by physical seed data

[Proof and replay](fixed-caps/README.md), [independent audit and provenance correction](audit/FIXED_CAPS_AUDIT.md).

Let F_T(j)=j! times the number of spanning covers by j unordered directed paths, and let Q_T(t)=sum_j binom(t,j) F_T(j). For fixed nonempty cyclic caps A,C, the inherited directed-triangle kernel gives explicit integers beta_j such that

    H(C3[A,D_1=>...=>D_r,C])
      =sum_(j=1)^(min(|A|,|C|)+1) beta_j product_i Q_Di(j).

Thus a fixed finite bank of arbitrary tournament seeds has a bounded-register exact compiler. Repeated seeds give constant-coefficient recurrences. The beta_j may be negative; the underlying path-cover formula has nonnegative coefficients and supplies the physical interpretation.

An explicit five-vertex seed has H=11 and Q(2)=100. Singleton caps around r ordered copies give a **strong tournament on 5r+2 vertices with H=100^r+11^r for every r>=1**.

This formula is a concrete connection between the two composition lanes: on the two-block cyclic specialization, singleton seeds give Q(j)=j and reproduce the block compiler's integer modes, including exactly (87,-189,72,36). Arbitrary seeds replace those modes by Q_D(j). No shared vocabulary or numerical resemblance is needed; both counts come from the same maximal-block-run bijection.

Symmetric caps C3[TT_m,T,TT_m] also give triangular path-cover observations with diagonal (m!)^2. These are an alternative explicit probe family. **Physical-context profile recovery and compressed known-order recovery already appear in Session176 and the later proof-principles packet.** Their prior results are credited; no generic new tomography theorem is claimed. Likewise, fixed-depth recurrences and the universal substitution kernel are inherited.

## What composes, and what still does not

- Arbitrary mask blocks compose using run boundary states. A scalar H cannot replace them.
- Ordinal seed words compose using the full path-cover polynomial, or the finite evaluations needed by specified fixed caps. This deliberately forgets the order of the seed blocks, not their path multiplicities.
- The new paired sources compose under O/Z using equal marked four-states. An arbitrary cyclic substitution need not preserve those marked equalities; unmarked profile composition is not automatically a paired-source operation.
- Fixed caps keep a bounded jet. Letting all three cyclic factors grow reintroduces the inherited moving-jet obstruction; these recurrence formulas do not evade it.

The next missing step is an all-target inverse construction or source-selection theorem. These exact infinite families and richer composition laws do not prove eventual coverage, F(8,137), or the full all-odd spectrum. No inherited unrestricted finite-coverage range is extended here.

## Verification and research lineage

All three mathematical packets have signed independent audits, source pins, reproducible programs and exact receipts. The all-k compiler has 7543 subset-DP controls, 43670 operator-algebra checks, 1200 EGF comparisons and 75 highest-mode controls, plus independent permutation tests and mode checks through k=6. Fixed caps have 160 physical cyclic-composition checks and independent profile/recovery controls. New carriers have independent full nine-coordinate and literal compiled-graph counts; the discovery search replays 12893 filtered roots, 341727 stored keys and 20 selected paired keys without assuming bank completeness.

Run each packet's README commands and the independent scripts listed in its audit. Normal and optimized Python receipts agree. Root SHA256SUMS covers this complete release except itself. RELEASE.json records scope and the earlier verified checkpoints; final publication is verified after the containing commit exists.

The research move was to retain the exact continuation information each operation needs, change the source class when its invariant blocked growth, and identify the same positive count in both the gap-assignment and cyclic-substitution descriptions. In particular, the sink's loss of one coordinate becomes useful only after the other three have been matched by actual graphs.
