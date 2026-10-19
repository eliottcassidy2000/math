# Independent review: saturated barriers and two free ports

19 October 2026. Verdict: **PASS for the explicitly stated induced-seed and two-free-port families.** No substantive proof defect was found. This review does not endorse a general seven-thin Gold theorem, the general balance conjecture, or any historical novelty claim.

Reviewed source: the five-file research package identified by the original hashes in `SOURCE_PROVENANCE.json`. The proof received a status-only update after this review.

## 1. Structural proof audit

The saturated-barrier lemma is valid for arbitrary finite induced exteriors, not only exterior chains.

- A degree-saturated core vertex has no exterior incomparable vertex in a k-thin ambient poset.
- If an exterior vertex lies on opposite sides of two incomparable members of S, transitivity contradicts their incomparability. Connectivity therefore forces one orientation across all S.
- Nonempty S gives A < B immediately.
- Domination is exactly what rules out a reverse relation q < a for a in A: choose s in S incomparable with q, and q < a < s is impossible. The upper argument is dual.

The same proof works when outside comparability with every S vertex is an explicit hypothesis instead of a consequence of saturation. Both hostile examples in the packet correctly identify necessary hypotheses of this particular structural lemma.

For the displayed Q, the exact degree tuple and S are correct. The remaining lower ports are minimal and lie below members of S; the upper ports are maximal and lie above members of S. Thus A is below the upper ports and B is above the lower ports. In particular, the later restriction to three coordinates on each side is justified structurally.

As a separate finite control, the reviewer enumerated all 46,656 possible one-point relation specifications comparable with S, without first imposing the proposed A/B orientation. Exactly 16 specifications extend Q to a poset; every one has the required one-sided orientation, and all 16 are seven-thin. This is a control, not the proof of the arbitrary-exterior lemma.

## 2. Compression and transfer audit

Condition on an extension order of A and an extension order of B. Because all cross relations have the directions established by the barrier, there are no omitted reverse cross constraints.

For a fixed core word, the number of A letters following each core vertex before the first S vertex is weakly decreasing and bounded by the terminal predecessor gaps. It determines a unique merge with the fixed A word. The upper side is the reversed dual. Letters preceding every permitted lower gap, or following every permitted upper gap, have a forced position and contribute a factor of one. Since A < B and S is nonempty, the two gap assignments do not interact. This validates the product formula and the event-count identity for arbitrary exterior posets.

The degree budgets are also valid: every A letter after q's last predecessor is incomparable with q, and dually for B. Every positive lower profile coordinate is incomparable with the last A letter, which bounds the lower support by k; the first B letter gives the upper support bound. Therefore the replacement chains satisfy the ambient degree bound.

Crucially, the SAME four cell inequalities are certified for every generator. Their transfer uses linearity of unnormalized cell counts, with extension-total weighting after normalization. It does not use the invalid inference that individual generators having potentially different successful trees implies one successful mixture tree.

The reviewer independently checked the full four-cell conditioned-mixture identity on 48 actual non-chain exterior posets, representing 924 pairs of exterior extension orders. Every identity held. These auxiliary examples test the general compression identity; not all 48 are asserted to satisfy the extra degree hypotheses of Theorem A or B.

## 3. Independent exact bounded certificate

The review implementation does not import any packet code. It constructs relations as sets of ordered pairs, takes transitive closure, and uses one four-vector dynamic program that records comparison orientations while growing ideals. This differs from the packet's four separate constrained-extension counts. Explicit failures remain active under Python optimization.

All 240 bounded atoms were reconstructed and recounted. Every four-cell count agrees with the packet, certified by agreement of the entire serialized count-table hash:

`b1c8c67270896d932f81472a78092beab6d6f6e39c1cd09ecdeeb342d8e6423a`

The independently reproduced cell-margin minima, in TT, TF, FT, FF order, are:

- 78113/453252
- 25363/180473
- 50695/453252
- 4768/24065

The first-comparison endpoints are exactly 66465/151084 and 84619/151084.

The stated FT witness has E = 453252 and cells (102269, 97126, 148700, 105157), hence margin 50695. It is an actual seven-thin induced-Q model. Thus 50695/453252 is attained and is sharp for this fixed tree in the induced-Q seven-thin family. Since 9 × 50695 = 456255 > 453252, the strict E/9 conclusion follows.

An additional 16 actual non-chain ambient posets satisfying all Theorem A hypotheses were counted directly. Each satisfied the asserted sharp bound.

## 4. Free-port proof and certificate

The affine-tail formula is exact. If the free lower vertex is not first in the boundary word, an earlier cap at most r bounds its coordinate, so increasing the free cap above r has no effect. If it is first, each new allowed coordinate value has the same bounded suffix choices. This is a positive affine expansion, not interpolation or a numerical extrapolation.

Minimality identifies the coefficient with extensions of the core after deleting that initial vertex. Maximality supplies the dual terminal deletion. The compared vertices 7, 9, 4, 5 avoid the deleted ports 2 and 10, so all four comparison cells transfer. Multiplication of the two expansions gives nonnegative coefficients for the four base/deletion terms.

All 144 generators were independently reconstructed and counted, including the double deletions. Their entire count-table hash matches:

`181cc53fb0202d38f5926c5b77217ba16e5ae7be145c384b96c1bb613f4b8e64`

The reproduced TT, TF, FT, FF minima are respectively 771/6245, 3231/24464, 386/6245, and 3116/16323. The reviewer also counted 112 larger-cap models directly, spanning all sixteen bounded coordinate choices and seven asymmetric or simultaneous free-cap choices. Every four-cell positive-cone identity held.

For every finite model, the positive base term retaining both free ports has all margins strictly above the proposed infimum, while every deletion term has margins at least the infimum. Therefore the final inequality is strict even though a deletion generator attains the cone minimum. Summation over finitely many exterior extension orders preserves strictness.

An additional 16 non-chain examples satisfying Theorem B and lying outside seven-thinness were counted directly. Their maximum incomparability degrees were 10 or 11, and all satisfied the strict 386/6245 bound.

## 5. Symbolic sharpness check

The reviewer derived the sharpness polynomials from the actual base, single-deletion, and double-deletion count vectors. In descending powers of n, the four cell polynomials have coefficients:

- TT: (3091, 18254, 26113)
- TF: (2281, 16367, 29046)
- FT: (4600, 26408, 36916)
- FF: (2518, 17755, 30813)

Consequently, for every integer n ≥ 1,

E(n) = 12490 n² + 78784 n + 122888,

m_FT(n) = 772 n² + 8213 n + 18243.

There is also an exact symbolic strictness check:

6245 m_FT(n) − 386 E(n) = 20879561 n + 66492767 > 0.

The leading-coefficient ratio is 772/12490 = 386/6245. The models satisfy the explicit Theorem B conditions and have unbounded degrees at the free ports. Thus the sharp infimum is justified by genuine finite posets approaching it. No claim of optimality across all possible comparison trees follows or is needed.

## 6. Scope and publication assessment

The two results are all-exterior theorems under structural hypotheses on a single induced seed, with exact fixed-tree constants. They are stronger than checks restricted to the unweighted seed, but remain conditional on containing the specified Q and on the ambient hypotheses in the relevant theorem.

The review supports:

1. The general saturated-barrier decomposition lemma as stated.
2. The exact compression/event-transfer argument.
3. The fixed-tree induced-Q seven-thin theorem with sharp constant 50695/453252.
4. The two-unbounded-port theorem with sharp strict infimum 386/6245.
5. The displayed extension and margin polynomials.

It does not support any assertion that all seven-thin posets contain Q, that all other core morphologies are handled, that arbitrary positive mixtures are realizable, or that the general seven-thin Gold/balance conjectures are solved. The packet already states these limitations appropriately. Its exploratory ten-point refinement remains outside the proved scope.

No mathematical revision is required by this review. An optional clarity improvement is to state explicitly near each sharpness claim that it is sharp for the fixed tree (1); this qualification is already present elsewhere in the packet.

## 7. Reproduction and artifacts

- `independent_controls.py`: standalone independent counting and structural/mixture controls.
- `independent-results.json`: full results, including all non-chain examples and end-to-end cases.
- `independent-count-tables.json`: complete independently reconstructed 240-atom and 144-generator tables.
- `independent-stdout.json`, `independent-optimized-stdout.json`: concise independent receipts.
- `original-normal.json`, `original-optimized.json`: fresh executions of the packet verifier.

The packet verifier passed 69,053 active guards in both normal and optimized execution; its receipts are byte-identical to each other and to the supplied verification receipt. The independent verifier passed 3,488 active checks in normal and optimized execution, with byte-identical concise receipts.
