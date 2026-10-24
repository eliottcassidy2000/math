# Independent audit: Session52 Fortunate / square-sum runs

**Status: PASS after one verifier hardening.** Mathematical claims are correct at the stated scope. The precise 17-bit identity is finite, the next comparison fails, and the ordinary/componentwise distinction is indispensable. No publishing was performed by this auditor.

## Frozen target and replay

Audited candidate manifest SHA256:

`610765931889818228afd6d25edd8ef1b4bf6f6befd35a0bdde26fa92398fae3`

Candidate results SHA256:

`7f3ab3cde28031a051d647814c44d0db2979a58cdd54ace0635af98d4cf14c32`

The exact delivered files were copied into `frozen_candidate/`. All twelve listed file hashes and sizes match its manifest. This is an actual relocated copy, not a claim about relocation based on the original directory. From `/tmp`, both `python3 .../frozen_candidate/check.py` and `python3 -O .../frozen_candidate/check.py` pass with identical stdout. Both versions of the supplied rejection suite pass, rejecting ten controls. Their saved output files are retained here.

A release review found that the initial outer audit manifest omitted the frozen candidate dependencies, even though the earlier local copied replay was genuine. This packaging defect is repaired: the final outer manifest includes the complete frozen candidate manifest and every one of its twelve payloads, as well as the independent checker, outputs, source-pin audit, and replay driver. No mathematical result or frozen candidate file changed.

`python3 verify_replay.py` verifies both manifest layers, copies exactly the outer root manifest and its listed payloads into a fresh temporary package, then runs the independent and candidate checkers plus the candidate hostile suite in normal and optimized modes from an unrelated empty working directory. The driver uses `-B` to suppress unlisted cache artifacts and verifies all payload hashes both before and after replay. `closure_replay.out` records the final manifest-only copy test. The thirteen independent semantic mutations are embedded in `independent_check.py`; no separate mutation script is required or claimed.

The independent script imports no producer code. Run:

    python3 independent_check.py
    python3 -O independent_check.py

Both pass and produce identical `independent_results.json` and stdout. Results SHA256 is `3eb06e3bb21c3ebd8926f64edbc175265b8ce374991bfd07081d0233ee25e341`.

## Independent arithmetic verification

All 200 inherited Lucas certificates were verified in increasing integer order, checking recursively certified complete factorizations of N−1, base domains, Fermat congruences, and the order/gcd conditions. This independently certifies the arithmetic, without accepting probable-prime passes or importing the producer's expected Fortunate values. All 860 positive-offset and 874 negative-offset endpoint decisions were checked, including every smaller offset before each first success. The negative side is an additional audit of the shared inherited certificate, not part of the claimed 17-bit comparison.

The resulting equality indices through19 are 1,2,3,6,7,8,14,16,17. Thus the positive equality-bit run lengths through17 are 3,2,3,5,1,1,2. The next-index test uses F_n=p_(n+1), not literal F_n=p_n+1; the latter holds only initially. Raw arithmetic gaps, gaps between equality indices, and sorted distinct values are different data.

The first extension is indeed n=18: F_18=109>67=p_19, while

    P_18+67=117288381359406970983337
           =4129×28406001782370300553.

The displayed factorization was checked exactly. The Fortunate lower-bound saturation equivalence follows directly from the primorial wheel and needs no Fortune conjecture.

### Verifier correction found and repaired before freeze

The initial new `mr_failure` helper lacked a base-domain guard, so `mr_failure(5,5)` incorrectly reported a compositeness witness. A failed strong Miller–Rabin calculation is sound only for an admissible base. The producer added `1<a<n` and a tenth hostile control before the frozen snapshot. All actual inherited bases are 2, so this fix changes no number or certificate conclusion. Importantly, the Session51 checker already checked the domain at its caller; this was a new helper/caller hardening issue, not a demonstrated error in Session51's certified results.

## Independent graph algorithm and controls

The audit uses a different mathematical state space from the producer's vertex-path DFS and visited-vertex bitmask counting. It enumerates spanning **linear forests as selected edge sets** by vertex-elimination dynamic programming:

1. Generate square-sum edges as unions of reflection matchings x↔s²−x.
2. At a vertex's elimination, choose zero, one, or two edges to later vertices, respecting the remaining degree capacity.
3. Keep connectivity labels through already eliminated vertices. Reject an edge joining vertices already in the same selected component, preventing cycles.
4. Maximize the number of selected edges and count all maximizers.

Every edge is decided exactly once. Every surviving set is acyclic with maximum degree two, hence a disjoint union of paths including singletons. Conversely every such forest follows exactly one branch. Therefore the maximum selected-edge count m gives pc=N−m, and the number of maximizers counts minimum unordered path covers (each undirected path counted without reversal). If pc=1, this is precisely the Hamiltonian-path count.

The implementation is cross-checked against literal edge-mask enumeration on **all 1,100 simple labelled graphs through five vertices**, including the empty graph. The general endpoint-aware birth iff and deficit identity are checked on **all 1,099 old-graph/new-neighborhood instances through four old vertices**. These finite controls supplement, not replace, the elementary general proof, which is valid including singleton and empty-graph boundary cases.

Through28 the independent optimum, positive witnesses, component partition, birth neighbors, touched component counts, cycle ranks, leaves, branching vertices, component path counts, and global path counts agree with the candidate. Every supplied Hamiltonian path is checked as a square-adjacent permutation, and reversal classes are distinct. Matching the independent total establishes completeness of the supplied path lists.

### Census and independently generated witnesses

The deficit-positive orders through28 are exactly

    13,14,18,19,20,21,22,24,

and every such deficit is1. All other orders have deficit0. Global Hamiltonian-path counts at positive orders are

    N: 1,15,16,17,23,25,26,27,28
    H: 1, 1, 1, 1, 3,10,12,35,52.

The independent minimum-cover counts at deficient orders are

    N: 13,14,18,19,20,21,22,24
    # optimal covers: 3,3,3,12,32,50,24,192.

`independent_results.json` contains an independently constructed minimum cover at every order and the **complete 192 minimum two-path covers at24**. This is a bounded supplementary census, not a new general theorem.

The negative mechanisms were also reviewed independently: at13,14,18 an offending component has at least three leaves; the displayed cuts refute19–22; and24 has the endpoint-aware forced-cycle argument already repaired in the repository. The new propagation respects marked endpoints and does not repeat the retracted unconditional degree-two rule. The triangle hostile gives three Hamiltonian paths, so the independent algorithm does not silently force both edges of a degree-two endpoint.

## Interpretation, window boundaries, and general proof

The exact graph predicate is pc(G_N)=c(G_N), equivalently a Hamiltonian path in **each component**. At10,11,12 there are three path components and no global Hamiltonian path. Replacing this with ordinary Hamiltonicity changes the initial five bits to zero and yields run lengths5,3,5,1,1,2, not the requested word.

Both graph-window boundaries are clipped. Componentwise traceability is true on all1..12, so the first displayed run of3 at10..12 is part of a run of12. The final run25,26 continues through27 and28, whereas the prime equality run16,17 ends at18. The match fails at the very first out-of-window index27, which has35 independently counted Hamiltonian paths. Thus no shared seven-run global process or fixed-shift universal identity follows.

The general vertex-birth law is sound:

    pc(H)−1 ≤ pc(H+v) ≤ pc(H)+1.

The lower value is attained iff some minimum old cover exposes neighbors of v as endpoints of two **distinct** paths. Deletion of v proves necessity; gluing proves sufficiency. The component recurrence c'=c+1−r gives the stated deficit update and lower bound. The theorem retains endpoint and component ownership, unlike the earlier failed forcing rule.

The note appropriately calls the graph predicate a natural reconstruction rather than claiming to know the user's unspecified intended observable. It declares no probability model, no p-value, no universal prime/graph map, no Fortune-conjecture proof, and no new literature priority. These limitations are mathematically necessary and correctly stated.

## Semantic hostile tests

The independent suite rejects13 meaningful mutations: ordinary/componentwise substitution; an invented path at24; a false minimum cover at24; suppressed extension failure; duplicate vertex in the27-path; raw gaps replacing equality bits; a false candidate factor; a forged wheel equality; invalid Lucas base; incomplete Lucas factorization; missing minimality interval; shifted primorial; and an invalid Miller–Rabin base equal to its endpoint. Normal and optimized Python reject the same mutations. The supplied producer suite separately rejects its ten controls.

## Inheritance and primary sources

Every one of the12 repository source pins was verified against the actual immutable Git blobs at commit `f50a8c2f3a307869ec51650851dbc43eaedfc176`, not against the stale shared `math/` checkout. See `source_pin_audit.json`.

Closest arithmetic inheritance is Session51's `prime-offsets/THEOREM.md`: the wheel saturation/cutoff mechanism and the certificate pack are inherited, not newly discovered. Closest graph inheritance is Session40's endpoint-marked path-cover/birth polynomial, together with `THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md`. The canonical hostile and repaired near miss are the unnumbered MISTAKES.md section **2026-09-25 Square-sum degree-two forcing forgot the endpoint coordinate** and `decoder_prime_square_20260925.md` §4. The current candidate credits these correctly. Session32's missing old published path is treated as provenance only, not a live dependency.

Primary web sources independently checked on24 October2026:

- [OEIS A035346](https://oeis.org/A035346): equality-index definition and the initial indices.
- [OEIS A060881](https://oeis.org/A060881): candidate sequence P_n+p_(n+1), with offset0 explicitly distinguished from the n≥1 comparison.
- [OEIS A065315](https://oeis.org/A065315/internal): least prime divisors of that candidate sequence; the audit only needs proper divisors and checks them directly.
- [OEIS A005235](https://oeis.org/A005235): the primorial-indexed Fortunate definition and conjectural universal primality.
- [OEIS A071983](https://oeis.org/A071983): square-chain reversal convention and the counts through28.
- [Čejchan–Křížek–Somer, JIS25 (2022), §3](https://cs.uwaterloo.ca/journals/JIS/VOL25/Krizek/krizek3.pdf): primorial definition, finite-wheel mechanism, and distinction from Fortune's conjecture.

No unverified all-order square-chain theorem, sequence-database probable-prime calculation, or statistical significance claim is needed for any audited conclusion.
