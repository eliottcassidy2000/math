# Independent order-8/9 census and physical weighted-cycle audit

This addendum supersedes REPORT.md's exclusions of the order-8/9 census and its statement that a nonzero compatible physical loop remained unproved. The non-ordinal extension census remains outside this audit.

## Completed independent completeness certificate

`orbit_mass.cpp` performs exhaustive degree-compatible canonicalization independently of the producer's individualization/refinement implementation. It verifies pairwise nonisomorphism, checks EVERY saved automorphism count, and computes exact labeled orbit sums:

- Order 8: 6,880 classes, 403,193 examined degree-compatible labelings, orbit mass 268,435,456 = 2^28.
- Order 9: 191,536 classes, 27,282,165 examined degree-compatible labelings, orbit mass 68,719,476,736 = 2^36.

The nonisomorphism certificate and exact orbit mass together establish completeness of the supplied representative banks; this is stronger than trusting the mass alone. No independent regeneration of the class list was needed.

The producer's classgen.c is byte-identical to inherited Session230 `primes/classgen.c`, SHA256 cb78b6418a7754d656f628eede3cb3628a5100ad93af14d947b25e48223f562e. Its inductive generation is also sound: every n-vertex tournament has an (n−1)-vertex deletion isomorphic to one previous representative, and extension over every incoming/outgoing neighborhood supplies an isomorphic copy. Canonical deduplication never removes all copies of a class. The independent orbit certificate avoids relying on that implementation argument alone.

## Independent complete census replay

`replay_census.cpp` is a separate implementation. It builds forward Hamilton-path endpoint tables (and forward tables in the reversed graph for starts), computes all cut coordinates by incremental subset recurrences, writes a record for every graph/cut/tail-bit combination, and SORTS records rather than using the producer's hash map. Runs are checked separately in each parity, so the replay does not silently assume homogeneous parity.

Results:

- Order 8: 6,880 classes; 3,522,560 records; 3,385,144 nonempty erased interfaces; maximal step-two run 3; maximal unsaturated run 2.
- Order 9: 191,536 classes; 196,132,864 records; 193,541,720 nonempty erased interfaces; maximal step-two run 3; maximal unsaturated run 2.
- Both universes have zero mixed-parity fibers. This independently validates the producer's sorted-adjacent step-two run test.
- Every one of the 64 order-9 bucket class counts, interface counts, overall maxima and unsaturated maxima agrees exactly with ERASE9.stats.

The all-topology grouping is correct. With tail bit b and selected leading singleton, the padded state is

    h' = h, q' = 4q,
    A' = 2(q+A)+hb,
    G' = h+2G+2qb,
    S' = h, E' = hb,
    T' = G+b(q+A)+(T+S+bE).

Thus the key (h,q,A,G,b) is in bijection with the full padded interface (h',q',A',G',S',E'); the local varying quantity is exactly T+S+bE. Every mask is included, and equal keys from distinct topologies are combined. Since h is positive and odd, b is recoverable from E', and gcd(h,4q)=gcd(h,q). Grouping by (savedH/2) modulo the bucket count cannot split an interface because H is part of the key, and every savedH is explicitly checked against an independent count.

Arithmetic bounds are safe: a nine-vertex graph code needs 36 bits, handled in unsigned 64-bit storage. Hamilton paths are at most 9!, and the number of ordered covering path pairs is at most 10·9! = 3,628,800. All seven-coordinate computations and hash keys use 64-bit arithmetic in the producer. The replay separately checks 32-bit record conversion rather than assuming it.

The earlier +1 interface artifact in the order-7/8 map census was the final insertion of an absent target key. Correct nonempty counts are 104,432 and 3,385,144. The bucketed order-9 implementation has no corresponding artifact.

These finite censuses bound original full-interface fibers at each enumerated order as well: fixed leading/trailing singleton padding preserves their T differences and maps them into a single erased fiber. They do not bound arbitrary orders or arbitrary padding lengths uniformly.

## Physical nonzero weighted cycle

`weighted.py` imports no producer counting or composition implementation. It independently verifies the new physical core

    adj = [34,20,9,35,13,86,31], mask = 17,
    (h,q,A,G,T,S,E) = (95,1276,270,372,40,20,30).

The core state is computed by literal permutation-and-split enumeration. Every bit-labelled transition in PHYSICAL_CYCLE.json is checked with exact integer arithmetic, including integrality and the strict finite-state bounds. The 312-step walk returns to the augmented state (270,372,0), has total cocycle gain +190, and is independently verified by a GLOBAL composition identity rather than merely summing the same local weight formula.

Forward search followed by reverse reachability from this physical state reconstructs its exact 2,232-node SCC, equal as a set to the producer's NONZERO_SCC.json. The explicit walk suffices to prove nonzero drift; no minimal-length or globally shortest nonzero-cycle claim is made.

### Boundary cancellation and proof for every repetition count

Let X,Y be the singleton words on the two columns of the closed walk. Write their states as (1,2^d,A_X,G_X,T_X,S_X,E_X) and (1,2^d,A_Y,G_Y,T_Y,S_Y,E_Y), with d=312. The closed A/G equations imply that C Y and X C have equal h,q,A,G. Put

    u = E_X = last(X), v = S_Y = first(Y).

Then the T difference between u C Y v and u X C v is

    Δ = A G_Y − A_X G + h(T_Y−T_X) − h u S_X + h v E_Y.

To see the cancellation, the unpadded difference is

    ΔT = A G_Y − A_X G + h(T_Y−T_X) + E S_Y − E_X S.

The singleton boundaries contribute u(S−hS_X)+v(hE_Y−E), canceling the physical core's S,E terms exactly. This formula is independently evaluated from the two complete word states and gives Δ=190. In the actual receipt all four endpoint bits first(X),last(X),first(Y),last(Y) are zero.

For every integer k >= 1, consider the k+1 deletion objects

    s1 ⇒ s0 ⇒ X^j ⇒ C ⇒ Y^(k−j) ⇒ s0,   0 <= j <= k.

Their h,q,A,G agree by successive replacement of C Y with X C. Their outside singleton boundaries force S=h and E=0. At each adjacent replacement, the common left context ends in u=0 and the common right context starts in v=0. The exact boundary calculation therefore gives the SAME T decrement 190 for every j, independent of j and k. Common singleton contexts have Hamilton count 1, so they do not scale this difference.

Hence this is an all-k theorem, not an extrapolation from six examples: an unsaturated full-interface arithmetic progression of k+1 distinct T values, with step 190, on 312k+10 deletion vertices. Strong star closure gives k+1 counts with the same step on 312k+12 vertices. The common h is 95 and common q is 1276·2^(312k+3), which are coprime. The leading selected singleton proves strongness.

Independent exact tests for k=1,2,3,6,10 agree, but the proof above establishes all k. This proves unbounded unsaturated fiber CARDINALITY. It does not provide a consecutive odd interval: the actual count spacing is 190, not 2.

## Other physical entries and smaller moving witnesses

An independent forward-path search through ALL H=95 classes in the complete banks through order9 reproduces PHYSICAL_BRIDGE.json. The five H95 classes at order7 all have q1276. Exactly four masks on two of those classes intersect the nonzero SCC:

    (code,mask,A,G)
    (1116497,17,270,372), (1116497,110,911,809),
    (1452353,24,372,270), (1452353,103,809,911).

The 20 H95 classes at order8 and 37 at order9 have no q1276 examples. Thus there are no further entries through order9 at this exact h,q. This is not a classification of larger orders.

The newly written Section 7 small single-core constructions were also checked on literal deletion graphs and direct forward subset Hamilton counting of the closures:

- C5 moving pair: common (h,q,A,G,S,E)=(5,1024,716,209,5,0), T=142,140; strong star counts1295,1297.
- C9 moving triple: common (h,q,A,G,S,E)=(9,1344,1090,427,9,0), T=334,338,336; strong star counts1855,1851,1853.

The latter is a DIFFERENT triple from the order7 cross-topology census example in REPORT.md. Both are valid; their interfaces should not be conflated.

## Reproduction

Compile orbit_mass.cpp and replay_census.cpp with g++ -O3 -std=c++17, then:

    orbit_mass 8 ../session231-dense-fibers/classes9_8.txt
    orbit_mass 9 ../session231-dense-fibers/classes9_9.txt
    replay_census 8 ../session231-dense-fibers/classes9_8.txt 16
    replay_census 9 ../session231-dense-fibers/classes9_9.txt 64
    python3 -B weighted.py

Run those commands from this audit directory. Exact receipts are MASS8.txt, MASS9.txt, REPLAY8.txt, REPLAY9.txt, WEIGHTED.json, and WEIGHTED.log. The weighted checker currently uses an explicit workspace path to the producer's receipt, as documented in its source. No producer file was edited and no commit was made by this audit worker.
