# A sharp feedback-number separation at H=533

Session197. **All-order classification within the strong C<=2 class, by a published strong-order bound and an independently reviewed complete finite bank.** No claim about all H533 tournaments or the unrestricted scalar spectrum is made.

Let H(T) count Hamilton paths as directed vertex words, C(T) count Hamilton cycles modulo rotation, and f(T) be the minimum number of vertices whose deletion leaves a transitive tournament.

**Theorem.** Up to isomorphism, exactly six strong tournaments satisfy H=533 and C<=2. All six have order11, exactly two Hamilton cycles, and feedback vertex number exactly3. In particular H533 is attainable in the full strong C2 class but is absent from the strong C<=2 class with feedback number at most2, at every order.

## 1. The quantitative hypothesis is strongness

We use Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006),N3, Theorem1, DOI10.37236/1141. The exact named primary source is

https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf

Its independently checked consequence, already used in the Session197 H67 classification `PROOF.md`, is

    H(T)^3 >= 5^(n-1)  for a strong tournament of order n>=3.

Here

    533^3 = 151419437 < 244140625 = 5^12.

Thus n>=13 is impossible and n<=12. This is the sharper base5 consequence of Busch's strong-tournament theorem, not the weaker base3 estimate floor(3 log_3 H). No cycle-count hypothesis is used in this order reduction. The inequality must not be applied to the total order of a weak graph with arbitrarily many neutral singleton SCCs.

## 2. Complete finite classification

`C3_COMPLETE.md`, Section1, proves that every strong tournament with C<=3 and order at least4 is generated from the unique strong four-vertex seed by strong singleton extensions that retain C<=3. Its exact canonical census through12 is independently reviewed, with all27,945 graph classes matched by an unrelated endpoint-matrix/isomorphism enumeration. This completeness proof, rather than a stabilization assumption, supplies the finite bank here.

The verifier regenerates that frozen bank from `verify_c3_complete.cpp` and checks its full aggregate SHA256 against `c3_complete_certificate.json`. It independently recounts H for all12,077 C1/C2 representatives through order12. Capping each nonnegative endpoint-DP entry at534 preserves exactly the distinction H=533 versus any other value relevant to the comparison.

There are no C1 hits. There are exactly six C2 hits, all at order11, with outgoing rows:

    [2,4,9,19,39,335,159,319,607,1279,511]
    [2,4,9,19,1319,1487,1951,1823,1551,1087,15]
    [2,4,521,531,1575,1615,1695,1855,1663,1027,15]
    [2,260,1289,1939,1959,1999,1951,1799,1537,1031,3]
    [2,772,1945,2003,2019,1999,1927,1795,1537,1025,3]
    [2,1028,521,1843,2023,1991,1935,1807,1543,1027,5]

For each graph a separate literal directed-word DFS gives H=533 and C=2. The canonical bank already supplies pairwise nonisomorphism and completeness. Orders at most3 have H at most3 and cannot contribute.

## 3. Exact feedback number

For each of the six graphs and every pair of deleted vertices, the remaining graph contains a directed triangle. The certificate records one literal triangle for every such pair, for330 independently checkable lower-bound witnesses. Deleting fewer than two vertices cannot make it transitive either: any such transitive remainder could be further reduced by a second deletion.

Conversely, respectively delete

    {2,5,8}, {1,4,5}, {2,3,6}, {1,2,4}, {1,2,3}, {0,2,3}.

Each remainder is transitive. The certificate supplies its complete forward order and verifies every arc. Hence all six feedback numbers are exactly3, proving the theorem.

## 4. Why this obstructs fixed ports but not changing ports

Start from two exposed vertices and append vertices in birth order. Each newborn may point arbitrarily to the two exposed vertices, while every other old vertex points to the newborn. If the same two original ports are retained forever, deleting them leaves the other vertices in their transitive birth order. Such a graph has feedback number at most2.

Conversely, a tournament with a feedback set of size at most2 can be built with two fixed ports: expose that set, adding an extra initial vertex to it if necessary, then append the transitive remainder in order and choose its exact incidences to the ports. Thus for the orders considered here, this fixed-two-port graph class is precisely feedback number at most2.

The theorem therefore proves an all-order obstruction for fixed ports in the strong C<=2 class. It does not obstruct changing which two vertices are exposed. For the first displayed graph, an actual changing-port construction starts with the forward edge on original labels8,9. The remaining birth order is10,7,6,5,4,3,2,1,0, with (newborn mask, dropped-port index) moves

    (1,1),(1,1),(2,1),(3,0),(2,0),(2,0),(2,0),(2,0),(2,0).

Masks are read against the current ordered pair; append the newborn and drop index0,1 or2. The verifier reconstructs all55 graph arcs and ends with original labels1,0 exposed. This is a physical construction, not an equal-count substitution.

The distinction is narrower than a claim about all feedback-two tournaments: graphs with C>=3 are outside the theorem. It is also narrower than a claim about weak H533 graphs or about all odd counts. No least-new-omission assertion is required here.

## 5. Reproduction and dependencies

Run `python3 verify_feedback533.py` and `python3 -O verify_feedback533.py`. The script compiles the frozen complete census and the adjacent independent path-count filter `verify_feedback533.cpp`, then performs the literal path/cycle, feedback and physical-word checks. It writes `feedback533_certificate.json`; normal and optimized outputs agree.

Named inherited inputs are `C3_COMPLETE.md`, `verify_c3_complete.cpp`, `c3_complete_certificate.json`, and the strong-order theorem sourced in the main `PROOF.md`. Exact hashes are recorded in the receipt. No generic repository history, restricted source, feedback-two census assumption, or unproved family completeness is used.
