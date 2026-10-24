# The 3,2,3,5,1,1,2 word: Fortunate saturation and square-sum path-cover deficit

Session52, 24 October 2026.

**PROVED (elementary, no priority claim):** the Fortunate lower-bound saturation test; vertex-birth laws for minimum path covers and their excess above component count. **FINITE-EXACT:** the 17-bit comparison, its first failure, and the complete graph census through28. **IDENTIFIED EXISTING SEQUENCE:** the prime-side equality indices are OEIS A035346. **NOT CLAIMED:** a universal arithmetic/graph correspondence, a statistical significance estimate, a proof of Fortune's conjecture, or an identification of the user's intended graph predicate.

## 1. Answer, with the two necessary distinctions

Let p_n be the nth prime (p_1=2), P_n=product_{i=1}^n p_i, and

    F_n = min{d>1 : P_n+d is prime}.

The user's word is exactly the run-length sequence of equality versus strict inequality between **F_n and p_(n+1)** for n=1,...,17. It is not a list of numerical differences, and p_(n+1), the next prime in index, must not be confused with p_n+1, one more than a prime's value.

Let G_N be the simple undirected graph on {1,...,N}, with distinct x,y adjacent when x+y is a square. Let pc(G) be the minimum number of vertex-disjoint paths whose vertex sets partition V(G), allowing singleton paths, and let c(G) be its number of connected components. Define

    delta(G) = pc(G)-c(G) >= 0.

Then delta(G)=0 iff **every connected component has a Hamiltonian path**. This is componentwise traceability, not ordinary Hamiltonian-path existence when G is disconnected.

The exact comparison is

    [F_n = p_(n+1)] = [delta(G_(n+9)) = 0],  1 <= n <= 17.

Both Boolean words are

    111 00 111 00000 1 0 11,

so their run lengths are 3,2,3,5,1,1,2. This is a precise natural reconstruction of the reported graph pattern. It is not evidence that this is the user's original definition; an unspecified graph property cannot be uniquely recovered from one finite word.

Under **ordinary** Hamiltonian-path existence, G_10,G_11,G_12 fail because each has three components. The graph word is instead

    00000 111 00000 1 0 11,

with run lengths 5,3,5,1,1,2. Thus ordinary traceability matches the proposed arithmetic bits only on n=4,...,17 (N=13,...,26), not on the first three positions.

### The extension fails immediately

At the next paired index,

    n=18, N=27,
    F_18=109 > p_19=67,
    P_18+67 = 117288381359406970983337
            = 4129 * 28406001782370300553.

Nevertheless G_27 has 35 Hamiltonian paths up to reversal. One certificate is

    18,7,2,14,11,5,20,16,9,27,22,3,13,23,26,10,6,19,17,8,1,15,21,4,12,24,25.

Every adjacent sum is a square and every vertex1..27 occurs once. Therefore the fixed-shift equality already fails at its first out-of-window test. No all-large-N theorem is required for this conclusion.

There is also a window-selection warning: componentwise traceability is true for **every G_1,...,G_12**. The first length3 is the last three terms of a length12 initial run, selected by starting at10. The last length2 at25,26 is truncated too: the run continues at27 and28. In contrast, the prime equality run16,17 really ends at18. The seven displayed runs are not seven shared complete global runs.

## 2. Inheritance and the concept board

The closest proved arithmetic mechanism is the primorial wheel and cutoff in Session51's `session51_structural_errors_relational_packets_20261024/prime-offsets/THEOREM.md`, Sections2–3. Its canonical hostile is a wheel survivor with a composite endpoint: here the first one is 210+11=13*17. The corrected near miss is to mistake sieve survival for primality or to sort the Fortunate values and silently preserve their original indices.

The graph side inherits `THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md`, the square-sum wave-six and decoder notes, and Session40's `session40_marked_ports_square_sum_20261023/square-sum/README.md`. The canonical hostile is a graph-degree2 vertex that is a path endpoint: degree2 alone does not force both edges. The correction is explicitly recorded under the unnumbered `01-canon/MISTAKES.md` heading “2026-09-25 Square-sum degree-two forcing forgot the endpoint coordinate”; no unrelated numerical MISTAKE ID is assigned to it. The least-used useful sidecar for this question is the **component label attached to every path**: it distinguishes path-cover number from its excess above the component lower bound. Session32's local graph certificate was inspected for provenance, but the shipped verifier does not import it or depend on its missing historical repository link.

Board: (i) primorial index versus sorted value; (ii) numeric gap versus equality bit versus run lengths; (iii) component lower bound versus path-cover excess; (iv) endpoint-marked birth operations; (v) selected finite windows versus continued processes. Each is retained explicitly below. The research moves used were “search the statement before the method” (which finds A035346), “correct the object before sharpening the technique” (ordinary versus componentwise traceability), and endpoint-aware hostile testing. This is additional evidence for those existing moves, not a new general research-method claim.

## 3. Why the prime-side bit has an exact primality interpretation

### Proposition 1. First wheel-survivor saturation

For every n>=1,

    F_n >= p_(n+1),
    F_n = p_(n+1)  iff  P_n+p_(n+1) is prime.

**Proof.** If 2<=d<p_(n+1), some prime factor r of d is at most p_n and divides P_n. Thus r divides P_n+d, and r<P_n+d, so that endpoint is composite. Hence every offset before p_(n+1) is rejected. The equality is now precisely primality of the first non-rejected candidate. This does not assume that F_n itself is prime for all n. QED.

The candidate sequence P_n+p_(n+1) is [OEIS A060881](https://oeis.org/A060881); the equality indices are already [OEIS A035346](https://oeis.org/A035346). [A065315](https://oeis.org/A065315) records the candidates' least prime factors. The standard Fortunate definition is [A005235](https://oeis.org/A005235). The wheel mechanism and its distinction from the general conjecture are also treated in [Čejchan–Křížek–Somer, JIS25 (2022), Section3](https://cs.uwaterloo.ca/journals/JIS/VOL25/Krizek/krizek3.pdf). Identification and elementary reproof here are not novelty claims.

In particular, define a_n=F_n-p_(n+1). The word encodes zero/nonzero runs of a_n, not the values of a_n. It discards how far the failed candidate lies below the first prime endpoint. The first failures are

    P_4+11 = 221 = 13*17,
    P_5+13 = 2323 = 23*101.

The exact data follow. Every F_n is independently checked from its inherited certificates, including compositeness of every earlier offset; every failed candidate has a checked proper factor. All prime certificates use a full recursively certified factorization of N-1 and the Lucas order criterion, not a probable-prime pass.

| n | p_n | p_(n+1) | F_n | F_n-p_(n+1) | equality |
|---:|---:|---:|---:|---:|:---:|
| 1 | 2 | 3 | 3 | 0 | yes |
| 2 | 3 | 5 | 5 | 0 | yes |
| 3 | 5 | 7 | 7 | 0 | yes |
| 4 | 7 | 11 | 13 | 2 | no |
| 5 | 11 | 13 | 23 | 10 | no |
| 6 | 13 | 17 | 17 | 0 | yes |
| 7 | 17 | 19 | 19 | 0 | yes |
| 8 | 19 | 23 | 23 | 0 | yes |
| 9 | 23 | 29 | 37 | 8 | no |
| 10 | 29 | 31 | 61 | 30 | no |
| 11 | 31 | 37 | 67 | 30 | no |
| 12 | 37 | 41 | 61 | 20 | no |
| 13 | 41 | 43 | 71 | 28 | no |
| 14 | 43 | 47 | 47 | 0 | yes |
| 15 | 47 | 53 | 107 | 54 | no |
| 16 | 53 | 59 | 59 | 0 | yes |
| 17 | 59 | 61 | 61 | 0 | yes |
| 18 | 61 | 67 | 109 | 42 | no |
| 19 | 67 | 71 | 89 | 18 | no |

### Failed arithmetic interpretations retained

1. **Literal F_n-(p_n+1):** the first17 values are 0,1,1,5,11,3,1,3,13,31,35,23,29,3,59,5,1. Equality occurs only at n=1, giving run lengths1,16. This fails universally after the first term: for n>=2, p_n is odd, p_n+1 is even, whereas F_n is odd; the wheel bound gives strict positivity of the difference.
2. **Raw F_n-p_(n+1):** 0,0,0,2,10,0,0,0,8,30,30,20,28,0,54,0,0. All are even, so they cannot themselves equal the mixed-parity proposed word.
3. **Gaps between equality indices:** the indices1,2,3,6,7,8,14,16,17 have successive gaps1,1,3,1,1,6,2,1. Those are not the run lengths.
4. **Sorted distinct Fortunate values:** through71 the complete set in increasing order is3,5,7,13,17,19,23,37,47,59,61,67,71. Subtracting the correspondingly indexed next primes gives0,0,0,2,4,2,4,14,18,28,24,26,28. Sorting has destroyed the repeated/indexed returns responsible for the word.
5. **Prime ranks rather than prime values:** in this finite table each F_n is prime, so testing whether its prime rank is n+1 gives the same equality bit. It is an equivalent zero test here, not a second independent occurrence or a proof that all Fortunate values are prime.

The supplied checker records the full surviving and failed arrays; it does not search arbitrary transformations until a matching word appears.

## 4. Graph census and certificates

Paths are undirected and counted up to reversal. The trivial one-vertex graph has one Hamiltonian path; a Hamiltonian cycle is a simple cycle of length at least3. The definition is the one used by [OEIS A071983](https://oeis.org/A071983). Through28, every connected nontrivial G_N has a leaf, so none has a Hamiltonian cycle. No cycle/path ambiguity can produce the proposed word.

For a new vertex N the neighbors are s^2-N with N<s^2<=2N-1, hence

    d_new(N) = floor(sqrt(2N-1))-floor(sqrt(N)).

If it touches r old connected components, then c_new=c_old+1-r. With beta=|E|-|V|+c, beta_new-beta_old=d_new-r. The following table is exact; pc is certified by an explicit covering family and a matching lower bound.

| N | new degree | edges | components | cycle rank | Hamilton paths | pc | delta |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 1 | 7 | 3 | 0 | 0 | 3 | 0 |
| 11 | 1 | 8 | 3 | 0 | 0 | 3 | 0 |
| 12 | 1 | 9 | 3 | 0 | 0 | 3 | 0 |
| 13 | 2 | 11 | 2 | 0 | 0 | 3 | 1 |
| 14 | 2 | 13 | 1 | 0 | 0 | 2 | 1 |
| 15 | 2 | 15 | 1 | 1 | 1 | 1 | 0 |
| 16 | 1 | 16 | 1 | 1 | 1 | 1 | 0 |
| 17 | 1 | 17 | 1 | 1 | 1 | 1 | 0 |
| 18 | 1 | 18 | 1 | 1 | 0 | 2 | 1 |
| 19 | 2 | 20 | 1 | 2 | 0 | 2 | 1 |
| 20 | 2 | 22 | 1 | 3 | 0 | 2 | 1 |
| 21 | 2 | 24 | 1 | 4 | 0 | 2 | 1 |
| 22 | 2 | 26 | 1 | 5 | 0 | 2 | 1 |
| 23 | 2 | 28 | 1 | 6 | 3 | 1 | 0 |
| 24 | 2 | 30 | 1 | 7 | 0 | 2 | 1 |
| 25 | 2 | 32 | 1 | 8 | 10 | 1 | 0 |
| 26 | 2 | 34 | 1 | 9 | 12 | 1 | 0 |
| 27 | 2 | 36 | 1 | 10 | 35 | 1 | 0 |
| 28 | 2 | 38 | 1 | 11 | 52 | 1 | 0 |

At N=12 the components are the paths

    2–7–9; 12–4–5–11; 10–6–3–1–8.

Earlier G_10 and G_11 restrict these paths. Thus their ordinary Hamilton count is zero but pc=c=3. Vertex13 joins two components via3 and12, creating a three-leaf tree component: its leaves are8,10,11. Vertex14 joins the last components, yielding a connected tree with leaves8,9,10. Three leaves cannot fit in a single path, proving positive deficit at13 and14. At18 the three leaves16,17,18 give the same obstruction.

For N=18,19, removing S={7} leaves3 components; for N=20,21, removing {5,7} leaves4; for N=22, removing {5,7,14} leaves5. Each violates the necessary condition c(G-S)<=|S|+1 for a spanning path. Exact remaining components are shipped in results.json.

For N=24, the leaf18 must be one endpoint. For each candidate other endpoint e, set target degree1 at18,e and2 elsewhere. Repeatedly apply only these sound rules:

- too many selected edges or too few available edges at a vertex is a contradiction;
- if its target is reached, delete all other incident edges;
- if exactly its target number remain, force all of them;
- a forced cycle is a contradiction to being a spanning path.

Nineteen of the23 possibilities produce a degree contradiction. For e=2,9,22, the forced proper cycle is

    1–8–17–19–6–10–15–21–4–12–24–1.

For e=11, the engine already forces the cycle

    2–14–22–3–13–23–2.

Full per-endpoint forcing traces are included. Endpoints of graph degree2 are marked before forcing, preserving the earlier correction. This is a sound finite proof procedure, not a claim that degree propagation decides Hamiltonicity in arbitrary graphs.

All eight deficient orders have explicit minimum covers. Here are witnesses for the less immediate ones:

- N=13: (8,1,3,6,10); (11,5,4,12,13); (2,7,9).
- N=14: (8,1,3,6,10); (9,7,2,14,11,5,4,12,13).
- N=18: (16,9,7,18); (2,14,11,5,4,12,13,3,6,10,15,1,8,17).
- N=19: (16,9,7,18); (2,14,11,5,4,12,13,3,1,8,17,19,6,10,15).
- N=20: (18,7,2,14,11,5,4,12,13,3,1,8,17,19,6,10,15); (9,16,20).
- N=21: (18,7,2,14,11,5,20,16,9); (4,12,13,3,1,8,17,19,6,10,15,21).
- N=22: (18,7,2,14,11,5,20,16,9); (22,3,13,12,4,21,15,1,8,17,19,6,10).
- N=24: (18,7,9,16,20,5,11,14,2,23,13,12,4,21,15,10,6,19,17,8,1,24); (3,22).

Each is checked as a vertex partition with square adjacencies. Their sizes exceed c by exactly1; the no-path certificates in the one offending component prove optimality. At all remaining orders through28, a path per component attains the trivial lower bound c.

### Alternative graph interpretations, not hidden

On the fixed inclusive window10..26:

- ordinary Hamilton-path existence: run lengths5,3,5,1,1,2;
- componentwise traceability, or equivalently zero path-cover deficit:3,2,3,5,1,1,2;
- birth degrees:1^3,2^3,1^3,2^8, hence lengths3,3,3,8;
- numbers of touched old components:1^3,2^2,1^12;
- component counts:3^3,2^1,1^13;
- connectivity: false^4,true^13;
- forest status: true^5,false^12;
- maximum degree at most2: true^3,false^14;
- raw minimum path-cover number:3^4,2^1,1^3,2^5,1^1,2^1,1^2;
- Hamiltonian-cycle existence: false at every position.

Edge counts, cycle ranks, leaf counts and branching-vertex counts are also retained in results.json and do not give the proposed run word. Shifting the physical window to9..25 or11..27 changes the boundary run lengths and is recorded, rather than choosing a successful shift silently.

A genuinely different convention that includes vertex0 also fails: on vertices0..N-1, both global and componentwise traceability give lengths5,4,2,6; on vertices0..N they give4,4,2,7. Those checks keep their stated physical vertices and are not confused with merely renaming the index of G_N.

## 5. A genuine general relation: path-cover deficit under vertex birth

### Proposition 2. Endpoint-aware minimum-cover birth law

For any finite simple graph H, let H'=H+v adjoin one vertex with any specified neighborhood. Then

    pc(H)-1 <= pc(H') <= pc(H)+1.

Moreover, pc(H')=pc(H)-1 iff **some minimum path cover of H contains two distinct paths with respective endpoints a,b in N(v)**. Singletons are allowed as paths and their sole vertex is an endpoint.

**Proof.** Adjoining {v} to a minimum cover gives the upper bound. Removing v from any cover of H' increases its path count by at most1, which gives the lower bound. If the lower bound is attained, v must be interior in its selected path: deleting it produces exactly pc(H) paths and exposes two different path pieces with endpoints in N(v). Conversely two such endpoints in distinct paths can be joined through v, lowering the number of paths by1. Joining the endpoints of the same path creates a cycle, not an allowed path, so the distinctness sidecar is essential. QED.

If v has no neighbors, the increase is exactly1. If it has only one neighbor, a decrease is impossible. In particular, a leaf attaches without increasing pc precisely when some minimum old cover exposes its neighbor as an endpoint; graph degree alone does not answer that question.

### Corollary 3. Separate component mergers from path repairs

If v touches r components of H, then

    delta(H')-delta(H) = pc(H')-pc(H)-1+r,
    delta(H') >= delta(H)+r-2.

The exact update follows from c(H')=c(H)+1-r and the lower bound follows from Proposition2. A vertex touching at least3 old components necessarily creates extra path-cover deficit, since one selected path can use only two of the new incident edges. For a degree2 arrival, the distinction is not the number2 but the locations and endpoint roles of its two neighbors.

In the present graph sequence:

- At13, r=2 but pc does not drop: delta increases from0 to1.
- At14, r=2 and pc drops3→2, but c also drops2→1: delta stays1.
- At15, r=1 and pc drops2→1: delta drops1→0.
- At16 and17, the new leaf extends an exposed endpoint of the unique path.
- At18, the new leaf's neighbor7 is internal in that path: pc rises1→2 and delta rises0→1.
- At23 and25, r=1 and the new vertex joins endpoints of a minimum two-path cover: delta drops1→0.

The checker stores minimum-old-cover and new-cover witnesses for14,15,23,25. This is an exact graph mechanism extending the inherited endpoint-marked two-path-cover viewpoint; it does not involve primes.

## 6. Connection contract and what remains

The only proved arithmetic-to-graph map here is the **finite index assignment** n↦N=n+9 for1<=n<=17, followed by the respective lower-bound saturation tests. It preserves one Boolean value. It loses the magnitude of F_n-p_(n+1), the primorial and its factorization, graph components and endpoints, path multiplicity, and every continuation claim. The repair sidecars are the full arithmetic row and full marked graph cover. At n=5, for example, the arithmetic excess is10 while graph delta at14 is1, so there is not even equality of the numerical defects inside the matching window.

The decisive continuation is n=18↔N=27, which fails. The strongest survivor is therefore the certified finite-word identity, plus two separate exact mechanisms explaining each bit. There is no declared probability model or account of all possible searched transformations, so no p-value or “unlikely coincidence” claim is warranted. Finite identification does not establish a common cause, recurrence, bijection of witnesses, asymptotic law, supersingular-prime interpretation, MRB-constant relation, or prediction for A000568. Other lanes may connect those objects by actual maps; this word alone supplies none.

A useful next question, if the user had a different graph observable in mind, is to specify that observable before further sequence searches. Under the explicit componentwise interpretation the question is already settled through the first hostile extension, and no longer computation is needed to decide the claimed fixed-shift universal equality.

## 7. Reproduction and scope

Run from this directory:

    python3 check.py
    python3 -O check.py
    python3 test_rejections.py
    python3 -O test_rejections.py

The normal and optimized check outputs agree byte-for-byte. No Python assertion is used as a validity gate. The arithmetic checker revalidates860 positive-offset minimality decisions and109 recursively used Lucas prime certificates from the pinned Session51 certificate. Strong Miller–Rabin failures are used only as exact compositeness witnesses; passes are never accepted as primality proofs.

The graph checker builds each graph independently from integer square tests, enumerates Hamilton paths in a vertex-set DFS, and independently counts them using a bitmask dynamic program with all possible starts. Their exact counts agree for every1<=N<=28. It checks minimum path-cover witnesses, all required endpoint-propagation negatives, the displayed cuts, birth laws, alternative conventions and the first extension failure. The ten hostile controls corrupt indexing, primality/minimality data, path partitions, adjacency and interpretation; all are rejected in normal and optimized modes.

All dependencies are local files plus the Python standard library. `--write` regenerates results.json; normal verification compares it against a fresh deterministic recomputation. Source hashes and immutable repository paths are in SOURCES.json; all deliverable hashes are in MANIFEST.json. These artifacts were prepared for independent review and not published by this worker.
