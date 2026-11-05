# Induced pruning, finite target reachability, and the blocked +4 child of the ternary clock

## 1. Physical language and exact scope

Start with any fixed nonempty tournament T0 and k distinct exposed vertices. A right birth adds a vertex beating any chosen subset of the current ports and losing to every other old vertex. A left birth adds a vertex losing to any chosen subset of the current ports and beating every other old vertex. After either birth, forget any one vertex from the k old ports and the newborn, keeping exactly k distinct ports. Forgetting exposure does not delete a graph vertex. An old vertex which has ceased to be a port cannot be re-exposed.

All subsets and all forgetting choices are permitted. Left and right births may be mixed arbitrarily and chosen adaptively. Port ordering only indexes these choices; no old edges are changed. This is the full k-port birth language considered here.

**Induced pruning theorem.** If a physical word constructs T from T0, and K is any vertex set containing every vertex of T0, then T[K] is constructible from the same initial graph and the same initial physical ports, using exactly the births of vertices in K outside T0. No other births are needed. The final exposed vertices need not coincide with ports that were deleted from T; every surviving original final port is still exposed. The final port-list order may differ, and must not be assumed to agree.

This is a theorem about induced graph construction. It does not claim that deleting arbitrary vertices preserves Hamilton counts.

## 2. Step-by-step simulation with inactive port representatives

Follow the original word, while constructing only its K vertices. Maintain a bijection f from the set of current original ports to the set of current simulated ports, with this invariant:

- If an original port p lies in K, then f(p)=p.
- Every simulated port is an already born vertex of K, and all k simulated ports are distinct.
- The simulated graph is exactly the induced graph on the already born K vertices.

If p is an original port outside K, its representative f(p) is an inactive kept vertex. In particular f(p) is NOT any current original port in K: such a port q would have f(q)=q=f(p), contradicting injectivity. It is therefore an ordinary nonport in the original run. This last fact prevents accidental exceptional arcs.

Initially every seed vertex is kept, so f is the identity on the initial ports. At each original birth z, there are two cases.

### A. Skip a vertex outside K

Do not perform a simulated birth.

- If the original move immediately forgets z, leave f unchanged.
- If the original move forgets old port p and exposes z, remove p from the domain of f and assign f(z) the previous representative f(p).

The simulated physical port set stays unchanged and still has k distinct members. If p was kept, it has just ceased to be an original port and its retained physical vertex can now act as an inactive representative. No old vertex has been re-exposed in the simulation: that physical vertex was continuously retained there.

### B. Perform a kept birth

Perform a birth of z on the same side (right or left). In its simulated mask select exactly the representatives of the original selected ports which lie in K. Select NONE of the inactive representatives of deleted original ports.

This gives every edge to an earlier kept vertex correctly. A kept original selected port receives its intended exceptional incidence. An inactive representative is an original nonport, so it must receive the default incidence; ignoring it in the simulated mask does exactly that. Every other kept vertex is also treated as in the original birth.

If the original move forgets z, forget the simulated newborn. Otherwise, if it forgets original old port p, forget the simulated vertex f(p), remove p from the domain of f, and set f(z)=z. Every other mapping is unchanged. The bijection and distinctness invariant are preserved.

Masks and forgetting positions must be translated through f and the CURRENT simulated port list. The lists need not stay aligned after skipped births. Since every subset and every forgetting position is allowed, this translation is always a legal move. The supplied implementation performs that translation explicitly.

Induction proves the theorem. Relabeling retained vertices in their original birth order makes this a usual consecutively labelled construction, without changing any edges or counts. The case k=0 also works: the bijection is empty and births are simply universal sources or sinks.

## 3. A finite, exact scalar reachability algorithm

Let the target N>1 be odd. Write Omega(N) for its number of prime factors with multiplicity, and let

`b(N) = max { t >= 0 : 5^t <= N^3 }`,

`B(N) = Omega(N) + b(N)`.

**Decision theorem.** If any word from T0 has Hamilton count N, there is such a word with at most B(N) births. If T0 is nontrivial and strong, of order n0, at most `B(N)-n0` births suffice.

To prove this, take any successful final tournament. Its strong components form an ordered chain, and its Hamilton count is the product of their counts. Every nontrivial strong component has odd count at least three. If there are s such components, then `s <= Omega(N)`.

For a nontrivial strong component of order n_i and count H_i, Busch's theorem gives `5^(n_i-1) <= H_i^3`. Multiplying over the nontrivial components yields

`5^(sum_i(n_i-1)) <= N^3`,

so their total order is at most `s+b(N) <= B(N)`.

Retain every vertex of T0 and every vertex of every nontrivial final strong component. Delete only NEW singleton strong components. This changes no Hamilton count: those components contribute factor one, and the other component graphs are untouched. The induced-pruning theorem constructs the retained graph with at most B(N) new vertices. If T0 is nontrivial and strong, all of its n0 vertices belong to one nontrivial final component, giving the subtraction by n0.

Thus enumerate every word up to the stated depth and count its actual graph. This is an exact all-order decision procedure, despite the possibly infinite dimension of a linear continuation representation. Reject branches once H>N, using induced-count monotonicity. Even targets are impossible by tournament parity. For N=1, reachability holds exactly when H(T0)=1, since no births are needed in that case and induced monotonicity excludes the other case.

Factoring N is optional for a looser bound: replace Omega(N) by `floor(log_3 N)`. For fixed k, the alphabet has `2^(k+1)(k+1)` operations and the depth is O(log N). Hence the word enumeration is pseudo-polynomial in the NUMERIC target N, not polynomial in its bit length. It can also be made polynomial in the input seed order: compute strong components first; a nontrivial component exceeding `1+b(N)` vertices already forces H>N, while every remaining component has O(log N) vertices and permits capped subset DP in N^O(1) time. The exponent depends on fixed k. No practical efficiency or unrestricted spectrum completeness is asserted.

The external bound used here is Theorem 1 of Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, DOI 10.37236/1141:
https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf
Its stated piecewise bound is stronger than the safe inequality used above.

## 4. An all-order missing child from the actual three-port clock state

Take the actual strong seed

`T0 = [2,12,145,149,99,143,175,19]`,

with exposed ordered ports `[0,1,4]`. It has H=179. The right birth selecting all three ports and forgetting the newborn gives H=539=3*179+2, and returns the same exposed ports. This is the first step of the independently audited `60*3^m-1` clock.

In contrast:

**No word of any length in the full mixed left/right three-port language from this same seed and these same initial ports has H=541=3*179+4.**

The target 541 is prime. Therefore any successful graph could have only one nontrivial strong component, necessarily containing the entire strong seed. Moreover

`5^11 <= 541^3 < 5^12`,

so its order would be at most twelve. The pruning theorem would give a successful word with at most four new vertices.

The complete capped census checks all 64 possible left/right mask/forget choices at every step through depth four, starting with exactly the specified seed and ports. The counts of accepted WORD PREFIXES at depths zero through four are

`1,36,440,5716,79768`.

These counts include different words yielding the same graph or interface; they are not claimed to be distinct states. No accepted prefix has count 541. The census is complete for the required finite depth and therefore, by the preceding theorem, proves the all-order exclusion.

This does not make 541 a missing tournament count, and it does not exclude a different graph with H179, different initial exposed vertices, arbitrary re-exposure, changing old edges, or a different compiler. It proves that the newly found returning +2 state cannot also supply its +4 child by this exact physical language, even with target-dependent words of unbounded length.

## 5. Verification and reproduction

Run `python3 verify.py` and `python3 -O verify.py` in this directory. The code is standard-library Python plus a C++17 compiler. Temporary census banks are generated under a temporary directory and removed; their deterministic SHA-256 is recorded in the receipt. No precomputed large bank is required.

The verifier:

- replays the shadow-port simulation on every retained subsequence of exhaustive short mixed words for k=0,1,2,3, checking the graph and bijection AFTER EVERY original step;
- tests longer mixed words on the actual clock seed, including arbitrary retained subsets;
- separately checks SCC pruning and the finite order bound on literal graphs;
- independently enumerates the first two census levels in Python;
- performs the full four-level census in C++, comparing forward and backward endpoint recurrences on every tested graph and using only the justified H>541 cap;
- verifies the strong seed, prime target, exact integer order bound, and positive H539 move.

All validation uses explicit exceptions or runtime errors, not removable assertions. The all-order claims come from the proved pruning and strong-order reductions, rather than extrapolating the finite census.
