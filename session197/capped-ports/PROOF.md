# Every bounded-count moving-port language is regular

Session 197. This is a separate all-order consequence of the full moving-port birth language. It constructs a finite deterministic machine preserving every capped Hamilton count under every continuation. It does not assert finite unrestricted linear rank, nor that its compressed states are physical descendants of the initial seed. No literature-priority claim is made.

## 1. Language, conventions and statement

Let H(T) count Hamilton vertex words of a nonempty tournament. Fix k distinct ordered exposed vertices, called ports. A right birth adds z with z beating precisely a chosen subset of the current ports and losing to all other old vertices. A left birth reverses that convention: z loses precisely to the chosen subset and beats all other old vertices. Then forget one of the k old ports and z, retaining exactly k ports in their remaining order. An old unexposed vertex can never be re-exposed. All subsets, sides and forgetting positions are allowed. The alphabet has 2^(k+1)(k+1) letters. The k=0 alphabet consists of a universal source or sink birth.

Fix an integer cap N>=1. Define the integers

    a = max{t>=0: 3^t <= N},
    b = max{t>=0: 5^t <= N^3},
    M = b+1,
    V = a+b+k+(a+k+1)(b+1).

The output of a graph is min(H(T),N+1).

**Theorem.** For every fixed initial tournament and ordered k-port list, there is an explicit finite deterministic output machine reading exactly this alphabet and giving min(H(T_w),N+1) after every word w. Every nonrejecting state is an actual tournament with at most V vertices and k ordered ports. Consequently each exact-count language {w:H(T_w)=h}, 1<=h<=N, and the bounded-count language {w:H(T_w)<=N} are regular.

One upper bound on the number of states, including one rejecting sink, is

    S(N,k) = 1 + sum_{n=max(1,k)}^V 2^(n(n-1)/2) * n!/(n-k)!.

The bound deliberately counts all labelled tournaments and port lists, not just reachable states. Any nonempty accepted exact-count language has a witness of length at most S(N,k)-1. The separate induced-pruning theorem gives a much smaller witness bound; the point here is exact finite-state behavior under all continuations, not an improved search-depth estimate.

## 2. A permanently unexposed transitive block

The strong components of a tournament have a unique total order, with all arcs from each earlier component to each later component. Consider a consecutive nonempty run B of singleton strong components, none of whose vertices is currently exposed. Its induced tournament is transitive. Every exterior vertex has the same incidence to every vertex of B, because its component lies wholly before or wholly after this consecutive run. Thus B is a module: exterior vertices see it uniformly.

Under every future allowed birth, every vertex of B receives the same default incidence. For a right birth all of B beats the newborn; for a left birth the newborn beats all of B. No vertex of B can ever be selected as an exceptional port because it is already unexposed and cannot be re-exposed. Changes to neighboring vertices' exposure do not affect this fact. Old arcs never change. Therefore B remains a transitive module, with identical incidences to all new exterior vertices, after every possible continuation.

Collapse B to one quotient vertex q. Every continuation yields the same exterior graph and the same quotient graph regardless of the positive number of ordered clones used for B. The ordered ports correspond exactly, and every prescribed letter selects and forgets corresponding ports.

There are precisely two possibilities in the quotient SCC decomposition:

- q is a singleton SCC. Then every vertex of B is itself a singleton SCC in the expanded graph. The block contributes the factor H(B)=1; all other SCCs and their counts are unchanged by its length.
- q lies in a nontrivial SCC. Then every vertex of B lies in the same nontrivial SCC as that component's exterior vertices. To see that all clones join, lift routes from q to an exterior vertex and back: their first and last arcs have the same incidence for every clone. Exterior-to-exterior reachability is likewise unchanged. In particular, a partial merge of only some clones is impossible.

This dichotomy is independent of the positive block length. It is also independent of future changes to exposure outside B.

## 3. Exact capped-continuation equivalence

Suppose the run has length l>M. Replace it by a transitive run of exactly M vertices with the same exterior incidences, preserving every current port. Call the resulting graph T'. For every finite continuation word w,

    min(H(T_w),N+1) = min(H(T'_w),N+1).             (1)

Indeed, in the singleton quotient case the entire Hamilton counts are equal, by SCC multiplication and the block's unit contribution. In the nontrivial quotient case, each expanded SCC contains at least M block vertices and at least one exterior vertex. Its order is at least M+1=b+2.

Busch's lower bound for a nontrivial strong tournament C is

    5^(|C|-1) <= H(C)^3.

Since 5^(b+1)>N^3 by the definition of b, both expanded strong components have H(C)>N. The total graph count, a product of positive SCC counts, is also greater than N in both graphs. Thus both capped outputs are N+1. The strict inequality is essential: the chosen M is b+1, not an equality-bound endpoint.

In fact a nontrivial quotient SCC contains at least two exterior vertices, allowing a slightly smaller cutoff, but that improvement is not needed or used here.

Equation (1) is an equivalence for every continuation, not just for the current graph, selected test words, or eventual outputs. Truncating several such runs sequentially preserves it by transitivity. At the current state every truncation preserves H exactly because the runs currently consist of singleton SCCs.

## 4. Why only finitely many vertices remain

Assume H(T)<=N. Let s be the number of nontrivial SCCs. Each has H at least three: a directed triangle has three Hamilton paths, and vertex insertion injects paths of an induced tournament into paths of the whole tournament. Thus 3^s<=H(T)<=N and s<=a.

For the orders n_i and counts h_i of those SCCs, multiply Busch's inequalities:

    5^(sum_i(n_i-1)) <= product_i h_i^3 = H(T)^3 <= N^3.

Their total order is therefore at most s+b<=a+b. This bound also covers s=0.

At most k singleton SCCs are exposed ports. Treat the s nontrivial SCCs and those exposed singleton vertices as separators in the total SCC order. There are at most s+k+1<=a+k+1 runs of unexposed singleton SCCs. Truncating each to length at most M leaves at most

    (a+b) + k + (a+k+1)M = V

vertices. Ports inside nontrivial SCCs have already been counted in a+b; counting all k again is only a harmless upper bound.

For these intermediate bounded-H states, a=floor(log_3 N) is required. Omega(N) cannot replace it merely because N is the final cap: an intermediate Hamilton count need not divide N. The number of runs is O(k+log N).

For N=1, a=b=0 and M=1. Every bounded state is transitive, and the bound becomes V=2k+1. The same argument and machine apply. For k=0 the port list is empty throughout; no special exception is needed.

## 5. Explicit deterministic machine and correctness

Use adjacency row masks on labels 0,...,n-1 together with an ordered list of k labels. No isomorphism testing is required. At any bounded graph, compute its ordered SCCs. In every unexposed singleton run longer than M, retain the first M singleton vertices in SCC order. Delete the others, then relabel the retained vertices in increasing old-label order and transport the ordered ports. This is a deterministic compression operation C.

Initialize as follows. If the seed has H>N, use the rejecting sink. Otherwise use its compressed representation C(seed). From a nonrejecting state, perform the indicated physical birth and forgetting operation. If the resulting graph has H>N, enter the rejecting sink. Otherwise compress it and use the resulting labelled graph and ordered ports as the next state. The sink loops on every letter and has output N+1.

Hamilton counts are monotone under induced extension: insert a new vertex before the first path vertex it beats, or append it if it beats none. The old path is recovered by removing the new vertex, proving injectivity. Therefore an actual graph with H>N can never return below the cap.

Maintain the invariant that the actual physical graph and the stored representation have equal capped outputs for every continuation word. It holds initially by Section 3. If it holds before a letter, applying that letter preserves the equivalence for every suffix. A sink transition is sound by monotonicity. Otherwise compression preserves the equivalence again by (1). Induction proves the theorem for every word. Each bounded stored graph has at most V vertices by Section 4, so the machine has the stated finite state bound. Outputs are its stored exact H, or N+1 at the sink.

Compression is a behavioral representation, not an authorized physical graph operation or an assertion that the shortened state retains the whole initial seed. An accepted word still denotes its original legal physical construction from that seed. This distinction separates the theorem from seed-preserving induced pruning.

## 6. Sources, scope and verification

The moving-port language and separate witness-shortening theorem are specified in `../../session197-hybrid-constructions/induced-pruning/PROOF.md`, frozen SHA256 `a5ec2d6842289f75e833881e4bb1c5e0ab715c38446a3486b5119b06c64bdee0`. That packet proves finite target reachability by shortening a successful word. The present argument instead identifies states with exactly equal capped behavior under every continuation and obtains regular languages. It does not replace that packet's smaller witness bound.

The only quantitative external theorem is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, Theorem 1, DOI10.37236/1141:
https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf

The full fixed-k grammar is essential. More importantly, permanent non-exposure is essential to the module argument. Arbitrary later re-exposure of old vertices is outside the theorem. Finite capped behavior does not imply finite unrestricted linear rank or unrestricted coverage of all Hamilton counts.

`verify.py` supplies independent physical controls on SCC ordering, run compression, the exact threshold, bounded vertex counts and continuation outputs. Its finite checks do not replace the all-order module and induction proofs above.
