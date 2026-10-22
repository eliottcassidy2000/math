# Full labelled attachment response is faithful: alternating-path Fourier isolation

Status: PROVED by the self-contained cancellation argument below, with independent exact controls in `faithfulness.py`. This extends the earlier Session27 response freeze. In particular it supersedes the earlier statement that full-response injectivity remained open; the universal TWO-base injectivity question remains open here. Root supplied the proof idea; this worker independently checked the combinatorics and implemented the direct-DP audit. No novelty or priority claim is made for the elementary coefficient-isolation mechanism.

## Theorem

Let Q and Q' be tournaments on the same fixed labelled k-vertex core. For each attachment word w, F_Q(w) counts Hamiltonian paths of the extension by an ordered transitive base, with bit u=1 meaning core u points to that base vertex.

For k≥2, the responses using exactly k−1 base vertices determine every labelled core arc. In particular the following are equivalent:

1. Q=Q' as labelled tournaments.
2. F_Q(w)=F_Q'(w) for every attachment word with exactly k−1 letters.
3. F_Q(w)=F_Q'(w) for every finite attachment word.

A fixed battery of at most C(k,2)·4^(k−1) response queries suffices. The query family is specified below and does not depend on the unknown core orientations. For k=0 or k=1 there is only one labelled core, so no measurement is needed. This is a reconstruction theorem, not merely equality decidability. It is not a proof that two base vertices suffice for all k, and it is not an optimal query/base-size bound.

## Experiment and exact formula

Fix distinct core vertices u,v, and list the other k−2 core vertices in any fixed order w_1,...,w_{k−2}. Take m=k−1 base vertices b_1,...,b_m, with their prescribed transitive orientation. Form the alternating vertex path

R=(v,b_1,w_1,b_2,...,w_{k−2},b_{k−1}).

This path visits every vertex of the extension except u. It has d=2k−3 edges, all cross edges between core and base. In particular d is odd. Orientations of these d cross edges are independently varied, and x_e is +1 when the orientation agrees with the displayed direction of R, −1 otherwise. Also independently vary the extra edge between u and b_{k−1}; call its bit t. This extra edge is distinct from every edge of R. Fix all remaining core/base attachment orientations arbitrarily, once for the experiment. Let T(x,t) be the resulting tournament, and write χ(x)=Π_{e∈E(R)}x_e.

Then

Σ_{x∈{−1,+1}^d} Σ_{t∈{0,1}} H(T(x,t)) χ(x) = 2(2Q_uv−1).

Equivalently, under independent uniform random orientations of these d+1 edges,

E[H(T(x,t))χ(x)]=(2Q_uv−1)·2^(−d).

Here Q_uv=1 iff u→v. Thus the exact weighted sum is +2 or −2 and directly returns the unknown labelled arc. The formula holds for every choice of the frozen attachment orientations. All tournaments queried have exactly k−1 base vertices.

## Proof by cancellation and four surviving vertex words

Expand H(T) as a sum over Hamiltonian vertex words P. Each word contributes its indicator that every consecutive arc is present. Exchange the finite sums over P and the varied orientations.

If P omits some underlying undirected edge of R, its indicator is independent of that edge's orientation sign x_e. Summing against the factor x_e cancels it exactly. Therefore only words containing every underlying edge of R survive.

In a Hamiltonian path word, every internal vertex of R already has its two incident path edges occupied by R. Thus R must occur as a single contiguous block, either in its displayed direction or in reverse, denoted R^rev. Since R includes every vertex except u, the only surviving candidate words are

uR, R^rev u, Ru, uR^rev.

This is also immediate from edge counts: R has d edges, while a full Hamiltonian path has d+1, leaving precisely one edge to attach u at an endpoint of R.

### Attachment at the core endpoint v

The word uR is valid exactly when u→v in Q and all R signs are +1. Its character is +1. It does not use the varied edge u–b_m, so both values of t contribute: total +2 when Q_uv=1.

The word R^rev u is valid exactly when v→u in Q and all R signs are −1. Since d is odd, its character is −1. Again both values of t contribute: total −2 when Q_uv=0.

Exactly one of these two possibilities is available, giving 2(2Q_uv−1).

### Attachment at the base endpoint b_m

The word Ru requires all R signs +1 and the extra orientation b_m→u. It contributes +1, for exactly one value of t.

The word uR^rev requires all R signs −1 and the extra orientation u→b_m. It contributes (−1)^d=−1, for exactly the opposite value of t.

These two contributions cancel. No frozen attachment edge or base/base edge appears in any of the four surviving words. Therefore the result is independent of all such fixed orientations, and the stated identity follows.

## Reconstruction, parity, and limits

Apply the experiment to each unordered core pair, using a fixed convention u<v. Each pair uses 2^(d+1)=2^(2k−2)=4^(k−1) assignments. The signed sum determines Q_uv, hence determines Q. This proves faithfulness and the query bound. Shared queries between different pairs may reduce the actual number; the theorem only claims the displayed upper bound.

Odd parity is essential to this particular construction: reversal of all d R edges changes χ by (−1)^d. Here d=2k−3 is odd, so it makes the two core-end orientations have opposite signs and cancels the two base-end contributions. An even-d path would not give the same orientation detector by this argument. This explains the parity mechanism, but does not prove any optimality of the number of base vertices or exclude other even-degree detectors.

The transitive base causes no obstruction: all surviving R edges are cross edges, and even R^rev has no consecutive pair of base vertices. In fact the proof works with any fixed orientation of the base edges. The theorem is stated using the requested transitive-base family, which is already sufficient.

These are exact noiseless counts and weighted sums. The normalized correlation has magnitude 2^(−(2k−3)); no efficient sampling/noise-robustness claim is made. The exponential query count and potentially large Hamiltonian counts are explicit. Labels, the distinguished u,v, the order of the remaining core vertices, and the chosen attachment ports are retained throughout. Forgetting them is not covered by the theorem.

## Relation to the earlier finite-automaton bound

The earlier common-slot realization and exact rational equivalence algorithm remain correct. Their generic bound 2k!3^k−1 does not become false; it is superseded for DISTINCT labelled tournament cores by the much smaller k−1 distinguishing-length upper bound proved here. The new argument uses a special combinatorial character of Hamiltonian paths that generic weighted-automaton reachability does not exploit.

Likewise, the quadratic-size two-port tomography theorem for (H,S,E,C) remains valid. It neither implies nor is implied here by injectivity of that sidecar. The full-response theorem recovers individual core arcs by a higher-order signed experiment rather than inferring orientations from the endpoint-cover current. In particular the transitive-three hostile showing current opposite to an actual arc remains relevant and unretracted.

## Exact controls and reproducibility

`faithfulness.py` builds each full extension and counts Hamiltonian paths by subset endpoint DP. It does not call the slot automaton, cover-count reconstruction, or an analytic formula for the response. It checks every unordered core pair for every labelled tournament of orders two, three, and four, using deterministic nonuniform frozen attachments. Additional controls use five distinct labelled order-five cores and both increasing/decreasing endpoint pairs; reversed-base controls confirm that fixed base orientation does not enter the survivor calculation.

Run:

    python3 faithfulness.py
    python3 -O faithfulness.py

Every numerical check raises explicitly on failure and survives optimized Python. Final counts and completion status are recorded in `faithfulness.out`; the optimized output must agree byte-for-byte. The proof does not depend on the finite controls, which serve to audit signs, conventions, and implementation.

## Relabelling and converse quotient consequences

Let π be a permutation of the core labels and let πa denote the same attachment mask with its bit coordinates permuted consistently. Relabelling core vertices in every Hamiltonian word gives

F_{πQ}(πa_1,...,πa_m)=F_Q(a_1,...,a_m).

Therefore the orbit of the marked response under simultaneous core/port-coordinate permutations is a faithful invariant of the tournament's isomorphism class. More precisely, two such response orbits agree iff the cores are isomorphic: one direction is this equivariance, and the other follows by applying labelled faithfulness after aligning the coordinates. This is not invariance under independently forgetting or arbitrarily reassigning port labels.

For a word w=(a_1,...,a_m), let τw=(complement(a_m),...,complement(a_1)), where every complement is on exactly k bits. Reverse every arc of the full extension and reverse the order of its base labels to restore the required increasing transitive base order. Reversing each Hamiltonian word preserves its count and gives

F_{Q^op}(w)=F_Q(τw).

The operations τ and simultaneous core-coordinate permutation commute. Consequently response orbits under their combined action agree iff the cores are isomorphic or converse-isomorphic. In a metagraph whose nodes merge exactly these two relations, this gives a faithful invariant of the nodes. It does not construct the metagraph's adjacency/wiring efficiently, count its nodes efficiently, or give a practical speedup for tournament enumeration such as A000568.

The stated query count counts calls to an exact Hamiltonian-response oracle. It excludes the computational and bit complexity of evaluating those calls, storing large counts, deciding graph isomorphism, or enumerating quotient classes. Those costs are not silently removed by the reconstruction theorem.

## Optional finite recognition of a proposed rational response series

There is a further exact decision consequence, not an implementation claim. Suppose k is fixed and a rational weighted automaton A of dimension n is explicitly supplied over the attachment alphabet {0,1}^k. We can decide whether its series equals F_Q for some labelled k-core tournament Q:

1. For every unordered core pair, evaluate A on the same fixed Fourier battery above. If any signed sum is not +2 or −2, reject. Otherwise use its sign to construct the unique candidate labelled tournament Q. If both ordered directions are tested instead, their signs must be opposite; this is an additional necessary consistency check, not an extra source of orientations.
2. Construct the inherited D=2k!3^k-state slot realization for Q, and test exact rational weighted-automaton equivalence between it and A. Accept iff equivalent. A direct-sum difference realization has dimension at most n+D, so reachability gives a finite test and, on inequivalence, a distinguishing word of length at most n+D−1.

Correctness: if A is a genuine response series, faithfulness forces step 1 to extract its unique Q, and equivalence passes. Conversely, acceptance in step 2 certifies equality to an actual Q response. Passing only the Fourier sign tests does not suffice. For k=0,1 the sole core can be compared directly, without extraction. This is a checkable recognition procedure for an explicitly represented rational series, not a claim that arbitrary degree/count data determine a tournament or its metagraph wiring.
