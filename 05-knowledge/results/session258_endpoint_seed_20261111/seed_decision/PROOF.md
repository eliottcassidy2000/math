# A global finite decision procedure for each exact ray multiplier

**PROVED reduction and constructor. FINITE-EXACT all-order classification for q=1,...,13.** The strong-order estimate is an inherited cited input. This does not construct seeds for every odd slope c and does not prove Chan–Pak.

**Later continuation:** [the strong-core classification](CORE_CLASSIFICATION.md) sharpens the general decision below to a finite strong-core search with an O(log q) conditional witness bound. The nine-state compiler remains an independent positive-state model and supports the backward Zeckendorf guards.

## 1. Target, source, and positive size bound

Use the inherited selected-set seed definition: B is a tournament, O is a selected vertex subset, H=H(B), W counts ordered pairs of nonempty spanning paths whose first endpoint and second starting point both belong to O, and U counts the same ordered covers with only the second start required in O, permitting the first path to be empty. A unit-offset seed of exact multiplier q satisfies

    H−W=1, U=q.

Then the actual tail construction B→TT_t→v, v→O, gives H=1+2^t q for every t≥1. Write q=c2^j, c odd. A slope-c seed is the existential statement that some j≥0 succeeds. This note decides a **fixed q**, not that unbounded existential quantifier.

**Lemma 1.** Every unit-offset seed of multiplier q≥1 has H(B)≤q+1.

Proof. Every W-cover is a U-cover, so 0≤W≤U=q. Thus H=W+1≤q+1. This uses the actual positive cover definitions, not a cancellation estimate. ∎

The closest inherited source is THM-5126, `01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md`, and `odd_c_seed_symmetry_20261110/EXTENSIONS.md`. The new point here is to use this positive size bound before capping the signed offset. The sibling [two-positive-probe compiler](../grammar/PROOF.md) is an independent way to avoid that cancellation.

## 2. A finite complete alphabet, despite arbitrary padding

Every tournament decomposes uniquely into an ordinal sequence of strongly connected components. Hamilton paths concatenate in that sequence, so H(B) is the product of the component H-values. Each component therefore has H-value at most q+1.

We use the **CITED** exact strong-order floor, recovered in `05-knowledge/results/nonzero_preparation_20261103/surgery/ORDER_OBSTRUCTION.md`, §5, from Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13 (2006), N3, Theorem1/Corollary1:

    m(3a)=3·5^(a−1), m(3a+1)=5^a, m(3a+2)=9·5^(a−1), a≥1.

These values increase strictly. Let R(q) be the largest n≥3 with m(n)≤q+1, or1 if none. Every nontrivial SCC of a seed has order at most R(q). Singleton SCCs may occur without bound; simply discarding them would lose the selected endpoint responses.

Define a finite alphabet A_q consisting of every labelled strong tournament of orders3,...,R(q) with H≤q+1, with every possible selected subset, together with the selected and unselected singleton. Different labels, isomorphic copies, or identical transfer vectors may be merged provided one literal representative is kept. **Every marked seed (B,O) is an ordinal word over A_q.** Conversely every word over A_q decodes to an actual marked tournament. Strongness is required for the alphabet-completeness argument, not for the resulting B.

An elementary weaker size input also suffices: a strong n-tournament has a Hamilton cycle (Camion), whose n cuts give H≥n. That yields a finite alphabet with n≤q+1. The sharper Busch input makes the small complete computations practical. Neither classical result is reproved or claimed as new here.

## 3. Only nine positive continuation coordinates

For an ordinal prefix B, partition all its vertices into an ordered pair of directed paths P1,P2, each allowed to be empty. Assign P1 one of three states:

    0 = empty; 1 = nonempty, last vertex unselected;
    2 = nonempty, last vertex selected.

Assign P2 three states:

    0 = empty; 1 = nonempty, first vertex unselected;
    2 = nonempty, first vertex selected.

Let z_ij count these ordered covers. The empty tournament starts at z_00=1 and all other entries0. For a nonempty B,

    H = z_01+z_02,
    W = z_22,
    U = z_02+z_12+z_22.                         (1)

All are nonnegative counts of specified disjoint path pairs. No graph-isomorphism quotient or sign cancellation is hidden in these coordinates.

For a next ordinal block C, compute its analogous vector L_ab exactly. Every prefix path and block path concatenate whenever both are nonempty, since all old vertices dominate all block vertices. A unique global ordered cover splits at that ordinal boundary into its old and local parts. The state transition is therefore

    z'_(a if a≠0 else i, j if j≠0 else b) += z_ij L_ab.            (2)

The first path's last mark is replaced by a nonempty local first path. The second path's first mark is retained once that path has started. These are the asymmetric endpoint roles that an unlabelled two-cover total would lose.

Equation(2) is a nonnegative integer linear map in z for each fixed block. Thus it retains everything needed for all later ordinal blocks and for the seed test(1). It is not a general nine-coordinate compiler for arbitrary vertex births.

## 4. Exact saturation, termination, and emitted seeds

Fix K=q+2 and cap each coordinate by pi(x)=min(K,x). Capped addition and multiplication make pi a semiring homomorphism on nonnegative integers, including pi(0·x)=0. Equation(2) therefore commutes with coordinatewise capping.

Accept a capped state precisely when the ordinary sums in(1) satisfy

    U=q, H=W+1, H≤q+1.                          (3)

**Lemma 2.** Acceptance(3) is equivalent to the actual uncapped seed equations.

Proof. If the uncapped state is a seed, Lemma1 bounds H, W, U, so every coordinate appearing in(1) is below K and the observation is unchanged. Conversely U=q implies all of its contributing coordinates are below K and exact. Hence W≤q is exact. The test H=W+1≤q+1 similarly prevents any H-contributing coordinate from being saturated. Thus(3) is the exact original pair of equations. Large unused channels are harmless and may remain capped. ∎

There are at most (q+3)^9 capped states. Breadth-first search from the empty state over the finite alphabet A_q terminates. If it reaches(3), retain predecessor letters and decode their ordinal sum to a literal (B,O). If it exhausts the graph, no seed of multiplier q exists in **any tournament of any order**: every possible seed had an alphabet word by §2, and its acceptance was preserved by Lemma2.

**Theorem 3 (complete fixed-multiplier compiler).** Unit-offset selected-set seed existence is decidable for every fixed positive integer q. If it exists, a seed with fewer than (q+3)^9 SCC blocks exists, and one can construct such a seed of order at most

    R(q)((q+3)^9−1).

Proof. A shortest accepted word cannot repeat a capped state; deleting the intervening segment would give a shorter accepted word. Each block has at most R(q) vertices. ∎

This is a computable all-order witness bound, not a practical complexity claim. The alphabet itself can be enormous. BFS minimizes SCC-block count, not vertex count; weighting each letter by its number of vertices and using Dijkstra gives minimum order inside this complete representation.

**Safe pruning used by the code.** H and U cannot decrease when an ordinal block is appended. H multiplies by the block's positive H. To inject U-covers of B into covers of B→C, append one fixed Hamilton path of C to P1, keeping P2 and its selected start unchanged; this works also if P1 was empty. Therefore a state with observed H>q+1 or U>q can be discarded. These tests are exact under the cap because their thresholds are below K. We do **not** discard a state merely because an unrelated coordinate saturated.

## 5. Exact small multiplier classification

The code enumerates all labelled marked strong blocks permitted by q≤13. The maximum strong-block order is5 because m(6)=15>14. There are15,762 marked presentations, reduced to146 distinct nine-state transfer vectors with literal representatives. Labelled strong blocks retained by H≤14 number1,2,24,480 at orders1,3,4,5 respectively. These are filtered labelled counts, not ordinary isomorphism counts.

The finite saturated closures decide, without any bound on the total number of singleton vertices:

    Seeds exist: q∈{1,2,4,6,8,11,12}.
    No seed exists: q∈{3,5,7,9,10,13}.

Every positive result includes explicit adjacency rows, selected mask, exact H,W,U, and independently counted one-/two-tail totals. Every negative result is the exhausted finite closure of the complete alphabet, not failure of a bounded graph search. Reproduction takes seconds in this environment.

For example q6 has base rows(4,1,2), O={0,1}, H=3,W=2,U=6 and actual tail totals13,25. It supplies slope c3 at j1. The q11 seed is retained in the receipt.

The negative q7 result forbids a single unit-offset base producing the paired totals15 and29, even though each total individually belongs to the inherited spectrum. Likewise q9 excludes paired19,37, and q13 excludes27,53. This is an exact simultaneous-endpoint obstruction that neither parity nor scalar count coverage detects. The q3 and q10 negatives also have immediate inherited forbidden-total controls7 and21; q5 would produce forbidden21 at its second tail point.

**Mechanism firewall.** The inherited line bank mixes selected-set tail roots with transitive-module doubling vertices. A doubling vertex starting at19 can certify slope9 at exponent1 without being a tail seed of multiplier9. Our absence statements concern the exact (B,O) tail equations, not every possible ray construction, and they do not contradict that bank.

## 6. The remaining universal step

Appending an unselected sink to B preserves H and W and doubles U: assign the sink to the second path to preserve a W-cover; a U-cover has exactly two extensions, by assigning the sink to either ordered path. Conversely deleting the sink recovers the original cover in each case. Thus q-seed existence implies2q-seed existence. Successful exponents on a fixed odd slope form an upper ray.

The compiler therefore isolates a precise next obligation:

    For every odd c, prove that the finite seed machine accepts
    at q=c2^j for some j, preferably with a computable uniform bound J(c).

An exponent bound would turn the unbounded existential statement into a finite family of complete searches. A universal nonemptiness invariant would supply the seeds. Neither is proved here. Even after that step, the finitely many below-start values on each slope must be handled uniformly to establish full Chan–Pak coverage.

The connection to [least-resistance recursion, Zeckendorf, and machines](../grammar/PROOF.md) is literal: cap states retain the two endpoint roles; BFS distance is a terminating rank; the recovered word constructs the graph; a Zeckendorf address names q but does not prove acceptance. Session257's loss ledger is respected: scalar agreement is not permission to splice two graph representatives. This constructor emits a new whole tournament, not a retained-name surgery.

## Reproduction and scope

    python3 check.py
    python3 -O check.py

`check.py` has no external imports and uses explicit RuntimeError checks. It independently checks148 ordinal-composition identities and896 saturation identities. Positive seeds are recounted by subset Hamilton DP on their actual tail graphs. The normal and optimized receipts are compared during session close-out. The proof supplies the general algorithm; only q1..13 is executed and claimed FINITE-EXACT here. No novelty claim for the general semiring saturation trick, SCC decomposition, or classical strong-order floor is made.
