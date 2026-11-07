# Every nonnegative integer is a tournament Hamiltonian-cycle count

**Status: elementary all-order proof; independent proof audit PASS; finite exact controls in `RESULT.json`. No literature-priority claim.** Hamiltonian cycles are directed spanning cycles counted modulo cyclic rotation. Reverse traversal is not a second orientation of a directed cycle. The ambient tournaments have at least three vertices.

## 1. Inheritance and the question being answered

Source snapshot: `c137c698e632cb76e232a79b48133987755e080a` on `eliottcassidy2000/math`. This packet combines independent binary-comb and rooted-endpoint constructions developed in this session. They are the same marked family after swapping the two apex labels, with the opposite apex arc exposing one extra cycle. The companion `ENDPOINT_BRIDGE.md` gives the exact size-indexed endpoint-spectrum dictionary and an independently derived endpoint proof. `SOURCES.json` pins the inherited results; neither the gate bijection nor the sole-exit closure is being rebranded as new.

Session210's recognizable closure `Gamma(T): a -> b -> T -> a` satisfies
`HC(Gamma(T))=H(T)`, by deleting the forced gate arc. Its degree-defined gate
is an arc `a -> b` with `d+(a)=d-(b)=1`. Session211 shows that deleting the
mark, or even retaining the complete labeled score vector, loses this
faithful path carrier. Its examples realize seven and twenty-one cycles,
although seven and twenty-one are forbidden **unconstrained Hamiltonian-path
counts** of tournaments.

The closest inherited mechanism is the single-exit attachment law of
THM-4862, `HC(U+a_exit)=s_r(U)`, where `a_exit` beats only `r` and `s_r`
counts Hamiltonian paths starting at `r`. One fixed starting vertex changes
the counting problem. The present family makes those endpoint constraints
explicit and realizes all integers by a binary grammar.

The concept board is: recognizable gates / fixed-endpoint paths / transitive
chains / binary subsets / sole-exit vertices / strong connectivity. The
anchor is the unbounded cycle spectrum; the niche is the restricted
nongated spectrum; the wildcard is a genuine doubling operation retaining
all marked tournament data. The hostile boundaries are zero, tiny orders,
and confusing a sole exit with a full degree-defined gate.

The connection contract is exact: binary words map to tournaments with two
named extra vertices and a named transitive chain; the preserved statistic
is the cycle count. Forgetting the names loses this particular construction
and its local recursion, but does not change the cycle count. This is not a
count-preserving transport to unrestricted tournament Hamiltonian paths and
has no consequence for the open Chan–Pak path-spectrum conjecture.

## 2. A binary comb and its exact two-apex response

Let `m>=1`, `q=sum_{i=0}^{m-1} epsilon_i 2^i`, with each `epsilon_i` zero or
one, and `tau in {0,1}`. Begin with the transitive chain

    x0 -> x1 -> ... -> x(m+1),

where every earlier chain vertex beats every later chain vertex. Add vertices
`a,b` with the following orientations:

- `a` beats `x1` and loses to all other chain vertices;
- `b -> x0`;
- `x1 -> b` exactly when `tau=1`;
- `x(i+2) -> b` exactly when `epsilon_i=1`.

There are two tournaments: `B+(m,q,tau)` has `a -> b`, and `B-(m,q,tau)`
has `b -> a`. These rules specify every unordered pair exactly once.
Both tournaments have `m+4` vertices.

**Theorem 1.** For every allowed parameter choice,

    HC(B+) = 1 + tau + q,
    HC(B-) = tau + q.                                      (1)

**Proof.** Normalize each cycle to start at `a`.

In `B+`, the cycles beginning `a,b` must continue through the entire chain
in increasing order. Exactly one exists:

    a, b, x0, x1, ..., x(m+1), a.

It closes because `m>=1` makes the final chain vertex different from `x1`.
This case is unavailable in `B-`.

Every other cycle starts `a,x1`. Before visiting `b`, it follows an
increasing chain subsequence, necessarily starting with `x1`. After `b`,
it follows the increasing complementary subsequence. This complement
contains `x0`, so it is nonempty. Its last vertex is different from `x1`,
so its last arc back to `a` always exists. Its first vertex is `x0`, and
`b -> x0` always exists.

If the subsequence before `b` ends immediately at `x1`, this contributes
one cycle exactly when `tau=1`. Otherwise its last vertex is `x(i+2)` for
some `i` with `epsilon_i=1`. Any subset of `x2,...,x(i+1)` may occur between
`x1` and `x(i+2)`, in increasing order. These `i` vertices yield exactly
`2^i` choices. No other vertices can precede this final vertex in the
increasing subsequence. Conversely, every listed choice gives a valid
spanning cycle. The alternatives are disjoint and exhaustive, proving (1).

Thus reversing this one specified apex arc subtracts exactly one cycle.
This is a family-specific identity, not a monotonicity law for arbitrary
arc reversal. QED.

## 3. Exact support, order bound, zero, and parity

**Corollary 2.** At each fixed order `N=m+4>=5`, the forward family `B+`
realizes every integer in

    [1, 2^(N-4)+1].                                         (2)

Indeed `tau+q` takes every value from zero through `2^m`. The full
Hamiltonian-cycle spectrum over tournaments of arbitrary order is exactly
`N_0`, and the spectrum over strong tournaments is exactly `N_{>0}`.

For the existence statement, `B+` supplies every positive integer and every
such tournament is strong because it has the explicitly exhibited spanning
cycle. A transitive three-vertex tournament supplies zero. Conversely,
a nonstrong tournament has no spanning directed cycle. The assertion that
a strong tournament of order at least three always has a Hamiltonian cycle
is the inherited classical Camion theorem; it is not needed to verify any
positive witness constructed here.

A concrete target compiler for `k>=2` chooses `tau=1`, `q=k-2`,
`m=max(1,bit_length(q))`. Therefore

    N <= max(5, 4+ceil(log2(k-1)))   for k>=2.                (3)

For `k=1`, the directed triangle has order three. This is an existence
upper bound, not a minimum-order or asymptotic-optimality claim. The
fixed-order interval (2) need not be the full spectrum at that order.

There is consequently no parity obstruction and no missing positive
integer for unrestricted Hamiltonian-cycle counts. In particular, the
path-spectrum exclusions seven and twenty-one do not survive this change
of counting problem.

**Boundary warning.** The forward formula fails at `m=0`: its alleged
direct cycle would need `x1 -> a`, whereas the definition has `a -> x1`.
We use `m>=1` throughout Theorem 1. The reversed variant does satisfy the
corresponding zero-bit formula at `m=0`, but that optional order-four
boundary is not used in the main construction or finite-size claim.

## 4. A complete positive spectrum without any degree-defined gate

**Theorem 3.** Every `B-(m,q,tau)` with `m>=1` is gate-free. Its sole-exit
vertex `a` does not make it a degree-defined gate carrier. Consequently,
every positive integer occurs as the cycle count of a **strong gate-free**
tournament.

**Proof.** A gate requires a vertex of outdegree one pointing to a vertex
of indegree one. Check all possible tails.

- The tail `a` has its sole exit to `x1`, whose incoming neighbors include
  both `x0` and `a`. Thus `a -> x1` is not a gate.
- The vertex `b` has at least two exits, to `a` and `x0`.
- Every chain vertex other than `x1` and the last vertex has at least two
  exits: a later chain vertex and `a`.
- The last chain vertex may have only the exit to `a`, but `a` has
  indegree `m+2`, at least three.
- The vertex `x1` has all `m` later chain vertices as exits, and also `b`
  if `tau=1`. It can be a sole exit only when `m=1,tau=0`, when its unique
  target `x2` receives arcs from both `x0` and `x1`.

No case yields a gate. Given `k>=1`, take `tau=1`, `q=k-1`, and
`m=max(1,bit_length(k-1))`. Equation (1) gives exactly `k` cycles, at every
order `N>=max(5,4+ceil(log2(k)))` by adding leading zero digits. Alternatively
`tau=0`, `q=k`, `m=bit_length(k)` gives the especially simple grammar below. Since `k>0`, a
spanning cycle exists and the tournament is strong. For `q=tau=0`,
`B-` instead has zero cycles and `b` is a universal source, giving the
nonstrong boundary directly. QED.

The recognizable-gate subclass is different: deleting its two gate vertices
recovers a core whose unrestricted Hamiltonian paths biject with the
cycles. Thus its spectrum inherits oddness and the exclusions seven and
twenty-one. Its all-order completeness is the original open path-spectrum
problem. The present theorem proves complete cycle spectra for the ambient
and gate-free classes, without resolving that restricted problem.

**Forward-family gate boundary.** The forward family is not always gate-free. In `B+`, if `tau=1` and all binary digits are one, the coding apex b has the unique exit x0 and x0 has the unique entry b, so b->x0 is a gate. This supplies the top value `2^m+1` of the interval. All reversed-family instances are gate-free, including its zero boundary, as proved above. The wider fixed-order interval is therefore an unrestricted statement; the guaranteed nongated positive interval is `[1,2^m]`.

**Sharp feedback-vertex refinement (connection identified with the classification lane).** Every member of this two-apex family becomes transitive after deleting a and b, so its feedback vertex number is at most two. Any tournament with feedback vertex number at most one has at most one Hamiltonian cycle: delete a vertex leaving a transitive tournament, whose unique Hamiltonian path is the only possible remainder of a cycle. Hence every constructed witness with HC>=2 has feedback vertex number exactly two. In particular every integer k>=2 occurs jointly with strong connectivity, no degree gate, a sole-exit vertex, and feedback vertex number exactly two. This refines the construction, not the separate classification of all tournaments with that feedback number.

## 5. A genuine marked binary growth operation

Keep the marks and chain order of `B-(m,q,tau)`. Insert one new chain vertex
`y` immediately after `x1` and before all the old bit vertices. Orient
`y -> a`; give it the chain comparisons dictated by that position; and
choose `y -> b` when `epsilon=1`, or `b -> y` when `epsilon=0`.
All old comparisons are retained. The resulting marked tournament is

    B-(m+1, 2q+epsilon, tau).

Hence its cycle count satisfies

    C_new = 2C - tau + epsilon.                              (4)

For `tau=0`, the two operations are literally `C -> 2C` and `C -> 2C+1`.
Starting from `B-(1,1,0)`, which has one cycle, the remaining binary digits
of any positive target build that target. Every intermediate is a strong,
gate-free tournament with the same named sole-exit tail and root. For
`tau=1`, the two maps are `C -> 2C-1` and `C -> 2C`.

This is an all-order count-controlling tournament grammar with a proof;
it is not extrapolation from a numerical initial interval. It succeeds
because rooted path counts can be arbitrary integers, while the complete
unrooted path count is subject to Rédei parity and the two forbidden values.

## 6. Verification and scope

Run from this directory:

    python3 verify_binary_comb.py
    python3 -O verify_binary_comb.py

The verifier uses explicit exceptions rather than assertions, so optimized
Python retains the checks. It checks:

- all 2,040 combinations of `m=1,...,8`, every `q`, both `tau`, and both
  apex orientations, through order twelve;
- a separately implemented subset dynamic program counting cycles rooted
  at vertex zero, so there is no unaccounted rotational factor;
- the literal cycle list in the proof, including distinctness and every arc;
- all permutations of the remaining vertices for every family case through
  order seven, checking equality of complete cycle sets;
- the fixed-order support interval, strongness, and the reversed-family
  gate-free claim;
- binary growth, including unchanged old comparisons and both recurrence
  choices;
- explicit witnesses including seven and twenty-one, and exact path counts
  to keep the two observables visibly distinct;
- the inherited gate law on all 75 labeled cores of orders one through four,
  and the transitive zero control.

Finite checks verify the implementation and expose boundary mistakes. The
all-order conclusions follow from the explicit cycle bijection, not from
these finite tests. An independent implementation and adversarial proof
review are supplied by the session's separate audit packet.

## 7. Limits and provenance

This is a complete classification of the **set of possible numerical cycle counts over all orders**, not a classification of realizing tournaments or of the full fixed-order spectra. The explicit order bound is not claimed minimal. For example Session211 already has an eight-vertex witness with 21 cycles, whereas this uniform binary construction uses nine vertices for 21. Seven and twenty-one remain impossible unrestricted Hamiltonian-path counts. The gated cycle spectrum remains equivalent to that open path-spectrum problem.

The independent endpoint verifier is `verify_endpoint.py`; its deterministic receipt is `ENDPOINT_RESULT.json`. It uses a different apex labeling, a separate dense subset-DP implementation, a separately assembled literal cycle family, full permutation checks through order eight on its tau=1 slice, and an explicit m=0 boundary. The audit packet supplies a third implementation. These strengthen reproducibility without elevating finite checks to all-order proofs.
