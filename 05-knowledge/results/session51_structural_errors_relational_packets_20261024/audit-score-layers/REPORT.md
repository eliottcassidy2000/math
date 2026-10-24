# Independent audit: score-layer separation and loop obstruction

24 October 2026. **PASS. No mathematical correction required.**

Frozen candidate: `session51-score-layer-separation`, manifest SHA256
`e1a91122cbc42f9de3edfc8e650959a1d79f55158e796c007b8dc34fd31432d9`.
The all-size proof is valid independently of the finite checks below. Its
transport and leg-law dependencies are inherited and properly attributed.
The two earlier blue-wiring audit packets were not modified.

## Independent all-size proof review

Let n=2p, with p vertices of each of two scores a<b. The total degree sum
forces a+b=2p-1, so k=b-a is positive and odd. Counting arcs inside the p
low-score vertices gives pa >= p(p-1)/2, hence

    ceil((p-1)/2) <= a <= p-1.

Thus p=1 and p=2 have only k=1. No existence claim for every arithmetically
allowed layer is needed.

The inherited actual Hamiltonian-path transport, for distinct endpoints s,t,
is

    d_K(v) = 2p-1-d_T(v)+1_(v=s)-1_(v=t).

Squaring, summing, and using the fixed total degree sum gives

    sum d_K^2 - sum d_T^2 = 2(d_T(t)-d_T(s)+1).

Since each transitive triple has a unique source, c3(T) equals binom(n,3)
minus the sum of binom(d_T(v),2). Consequently

    c3(K)-c3(T) = d_T(s)-d_T(t)-1.

This gives the unrestricted necessary jumps {-1,-k-1,k-1}. On a blue
presentation, reversal is an anti-automorphism and takes s to t, so their
scores sum to 2p-1 and must be opposite. Only {-k-1,k-1} remain. The converse
implication is not asserted: opposite endpoint scores do not imply blue.

For an independent algebraic check, the source triangle count can be written

    c_p(k) = p(4p^2-1-3k^2)/12.

If a increases by d>=1, then k decreases by 2d and

    c(a+d)-c(a) = p*d*(k-d).

Writing I_a=[c(a)-k-1,c(a)+k-1], the separation between the later interval's
lower endpoint and the earlier interval's upper endpoint is

    (k-d)(p*d-2).

For any distinct feasible layers, p>=3 and k>=2d+1, so this is strictly
positive. In particular d=1 recovers exactly (p-2)(k-1). This directly verifies
both adjacent and nonadjacent interval separation, including skipped or
unrealized layers.

A target class has a single converse-invariant triangle count, so distinct
layers cannot share any native-line neighbor. Also c(a) belongs to I_a, so a
node in one layer cannot itself be a neighbor of a node in a different layer.
Within a layer, preserving c3 forces d_T(s)-d_T(t)=1, which in a two-score
layer requires k=1 and start score b, end score a. These are necessary
conditions only, not sufficient adjacency tests.

The argument applies to both line colors. It requires no primality, cyclic
cell structure, refinement partition, or self-converse hypothesis for the
unrestricted path statement.

## Distance and loops

Different layers have neither an edge nor a common neighbor. Their distance,
if finite, is therefore at least three in the full converse-merged native-line
metagraph, including both colors. Longer paths through other score types are
not ruled out. This is not a global-disconnection theorem.

At a merged loop, an actual partner is isomorphic to T or to its converse.
Both preserve c3. The same endpoint identity therefore forces

    d_T(start)=d_T(end)+1.

Thus a degree multiset without consecutive values forbids such loops;
unbalanced two-score layers are a special case. The general leg law is
explicitly inherited from THM-4731 statement 4, which covers ordinary and
reversal-neutral self-lines. Nothing here treats adjacent scores as sufficient
for a loop or claims that inherited obstruction as new.

The n=2 case is handled correctly: C is the identity because there are no
non-path arcs. The path identities remain true, but no nondegenerate edge or
loop convention is imposed. The independent receipt calls the algebraic
self-match `self_or_converse_partner`, including this degenerate case.

## Independent finite controls

`independent.py` imports no producer code. Instead of enumerating labelled
sources and traversing their Hamiltonian paths, it enumerates the complete
base-path presentation cube: fix 0->1->...->n-1 and vary every other arc.
It filters the actual two-score sources, reverses the non-path arcs literally,
and directly counts cyclic triples in the source and partner matrices.

There is a bijection between (a normalized presentation, a permutation of
vertex labels) and (a labelled tournament, an actual Hamiltonian path).
Hence each normalized presentation represents exactly n! labelled actual
ports, regardless of automorphisms. Blue membership is invariant under this
relabeling. The independent normalized counts and weighted totals are:

- n=2: one presentation, one blue; 2 actual ports and 2 blue ports.
- n=4: five presentations, three blue; 120 actual ports and 72 blue ports.
- n=6: 158 presentations, 22 blue; 113760 actual ports and 15840 blue ports.

A separate direct enumeration of all labelled tournaments gives 2, 24, and
2720 two-score sources respectively. All labelled-source counts, all 113882
actual ports, and the entire joint histogram of order, layer, c3 jump, and
blue predicate match the frozen candidate exactly.

Full degree transport, the exact blue anti-reflection predicate, actual
layer-to-layer adjacency restrictions, and literal C involution are checked
on every retained presentation. Exhaustive full vertex-permutation tests also
identify whether each actual partner is isomorphic to its source or converse;
all such self-matches obey the endpoint leg law.

Arithmetic checks independently use the closed quadratic in k on all 62750
arithmetically allowed layers through p=500. All adjacent intervals are
checked, and every distinct pair of layers through p=100 is also checked
against the nonadjacent formula above. These controls are not substituted
for the proof.

The order-ten application is exactly:

    scores 3^5,6^5: c3=30; blue neighbor c3 in {26,32};
    scores 4^5,5^5: c3=40; blue neighbor c3 in {38,40}.

For arbitrary paths the additional necessary values are 29 and 39. Actual
same-score black endpoint examples are retained in the independent receipt,
so applying the two-value blue restriction to black paths is explicitly
rejected.

At n=6, the 1^3,4^3 source has c3=2 and a formal necessary entry -2, which
cannot be attained. The checker preserves this negative-value hostile rather
than interpreting the necessary list as an attainment claim.

## Replays, mutation controls, and verdict

Normal and -O independent receipts are byte-identical. Actual copied
directories, not symlinks, reproduce both the independent and producer
receipts normally and under -O from an unrelated working directory. Strict
candidate manifests and the original frozen file inventory are verified.

Normal and optimized modes reject: corrupted inherited source text, removing
the black same-score endpoint value, the wrong interval-gap coefficient,
promoting the negative formal value to a nonnegative attainment assertion,
replacing C with the identity, and corrupting the candidate theorem.
All runtime checks use explicit exceptions rather than removable asserts.

The candidate's proof, finite controls, scope boundaries, inherited attribution,
and necessary-versus-sufficient distinctions pass without correction.
No source packet was changed and no publication was performed.
