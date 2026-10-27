# Session75 Phase2: the black cut bound fails, the core bound survives this family

**Status: proof candidates and FINITE-EXACT certificates; independent external
audit pending.** No reserved theorem namespace is used or promoted.
Research input pin: `96e35642110104cc547516a63f37eebcf930f707`.

## 1. All-order statement, separate from the finite minimum-order gate

For **every n≥10**, an explicit strong, black-active tournament T has

    λ(T)=1, δ(T)=2,
    d_black(T, nonstrong)=2,
    d_black(T, singleton + strong core)=2.

Thus proposed statement A, black distance to nonstrong equalsλ, is REFUTED.
Statement B, black distance to singleton plus a strong core equalsδ for every
strong black-active tournament, remains OPEN. It holds uniformly on these very
same counterexample families; their λ/δ distinction cannot be erased.

The mechanism is a **rooted Hamilton-word orbit**, not a pure-blue strong atom.
Let U=C3(T_a,C3,1), root its lastsingleton, and let V=U^op. The even tournament
is V→U with the root-root arc reversed. It has a unique minimum outgoing shore
U, and its complete one-step exit cone consists of exactlynine blue words.
At odd order insert a middlevertex w and replace that bridge by
root_U→w→root_V, with other arcs V→w→U. There are exactlytwo minimum shores,
U and U+w, and both exit cones are the same nine blue words. Literal two-black
routes give the upper bound; the separate core supplement supplies a two-black
route ending at one named source with a genuinely strong deleted core.

Proofs and self-contained control scripts are in `family/`; uniform core
routes and the repaired unique-unbalanced-cut theorem are in `survivors/`.
The original family packet is unchanged after its frozen handoff. The general
cut/endpoint current and both complete search cones are proved in
`CONE_AND_LOWER_BOUNDS.md`.

## 2. Exact first member and meaningful hostile

At n10 one labelled member has rows

    (1008,1001,995,997,494,448,640,768,576,48).

Its unique minimum shore is {5,6,7,8,9}. Exhausting ALL657 directed Hamilton
words finds558 black words, and reversing each one leaves a strong tournament.
The nine minimum-shore-decreasing words are all actually blue. The independent
set-of-arcs checker also gives a distinct two-black route to the stronger core
target:

    (0,5,6,7,8,9,4,3,2,1),
    (2,3,7,9,8,6,5,0,4,1).

This is a complete one-step exclusion, not an unsuccessful bounded witness
search. The n11 family head is independently exhausted in the family checker
(1,215 Hamilton words, 1,116 black, all black neighbors strong).

## 3. Minimum order10: complete ORDINARY class authentication through9

The finite gate proves BOTH exact equalities for every strong black-active
ordinary tournament class of orders4 through9. It is not a converse-merged
class count. Pure-blue R5 is explicitly excluded at order5; strong C3 is
separately exhaustively checked as the order3 boundary. There is no non-pure-
blue strong source at orders1–3 within the question's n≥4 domain.

| n | All ordinary classes | Strong black-active classes | Pure-blue strong exclusions |
|---|---:|---:|---:|
|4|4|1|0|
|5|12|5|1 (R5)|
|6|56|35|0|
|7|456|353|0|
|8|6,880|6,008|0|
|9|191,536|178,133|0|

Every source index has one compact record in `finite/witnesses.jsonl.gz`.
For each retained strong source it records a literal named core route of length
δ. The independent replayer recomputes λ=δ, every R/K edge, full path colour at
both ends, actual intermediate strongness, and the final strong core. Therefore
that same route proves the nonstrong upper boundλ. Both lower bounds are
independent endpoint/cut currents. The n7 gate covers ALL353 strong classes,
not just its three regular classes.

Source-universe completeness is independently reauthenticated here, rather
than inferred solely from an inherited class count. For every source,
`finite/authenticate_universe.cpp` enumerates every permutation inside each
outdegree block, recomputes its canonical key and automorphism order, rejects
isomorphic duplicates, and sums n!/|Aut|. These disjoint orbit masses equal
2^(n choose2) at every n4..9, including 2^36=68,719,476,736 at n9. Canonicalization
examines27,282,165 score-block permutations at n9. Since every tournament
isomorphism preserves outdegrees, this tests all possible isomorphisms.
Distinctness plus full labelled mass proves complete ordinary class coverage.
The source inputs and their exact repository paths/hashes are in SOURCE_BINDINGS.

Hence the all-order family's first member at n10 is a **minimum-order**
counterexample to A, not merely the smallest construction tested.

## 4. Sharp surviving all-order repair

If a strong tournament has a unique minimum shore A with |A|≠n/2, then its
black nonstrong distance isλ. Along an A-decreasing geodesic that shore stays
uniquely minimum. An anti-automorphism would map it to the equally minimizing
complement of its image, forcing |A|=n/2. Thus every intermediate source is
non-self-converse and every chosen Hamilton word is black. The even obstruction
has a balanced unique shore; the odd obstruction needs two shores.

## 5. Reproduction and scope boundaries

Run `bash run_verification.sh` from this directory. It builds fresh C++
authenticators, verifies the frozen source hashes, replays all compact paths,
and compares normal/optimized controls. It does not need the shared checkout,
network, a graph library, or the large exploratory JSONL dumps.

The all-order proof, finite minimum-order authentication, and B evidence are
separate. Additional bounded probes in the parent scratch directory are NOT
part of this frozen milestone and are not used to claim minimum order or an
all-order core theorem. No all-order black connectivity assertion is made.
