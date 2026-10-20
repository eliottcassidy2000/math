# Independent audit: closure, orientation, asymptotics, and reversal criterion

Verdict: **PASS** for the four requested deductions in the reconciliation theorem and audit notes, subject only to the scope safeguards below. No mathematical blocker was found. This receipt does not independently recertify the entire incoming THM-4730 or the growth packet.

## 1. Fully closed tournaments

**PROVED from an already audited CITED theorem plus FINITE-EXACT small cases.** For a tournament T of order n >= 3, write H(T) for directed Hamiltonian vertex words and HC(T) for directed Hamiltonian cycles modulo rotation. Exactly n HC(T) path words close, so fully closed means H(T) = n HC(T).

The previously source-audited Benhocine–Wojda input, reproduced in `../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md`, section 3, guarantees a Hamiltonian bypass at every n >= 7. Its first vertex beats its last, so its path is open. No fully closed tournament exists at these orders.

At even orders, H(T) is odd by Rédei's parity theorem, whereas n HC(T) is even. Thus equality is impossible. This also excludes n = 4,6, with no census assumption. The parity input is current canon in THM-466, `01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md`, part (i), m = 1.

The independent script exhausts every labelled tournament at n = 3,4,5, uses bitset DFS to count Hamiltonian paths, and canonizes only fully closed candidates by all relabellings. It finds:

- n = 3: 2 fully closed labelled tournaments, one class C3; H = 3, HC = 1, Aut = 3; scores (1,1,1).
- n = 4: none among 64 labelled tournaments.
- n = 5: 40 fully closed labelled tournaments, one class Q5; H = 15, HC = 3, Aut = 3; scores (1,2,2,2,3).

The explicit Q5 model `a -> C3 -> b -> a` has the same canonical code as the sole n = 5 class. The regular cyclic C5 is a negative control: H = 15, HC = 2, with 5 open path words. Hence the classification is precisely C3 and nonregular Q5, not C3 and regular C5. Converse merging does not change this classification because closure is converse invariant.

This settles HYP-9345(a) outright. It does not rely on a census at n = 6,7,8,9, or on an unproved statement beyond order 62.

## 2. Open-to-closed orientation and cuts

**PROVED.** For n >= 3, the endpoint edge of the fixed Hamiltonian path is a nonpath edge. Complementing the nonpath edges reverses it. Every literal complementary pair therefore has one open and one closed end. The colour is unchanged by complement.

Rooting a Hamiltonian cycle gives n closing path words. The free action of Aut(T) on path words gives ordinary class masses

    C(T) = n HC(T) / |Aut(T)|,
    O(T) = (H(T) - n HC(T)) / |Aut(T)|.

Converse-pair merged nodes have twice these masses; self-converse nodes do not receive an extra factor two.

For a node shore S, let I be the number of internal line instances, in the number of lines entering S, and out the number leaving. Counting half-edges gives C(S) = I + in and O(S) = I + out. Consequently

    in - out = C(S) - O(S),
    d(S) = in + out >= |C(S) - O(S)|,
    d(S) = C(S) + O(S) - 2 I.

The same proof works after restricting all three counts to blue lines or to black lines.

Camion gives C(T) > 0 exactly for strong T at n >= 3. Nonstrong nodes are sources, and every directed line ends at a strong node. The classification above gives O(T) > 0 for every strong T at n >= 6. Therefore the induced directed strong core has positive outdegree everywhere. A finite directed walk staying in that core must repeat a vertex, so every strong node reaches a directed cycle.

Scope safeguards:

- This is a bipartition of half-edges, not of nodes. A strong node can carry both kinds of half-edge.
- Directed cycles may be self-loops. Positive outdegree alone does not prove a cycle with two or more distinct nodes.
- A loop contributes one internal line, two half-edges, and one each to its vertex's indegree and outdegree.
- Positive incoming degree is asserted in the full directed multigraph. Positive outgoing degree is also present in its induced strong core. No induced-core positive-indegree claim is needed.
- No connectivity, strong connectivity, unique arrangement, or sufficiency of the cut constraints follows.
- This is not the cyclic-triangle-gradient orientation.

Independent finite controls construct every fixed-path tiling at n = 3..6, form the converse-merged class quotient, and check orientation and colour directly. All shores are checked at n = 3..5, independently for all lines, blue, and black. At n = 6 all singleton shores, the strong shore, the nonstrong shore, the empty shore and the full shore are checked. The n = 6 graph has 1,024 tilings, 512 literal lines, 34 merged nodes and 22 strong nodes; every strong node has positive full indegree and outdegree.

## 3. The nonstrong asymptotic

**PROVED.** Let T_n = 2^binom(n-1,2), n >= 1, and let a_n and NS_n denote strong and nonstrong fixed-path tilings. The ordered strong-component decomposition and the fixed path force the first component to be an initial interval. Its unique size gives

    NS_n = sum_{k=1}^{n-1} a_k T_{n-k},
    a_n = T_n - NS_n,
    a_1 = 1.

The exact exponent identity is

    T_k T_{n-k} / T_{n-1} = 2^{-(k-1)(n-k-1)}.

At size n-1, bounding every a_k by T_k and each product by T_{n-2} gives

    NS_{n-1}/T_{n-1} <= (n-2) 2^{-(n-3)}.

For n >= 4, separate k = 1 and k = n-1. They contribute T_{n-1} and a_{n-1}; the n-3 middle terms are nonnegative and are each at most T_{n-1} 2^{-(n-3)}. Therefore

    2 - (n-2) 2^{-(n-3)}
       <= NS_n/T_{n-1}
       <= 2 + (n-3) 2^{-(n-3)}.

In particular

    NS_n = 2 * 2^binom(n-2,2) * (1 + O(n 2^{-n})).

Both endpoints are asymptotically equal contributors. The coefficient-one proposal in HYP-9345(f) is false; the corrected factor two is proved. Neither the full asymptotic expansion nor the separate transform question is settled here.

Independent integer arithmetic checks the recurrence, exponent identity and both inequalities through n = 100, with all inequality denominators cleared. The quoted ratios at n = 8,16,40 are reproduced exactly before decimal formatting. No proof step depends on floating-point arithmetic.

## 4. Nonconstant pair sums are not a valid anti-automorphism test

**EXACT COUNTEREXAMPLE to the general criterion, not a counterexample to the atom classification or to the incoming explicit path family.** In Circ(Z/9Z, {1,3,4,7}), take

    P = (3,1,2,6,4,5,0,7,8).

It is a directed Hamiltonian path. Its positional reversal is the permutation, indexed by vertices 0..8,

    r = (2,7,0,8,4,6,5,1,3).

The independent audit checks every one of the tournament's 36 directed arcs u -> v and verifies r(v) -> r(u). Yet its reflected path-label sums are

    (2,8,2,2,8,2,2,8,2) modulo 9.

Structural explanation: the three residue classes modulo 3 are directed-triangle modules in a cyclic order. Shifting the residue-zero module by +3, while fixing the other modules pointwise, is an automorphism sigma. Applying sigma to the standard path (0,1,...,8) produces P. Standard reversal x -> 8-x is an anti-automorphism; its conjugate by sigma is r. The conjugate need not be an affine reflection of the original cyclic labels.

Thus nonconstant reflected label sums exclude only the affine reflection form r(x) = c-x. They do not exclude all anti-automorphisms. The incoming proof must supply a restricted-family argument that really rules out non-affine reversals, or exhibit an actual violated arc. The present hostile establishes the unsupported general implication, not that an incoming three-interval or descending-pair witness fails.

Independent direct-arc retesting of all normalized skew connection sets at odd n = 7..23 gives 2,044 connection sets and 41,992 candidate paths. Every candidate has an actual violated arc; none is anti-reversible. This is FINITE-EXACT evidence at the stated orders, not an all-order repair. The existing independently audited prefix/violated-arc proof is separate and is not re-certified by this receipt.

## 5. Reproduction and audit boundaries

From this packet directory:

    python3 independent_audit.py

The output is retained as `independent_audit.json`. The reviewed source digest is retained as fixed provenance metadata in the JSON. The independent script reads no files and can be copied elsewhere without changing its output. The implementation does not import or execute `audit_incoming.py`; it uses adjacency bitsets, DFS Hamiltonian counts, independently assembled fixed-path cubes, and explicit directed-arc tests. Failure checks use explicit exceptions, not removable assertions.

The source-audited bypass theorem is a stipulated input to this audit, not newly rechecked literature. The script and this receipt make no changes to the math checkout, incoming snapshots, other audit files, or frozen manifests. Publication freshness, exact-prior fingerprints and the wider pure-blue proof remain separately checked prerequisites.

## Final text review

The final theorem and audit notes and both narrow canon overlays passed an additional wording/scope review. One requested proof clarification was applied: colour preservation follows because complement commutes with reflection, not merely because reflection fixes the endpoint coordinate. No mathematical blocker remains.

The final targeted refresh at eb6c3e9818e was also spot-checked: both overlays preserve every incoming body byte and warning, the current auxiliary-ML scope is correct, and the new source-reported self-line data are kept separate from merged loops and path rigidity. The newer remote hostile was independently checked on all 36 directed arcs. No additional proof change was required.
