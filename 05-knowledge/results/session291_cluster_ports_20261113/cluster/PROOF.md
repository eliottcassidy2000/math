# Tournament mutation: the frozen scalar corridor and a lossless endpoint exchange atlas

**Status:** PROVED by elementary arguments below; finite controls in `verify.py`; session-level independently reviewed (see ../audit/REVIEW.md).

## 0. Inheritance and type ledger

- The March 7 S34 experiments `04-computation/tournament_quiver_mutation.py`, `quiver_mutation_h_invariance.py`, and `quiver_mutation_mechanism.py` already observed the preservation of H when a single Fomin–Zelevinsky (FZ) mutation stays a tournament. `07-reflections/breakthrough-candidate-map.md` §12 proposed using its exchange graph. We explain the entire all-order tournament-only component, not claim the original observation as new.
- [THM-4602](../../../../01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md) identifies real rank-two Plücker sign patterns precisely with locally transitive tournaments, and explicitly separates positive friezes from arbitrary tournaments. Thus a generic tournament is not automatically a type-A Grassmannian point with those signs.
- [THM-4846(F)](../../../../01-canon/theorems/THM-4846-run-observers-of-tournaments-the-nauru-orientation-the-permanental-compound-t5-stars-and-the-ghost-lemma-behind-its-four-spurious-z6-factors.md) already supplies the squarefree path matrix `P=X(I−AX)^−1` and Jacobi's complementary-minor formula for its sign sector. We extend that carrier to an explicit endpoint-preserving pivot atlas; the path matrix is inherited, not newly discovered.
- Rank-two pentagon periodicity, Schur complementation, and weighted-automaton state elimination are standard mechanisms. Their concrete tournament scope and the endpoint loss ledger are the application here. No new general cluster-algebra theorem is asserted.

**Concept board:** tournament quiver / legal FZ mutations; endpoint-labeled path matrix / nilpotent vertex occupancy; pivot order / exchange relation; port equivalence / continuation tests. The useful bridge is that endpoint-preserving formal elimination has exact exchange laws even where physical tournament mutation is frozen.

## 1. The tournament-only FZ graph is completely classified

Let B=A−Aᵀ be the skew adjacency matrix of a finite tournament. Use ordinary FZ mutation, with integer arc multiplicities and cancellation of opposite arrows.

**Theorem 1.** Mutation at v is again a tournament iff v is a source or sink. It then preserves H.

**Proof.** If i→v→j exists, the new i,j arrow multiplicity is either two (if i→j) or zero (if j→i). A tournament cannot have either. Therefore one of the incoming/outgoing shores is empty. Conversely a source/sink has no through-two-path, so only its incident arcs reverse. A source is forced first in every Hamiltonian path, and after mutation is forced last. Removing/reappending it gives a bijection. ∎

**Theorem 2 (all components, labeled).** If T is transitive of order n≥1, its tournament-only component consists of its n cyclic orders (n=1 is one vertex, n=2 one undirected edge). Otherwise write its unique ordered SCC decomposition as

`a_1 → ... → a_p → Q → b_1 → ... → b_q`,

where Q starts and ends in nontrivial SCCs; its internal SCC word is retained. Put m=p+q. The component is a path on m+1 labeled tournaments. No nontrivial SCC or singleton inside Q changes.

**Proof.** A legal mutation can only remove the first exterior singleton and append it, or remove the last exterior singleton and prepend it. The fixed linear list relative to the core is the concatenation `b_1,...,b_q,a_1,...,a_p` cyclically cut by Q; equivalently choose the position of Q among the ordered exterior list determined by the initial state. Forward source moves end after p moves when Q becomes first; backward sink moves end after q moves when Q becomes last. These m+1 states are distinct and are the entire path. If all SCCs are singletons, every move rotates the entire transitive order instead, giving the n states. ∎

Passing to isomorphism classes makes the transitive component one vertex. In the nontransitive case the m+1 states remain distinct because the number of singleton SCCs before the first nontrivial SCC differs. H is constant throughout each component. Thus a sequence required to stay within tournaments cannot enlarge the scalar spectrum at all.

This is not a statement about mutation sequences allowed to leave the tournament category.

### Endpoint data is not preserved by the scalar bijection

For U with endpoint matrix E_U, the tournament v→U has

`E_T(v,j)=Σ_i E_U(i,j)` and all other rows zero.

After the mutation U→v,

`E_T'(i,v)=Σ_j E_U(i,j)` and all other columns zero.

H agrees, but a fixed named endpoint observable need not. This is a small exact model for why the two ports of one base graph cannot be identified merely because the initial count agrees.

## 2. The ambient pentagon is exact port transport, not new count supply

**Theorem 3.** For any two distinct vertices u,v of a tournament, the ambient quiver word

`μ_u μ_v μ_u μ_v μ_u`

(where operations are applied in the displayed order) returns the original tournament with u and v relabeled. Intermediate quivers may contain missing/doubled arrows. Consequently it transports any marked endpoint by the same transposition and preserves H.

This is the familiar rank-two A₂ pentagon in this special setting. An elementary all-order proof only needs four vertices. Relabel so b_uv=1. For an exterior vertex with initial pair `(b_iu,b_iv)=(x,y)`, its successive pairs are:

| initial | after μu | μv | μu | μv | μu |
|---|---|---|---|---|---|
| −1,−1 | 1,−1 | 1,1 | −1,2 | 1,−2 | −1,−1 |
| −1,1 | 1,1 | 2,−1 | −2,1 | −1,−1 | 1,−1 |
| 1,−1 | −1,0 | −1,0 | 1,0 | 1,0 | −1,1 |
| 1,1 | −1,2 | 1,−2 | −1,−1 | −1,1 | 1,1 |

The final pair is `(y,x)`; the u,v arrow reverses. For any two exterior vertices, each change to their mutual entry is an integer determined only by their two pairs and the current pivot. Substitution of the four rows gives total change zero in all16 cases, regardless of the old exterior entry. That proves the assertion for every order. Opposite initial b_uv follows by exchanging u,v or negating B. `verify.py` checks the 16 local cases as part of its exhaustive order-four controls, as well as every pair in all labeled tournaments through five.

This supplies lawful **ambient chart transport**, but it is only a relabeling on the physical tournament endpoints. It does not implement a change between inequivalent port choices of the same graph. Nor does it prove that every ambient return is such a relabeling.

## 3. Positive state elimination retains every endpoint and vertex subset

For a loopless digraph on V, work over

`R_V=Z[x_v:v∈V]/(x_v²:v∈V)`.

The ideal of zero-constant polynomials is nilpotent. Therefore every `1−r`, r in that ideal, has the finite inverse `1+r+...+r^|V|`.

Build a weighted network with internal vertices V, one source s_i and sink t_i for each i∈V:

- `s_i→i` has weight x_i;
- `i→j` has weight a_ij x_j;
- `i→t_i` has weight1;
- no other edges incident to the sources/sinks.

A walk that repeats an internal vertex has weight zero. Thus the s_i-to-t_j walk polynomial is exactly the inherited `P_ij`: it records every simple directed path, including the singleton i=j, by its vertex set. The coefficient of `x_V` is the Hamiltonian endpoint count.

**Theorem 4 (lossless pivot compiler).** Eliminate any internal vertex k by

`Q'_ij=Q_ij + Q_ik(1−Q_kk)^−1 Q_kj` for i,j≠k.

For every elimination order, the final source/sink matrix is P. Every intermediate Q has nonnegative coefficients and Q_kk is nilpotent. Consequently all inverses exist and the procedure preserves every marked path subset coefficient, not merely H.

**Proof.** Every walk between surviving states has a unique decomposition into excursions through k, and the displayed geometric series sums the consecutive k-return excursions. Nilpotence deletes exactly those concatenations reusing an original internal vertex. All contributing weights are nonnegative. Internal self-return terms have positive degree, so their constant coefficient is zero. Eliminating all internal vertices produces exactly the original source-to-sink path series. Equivalently, these are Schur complements of I−Q, whose principal pivot blocks have identity constant term. ∎

The same argument applies if some real ports are retained rather than eliminated. Any later external attachment on disjoint new vertex variables that meets the old graph only at these retained ports may be performed before or after eliminating the old interior, with identical marked output; the original occupancy variables are retained. Attachments to an eliminated interior vertex are not covered unless that vertex was retained as an addressable interface. This is a genuine compositional interface. It is not small: the algebra has2^|V| basis monomials, and the theorem makes no compression or complexity improvement claim.

## 4. Exact two-pivot exchange and its cluster boundary

Suppose two surviving internal states u,v have loop weights a,d and cross weights b=Q_uv,c=Q_vu. Put

`D=(1−a)(1−d)−bc`.

Its constant term is1. For any remaining i,j, eliminating both gives

```
Q''_ij = Q_ij + [ Q_iu(1−d)Q_uj + Q_iu b Q_vj
                    + Q_iv c Q_uj + Q_iv(1−a)Q_vj ] / D.
```

The two orders are the two positive state-elimination charts of the same rational expression. Their pivot denominators satisfy the exact exchange

`(1−a)(1−d)=D+bc`.

Every permutation of elimination order is generated by adjacent swaps, so this two-pivot identity transports the whole endpoint-preserving compiler between all orders. Squares and braid hexagons in the permutation graph express consistency of these transports. This is an explicit bridge from a cluster-like exchange equation to literal tournament endpoint counts.

**Do not overstate the bridge.** R_V has zero divisors. We invert only the guaranteed units with constant1, not arbitrary path entries or cluster variables. The equation alone is not an embedding of a cluster algebra into R_V, nor a Laurent-positivity theorem for Hamiltonian counts. For C3, `det(I−AX)=1−x_0x_1x_2` has a negative coefficient even though all path counts are nonnegative. Positivity belongs to the geometric-series pivot charts, not every determinant coordinate.

Also a commutative-ring evaluation x_v↦1 is impossible: it would send the relation x_v²=0 to1=0. Taking the full-support coefficient is a linear observer, not a ring homomorphism; it cannot be passed blindly through products or exchange division. Forgetting occupancy would count walks instead of Hamiltonian paths.

## 5. Consequence for the constructive problem

The original tournament-only mutation corridor is H-frozen. Ambient A₂ pentagons restore named-port transport but generate no new counts. The squarefree exchange atlas does retain exactly the endpoint information a core-changing move needs, and supplies interchangeable elimination charts. The missing step is a **physical lift**: turn a useful formal pivot/interface replacement into another tournament with controlled absolute count and a renewable marked interface. No such lift for arbitrary targets is proved here. Chan–Pak remains OPEN.

A precise next hypothesis, not evidence: find a bounded collection of marked physical templates whose squarefree pivot interfaces agree on every continuation coordinate required by a given target-decreasing recursion. Each claimed replacement then has a finite coefficient identity to prove; scalar equality alone is insufficient.

## Reproduction

Run `python3 cluster/verify.py` from the session root (or its full path from the repo root); repeat with `python3 -O`. The verifier exhausts all labeled tournament vertices/pairs through n=5 for the FZ gate, component classification and pentagon; it independently compares DFS path enumeration to the pivot compiler on every labeled tournament through n=4 and deterministic samples at n=5,6, over593 elimination orders. It checks2050 two-pivot exchange identities. These finite checks supplement the all-order proofs, not establish their quantifiers.
