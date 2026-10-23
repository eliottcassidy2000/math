# Independent audit: native C3-template complexity

Status: PASS for the precise decision problems in the source THEOREM.md, not for ordinary tournament Hamilton-path existence or a complete-isomorphism transporter. Audit performed 2026-10-23. No literature-priority assertion is audited or made.

## Independent mathematical audit

1. Connected simple graph 3-colourability is a valid NP-hard source. Connecting components in a chain preserves feasibility because each new component can have its three colour names permuted to satisfy the new bridge. Connected inputs with at least two vertices have no isolated vertices, so incidence copies omit no original colour variable.
2. Incidence copies have exactly one inequality edge and at most two equality-chain edges. Equality chains force one phase per original vertex; hence source colouring and pre-Eulerization phase assignments are bijective. The auxiliary graph is simple and connected.
3. Pairing odd vertices with fresh two-edge unrestricted paths never creates a repeated cell pair. Degrees become even and at most four. Adding two leaves at one original cell raises its degree to at most six. The only odd vertices are the two leaves. An Euler trail uses every edge once; each internal degree-d cell occurs d/2 times and each leaf once. Thus every cell occurs at most three times.
4. Padding the first visit to length 4-r gives exactly three local vertices per cell. Chronological labels 0,1,2 make every padded edge a directed cyclic edge. A shift preserves these edges. Cross-cell boundaries occur once per unordered pair because the constraint graph is simple.
5. The actual full block D=E+(v-u) gives shifted arc condition (v+a_j)-(u+a_i) in D, exactly a_j-a_i in E. All three palettes in the reduction are negation invariant, so reversing a constraint traversal does not alter its meaning. Unused blocks may be filled uniformly without imposing new template constraints. Circulant blocks have equal row and column sums, so the C3 partition really is equitable.
6. Every variable phase tuple gives a distinct literal vertex order, even when the host has automorphisms. Each parity vertex and each leaf is independently free. The count is therefore exactly 3^(k+2) times the number of labelled proper source colourings. A global simultaneous shift acts freely and preserves every condition, so one original-cell indexed-owner pin divides the count by three and preserves feasibility.
7. The regular enlargement uses BOTH internal copies of T, not T and its converse. With diagonal cross entry zero and off-diagonal entry T[j,i], every A vertex has score d_i+1+(N-1-d_i)=N and every B vertex d_i+(N-d_i)=N. The center has score N. Matched C3 cross blocks are transposed directed-cycle matrices with zero diagonal; unmatched ones are transposes of circulant/uniform blocks. These are row- and column-regular, establishing the equitable partition including the center singleton.
8. A B-suffix Hamilton path is found by insertion in polynomial time. An indexed occurrence of vertex B_(i,0) required to stay B_(i,0) forces phase a_i=0 in that cell, regardless of where its other occurrences lie in the suffix. Thus these are genuine owner pins, not an assumed independent suffix orientation. Every A vertex points to the center, and the center points to every B vertex. The accepted regular-template orders have exactly the same free A-phase assignments as before regularization.
9. The explicit tournament has O(m) vertices and O(m^2) adjacency entries, with all preliminary transformations polynomial. A phase assignment and owner constraints are polynomially checkable. NP-hardness plus membership establishes NP-completeness for the stated two problem types.
10. The bounded-treewidth positive result is the usual exact finite-domain binary-CSP dynamic program with states carrying shared phases. Within-cell boundaries contribute either tautologies or impossibility, not independent phase variables. Additional guards require explicit state size and update semantics; the source correctly makes no arbitrary electrical-guard claim.

## Guard and scope audit

The theorem preserves supplied cell identities and actual indexed owners. It does not identify all isomorphic tournament paths, does not allow independently chosen phases at different occurrences of a cell, and does not assert that allowed rotations are global automorphisms. The pinned regular-host statement genuinely needs its stated pins and singleton; neither may be silently deleted.

Converse transport is available by reversing the entire path and tournament and transporting all cell labels and owner pins. This is a bijection, not an extra source of NP-hardness. No theorem about an additional nontrivial mod-8/electrical restriction follows automatically: imposing such a guard can remove all witnesses or otherwise change the language. The source's absence of that stronger assertion is appropriate. This packet audits the native-template complexity layer, not any unspecified guard interface.

## Independent computational checks

check.py imports no producer code. It uses a different incidence ordering and Euler-trail implementation, builds an actual tournament, and verifies every one of the nine local phase pairs for every cross boundary. These local identities certify equivalence for every global phase tuple, not merely sampled tuples. It checks literal phase counts on the smallest hosts, ordinary Hamilton paths, regular degrees, all equitable-cell row sums, actual owner pins, and automatic suffix joins.

The 47 source cases comprise every connected labelled graph on 2,3,4 vertices, plus K5, C5, K6, C6. K4/K5/K6 are infeasible prescribed families while their actual tournament hosts retain ordinary Hamilton paths. The independent checker and optimized (-O) execution produce identical PASS receipts. The producer separately reports 56 cases with generic native-CSP counts; that additional coverage is not represented as independently repeated here.

## Primary complexity source

M. R. Garey, D. S. Johnson, L. Stockmeyer, Some simplified NP-complete graph problems, Theoretical Computer Science 1(3), 237–267 (1976), DOI https://doi.org/10.1016/0304-3975(76)90059-1 . The publisher abstract explicitly identifies Graph 3-Colorability among the problems whose NP-completeness survives bounded-degree restrictions. Verified publisher entry: https://www.sciencedirect.com/science/article/pii/0304397576900591 . Accessible author-paper scan also located at https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/JohnsonDavid2.pdf . The audit uses only ordinary graph 3-colourability as an inherited NP-completeness premise; the connected-input reduction is supplied above.

## Limits

The audit proves no unrestricted regular-tournament path-selection hardness, no general native-port theorem with all historical guards attached, and no general counting complexity classification. Exact source-colouring multiplicity is a proved finite-fibre identity. Computational certificates corroborate the implementation; the all-size conclusions rest on the explicit reduction proof.
