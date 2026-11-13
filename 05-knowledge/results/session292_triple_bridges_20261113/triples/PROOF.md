# Occupancy before mutation: a polynomial pivot atlas and a Sierpiński precision law

**Status: PROVED by the elementary arguments below; independent session review pending.** No Chan–Pak completeness, physical interface replacement, or new general cluster theorem is asserted.

## Inheritance and live concepts

The closest mechanism is [Session291's marked squarefree elimination atlas](../../session291_cluster_ports_20261113/cluster/PROOF.md), itself based on THM-4846(F)'s path matrix. This note strengthens its concrete compiler rather than rediscovers that carrier. The hostile is the impossible ring substitution `x_v -> 1`: occupancy is not optional. The corrected near miss is that arbitrary formal Schur exchanges are not physical tournament changes. A least-used sidecar is the **destination ideal** `(x_j)` of each path entry.

The other recovered sources are [THM-4754](../../../../01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md), whose Pascal/Frobenius reduction supplies a genuine Sierpiński operation; [THM-4755](../../../../01-canon/theorems/THM-4755-blow-up-gains-skew-spectrum-lexicographic-law-ternary-tower.md), whose ternary symmetry tower has only finite-exact valuation statements; and [geometry567](../../geometry_567_collatz_observers_20261004.md), where `(2,3,5)` and `(2,3,7)` denote spherical/hyperbolic incidence, with local rotations and circulation as indispensable data. We do not infer equal structure from those triples. The divided-power product law is standard; the application to this endpoint compiler and its sharp tournament denominator precision bound are the contributions here.

Concept board: named occupancy / source-to-destination path entries / pivot-order exchange / aggregate cycle denominator / prime-power precision / digit carries. Two proof loops below change the first proposed viewpoint: rational pivots collapse further to polynomial ones; ordinary Frobenius truncation strengthens to an optimal prime-power factorial bound.

## 1. Destination ideals remove every individual pivot denominator

Let `R=Z[x_v:v in V]/(x_v^2:v in V)` and `A` be any loopless digraph adjacency matrix. Start with `Q_ij=a_ij x_j` for internal states. As in Session291, add named external source and sink states to retain the required endpoint matrix, with source-to-i weight `x_i` and i-to-sink weight 1. Eliminate internal vertices only. Retain all original variables throughout.

**Theorem 1.** At every stage, for each surviving internal vertex j and any surviving state i, `Q_ij` belongs to `(x_j)`. Consequently, for a pivot k,

`Q_kk^2=0`, `Q_ik Q_kk=0`, and

`Q'_ij = Q_ij + Q_ik Q_kj`.

This polynomial update equals the Session291 Schur update `Q_ij+Q_ik(1-Q_kk)^(-1)Q_kj` integrally, not merely modulo a prime.

*Proof.* Initially every edge entering j has factor `x_j`. If k is eliminated, both `Q_ij` and `Q_ik Q_kj` retain that factor when j is internal. Induction establishes the ideal statement. Both `Q_ik` and `Q_kk` contain `x_k`; their product is zero. Since `Q_kk^2=0`, its inverse factor is `1+Q_kk`, whose second term vanishes after multiplication by `Q_ik`. The path proof is equally direct: a simple path can visit k once, and any two segments sharing another internal vertex multiply to zero. ∎

The statement covers weighted edges over any commutative coefficient ring. It does not cover an interface that has forgotten or identified occupancy variables, and it does not license arbitrary attachments at an eliminated vertex.

**Theorem 2 (five-term two-pivot identity).** Eliminating distinct internal u,v in either order gives, for remaining i,j,

`Q''_ij = Q_ij + Q_iu Q_uj + Q_iv Q_vj + Q_iu Q_uv Q_vj + Q_iv Q_vu Q_uj`.

*Proof.* Expand two polynomial updates. The only extra term from the u-then-v order is `Q_iu Q_uv Q_vu Q_uj`, which contains `x_u` twice. It vanishes. The other order gives the same four added simple words. ∎

For r pivots, the universal template consists of all ordered lists of distinct pivots, with number of words `sum_{j=0}^r r!/(r-j)!` (including the empty word). These counts start 1,2,5,16,65. This is a meaningful occurrence of five: two distinguished pivots have five admissible words. It does not identify a five-element set with Platonic solids or triangle-group parameters. The template is independent of the original graph order, but its squarefree coefficients may still occupy exponentially many monomials.

This yields an exact **positive polynomial atlas** behind the rational cluster-like equation. It supplies a checkable candidate replacement identity without divisions; it does not supply a new tournament realizing that identity.

## 2. Aggregate inverses have an optimal factorial precision cutoff

For any zero-constant `r in R`, write `r=sum_{nonempty S} a_S x_S`.
Define `r^[k]` as the sum over unordered collections of k pairwise disjoint nonempty supports S of `prod a_S x_(union S)`; set `r^[0]=1`.

**Theorem 3 (integral divided powers).** `r^k=k! r^[k]`, and

`r^[a] r^[b] = binom(a+b,a) r^[a+b]`.

*Proof.* A surviving term of `r^k` has k pairwise disjoint, hence distinct, supports. Its k! orders have equal weight. To partition an unordered `(a+b)`-collection into groups of sizes a and b takes exactly the displayed binomial number. ∎

Thus the exact mod-p multiplication table of these coefficients has the Lucas/Sierpiński carry pattern: if adding a and b in base p has a carry, the binomial coefficient is zero modulo p. More precisely its p-adic valuation is the number of carries, by the factorial identity `v_p(k!)=(k-s_p(k))/(p-1)`; this is a standard Kummer/Lucas interpretation, not a new number-theoretic theorem. The characteristic-two Sierpiński triangle is therefore an actual multiplication law of marked occupancy packets, not only a picture resembling tournaments.

Let `K(p,s)=min{k>=1: sum_{j>=1} floor(k/p^j)>=s}` for a prime p and positive s.

**Theorem 4 (uniform precision bound).** In `R/(p^s)`, every zero-constant r satisfies

`(1-r)^(-1)=sum_{k=0}^{K(p,s)-1} r^k`.

The bound is independent of |V| and is sharp uniformly over such algebras. For s=1, K=p, recovering Frobenius. For p=2 and s=2, K=4; for p=2 and s=3, K=4, rather than the naive 2s=6.

*Proof.* Theorem3 makes `r^K` divisible by K!, hence by p^s; the geometric identity proves the inverse. Sharpness: take `r=x_1+...+x_(K-1)`. The coefficient of the full monomial in `r^(K-1)` is `(K-1)!`, not divisible by p^s. ∎

**Theorem 5 (sharpness already in tournament cycle denominators).** Put `D_T=det(I-AX)` and `r_T=1-D_T`. The bound in Theorem4 remains sharp uniformly over tournaments. For K>=2, take the ordinal sum of K−1 directed triangles. Let `t_i` be the product of the three variables in triangle i. Then

`D_T=prod_i(1-t_i)`;

in `r_T^(K-1)`, the coefficient of `prod_i t_i` is `(K-1)!`. Thus no smaller universal geometric-power cutoff works modulo p^s, even for tournament denominators.

*Proof.* The adjacency matrix is block upper triangular, and each triangle determinant is `1-t_i`. The lowest possible total t-degree in a power of r is its power; at degree K−1 only the linear terms `sum t_i` contribute. Apply the factorial count. ∎

## 3. Hostile control: fewer geometric powers does not mean fewer cycles

For the same ordinal sum, `D_T^(-1)=prod_i(1+t_i)` integrally. Every full-support coefficient is 1. Modulo2 the formula is also `1+r_T=D_T`, because all signs disappear. Although the geometric series needs only two terms modulo2, r_T already contains arbitrarily large unions of cycles. Hence the precision theorem is **not** a bounded-cycle-order algorithm, a small state grammar, or a proof of scalar odd-count coverage. It controls one algebraic dimension (geometric power), not support complexity.

Likewise applying the full-support observer before multiplying divided powers destroys the disjoint-support guard. The scalar count is not a multiplicative evaluation of the occupancy algebra.

## 4. A precise boundary for the triples (2,3,5) and (2,3,7)

Let `Aut±(T)` consist of vertex permutations that either preserve every arc or reverse every arc of a fixed tournament T, with at least two vertices. Its sign map to `{+1,-1}` is a homomorphism, with kernel `Aut(T)`.

**Theorem 6 (triangle-action obstruction).** For every odd r coprime to 3, every homomorphism

`<s,t | s²=t³=(st)^r=1> -> Aut±(T)`

is trivial. In particular, neither a (2,3,5) nor a (2,3,7) triangle action can act nontrivially by global symmetries/converse symmetries of one tournament.

*Proof.* An involutory tournament automorphism is the identity: a nonidentity involution exchanges some pair of vertices, reversing that pair's arc. The image of t has positive sign since its cube is the identity. Since r is odd, the sign of st must also be positive; hence s has positive sign. Thus s is an automorphism whose square is the identity, so s is the identity. The remaining t satisfies both t³=1 and t^r=1, so t is the identity because gcd(3,r)=1. The one-vertex case is immediate. ∎

This is an application of the standard no-involution observation, not a new theorem about triangle groups. Its role is to tell us **where a transfer must live**: in an ambient quiver exchange space, a family of different tournaments, or an auxiliary incidence/flow carrier, not as this action on one tournament. The Sierpiński tower's order-three symmetries survive; demanding the complete two-generator triangle relation collapses them. This matches the earlier geometry567 distinction between a graph, its rotations, and its surface map. It is a group-action obstruction, not a ban on all mathematical connections involving those triples.

If r is divisible by 3, the C3 quotient obtained by s=1 can survive. If r is even, the sign argument fails, as it should; a directed triangle's rotations and converse reflections give a nontrivial (2,3,2) action. These are sharp controls on both hypotheses.

## 5. What transfers to the constructive problem

Source: endpoint-weighted tournament network. Target: positive polynomial pivot charts and prime-power marked coefficient packets. Preserved: every endpoint/subset coefficient, and therefore H after full-support extraction. Lost only if one subsequently projects: the occupied vertex support and named ports. Needed sidecar: the full support mask, or a proved continuation-sufficient quotient.

A useful new task is now purely polynomial: find distinct physical marked tournament templates for which these five-term (or r-pivot) endpoint polynomials agree on a renewable collection of required interfaces, while their desired output readouts differ by a controlled amount. Modular tests can use the optimal K(p,s) cutoff, but must retain support. This is a sharper experimental formulation, not a claim that the templates exist for every target.

The `(2,3,5)/(2,3,7)` curvature comparison remains a different carrier. Ordinary FZ mutation and the Sierpiński divided-power law can inform which operation coordinates to retain; no map from a triangle group to tournament count coverage is supplied. Chan–Pak remains OPEN.

## Reproduction

Run `python3 triples/verify.py` from this session's directory, or by full path from the repository root; repeat with `python3 -O`. The script compares literal DFS with polynomial elimination, tests random pivot orders and the five-term identity, verifies factorial truncations and the ordinal-triangle sharpness witnesses, and writes `receipt.json`.

## 6. A failed physical proposal yields a reciprocal seam theorem

The first proposed physical synthesis target was a **single saturated two-traversal seam**: two interior tournament choices on I, with identical incidences to retained ports a,b, whose eliminated boundary matrices differ only by `Delta Q_ab=2 x_I x_b`. Global source s attaches to every original vertex with its occupancy variable, and every original vertex attaches to global sink t with weight1. Every entry involving s or t was required to agree. This is stronger than merely keeping H fixed on some smaller graph.

The cheap exhaustive test found no pair at orders3..6, even though many distinct interiors shared a signature with the requested entry omitted. Rather than enlarging the search, the transpose/complement identity gives an **all-order refutation** of that single-entry proposal.

### 6.1 Marked reciprocal identity

For an interior tournament on I, put `X=diag(x_i)`, `P(X)=X(I-AX)^(-1)`, `bar P=P(-X)`, and `bar F=1+1^T bar P 1`.

**Theorem 7.** In the squarefree occupancy ring,

`P(X)^T = -bar P + (bar P 1)(1^T bar P)/bar F`.

Summing every entry also gives `(1+1^TP(X)1) bar F=1`.

*Proof.* The tournament identity is `A^T=J-I-A`. Set `B=I+AX+X=(I+AX)(I+X)`, using `X^2=0` as a matrix. Then `I-A^TX=B-JX`. Since `B^(-1)=(I-X)(I+AX)^(-1)`, both `XB^(-1)=-bar P` and `1^TX B^(-1)=-1^Tbar P`. Apply the rank-one inverse identity to `B-1(1^TX)`, and multiply on the left by X. Its denominator is `1-1^TX B^(-1)1=bar F`, a unit with constant1. This proves the displayed identity. Its entry sum simplifies to `F=1/bar F`. ∎

This is a marked squarefree version of standard matrix inverse algebra, not a new Sherman–Morrison theorem. The search for antecedents included THM-3181's path-colour reciprocity and THM-3315's centered coronal switching compiler; neither is this particular marked seam statement, and no global priority is claimed.

### 6.2 Frozen observer margins force the reverse seam change

Let `r_i=1` iff `a->i` and `c_i=1` iff `i->b`; these incidence vectors are fixed between two interior choices. After eliminating I,

`Q_ab/x_b = a_ab + r^T P c`,

`Q_ba/x_a = a_ba + (1-c)^T P(1-r)`.

The quotients here mean removal of the known final-port variable, not inversion of a nilpotent variable. The global boundary entries encode

- `Q_st=1^TP1`;
- `Q_at=1+r^TP1`, `Q_bt=1+(1-c)^TP1`;
- `Q_sb/x_b=1+1^TPc`, `Q_sa/x_a=1+1^TP(1-r)`.

**Theorem 8 (reciprocal seam invoice).** If all five displayed observer entries agree between the two choices, then

`Delta(Q_ba/x_a)(X) = - Delta(Q_ab/x_b)(-X)`.

*Proof.* Frozen totals and margins imply `1^TDelta P1=0`, `r^TDelta P1=0`, `c^TDelta P1=0`, `1^TDelta Pc=0`, `1^TDelta Pr=0`. Therefore `(1-c)^TDelta P(1-r)=c^TDelta Pr`. Apply Theorem7 between r and c. The rank-one correction has identical numerator and denominator for the two interiors because its factors are precisely the frozen r-row margin, c-column margin and total, evaluated at -X. Thus `r^TDelta P^Tc=-r^TDelta P(-X)c`, which is the result. ∎

In particular, if `Delta Q_ab=d x_I x_b`, the reverse change is forced to be

`Delta Q_ba=d (-1)^(|I|+1) x_I x_a`.

A nonzero single-entry seam with those observer guards is impossible at **every** size. If |I| is odd, the two changes have the same sign; if |I| is even, their signs oppose. No claim is made from these guards about arbitrary other boundary entries; a physical replacement test must still compare them.

### 6.3 Repaired constructive target, with a sharper hostile test

The repaired candidate is a **paired saturated seam**: the only changed boundary entries are the two forced entries above. The retained boundary must include every exterior vertex incident to I (all exterior vertices for a tournament), not just the two specially named ports, together with the global observer terminals.

For nonempty I, the perturbation E of the boundary matrix satisfies `E R E=0` for every R, because each entry of E contains x_I. Hence for the exterior network resolvent `R=(1-Q)^(-1)` the exact correction is `R'-R=R E R`. This is a finite polynomial identity; the inverse exists since the zero-constant boundary network is acyclic and every return has positive occupancy.

The full-support count difference is therefore

`Delta H = d [D_ab + (-1)^(|I|+1) D_ba]`,

where each D is the nonnegative integer count of disjoint exterior prefix/suffix completions after its saturated seam is used. For odd |I|, d=1 and `D_ab+D_ba=2` would give the desired increment2. At even |I|, a signed difference, not a positive sum, must be controlled. These exterior readiness integers cannot be replaced by the mere existence of two ports.

The first repaired proposal was to find a physical paired saturated seam and then renew it from H=23. The same exhaustive search through order6 found no paired seam with nonzero changed coefficients either. The stronger theorem below **refutes that paired proposal as well whenever all the stated observer guards remain frozen**. Thus the right next constructive target must change a source/sink or interior subset-count response and compensate it elsewhere. Merely adding more port-to-port changed entries while freezing those observers cannot work.

The generation-to-proof loop here is substantive: a positive formal single-use perturbation looked sufficient; exact physical tests found no small witness; a tournament-specific inverse identity refuted it at all orders and specified the mandatory paired replacement. Chan–Pak remains OPEN.

### 6.4 Exact search scope and receipts

`seam_probe.py` and `paired_seam_probe.py` exhaust every labelled tournament at n=3,4,5,6, with ports0,1, interior2..n-1, fixed port/interior incidences and the port edge. The boundary retains both global observer terminals. Signatures omit respectively one or two full-interior port coefficients, and compare every other coefficient exactly. The signature counts are8,56,832,26728 out of8,64,1024,32768 tournaments; all collisions have identical omitted coefficients. These are finite controls, not the proof of Theorem8 or a general obstruction to multi-entry moves.

## 7. Full frozen-observer rigidity: even paired seams must pay an invoice

**Theorem 9 (boundary reconstruction from subset counts and degree-two port data).** Fix a labelled interior I and a labelled retained boundary B. Consider two tournaments on I with the same incidences between I and B and the same edges within B. Include global source/sink observer terminals as above. If

1. the eliminated `Q_st` polynomials agree, and
2. for every retained vertex p, the degree-two part in interior variables of `Q_pt` agrees,

then **every** eliminated boundary matrix entry agrees. No port-to-port change, single, paired, or otherwise, is possible under these guards.

*Proof, step1: detect the modules.* Let D be the undirected graph on I whose edges are precisely the interior pairs whose direction differs between the two tournaments. The coefficient of x_S in Q_st is H(T[S]) for every nonempty S. Equality on triples in particular means identical cyclic-triangle predicates, since a transitive triple has1 path and a cyclic triple has3. If uv is an edge of D and w lies outside its connected component C, both wu and wv have unchanged directions. If w saw u and v differently, changing uv would change whether the triple is cyclic, contradiction. Thus w sees the ends of every D-edge in C uniformly; connectivity makes w see every vertex of C uniformly. Every component C is a common tournament module. Edges between components are unchanged.

*Step2: ports also see the modules uniformly.* Write r_i=1 iff p->i. The coefficient of x_u x_v in Q_pt is `r_u A_uv+r_v A_vu`. On a flipped pair uv, equality of that coefficient forces r_u=r_v. Connectivity again makes every retained port uniform on every component C. This is why the small degree-two observer, not merely a large final count, is decisive.

*Step3: module visits are controlled by subset path counts.* For a subset S of C and k>=1, the number of ordered covers of S by k nonempty directed paths is

`sum_(S_1,...,S_k ordered nonempty partition of S) prod_j H(T[S_j])`.

All these numbers agree for the two interiors by condition1. Any boundary-to-boundary path decomposes uniquely into maximal consecutive runs through the D-components. Its sequence of visited components and retained endpoint incidences is determined by the unchanged quotient tournament. For each component, its ordered runs use disjoint labelled subsets and are counted by the displayed ordered-cover formula. Hence every coefficient of every boundary entry is identical. Closed internal excursions at a retained vertex use the same decomposition; direct edges are fixed. This proves the complete boundary equality. ∎

This is a scoped reconstruction theorem about **two fixed labelled interiors with fixed exterior incidences**, not a classification of all tournaments, a general reconstruction theorem from unlabelled decks, or a claim that arbitrary endpoint observations determine an interior tournament. Distinct interiors can remain indistinguishable (as the search collisions show); the conclusion is exactly continuation-interface equality.

**Constructive corollary.** Any useful physical replacement in this fixed-incidence setting that changes a boundary response must pay at least one of two invoices: some interior induced-subset Hamiltonian count changes, or a degree-two source response at a retained port changes. A globally invisible saturated correction cannot be hidden until the final full-support coefficient. This is stronger than simply requiring the reverse port change.

**Next target (OPEN).** Search for compensated moves on two disjoint interior regions where a deliberately changed degree-two port response or subset-count response cancels in the common exterior, while the full Hamiltonian readout changes by2 and a new usable interface survives. The required proof must retain the changed response explicitly through composition; Theorem9 specifies what cannot be frozen. No such renewable pair is claimed here. This is a new constrained synthesis problem, not a restatement that positive polynomials automatically have tournament realizations.

The method has now completed two hostile loops: one-entry positivity led to reciprocal pairing; paired positivity under fully frozen observers led to module reconstruction rigidity. Both failures identify an exact coordinate that a genuinely constructive move must change rather than hide.

### 7.1 A finite-degree complete observer, using the corrected C3-fibre theorem

**Corollary 10 (PROVED using a CITED structural input).** In Theorem9, replace condition1 by equality of the degree-three part of Q_st, and retain condition2 and the fixed incidence assumptions. The complete boundary matrix still agrees.

*Dependency and proof.* The corrected invariant statement in [THM-4803](../../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md) §3 says the full ordered path-cover profile, including H, is determined by the labelled cyclic-triple hypergraph. Its realization-fibre structure uses the cited Boussaïri–Ille–Lopez–Thomassé and Boussaïri–Chergui–Ille–Zaidi results; we do not reprove that external classification. Equality of cubic Q_st coefficients is precisely equality of the cyclic-triple hypergraphs. This remains true after restricting to every labelled subset S, so the inherited invariant gives H(T[S])=H(T'[S]) for every S, which is Theorem9's condition1. Apply Theorem9. The ordered-run substitution mechanism is the inherited [THM-3121](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md), also proved directly in the displayed ordered-partition count above. ∎

The source theorem's automorphism-order and ordinary component-isotype assertions were retracted; **neither** is used. Only its corrected path-cover invariant survives here. We are not claiming that cubic data reconstructs the interior tournament: it reconstructs this continuation boundary after the named incidence and port guards are included. The cubic inventory still has binomial(|I|,3) labelled entries; the conclusion is an exact observer-degree compression, not constant memory, a fast evaluation algorithm, or new count supply.

**Degree-two port control is necessary.** Let I be a directed triangle and compare it with its converse. All induced-subset H values agree. Add a retained port p that beats only vertex0, keeping these incidences fixed. The degree-one Q_pt values agree. On the pair {0,1}, the coefficient of x_0 x_1 in Q_pt is1 in orientation0->1 and0 after reversal. Thus replacing degree-two port observations by degree-one observations is invalid, already for one retained port and three interior vertices.

This supplies a positive answer to the structural question behind the failed seam: a fixed finite **degree** of correctly marked observations can determine arbitrarily large continuation kernels, even though a few scalar samples without the right labels cannot. Its precise domain is the fixed-incidence tournament fibre above.

**Corollary 11 (complete relative C3-fibre classification; same cited dependency).** Under fixed port/interior incidences and fixed edges between retained ports, two interior tournaments have identical complete eliminated boundary matrices (including s,t) **if and only if** their augmented tournaments on I union B have identical labelled cyclic-triple hypergraphs.

*Proof.* Full matrix equality gives the cubic interior inventory via Q_st. It gives every quadratic Q_pt coefficient as well. For an interior pair u,v and port p, the triangle p,u,v can change its cyclic status under reversal of uv exactly when p sees u and v differently. This is also exactly when the quadratic expression `r_u A_uv+r_v A_vu` changes. Thus equality of those coefficients is equivalent to preserving all cyclic triples with one retained port. Triples with two or three ports are already fixed by the incidence assumptions. This proves the forward implication. Conversely the same case split recovers the hypotheses of Corollary10, which gives full boundary equality. ∎

Thus the corrected C3-fibre theorem, combined with quadratic marked port observations, classifies this entire continuation-equivalence relation. It does not identify the actual tournament isomorphism class or create new scalar counts, but it replaces an exponential list of boundary coefficients by a precise finite-degree structural predicate. The statement concerns a relative fibre with the exterior incidence pattern fixed; removing that condition is a different problem.
