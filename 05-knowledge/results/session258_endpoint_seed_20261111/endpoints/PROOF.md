# Endpoint constraints as exact nonnegative cut capacities

**Status: PROVED reduction and conditional constructor; FINITE-EXACT controls. Chan–Pak remains OPEN.** No general novelty claim about weighted-cut representations is made. The new use here is a support-sensitive endpoint obstruction/constructor, replacing an integer-ideal test by the actual feasible cut image.

## Target and inheritance

The tournament target is every positive odd Hamilton-path total except **7 and21**. The primary Chan–Pak statement and original2016 provenance are recovered in `05-knowledge/results/odd_c_seed_symmetry_20261110/EXTENSIONS.md`; impossibility of7 and21 is inherited from THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. Write N for a total and c for an odd ray slope. Slopes7 and21 are allowed; the exact unit-offset slope problem is treated in §7. This note does not prove universal support or a universal logarithmic-order constructor.

Closest mechanism: the quadratic one-birth polynomial in `05-knowledge/results/berggren_target_return_20261111/fibres/EXTENSION_CONTENT.md`, and the general extension-content filtration C_k. Corrected near miss: gcd one gives an integer combination, not a feasible orientation. Hostile: a four-vertex strong H=5 source has gcd-one extensions but no H=7 extension. Underused sidecar: the nonnegative odd-cycle expansion, inherited via THM-5071, `01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md` (which credits Grinberg–Stanley and Irving–Omar).

Concept board: endpoint orientation; nonnegative odd-cycle packs; directed cuts; reachability; min-cut versus exact cut; target-aware heavy arcs. Anchor is actual odd-count support, niche is directed cut structure, wildcard is a finite-state constraint solver retaining implications rather than only congruences.

## 1. The exact capacity matrix

Let T be a labelled tournament with vertex set V, h=H(T), and H(empty)=1. For distinct a,b define

    W_ab(T) = sum over directed simple paths P from a to b
              having a positive even number of vertices, of H(T minus V(P)).

W_aa=0. These are nonnegative integers. In particular every arc a→b of T contributes the term H(T−{a,b})≥1. A nonarc can also have positive W_ab, through a longer odd-length path. **Do not discard these nonarc capacities.**

For X⊆V add a new vertex v beating exactly X. Let

    q_T(X) = sum_{a in X, b outside X} W_ab(T).

**Theorem 1 (exact cut image).**

    H(T+v_X) = h + 2 q_T(X).

Consequently a proposed one-birth target N is possible if and only if N≥h, N≡h mod 2, and (N−h)/2 is in the actual directed-cut spectrum {q_T(X): X⊆V}. Given such X, all orientations are explicitly prescribed, so this is a physical constructor, not an abstract state existence statement.

**Proof.** Use the inherited tournament odd-cycle formula

    H(T)=sum over collections C of vertex-disjoint directed odd cycles,
         each of length at least 3, of 2^{|C|}.

Cycle collections in T+v not containing v contribute h. Each remaining collection has a unique odd cycle containing v. Delete v from this cycle. The remainder is a directed simple path P on an even number of old vertices, starting at a∈X and ending at b∉X. Conversely every such path closes uniquely through v. The other cycles form an arbitrary odd-cycle collection on T−V(P), of total weight H(T−V(P)). The v-cycle contributes another factor 2. Summing these disjoint cases proves the formula. ∎

This also proves the exact quadratic identity

    q_T(X)=Σ_ab W_ab x_a(1−x_b).

Hence q is submodular. Its mixed second Boolean difference on a,b is −(W_ab+W_ba). This identifies the nonnegative interaction hidden in the inherited endpoint/two-path-cover polynomial. No assertion that a scalar endpoint count determines W is made.

## 2. Complete zero-gap and unit-gap classifications

**Corollary 2 (zero gaps).** q_T(X)=0 iff no original T-arc goes from X to its complement. Equivalently the complement dominates X; the newborn is a singleton SCC inserted between SCC blocks. If T has s strongly connected components, there are exactly s+1 zero-gap cuts, including empty and full.

Proof: every original arc has positive capacity. Conversely a positive W_ab comes from an actual T path; any path leaving X must use an original arc leaving X. In a tournament the SCC condensation is a total order, and sets with no outgoing crossing arc are precisely its suffixes. ∎

**Theorem 3 (unit-gap criterion and constructor).** A birth with H=h+2 exists iff some ordered pair (a,b) has W_ab=1 and, after deleting that one capacity arc from the support digraph G_W={u→v:W_uv>0}, there is **no directed path from a to b**.

For every such pair take X to be the vertices reachable from a in G_W−(a,b). This gives the required actual birth.

Proof: an integer nonnegative cut has capacity one iff exactly one positive capacity arc crosses and that arc has weight one. Such a cut blocks every alternative a-to-b path. Conversely the reachable set X contains a, excludes b, and is closed under all remaining capacity arcs. Thus its sole crossing capacity is W_ab=1. ∎

The condition “the original arc a→b is a strong bridge” alone is insufficient: other W-arcs may cross, and a nonarc may carry positive capacity. The theorem uses the entire capacity support, not T alone.

A weight-one arc itself has a local combinatorial classification: exactly one even-vertex directed path P from a to b exists, and T−V(P) is transitive. This follows because every path contributes a positive odd integer and H(S)=1 iff S is transitive (the odd-cycle expansion and absence of triangles). For an actual original arc a→b, this says T−{a,b} is transitive and no longer even-vertex a-to-b path exists. This condition must still be combined with the no-bypass test above.

## 3. Target-specific heavy-arc elimination

For a target half-gap r≥0 let G_{>r} have precisely the arcs with W_ab>r.

**Theorem 4 (exact implication reduction).** Every cut of capacity r is closed under every arc of G_{>r}: a∈X and W_ab>r imply b∈X. Therefore X is a union of SCCs of G_{>r}, closed under reachability in its condensation DAG. Conversely enumerate those closed unions and retain exactly the ones whose cut capacity is r; this returns every possible newborn orientation, with no omissions and no false positives.

Proof: one crossed heavy arc would already exceed r, since all capacities are nonnegative. Closure under arcs is equivalent to closure under reachability and to a successor-closed union of SCCs. The final exact capacity test is necessary and sufficient by Theorem 1. ∎

This is a support-aware solver: heavy arcs encode forced endpoint implications, SCCs are inseparable orientation blocks, and light capacities retain the exact budget. It does not claim a polynomial-time algorithm for exact cut values or an all-target coverage theorem. The minimum positive q is a useful coarse obstruction, but even values above it may be absent.

For a strong T all proper nonempty cuts are positive. Their minimum is the minimum directed s-t cut over all ordered distinct pairs s,t in the weighted graph W. Thus a polynomial-time min-cut computation *after W is supplied* gives a proved forbidden interval

    h < N < h + 2 lambda(T).

It does not fill the interval above the threshold. Computing W itself by the supplied exact path enumeration is exponential/factorial in general; no polynomial-time Hamilton-count algorithm is implied.

## 4. What this adds to C_k=1

At depth one,

    C_1(T)=gcd(h, {2q_T(X):X⊆V}).

Since h is odd, this equals gcd(h,{q_T(X)}). The gcd projection forgets all the positive constraints of Sections 2–3. A gcd-one source may have no cut of capacity one, and may have internal gaps in its cut spectrum. The H=5 control below realizes exactly this hostile.

At subsequent births one recomputes the capacity matrix of the *actual new tournament*. One may not retain only h or q or reuse W from the parent. A finite-state or cellular-automaton analogy is lawful only with a sidecar encoding the transition-relevant endpoint/path information. This is a precise reason why unrestricted integer numeration, a prime residue wheel, or Zeckendorf digits do not by themselves construct tournament seeds.

The exact remaining support obligation can now be posed without gcd language: for every allowed target N, produce a realizable predecessor T with H(T)≤N and a successor-closed union of its target-heavy SCCs whose light crossing weight is exactly (N−H(T))/2. This is equivalent to existence of a nontrivial one-birth realization when such a predecessor is furnished. It does **not** prove that predecessors with the needed property exist for all c.

## 5. Executable controls and scope

Run from this directory:

    python3 check.py > results.json
    python3 -O check.py > results_optimized.json
    cmp results.json results_optimized.json

The test universe is all labelled tournaments on 1 through 5 vertices and all 2^n newborn orientations for each. Independent paths are (i) literal subset-DP Hamilton counting of each extended graph and (ii) odd-length simple-path capacities plus cut summation. Every unit-gap criterion and returned reachable-set constructor is checked. Heavy-arc pruning is checked at every realized half-gap. Controls include C3 (positive +2) and a strong four-vertex H=5 parent (gcd one, forbidden +2). Checks use explicit RuntimeError guards, never Python assert, and remain active under -O.

The output records exact totals, unit-gap parent counts, representative witnesses, and the hostile extension spectrum. These are finite observations, not a complete all-order list. No arbitrary c seed construction or Chan–Pak proof is claimed.

## 6. Sharp all-order obstruction: +2 births stop after 5

**Theorem 5.** If a tournament U is obtained from T by adding one vertex and H(U)−H(T)=2, then (H(T),H(U)) is either (1,3) or (3,5). Both possibilities occur. This is not a statement about one-edge changes, two births, or arbitrary pairs of tournaments.

Proof. Put X=N^+(v), Y=V(T)−X. A unit capacity cut has exactly one crossing capacity, of weight one. Since every original arc has positive capacity, there is precisely one original arc a→b from X to Y. All other comparisons between the sides point Y→X. The direct path (a,b) contributes H(T−{a,b}) to the unit cut, so that complementary tournament is transitive, and no other even-vertex directed path may cross from X to Y.

There is no directed triangle within X: a triangle avoiding a would lie in the transitive complement; a triangle containing a can be written x→y→a→x, and x→y→a→b would be a forbidden four-vertex crossing path. Similarly a triangle in Y must contain b, and b→y→z→b would yield a→b→y→z. Hence X and Y are transitive.

Let p be the number of vertices preceding a in the transitive order of X, and q the number following b in Y. If p≥2 there is a path x1→x2→a→b. If q≥2 there is a path a→b→y1→y2. If p,q≥1 there is a path x→a→b→y. These are all forbidden. Therefore p+q≤1.

Every directed cycle in T must cross X→Y via the sole arc a→b. A cycle can cross this cut only once, so it has the form a→b, then a forward segment of Y after b, then a forward segment of X before a, then a. Since p+q≤1, there is no cycle when p=q=0, and exactly one directed triangle when p+q=1. The odd-cycle expansion gives H(T)=1 or 3 respectively. A transitive two-vertex base extended to C3 gives 1→3; a C3 extended with any nonuniform cut gives 3→5. ∎

Thus a proposed strategy repeatedly attaching one vertex to increase H by two fails universally after H=5. This sharp obstruction is compatible with arbitrary even-gap constructors and with all known multi-step/reroot operations.

## 7. The requested odd slope c: exact unit-offset seed equations

**Parameter correction.** In the user's seed problem c is the odd ray slope, not necessarily a Hamilton count. Slopes c=7 and c=21 are allowed. Write **N** for an actual desired Hamilton count. The inherited HYP-9658(A) asks for a unit-offset seed for every positive odd c; the certified ray bank already covers all odd c≤4097, using a mixture of selected-set tails and doubling-vertex mechanisms. Section 1's arbitrary count target should be read as N when discussing this ray problem.

Use the seed construction and definitions of
`05-knowledge/results/odd_c_seed_symmetry_20261110/README.md` and
THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md:
B→TT_t→v, v beats selected O⊆V(B), and

    H = H(B)−W_B(O)+2^t U_B(O), t≥1.

Here W_B is the ordered two-path-cover internal-selected sum; it is **not** our capacity matrix. To distinguish them write K_ab for the ordered cover count, e_a for endpoint counts in B, and A_ab for the matrix W_ab of Section 1. Put

    I(O)=Σ_{a,b in O} A_ab,
    V(O)=Σ_{a in O,b in V(B)} A_ab,
    E(O)=Σ_{a in O} e_a.

**Theorem 6 (capacity form of the exact seed constraints).**

    W_B(O)=2 I(O),
    U_B(O)=E(O)+2 V(O).

Thus (B,O) is a seed for odd slope c iff

    I(O)=(H(B)−1)/2,
    E(O)+2 V(O)=c·2^j  for some integer j≥0.

The resulting physical ray is N=1+c·2^{j+t}, t≥1. These are necessary and sufficient constraints, not merely necessary parity or divisibility conditions.

Proof. The inherited two-path-cover extension polynomial for a newborn beating x is

    H(B+v_x)=H(B)+Σ_j(s_j−e_j+Σ_i K_ij)x_j
                       −Σ_{i<j}(K_ij+K_ji)x_i x_j.

Compare its unique multilinear coefficients with Theorem 1:

    K_ij+K_ji=2(A_ij+A_ji),
    s_j−e_j+Σ_i K_ij=2Σ_b A_jb.

Summing the first identity over selected unordered pairs gives W_B(O)=2I(O). Since U_B(O)=Σ_{j∈O}(s_j+Σ_i K_ij), the second gives U_B(O)=E(O)+2V(O). Substitution into the inherited tail formula proves the equivalence. ∎

This exposes two different budgets that must be solved **simultaneously**: an exact induced-edge weight budget and an endpoint-corrected row-volume budget. Neither C_k=1 nor endpoint-count universality solves these coupled equations.

**Exact positive constraints.** Let r=(H(B)−1)/2 and d_ab=A_ab+A_ba for unordered pairs. Every d_ab≥1 since one of a→b,b→a is an original arc. Therefore:

- Any seed has binom(|O|,2)≤r.
- If d_ab>r, at most one of a,b may be selected. These are hard incompatibility edges, unlike the heavy-cut implication arcs of Section 3.
- The first budget is exactly Σ_{unordered {a,b}⊆O}d_ab=r; the second uses the explicitly attached vertex weights e_a+2Σ_b A_ab.
- Consequently seed search is an exact weighted induced-subgraph budget problem with a vertex-weight odd-part target. Enumerating the independent sets of the incompatibility graph, then enforcing both exact budgets, is complete for the supplied base B.

**Recovered inherited corollary, not a new obstruction:** HYP-9658 already records that transitive bases yield only slope 1. In the present capacity language, for a transitive base H(B)=1, the first budget forces |O|≤1. A singleton's U is a power of two: if a has m successors in the transitive order, its endpoint term is 1 when m=0 and 0 otherwise, and its doubled outgoing capacity sum is 2^m when m≥1. Hence transitive bases supply only slope c=1. This is an all-order seed obstruction, not a bounded-census observation. Nontransitive bases are indispensable for every c>1.

**Remaining gap.** No theorem here constructs bases and selected sets satisfying both budgets for every odd c. Universal seeds also need treatment of finitely many below-start values on each line to imply full count coverage. The reduction now identifies precisely what a Zeckendorf/finite-state/natural-recursion construction would have to preserve: internal selected capacity and the endpoint-adjusted volume, not only total H or a gcd.

### Recorded replay

Both runs completed successfully and their JSON outputs are byte-identical: 1,099 labelled base tournaments, 33,866 births, all selected subsets checked against independently enumerated ordered two-path covers. Only H=1 and H=3 parents admit a +2 birth in this finite universe, in agreement with the all-order Theorem 5. The selected H=5 hostile has extension spectrum {5,9,11,13,15}, gcd one, and minimum positive half-gap two. The actual all-order proofs do not depend on these census boundaries.
