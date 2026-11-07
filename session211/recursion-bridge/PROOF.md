# Two vertex recursions, one apex flip, and a folded permutohedron

**Status: PROVED elementary statements, independently audited; FINITE-EXACT evidence is separately identified.** No literature-priority claim. The unrestricted tournament Hamilton-path spectrum conjecture remains OPEN.

## 1. Inheritance and the variables that must not be conflated

Throughout, T is a nonempty finite tournament of order n. H(T) counts directed Hamiltonian vertex words; HC counts directed Hamiltonian cycles modulo cyclic rotation. For algebraic purposes H(empty)=1. P_T(q) counts ordered q-tuples of possibly empty directed paths partitioning the vertices. Thus P_T(1)=H(T), and

    P_T(2)=sum_(S subset V(T)) H(T[S]) H(T[V(T) minus S]).

The close inherited mechanisms, at repository revision `214586af6c8803d8b21959951e8f151e64bff36f`, are:

- THM-291, `01-canon/theorems/THM-291-mode-b-multilinear-recursion.md`: the **order** n to n+2 staircase extension. With zero boundary wiring, add a source s and sink t, obtaining Sigma(T)=s -> T -> t, with s -> t. H(Sigma(T))=H(T). Each boundary arm has n-1 free tiles, and the full boundary has 2n-1; the older n-2 arm count is corrected in that file.
- THM-447, `01-canon/theorems/THM-447-skew-sylvester-doubling.md`: the **order** n to 2n skew-Sylvester double D(T). Its dominance matrix is [[M,M+I],[M-I,-M]], where M is T's dominance matrix. This is not H to 2H. THM-450, `01-canon/theorems/THM-450-doubling-family-trichotomy.md`, classifies only its explicit span{M,I} block ansatz, not every imaginable doubling.
- Session210, `session210/graded-spectrum-gates/PROOF.md`: Gamma(T) has t -> s -> T -> t, with HC(Gamma(T))=H(T), H(Gamma(T))=H(T)+P_T(2), and cover operator Lf(q)=q^2 f(q)+q(f(q+1)-f(q-1)).
- Session210, `session210/coherent-permutohedron/PROOF.md`: tournament Hamilton paths are sinks of the parallel-coherent permutohedron orientation, with edge moves decreasing the all-pairs feedback energy by one.
- THM-002, `01-canon/theorems/THM-002-ocf.md` (pinned in SOURCES.json). The odd-cycle collection formula implies H(T)=1 iff T is transitive, and H(T)=3 iff there is exactly one directed odd cycle. The proof below uses only these elementary low-count consequences.
- THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, proves the omissions 7 and 21. Session195's `session195/combined-growth/PROOF.md` explains why a genuinely closed growth grammar needs marked endpoint/cover/deletion coordinates, not just H. Its finite support results are not eventual coverage.

The corrected near miss is to read the order recursions as affine count operations. Literal H to 2H cannot be a tournament count map, because tournament Hamilton-path counts are odd. A universal H to H+2 constructor is also impossible: H=5 is realized while H=7 is excluded by the inherited theorem. The canonical hostile is the pair C3 -> 1 and 1 -> C3: both have n=4 and H=3, but their skew doubles have H=333 and H=189. These numerical values are inherited in THM-450 and independently replayed here. The least-used sidecar is the endpoint distribution resolved by the two sheets.

Concept board: staircase boundary / gate cycles / coherent permutohedron / reflected sheet words / endpoint norms / count spectrum. Anchor: distinguish and connect the actual operations. Niche: sharp low gate strata. Wildcard: interpret the reflected doubling map inside a cubical face of the permutohedron.

## 2. The apex flip is an exact bridge, not a count +2 step

**Proposition RB1 (apex flip and sink creation).** Sigma(T) and Gamma(T) differ in exactly the arc on their two new vertices. Reversing s -> t into t -> s leaves all H(T) old Hamilton words s,P,t valid and creates exactly P_T(2) new words. No old Hamilton word is destroyed. In the coherent permutohedron, reversing this one parallel edge class therefore creates P_T(2) sinks and destroys none.

**Proof.** Any Sigma(T) Hamilton word must start with s and end with t, hence is s,P,t for a Hamilton word P of T. Since n>=1, its consecutive arcs never use s -> t, so it survives the reversal. Any genuinely new Gamma(T) Hamilton word must use the only changed arc t -> s. Cutting at this consecutive pair gives

    P_left, t, s, P_right,

where the two possibly empty directed paths partition T. Conversely every such ordered two-path cover gives a valid word, because T -> t and s -> T. This correspondence is bijective, including the two cases with an empty path. The permutohedron statement follows from its exact sink/Hamilton-word correspondence. QED.

The scalar identity is inherited; the point is its precise operation-level realization between the old staircase boundary and the new gate carrier. It is a special unsigned situation of the general arc-flip formula, not a new general flip law. The nonempty-core hypothesis matters: at n=0 one old word is lost and one gained.

**Connection contract.** Source: a fixed n-tournament with the source/sink padding. Target: a strong gated tournament, and its coherent permutohedron. Map: flip the unique apex arc. Preserved: old induced tournament and its H embedded Hamilton words. Added: ordered two-path covers as new sinks. Lost by retaining only H: P_T(2), hence the response. Required sidecar: the two-bin cover count (or the whole cover polynomial for repeated closure). Cheapest decisive test: the singleton core gives H=1 to H=3, while the two-vertex core gives H=1 to H=5.

## 3. The gate has sharp first and second path-count strata

**Theorem RB2 (sharp two-level gate bound).** For fixed n>=1,

    H(Gamma(T)) >= 1+2^n,

with equality exactly when T is transitive. For n>=3, among nontransitive T,

    H(Gamma(T)) >= 3+3*2^(n-1),

with equality exactly when T is a dominance sequence of singleton components and exactly one directed triangle component. In particular these are the exact smallest and second-smallest distinct H values among gated tournaments of order n+2. There is no second level for n=1,2.

**Proof of the minimum.** Every induced tournament has at least one Hamilton path, so each of the 2^n summands in P_T(2) is at least one. Also H(T)>=1. For a transitive T every summand is one, giving equality. Conversely equality forces H(T)=1, hence transitivity.

**Proof of the second bound.** Fix a directed triangle C in a nontransitive T. There are exactly 2^(n-2) subsets S for which all three vertices of C are contained in S or all are contained in its complement: choose the shore in two ways and place the other n-3 vertices freely. On each such partition the factor on the C shore is at least three, while the other factor is at least one. All remaining summands are at least one. Therefore

    P_T(2) >= 2^n + 2*2^(n-2) = 3*2^(n-1),

and H(T)>=3. Equality in the final H bound forces H(T)=3. By OCF there is exactly one odd cycle, necessarily a triangle. Conversely a tournament with exactly one directed triangle has that triangle as a module and no other nontrivial strong component: if an outside vertex had mixed orientations to the triangle, one of its three cyclic edges would form another triangle with it. All outside vertices therefore see the triangle uniformly, and the quotient has no triangle, hence is transitive. This is precisely a dominance sequence with one C3 and singleton components. Its H is three, and every induced subgraph has H=3 if it contains that C3, otherwise H=1. No partition puts C3 on both shores. The displayed bound is attained. QED.

The use of OCF can be confined to the standard low-count fact H=3 iff exactly one odd cycle. The partition estimate itself is a direct counting argument and does not invoke polynomial real-rootedness or any conjectural positivity.

**Triangle-resolved strengthening.** If c3(T) is the number of directed triangles, then

    H(Gamma(T)) >= 1+2^n+c3(T)*(2+2^(n-1)).

Indeed OCF gives H(T[S])>=1+2c3(T[S]). On multiplying the two shore bounds, discard only nonnegative terms. Every fixed triangle lies on one shore in 2^(n-2) partitions, so summing the linear triangle contributions gives P_T(2)>=2^n+c3(T)*2^(n-1). Add H(T)>=1+2c3(T).

**Exact small gated spectrum.** The Hamilton-path values at most 31 in the nonempty-core gate class are exactly {3,5,9,15,17,27}. For n>=5 the minimum is at least 33. At n=1,2 the cores are transitive; at n=3 the two core classes are transitive and C3. At n=4 the first two levels are 17 and 27, while every remaining core has at least two triangles and hence H(Gamma(T))>=37 by the triangle bound. Each of the six listed values is attained by the equality constructions. Thus 7 and 21 can be excluded in this local class without importing the global omissions theorem; this is not a proof of their unrestricted exclusion.

**Consequences for the spectrum question.** The gate-cycle observable faithfully transports the original graded spectrum: HC(Gamma(T))=H(T). The gate-path observable emphatically does not: at each fixed order its first two levels are exponentially separated. For example, at gate order six they are 17 and 27, so 19,21,23,25 all fail within that fixed-order gate family. This does not assert that 19,23,25 are absent from the unrestricted tournament spectrum; they are present. It is a structural reason that the canonical strong gate, used only as a path-count constructor, is too sparse to justify a scalar +2 induction.

## 4. Doubling has a folded permutohedron, not count multiplication

Write the two sheets as v and v'. Within the first sheet D(T) is T, within the second it is T^op; for distinct u,v, u -> v' iff u -> v in T, and every v -> v'.

For any permutation p=(p_1,...,p_n), set

    F(p)=(p_1,...,p_n,p_n',...,p_1').

Write E_T(p) for the number of backward all-pairs arcs in the word, not only backward adjacent pairs.

**Proposition RB3 (folded energy and squares).**

    E_D(T)(F(p)) = binom(n,2)+2 E_T(p).

F(p) is a Hamilton word of D(T) iff p is a Hamilton word of T. Each adjacent interchange in p becomes the opposite-corner move in a square face of the doubled permutohedron: swap the corresponding adjacent unprimed labels and the corresponding adjacent primed labels. Both constituent edges point in the direction of the original coherent interchange. A length-d path in the original permutohedron thus lifts to a length-2d path by either order of each pair of swaps.

**Proof.** The within-sheet contributions both equal E_T(p), since the second sheet reverses both order and arcs. All unprimed labels precede all primed labels. For each distinct pair {u,v}, exactly one of the two cross arcs points backward, regardless of p; diagonal twin arcs point forward. Thus cross energy is binom(n,2). The word is directed exactly when p is, because the second half reverses both order and arcs, while the central twin arc always points forward. For the square, the two adjacent swaps occupy disjoint positions and hence commute; their arc preferences agree after reversal of both the primed order and its tournament. Each gives the same energy decrement as the original interchange. QED.

**Actual affine section.** In label-position coordinates, the face that puts the entire first sheet before the entire second is P_(n-1) x (n+P_(n-1)). Intersect this product face with

    x_v+x_(v')=2n+1  for every v.

The section is exactly the graph x -> (x,2n+1-x) for x in P_(n-1), since (n+1)-P_(n-1)=P_(n-1). Its vertices are precisely F(p), and it is an affine permutohedron with Euclidean scale factor sqrt(2). Its edges are the square diagonals described above. For n>=2 this section is not a face of the product; at n=1 it is a point and the degenerate exception is harmless.

This recovers THM-447's reflected Hamilton-path injection and supplies its exact geometric/energy sidecar. The section is not a bijection onto all sinks of D(T). The cross-sheet adjacency sectors below already contain many additional words.

**Theorem RB4 (exact one-switch sectors and a lower bound).** Let e_v and s_v count Hamilton words of T ending and starting at v. Set h=H(T), E2=sum_v e_v^2, S2=sum_v s_v^2. Among Hamilton words of D(T), the number using the whole first sheet and then the whole second is

    B_01=(h^2+E2)/2.

The number using the whole second sheet and then the whole first is

    B_10=(h^2-S2)/2.

These sectors are disjoint, so

    H(D(T)) >= h^2+(E2-S2)/2 >= B_01 >= h(h+1)/2.

**Proof.** For first-then-second words, choose two Hamilton words of T and reverse/prime the second. If their endpoints are u,v, concatenation is legal iff u=v or u -> v. Hence the count is sum_u e_u^2+sum_(u->v)e_u e_v. Since each unordered pair contributes once, this is (h^2+E2)/2. For second-then-first, the corresponding starting vertices u,v require u -> v with u!=v, because no reverse twin arc exists. Thus the count is sum_(u->v)s_u s_v=(h^2-S2)/2. Finally E2>=sum_v e_v=h, since the e_v are nonnegative integers. QED.

The remainder R(T)=H(D(T))-B_01-B_10 counts paths with at least two sheet changes and is nonnegative. On the hostile C3 -> 1, (E2,S2)=(9,3), B_01+B_10=12 and R=321; on 1 -> C3 the norms are (3,9), B_01+B_10=6 and R=183. Thus even the endpoint-norm asymmetry accounts for only 6 of the total 144 count difference. No closure on the scalar h, nor any closure on h and the two norms, is asserted.

## 5. Equal size arithmetic does not make the recursions commute

**Proposition RB5 (physical noncommutation).** The following pairs have the same order 2n+4, but are never isomorphic:

    D(Sigma(T))  and  Sigma(Sigma(D(T))),
    D(Gamma(T))  and  Gamma(Gamma(D(T))).

**Proof.** For any tournament X of order N, a first-sheet vertex v of D(X) has outdegree 2d_X^+(v)+1, and every second-sheet vertex has outdegree N-1. For X=Sigma(T), N=n+2>=3, so every vertex of D(X) has positive outdegree: it has no sink. The doubly padded object has a sink.

For X=Gamma(T), X is strong, so d_X^+(v)>=1. All first-sheet vertices in D(X) have outdegree at least three; all second-sheet vertices have outdegree N-1>=2. Thus D(Gamma(T)) has no vertex of outdegree one, and hence no gate. Gamma(Gamma(D(T))) has its outer canonical gate. QED.

This is an exact obstruction to treating the arithmetic identity 2(n+2)=2n+4 as an operation identity. It neither disproves any proposed marked natural transformation nor says that all useful mixed recursion diagrams fail. It identifies the gate/degree sidecar such a diagram would have to change.

## 6. Scope, reproducibility, and the remaining question

`verify.py` uses exact integer subset dynamic programming, direct permutation checks, and independent explicit products of path lists. Its frozen universe is all 1,099 labelled tournaments of orders 1 through 5. It checks padding/gate counts, sharp minima/equality, folded energy/sink equivalence on all 124,469 permutation words, the endpoint sector lower bounds, and the gate noncommutation invariant. Explicit two-block enumeration, independent of the doubled subset-DP total, covers all 75 tournaments through order four. Higher n statements are proved above; finite tests are controls, not proofs of all-order claims.

Reproduce with:

    python3 session211/recursion-bridge/verify.py
    python3 -O session211/recursion-bridge/verify.py

The portfolio outcome is an operation dictionary, two sharp gate strata, a square-face doubling lift, and an endpoint-resolved count lower bound. It does not prove that every odd value except 7 and 21 is realized. The productive remaining problem is a growth interface that controls the newly created sinks while returning to its own marked input class; Session195's four-coordinate grammar is a concrete inherited example of that requirement. Merely composing the two order recursions does not supply the needed count-surjectivity theorem.
