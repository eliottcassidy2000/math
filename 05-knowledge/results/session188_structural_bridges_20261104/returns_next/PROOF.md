# Simultaneous module reversals: disjoint detection and overlap growth

2026-11-04. **PROVED:** arbitrary families of pairwise disjoint nontrivial modules cannot cancel in the full named two-port response; exactly two properly overlapping original modules strictly increase the scalar Hamilton count; a nested pair with one uniform cross orientation has a strict actual two-port separator. **OPEN here:** arbitrary nested pairs with mixed complement incidence and general interacting families. No cubic truncation is used for an ambient comparison.

## Conventions and incoming dependency

Vertices retain physical names. A module has uniform adjacency to every exterior vertex. Reversing a selected subset means reversing each internal arc. For simultaneous selections, an arc is reversed iff it belongs to an odd number of selected subsets. For overlapping selections these are modules of the original tournament; this definition does not assume the second subset remains a module after the first reversal.

F_T(a,c) is the Hamilton count after adding x->y, with a_u=[u->x] and c_u=[u->y]. It is the complete actual response, including every cubic and quartic term.

The incoming [76b3b1a parity theorem, E3/E4](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md) supplies:

1. Every named tournament M of order at least2 has a core query with F_M-F_(M^op)=2 modulo4.
2. Lift that query through a module by giving x and y the module's uniform exterior incidences. For a quotient Q, the actual ambient difference is congruent to H(Q) times the core difference modulo4, hence to2.

The second statement is a congruence, not exact integer factorization. The incoming theorem is used only for one isolated module. The additional arguments below explain why other reversals cannot cancel that probe in the stated scopes.

Let q_k(A) count ordered spanning covers of A by k nonempty directed paths. Reversal of every path, without changing the ordered slot labels, gives

    q_k(A^op)=q_k(A) for every k.

A module substitution satisfies H(Q[A])=sum_k b_k(Q)q_k(A), with b_k>=0 and b_1=H(Q). This is the exact occupied-run expansion. It is used over the integers, without dropping high run counts.

## D1. The full disjoint reversal cube is detected modulo four

Let M_1,...,M_r be pairwise disjoint modules of T, each with at least2 names. For S subset {1,...,r}, let T_S reverse precisely the selected modules. Then the 2^r complete named response functions F_(T_S) are pairwise distinct, already modulo4.

Proof. Compare two different states and select an index i at which their choices differ. The selected M_i is still a module in either state. Give the ports its E3 core separating query, and give both ports the uniform exterior incidences of M_i.

Every other selected M_j remains a module in the augmented graph. Indeed its old exterior incidences were uniform, and its incidence to either new port is uniform because M_i and M_j are disjoint modules. Reversing M_j therefore leaves the actual integer query count unchanged: its entire ordered path-cover profile is converse-invariant. This holds while any subset of the other disjoint modules has already been reversed.

Undo all these other differing reversals, keeping the query fixed. The difference of actual query counts is now exactly the difference for the single M_i reversal. The inherited ambient congruence makes that difference2 modulo4. Thus it was nonzero before the other changes were undone. QED.

This is an isolation argument on actual attachments. It neither asserts that the ambient contrast is cubic nor that the individual global polynomials add without interactions. Quartic terms may remain; equality of the full responses is contradicted by one actual query.

## D2. Ordered quotients minimize every ordered path-cover count

Let A_1,...,A_s be nonempty named tournaments, and let L be their ordered sum in some fixed transitive order. Let R be any tournament on the block labels. In R[A_1,...,A_s], optionally reverse any of the individual blocks internally. Then, for every k,

    q_k(R[blocks]) >= q_k(L).

Moreover

    H(R[blocks]) >= H(R) product_i H(A_i).

Thus H(R[blocks])>H(L) whenever R is nontransitive; since H(R)>=3 in that case, H(R[blocks])>=3H(L). If R is transitive, every q_k agrees with q_k(L), regardless of its named block order.

Proof. In each slot of an ordered cover of L, the occupied blocks occur at most once and in the original block order. Reverse the internal path segment in a block that was conversed. Choose once and for all a directed Hamilton order of every nonempty induced quotient R[J], which exists for a tournament. Rearrange the slot's block segments into that order. This is a valid target path. Repeat independently in every labelled slot.

The map is injective: physical block labels identify every segment, and the original path is recovered by undoing the chosen internal reversals and restoring the fixed original block order. Slot labels are unchanged. Each target slot in the image has one run per occupied block, so no ambiguity comes from possible target paths with repeated block visits.

For k=1, every block is occupied. Each of the H(R) quotient Hamilton orders and each internal Hamilton path choice gives a distinct target Hamilton path, proving the lower bound. If R is transitive, every target path also has one run per block, and restoring the original transitive order is an inverse on all covers. QED.

No claim that arbitrary strongly connected quotients are interchangeable follows. The source quotient must be transitive for this monotone injection.

## D3. Exactly two properly overlapping original modules strictly increase H

Let A and B be modules of T with all three sets

    P=A\B, I=A intersection B, R=B\A

nonempty. Put U=A union B. Let T' reverse A and B, with arcs in I reversed twice and therefore unchanged. Then

    H(T')-H(T) >= 2 H(Q) H(T[U]) > 0,                 (1)

where Q is the quotient obtained by collapsing U to one vertex. In particular no full named two-port collision is possible, and the obstruction is already scalar.

### Module and orientation checks

U is a module: an exterior vertex has one incidence to A and one to B, and their nonempty intersection forces these incidences to agree. I is a module, because a vertex outside I is outside A or outside B.

Choose p in P, i in I, r in R. Since p is outside B, p sees i and r identically. Since r is outside A, r sees p and i identically. Consequently the triple is p->i->r with p->r, or the reverse. Fixing one r and then one p shows the orientation is uniform across all choices. Hence, after possibly interchanging P and R, the quotient is

    P -> I -> R.

This also shows P and R are modules: vertices in the other two parts have the uniform incidences just proved, and exterior vertices see all of U uniformly.

### Reversal and injection checks

The two reversals converse P and R internally, preserve I internally, reverse the P-I and I-R cross cuts, and preserve P-R. The quotient therefore becomes the directed cycle

    I -> P -> R -> I.

Apply D2 to this three-block quotient. Every ordered q_k is nondecreasing, and q_1(T'[U])>=3q_1(T[U]). The unchanged exterior sees U uniformly in both graphs. In its exact run expansion, every coefficient b_k(Q) is nonnegative and b_1(Q)=H(Q)>0. Therefore

    H(T')-H(T)
      = sum_k b_k(Q)(q_k(T'[U])-q_k(T[U]))
      >= H(Q)*2H(T[U]).

This proves (1). The bound is sharp on three singleton blocks with no exterior: the transitive triangle becomes a directed triangle, and the count changes1->3.

Identical selections, disjoint selections, and nested selections do not satisfy the three-nonempty-parts hypothesis and are not silently included. In particular repeating the identical module twice produces no change at all.

## D4. A uniformly split nested pair has a strict two-port separator

Suppose A is a proper nonempty module inside a module B, R=B\A is nonempty, and all cross arcs between A and R have one direction. Let T' reverse both A and B. If A->R, add x with A->x->R; if R->A, add x between R and A in that original order. Give x the exterior incidences of B. Add y as a global sink, with x->y.

Before the reversals, B union {x} is a module with a transitive three-block quotient. Afterwards it has a directed three-cycle quotient; A retains its internal graph, R is conversed, and x is a singleton. D2-D3's profile comparison therefore applies to this augmented module. If Q collapses B in the original T, the actual query difference obeys

    F_(T')(a,all_ones)-F_T(a,all_ones)
         >= 2 H(Q) H(A) H(R) > 0.                    (2)

The second port y is forced last, so it does not change the one-port Hamilton count. This is an actual named two-port query, not a coefficient truncation.

A module A need not have the same incidence to every vertex of B\A: the complement can contain both in-neighbors and out-neighbors. D4 assumes one uniform cross direction and does not prove the general nested case.

## Larger-family boundary

For a transitive tournament on names1,...,n, every interval is a module. Let f_ij specify any desired edge-flip pattern for i<j. Selecting interval[a,b] with parity

    sigma_ab = f_ab + f_(a-1,b) + f_(a,b+1) + f_(a-1,b+1)  modulo2,

with out-of-range f terms zero, makes edge(i,j) flip by exactly f_ij. This is two-dimensional finite-difference inversion of the sum over intervals containing both i and j.

Therefore unrestricted families of overlapping original modules can change a transitive quotient into an arbitrary quotient; they are not merely independent converse bits. D2 handles the scalar growth when a transitive quotient becomes nontransitive. It does not prove named response separation for every residual nested/order interaction in an arbitrary modular decomposition.

The current positive statements are D1-D4 in exactly their stated scopes. General full two-port reconstruction and arbitrary interacting module families remain outside this packet's conclusion.

## D5. Every tournament module is pair-generated; incidence vectors are independent

For a nonempty set S, its module hull is the intersection of all modules containing S. Intersections containing S are again modules, so the hull exists and is the least such module. Equivalently, repeatedly add an exterior vertex that sees the current set in both directions. Every such vertex belongs to every module containing the current set; the process terminates at its hull.

**Pair-generation theorem.** Every nontrivial module M of a finite tournament is the module hull of a pair of distinct vertices.

Proof. Work inside the induced tournament on M. Its modules lift to modules of the original tournament because every exterior vertex sees all of M uniformly. Every vertex lies in a maximal proper module of M, by extending its singleton; these maximal proper modules cover M.

If the maximal proper modules are pairwise disjoint, choose p and r in two different ones. If a proper module contained both, extend it to a maximal proper module D. Since distinct maximal modules are disjoint in this case, D could not intersect both chosen maximal modules. This is a contradiction.

Otherwise choose two overlapping maximal proper modules A and B. Neither contains the other. Their union is a module strictly containing A, so maximality gives A union B=M. As established in D3, their nonempty parts P=A\B, I=A intersection B, R=B\A have a transitive quotient; choose p in P and r in R.

Suppose a proper module contains p and r, and extend it to a maximal proper module D. Then D overlaps both A and B: it meets A at p and has r outside A; it cannot contain A because A is maximal proper. The same argument applies to B. Therefore A union D=B union D=M. Taking complements inside M shows R subset D and P subset D. Every i in I sees p and r in opposite directions, by the transitive three-part quotient, so it cannot lie outside the module D. Thus I subset D also, forcing D=M, a contradiction. In both cases no proper module of M contains p,r, and their hull is M. QED.

**Incidence independence.** For each nontrivial module M choose one generating pair e_M. If N is any module, then

    e_M subset N  iff  M subset N.

Consequently the matrix with rows e_M and columns N, entries `[e_M subset N]`, is the inclusion-zeta matrix of the nontrivial modules. Order the modules by increasing size; it is triangular with diagonal1. Its determinant is1 over the integers. Therefore the vectors

    1_{unordered vertex pairs contained in M}

of all distinct nontrivial modules are linearly independent over every field, in particular over F2. Their integer span has a unimodular minor.

An equivalent direct F2 argument is especially useful for simultaneous reversals. From a nonempty selected family choose an inclusion-maximal M and its generating pair e_M. If another selected N contained that pair, it would contain M, contradicting maximality unless N=M. Thus exactly one selected module flips this pair. No nonempty family of distinct nontrivial original modules can compose to the unchanged labelled tournament.

The selected family can be recovered from its edge-flip vector by descending inclusion:

    sigma_M = flip(e_M) - sum_(N strictly contains M) sigma_N.

Over F2 subtraction is addition. Repeated copies of the same module, or singleton modules that change no arc, are excluded from this statement for the necessary reasons.

**Sharp module-count bound.** There are at most binom(n,2) nontrivial modules. A transitive tournament attains the bound: its modules are exactly its intervals, and each pair generates the interval between its endpoints. Conversely, if a directed triangle exists, each pair in that triangle forces the third vertex into its hull. Its three pairs therefore have the same hull. Because pair hulls surject onto all nontrivial modules, a nontransitive tournament has at most binom(n,2)-2 nontrivial modules. In particular equality characterizes transitive tournaments.

This rules out a trivial exact collision from an unchanged edge set. It does not by itself prove injectivity of the full named two-port response for the distinct tournaments created by general interacting families.

## A parity hostile for quotient-only lifting

One must not extend D1 by multiplying an E3 quotient separator through several inflated vertices. Let Q=C3 and put a transitive two-vertex block at each quotient vertex. Compare Q[TT2,TT2,TT2] with Q^op[TT2,TT2,TT2]. All64 queries whose attachment bits are constant on each block have contrast0 modulo4, although their integer contrasts are not all zero (maximum absolute contrast64). The analogous TT3-block example has the same modular invisibility on all64 coarse queries and maximum absolute integer contrast2508.

These are not full named response collisions. They show that the coarse uniform-port parity observer can lose every quotient separator. D1 avoids this loss by putting arbitrary named queries inside one actual module and making only the exterior homogeneous. General nested comparisons cannot reverse that argument without an additional theorem.
