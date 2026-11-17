# A two-run operator that re-arms itself at every order

**Proved.** Every nonempty tournament can be wrapped to give a normalized ready root, and that wrapper can be applied repeatedly to its own output, with no third visit to the newly inserted module. Its exact action on all ordered path-color counts is a polynomial shift operator. Starting at a literal ready count 511, an unbounded orbit remains on 7 times odd numbers coprime to 3, with **all six eligible quotient phases modulo 18** on one exact period-12 cycle. The quotients need not be prime. This is closure of a physical readiness mechanism, not arbitrary absolute placement.

## 1. A genuine re-arming construction

For any nonempty tournament B, take five fixed vertices w,a,b,r,s, with rows (18,20,17,7,8) in that order. Every member of B loses to b and beats w,a,r,s. Equivalently, after removing r the deletion is

    C3[B→w, a, b] → s,   cyclic direction (B→w)→a→b→(B→w).

Here C3[...,...,...] means that the three named blocks replace the corresponding cycle vertices. The root r beats exactly w,a,b; s beats r and loses to all its other old vertices. Call the resulting tournament X(B). The only exterior vertex entering B is b, so a whole Hamiltonian path uses at most two B runs (an initial run, then possibly one return). This is a property of the latest wrapper interface; deeper interfaces in a nested read may need additional cover coordinates.

Write u=H(B), R=P_B(2), the number of ordered two-slot path covers when either slot can be empty. Then the exact rooted profile of X(B) is

    H = 3u+14R,       g = u+2R,
    E = 2u+8R,       O = E+g,

so H−2E+g=0 and O−E=g. This holds for odd and even B, not just growing transitive twins. In particular the operation *creates* a fresh ready root; it does not require readiness of the input. “Ready” means the two established profile equalities, not that all three ternary digits already exist.

### Direct proof of the profile

Put M=B→w. Its ordinary count is u and its two-slot count is 2R. In C3[M,a,b], paths that visit M once contribute 3u; those that revisit it contribute P_M(2)−2u, corresponding to the two nonempty M runs on either side of a→b. Therefore the deletion count is u+2R.

Use the [positive readiness test](../ready_core_growth_and_mod18_20261117/README.md): K counts weighted even paths of the deletion with both endpoints in U={w,a,b}. The weighted optional subsets of a nonempty B have equal even and odd totals R/2 by odd-pack cancellation. For the unordered contact pairs {a,b}, {w,a}, {w,b} the contributions are respectively
- R: a→b, plus b→(nonempty even M path)→a;
- u+R/2: w→a, plus a→b→(odd B path)→w;
- R/2: b→(optional even B path)→w.

No other path has both endpoints in U, since w is the last vertex of its block. The guard s cannot occur on such a path. Thus K=g exactly, not just at leading scale.

The guard s makes O−E=g: a rooted non-singleton path without s can either append it or leave it to its complement, with equal weight and opposite parity. The root alone contributes g to O. It follows that E is the unsigned weighted start sum in the unguarded deletion from the three contacts. The totals starting at a,b,w are u+2R, 4R−u, and 2u+2R respectively. For instance, from a there is the singleton a or a→b with any optional M path; from b a nonempty M path may terminate in M or go on to a; from w one may stop, take a, or take a,b and any optional B path. In each case the complement weight gives the stated optional totals. Hence E=2u+8R. Combine K=g, normalization and H−2E+g=2(O−E−K) to get H.

As a contrast, insert a new ordinal block of count u anywhere in a retained-root deletion while leaving all old arcs intact. The old K contributions scale by u, and if there are any new root outcontacts as well as old ones an old/new pair contributes an additional positive term. Meanwhile g simply scales by u. A normalized ready root cannot stay ready under that active ordinal insertion. The nonordinal cycle around B and the fresh root in X are doing essential work.

## 2. The all-order path-color recursion

For integer q≥0, P_B(q) counts q ordered path slots, all allowed empty, partitioning the vertices of B. Its polynomial extension is given by the [odd-cycle packing identity](https://arxiv.org/abs/2307.05569): each pack with c odd directed cycles and s uncovered vertices contributes 2^c q^(c+s). In particular P_B(1)=H(B), P_B(0)=0 for B nonempty, and P_B(2)=R.

The exact operator is

    P_X(B)(q) = (q^5+2q^3) P_B(q)
              + (2q^4+3q^3+6q^2+3q) P_B(q+1)
              + (−2q^4+3q^3−6q^2+3q) P_B(q−1).             (1)

**Proof.** Packs that do not cross the module boundary contribute P_E(q)P_B(q), where E is the fixed exterior. A mixed cycle must use the unique ingress b, so at most one such cycle occurs. Its B segment is a nonempty path. Its exterior segment ends at b, has length at least two vertices and can start at any other exterior vertex. The complement exterior read weights those fixed paths. Direct five-vertex enumeration gives
- P_E(q)=q^5+8q^3+6q;
- weighted exterior paths ending at b of even size: A(q)=2q^3+6q;
- those of odd size ≥3: D(q)=3q²+3.

For a distinguished bonus slot in B, the sum over its contents S of H(B[S])P_(B−S)(q) equals P_B(q+1). With sign (−1)^|S| the sum instead equals P_B(q−1): each odd pack block or leftover singleton independently has q ordinary choices and one negative choice in the distinguished slot. Separate even and odd S, subtract the empty S from the even term, and multiply by 2q for the mixed odd cycle. The exterior even term gives q A(q)[P_B(q+1)−P_B(q−1)], and the exterior odd term gives q D(q)[P_B(q+1)+P_B(q−1)−2P_B(q)]. This is (1).

At q=1 this recovers the graph count formula, and at q=2 it gives 48P_B(2)+86P_B(3)−26H(B). Thus the recurrence does not close on just H and R. In fact for starting B=TT_m (so P_B(q)=q^m), the function m↦H(X^j(TT_m)) has leading exponential base j+1 with positive coefficient. Only the all-upward choice of color shift j times supplies that largest base. These functions for different j are linearly independent, so no fixed finite-dimensional linear profile of those transitive inputs recovers every iterated read. The whole path-color polynomial retains what the two-scalar profile forgets.

## 3. An unbounded ready orbit in the 7-channel

An ordinary near-sink move N at the retained root preserves readiness and sends H+g to 2(H+g). Starting with B=TT₃, X gives (H,g,E,O)=(115,17,66,83); two N's give (511,17,264,281), where 511=7·73. From any nonempty input of count h divisible by 7 but not by 3, X preserves these two properties. The first is immediate from H'=3h+14R. For the second, P_B(2)≡2^|B| H(B) mod 3, since the exponent c+s in every pack term has the parity of |B|. Hence H'≡(−1)^(|B|+1)h mod 3. All counts remain odd.

The seed 511 has path-color state modulo 9 (P(0),…,P(8))=(0,7,7,0,7,7,0,7,7), verified by literal subset convolution. Reduce (1) modulo 9, with q±1 interpreted cyclically. This is an exact linear finite-state recurrence regardless of the growing graph order. The seed state is on a cycle of length 12. Since H is always divisible by 7, multiply H mod 9 by 7^−1=4 and take the odd lift mod 18. The resulting quotient phase cycle for H/7 is

    1, 17, 17, 7, 1, 5, 17, 1, 1, 11, 17, 13.

Thus all six prime-eligible phases recur infinitely, on one physically re-arming orbit, not on unrelated test cores. Its first ordinary counts are 511, 357581, 532719551. The orbit still is sparse: n increases by five at each step, and P_B(2)≥(n+1)H(B) by cuts of full paths. This bounds growth below by a factorial scale; the n! bound gives the other direction. It has Θ(log X/loglog X) entries through X. Moreover N^6 preserves H modulo 63 and readiness, so arbitrary choices of X or N^6 also stay in the odd 7-channel coprime to 3, but branching that way is not claimed to give sufficient address density.

## 4. The next exact inverse this gives

After wrapping, t near-sinks give

    H_t = (2^(t+2)−1)u + 2(2^(t+3)−1)R.

In particular for prescribed 7p use t=2 and seek a *physical pair* (u,R) satisfying

    15u+62R=7p,  u≡17p (mod 62),  u≤7p/139.

The bound follows from R≥2u. This gives a strictly smaller required Hamiltonian count u and a fully determined two-cover count R once u is selected. The TT₃ example pays u=1,R=8 for p=73. At t=0 the inverse stays within the 7-channel itself, p=3(H(B)/7)+2P_B(2). Realizing a smaller count u alone does not ensure its required R. This is now an explicit marked two-cover selection problem, not a lack of a recurrent ready operation.

No universal selector, arbitrary fresh-digit closure, or infinitude of prime quotients is asserted. The progression from the [growing-deletion work](../growing_deletion_lanes_20261117/README.md) is that the growing output is now itself ready and can serve as the next module, without maintaining an exact +2 companion at each step.

Run python replay.py and python modular.py. The first compares every rooted profile with literal subset DP on all 75 labeled nonempty inputs of order at most four and recomputes the exterior polynomials. The second replays the mod-9 orbit, checks the 511 seed and verifies the shifted operator against independent q-slot subset calculations in small literal examples.
