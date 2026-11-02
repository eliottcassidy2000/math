# Exact weighted stability at every order n>=748

**Historical checkpoint: valid but SUPERSEDED as the current frontier by [n>=514](STABILITY_514.md).** Its original remaining-gap statements below record the pre-fetch state. Concurrent Session139 independently proved this threshold and then526.

**Status: PROVED; same-session proof review, exact rational verification and a separate combinatorial implementation. No independent-agent audit has been performed.** This file strengthens the sufficient threshold1500; it does not assert exact optimality at every smaller order. No external novelty claim.

Let A be a weighted tournament: diagonal zero, A_ij+A_ji=1, entries in[0,1]. Write H for its weighted directed Hamiltonian-path sum, V=sum_(i<j)q_ij with q_ij=A_ij A_ji, and Delta=H-det(I+A)-|det A|. Put R=floor(n/2) and

    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

**Theorem.** For every n>=748, Delta>=kappa_n V. The coefficient is optimal. Equality at V>0 consists exactly of an ordinal sum of R fair pairs and, at odd order, one singleton. At V=0 equality is the inherited deterministic transitive case.

## Inherited proved inputs

All inputs are from Session138's audited [asymptotic proof](../session138_stability_and_actual_contexts_20261102/asymptotic_stability/PROOF.md), [matching-support theorem](../session138_stability_and_actual_contexts_20261102/stability/MULTISTAR_STABILITY.md), [triangle variance lemma](../session138_stability_and_actual_contexts_20261102/stability/TRIANGLE_VARIANCE.md), and [H=1 classification](../session138_stability_and_actual_contexts_20261102/stability/H1_STABILITY.md). Let B=E Delta(T) under independent binary rounding and let W be the sum q_e q_f over unordered incident edge pairs. They prove:

1. B>=W and W>=2V^2/n-V/4.
2. Delta>=B+(1-exp(-4V))/4.
3. For any retained matching of m edges and weight M, Delta>=B+a_m M when n>2m, where a_m=(4/m)(1-(3/4)^m). In all cases the weaker Delta>=B+kappa_n M holds.
4. The descending-weight greedy matching has selected weight V-D. Its discarded edges are assigned to the selected edge that removes them. A group with selected weight w and discarded mass S contributes at least wS and at least S^2/4 to W. Contributions of different groups are disjoint. Hence W>=D^2/(2n).

The new argument retains the same B simultaneously in inputs2 and3, and later adds the disjoint heavy/light group costs before estimating them.

## 1. Coupled mass bound for a candidate-beating configuration

Suppose V>0 and Delta<=kappa V, with kappa=kappa_n. For n>=748,

    kappa <= 4/R <= 2/187,   kappa <=8/(n-1),   kappa R<=4.

The elementary exponential bound exp(4)>50, obtained from its Taylor sum through degree7, gives c0=(1-exp(-4))/4>49/200. If V<=1 then input2 gives Delta/V>=c0, impossible. Thus V>1 and kappa V>49/200, which also implies V>2. Since exp(8)>2500, input2 now yields

    kappa V > B+c,       c=2499/10000.                 (1)

Keep the greedy matching in input3. Then

    B<=kappa D,   D^2/(2n)<=W<=B.                     (2)

Combining(1) with the second inequality in(2) before discarding B gives

    V-D > c/kappa + D^2/(2n kappa)-D
        = c/kappa-n kappa/2+(D-n kappa)^2/(2n kappa)
        >= c/kappa-n kappa/2.                        (3)

This replaces the earlier separate loss bound D<=2n kappa. Its subtraction cost is n kappa/2, not2n kappa. Separately, W<=kappa V implies

    V<=n/8+n kappa/2<=n/8+4n/(n-1)<=n/7.             (4)

The last comparison needs only n-1>=224.

## 2. Count light selected edges without losing parity

Set eta=1/55. Suppose L>=1 selected greedy edges have weight at most eta, and delete all of them. The retained matching has size r_h<=R-L and mass M_h. By(3),

    M_h >=V-D-eta L
        >c/kappa-n kappa/2-eta n/2
        >=c(n-1)/8-4n/(n-1)-n/110
        >(n+2)/60.                                  (5)

The last strict inequality holds for every n>=748: the difference at748 is1736491/59760000, and its derivative is14467/2640000+4/(n-1)^2>0. Since every weight is at most1/4, r_h>=4M_h>(n+2)/15. In particular the retained matching is nonempty and n>2r_h.

Let epsilon=(3/4)^((n+2)/15). Then (3/4)^r_h<epsilon and

    epsilon <=7/(16000n),  n>=748.                   (6)

Indeed n(3/4)^((n+2)/15) decreases for n>=748, since log(4/3)>1/4. At748, the exponent is50 and the endpoint is the exact integer inequality16000*748*3^50<7*4^50.

The retained matching inequality and Delta<=kappa V give, after multiplying by R-L>=r_h,

    kappa V(R-L) >=4(1-epsilon)(V-D-eta L).

Using kappa R<=4 and rearranging yields

    L[kappa V-4eta(1-epsilon)] <=4epsilon V+4D.       (7)

By(4),(6), delta=4epsilon V<=1/4000. Since L>=1, absorb delta into delta L. By(1),

    L[c-4eta-1/4000]<4D.

Here c-4eta-1/4000>32/181, with exact surplus5063/39820000. Therefore

    L<(181/8)D,   D>0.                              (8)

## 3. One energy budget for the two group types

Split D=D_h+D_l according to whether the owning selected edge is heavy or light. The assigned groups have disjoint wedge contributions. Therefore

    B>=W>=eta D_h + sum_light S_i^2/4
        >=eta D_h + D_l^2/(4L)
        >=eta D_h + 2 D_l^2/(181D).                 (9)

The final weak inequality is sufficient; it follows from(8), including D_l=0. Set x=D_h/D in[0,1] and divide by D>0. Completing the square gives

    B/D >=eta x+(2/181)(1-x)^2
         >=eta-(181/8)eta^2
          =259/24200
          >2/187>=kappa.                            (10)

The strict rational surplus is3/411400. This contradicts B<=kappa D in(2). Notice that the heavy and light terms are added before using B; two separate inequalities bounded by the same B would lose this improvement.

Thus L=0: every selected positive weight exceeds eta>kappa. If D>0, at least one assigned group has S>0, so B>=W>=sum wS>eta D>kappa D, again impossible. Hence D=0 and B=0. The uncertainty graph is a matching, and the inherited matching theorem proves Delta>=kappa V.

## 4. Optimality, equality and remaining gap

Ordinal maximal fair pairs attain kappa, as calculated in the inherited H=1 minimization. The preceding argument applies to every equality configuration with V>0 and forces D=B=0. At n>=4 every supported binary completion with zero gap is transitive, so H=1. The inherited complete H=1 classification and strict scalar minimization give precisely the stated fair-pair equality locus.

The argument proves the exact coefficient for every n>=748, hence newly covers748..1499. Together with inherited orders3..9, the unresolved interval becomes10..747. It does not settle the unrestricted Hamiltonian-path value spectrum.

## 5. The limit of this two-bin budget

An idealized version with c=1/4 and epsilon=0 has L/D<=4/(c-4eta). The best coefficient obtainable from exactly this additive heavy/light estimate is

    max_eta [eta-4eta^2/(1/4-4eta)] = (3-2sqrt(2))/16.

The maximizing eta is(1-1/sqrt(2))/16. This is a limitation of this particular scalar relaxation, not a weighted-tournament counterexample. The near coincidence at R=373 is exact:

    1055^2-8*373^2=-7.

It implies4/373>(3-2sqrt(2))/16, whereas4/374 lies below it. The actual kappa_746 also exceeds this ceiling; the verifier bounds its exponentially small correction by the explicit positive gap. Thus another improvement below748 must use more than perfect constants in these same two bins. Retaining detailed selected weights, directed triangle information or a different local-minimum argument are genuine new routes.
