# Internal quartics recover a transitive module up to converse

Session 189. **Proved by an exact ordered-run count, subject to independent review.** A fixed physical transitive module of order at least five has its internal named order determined, up to whole converse, by its shared-name quartics in every fixed ambient. This gives a new necessary restriction on full-response collisions arising from original-module XOR: at any highest changed prime node, a transitive child of size at least five can only remain unchanged or be wholly conversed. Arbitrary internal permutations in that child cannot survive.

This does not exclude the resulting reversal cube, prove unrestricted response injectivity, or exhibit a collision. No literature-priority claim is made.

## 1. Physical observer and frozen input

Let F_T(a,c) be the number of Hamilton paths after adjoining x->y, with a_v=[v->x] and c_v=[v->y]. Names remain physical. We use the shared-name quartic identity from the reviewed Session188 normal form:

    [a_i a_j c_i c_k] F_T = C_jk(T-i)+C_kj(T-i).       (1)

Here i,j,k are distinct, and C_jk counts ordered spanning covers by two nonempty directed paths whose first path ends at j and second starts at k. The two path supports are disjoint and cover the complete indicated tournament; the gap arc j->k is not tested.

The source is [Session188 response/PROOF.md](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/response/PROOF.md), Section2, pinned at commit `8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7`. Its complete named cubic and quartic conventions are retained. We do not identify physical names or replace a joint cover by independent endpoint marginals.

## 2. The transitive two-cover distance polynomial

Give the n names of a transitive tournament their physical order. Fix distinct names j,k at order distance d, with 1<=d<=n-1. In every directed path-cover weight replace each potential internal edge p->q by

    [p precedes q] + t.

The t term means that an exterior excursion is to replace the transition; it is allowed even when the direct transitive arc is absent. Let f_(n,d)(t) be the total weight of two spanning paths, summed over the endpoint conditions j|k and k|j. Then

    f_(n,d)(t) = sum_(r=0)^(n-2) f_(n,d,r) t^r,

    f_(n,d,r) = (r+1) sum_(s=0)^r (-1)^(r-s) binom(r,s)
       (s+1)^(n-d-1) [s^(d-1)+(s+2)^(d-1)].             (2)

The convention 0^0=1 is used.

**Proof.** Expand a contribution with exactly r transitions using t. Cutting those transitions, as well as the distinguished two-cover boundary, produces r+2 increasing nonempty runs. The run ending at j and the run starting at k have distinguished adjacent slots; the r other runs occupy the remaining ordered slots. There are r+1 possible positions for the distinguished boundary. Fix the ordinary labels in their order along the remaining slots; their color assignments already distinguish these slots. Thus the factor is r+1, with no extra r! factor.

Color the old names by the r+2 run labels J,K,1,...,r. The names j and k belong to J and K respectively, so those two colors are automatically nonempty. Inclusion-exclusion requires every other color to be used. If j precedes k, a name strictly between them can use only one of the s surviving ordinary colors; names outside their interval have s+1 allowed colors. If k precedes j, an interior name can also use J and K, giving s+2 choices. The endpoints themselves are already fixed. This proves (2), including its independence from the position of the interval in the total order. Each coloring gives a unique increasing run in each color and a unique reconstructed ordered two-cover with its selected t transitions. QED.

The first coefficient is

    f_(n,1,0)=f_(n,2,0)=2,
    f_(n,d,0)=2^(d-1) for d>=2.                          (3)

The following stronger comparison is essential:

    f_(n,1)=f_(n,2),
    f_(n,d+1)-f_(n,d) is coefficientwise nonnegative
       and has strictly positive constant term, 2<=d<=n-2.  (4)

To prove (4), subtract (2) at consecutive distances. The coefficient difference, divided by r+1, is the r-th forward difference at zero of

    g(s)=(s+1)^(n-d-2) [(s+2)^(d-1)-s^(d-1)].           (5)

For d>=2 this polynomial has nonnegative ordinary coefficients and positive constant term. Every forward difference at zero of a polynomial with nonnegative coefficients is nonnegative: on a monomial s^q it equals r! S(q,r), the number of surjections onto r labeled colors. At r=0 the difference equals 2^(d-1)>0. When d=1 the difference is identically zero. This proves the full coefficientwise statement rather than only scalar monotonicity.

## 3. Exact transfer through every fixed exterior

Let A be a transitive module of T, of order m, and put R=T-A. The exterior may be arbitrary, including strongly connected ambients. Work over the squarefree algebra on the physical names of R. For a directed exterior path use the product of the variables of its vertices; products with repeated support vanish.

Define the nonnegative exterior path polynomials:

- alpha: the empty prefix, or an exterior path whose last vertex points into A;
- beta: the empty suffix, or an exterior path whose first vertex receives an arc from A;
- tau: an exterior path whose first vertex receives from A and whose last vertex points into A.

The prefix and suffix may be empty; each tau excursion is nonempty. Each exterior vertex separately has its fixed common incidence to A. A path may return to A several times. These definitions impose all necessary boundary arcs and no others.

For r>=0 put

    theta_r = [R] alpha beta tau^r >=0.                 (6)

All extraction is joint in the same squarefree algebra. In particular this is not multiplication of unrelated endpoint counts. The constant weight has the exact interpretation

    theta_0 = H(R+z)>0,                                (7)

where z is one singleton having A's exterior incidences: cut a Hamilton path at z to obtain precisely a valid prefix and suffix.

For distinct i,j,k in A, delete i and let d_i(j,k) be the distance between j,k in the remaining transitive order on m-1 names. Maximal exterior-run decomposition in the two paths from (1) gives

    L_ijk(T):=[a_i a_j c_i c_k]F_T
      = W_(d_i(j,k)),
    W_d = sum_(r=0)^(m-3) theta_r f_(m-1,d,r).          (8)

The same alpha,beta,tau apply to every i,j,k and every reordering of A. Indeed both paths in (1) meet A; one has an external prefix and the other an external suffix, while all other exterior runs are excursions between consecutive A names within a path. The marked two-cover gap is not an excursion and is not tested. Conversely these data reconstruct each cover uniquely. No exterior-only third component can occur in a two-cover whose two specified endpoints lie in A.

For m>=5, (4), (6) and (7) now give

    W_1=W_2<W_3<...<W_(m-2).                           (9)

The strict gaps survive every exterior because theta_0 is strictly positive. This is an exact integer ambient theorem, not an inferred generic sign lift.

## 4. Reconstruction of the physical order up to reversal

For each unordered pair of names j,k in A, sum all its shared-name quartics:

    S_jk = sum_(i in A-{j,k}) L_ijk.

If j,k have distance d in the original m-name order, deleting one of the d-1 interior names shortens their distance by one; deleting one of the m-d-1 exterior names does not. Hence

    S_d=(d-1)W_(d-1)+(m-d-1)W_d.                        (10)

At d=1 the first summand is absent; at d=m-1 the second is absent. Equation (9) implies

    S_1=S_2<S_3<...<S_(m-1).                            (11)

For 2<=d<=m-2 the consecutive difference is

    (d-1)(W_d-W_(d-1))
       +(m-d-2)(W_(d+1)-W_d).

At d=2 its second coefficient is m-4>0. At larger d its first term is strictly positive. The d=1 equality follows directly from W_1=W_2.

Thus the unique pair maximizing S_jk consists of the two physical endpoints. Choose one as the first endpoint. Its S-values distinguish every distance >=3, with only distances1 and2 tied. Those two names are distinguished by their distances from the opposite endpoint: they have distances m-2 and m-3, which are distinct under (11) when m>=5. All positions are therefore determined. Choosing the other endpoint reverses the complete order.

**Theorem.** If two transitive orders on the same module A of size m>=5 have the same exterior incidences and all the same internal shared-name quartics, they are equal or exact whole converses. Conversely these internal shared-name quartics are identical under whole converse of A, since all distances remain unchanged.

The threshold m>=5 is sharp for this shared-name battery: at m<=4, after deleting a shared name every possible distance is at most2, and all available entries coincide for every transitive order. This does not claim that the full two-port response fails to reconstruct smaller transitive modules.

The proof does not need the numerical values theta_r or W_d in advance. Only the ordering and equality relations among the observed S_jk are used.

There are m*binom(m-1,2) displayed named coordinates. Each is a signed finite difference of sixteen actual physical queries, with all other core bits zero and any prescribed ambient bits unchanged. Thus at most 16*m*binom(m-1,2) actual queries recover the unoriented order; duplicate queries may be removed. Every such coordinate has the complete port bidegree (2,2), so no omitted cubic or quartic alias can enter it.

## 5. Restriction on original-module collision candidates

Let U be obtained from a tournament T by simultaneous XOR of modules of the original T. Assume F_T=F_U. The reviewed Session188 hereditary monotonicity theorem supplies equality of all ordered-cover profiles on every fixed induced subset. Its highest-node reduction supplies a common original strong-module tree and says each highest changed node M has a prime quotient, changed to its whole converse. Other highest nodes are disjoint and can be restored in any query whose ports imitate M on the exterior.

Take a child A of M for which T[A] is transitive and |A|>=5. Hereditary equality forces U[A] to be transitive too: every fixed triple retains Hamilton count one. We claim

    U[A] is either T[A] or T[A]^op.                     (12)

Use the physical internal quartics L_ijk with i,j,k in A, while both ports imitate M at all exterior names. The coefficient is saturated in both port degrees, so nonzero prescribed exterior bits introduce no extra alias terms. Restore the other highest nodes by their equal complete profiles.

For completeness, the all-cover top-degree restoration argument is the same one already reviewed in `quartic_cover_transfer/PROOF.md` and `ordered_pair_ambient/PROOF.md`. Let P_h be the ordered h-cover response of the augmented core. A degree-(2,2) term has each port internal to one of the paths, and no port-port edge is used. Reversing each directed cover path proves invariance of this top degree under whole converse of the old core. Apply this to U[M]: its quotient becomes T's quotient and its A interior becomes U[A]^op. Every other child is uniform to both ports and has the same complete profile as its source child, so restore those children in every h-coordinate. Run substitution then gives the same equality in every fixed ambient, with A the only changed old interior.

Full-response equality therefore makes all the internal shared-name quartics of this local comparison agree. Section4 forces U[A]^op to be T[A] or T[A]^op, which is equivalent to (12).

This covers arbitrary crossing or nested descendant selections within A. It does not assume they form a legal sequential family. It uses the original-module hereditary theorem before profile restoration.

If every child of a highest changed prime node is transitive of size at least five, all surviving full-response candidates at that node have only binary child choices, unchanged or wholly conversed. The quotient remains conversed. This is a reduction to a reversal cube, not its exclusion. Fully covered quotient-only converses and the larger-prime obstruction remain open.

## 6. Reproducibility and scope

`verify.py` checks (2) against literal permutations with marked cuts, the positivity identities, physical-order reconstruction, and exact ambient capacity counts. `bank.cpp` computes complete named response coefficient vectors by literal augmented Hamilton subset DP followed by Boolean finite differences; the canonical basis includes every zero, every cubic, every quartic, and every lower coefficient. `receipt.json` reports separately which tests are finite. Checks use explicit exceptions, so optimized Python retains them.

The proofs in Sections2--5, rather than the finite examples, establish the unbounded claims. The packet makes no full-response collision assertion. A final independent review is required before promoting this result.
