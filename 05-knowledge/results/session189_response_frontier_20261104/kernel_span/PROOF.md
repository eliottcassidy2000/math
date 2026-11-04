# Full kernel span and matching-based local noncancellation

**New restricted-battery structural extension of the exact two-channel identity.** Generic detectability of a nontrivial local module converse is already supplied by incoming E4, pinned at [two_port_parity_separation_20261104/PROOF.md](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md). This note does not claim first local noncancellation. Its contribution is the full odd-transitive mixed-quartic kernel span, exact moment/jet description, and a constructive maximum-matching witness within that restricted explicit battery for every nonempty physical exterior once the stated size condition holds. The theorem permits arbitrary grouping of exterior marks provided each group sees the transitive module uniformly. It is a local-converse theorem. It does not by itself supply a quotient-only, source-relative root-family theorem in an arbitrary fixed ambient.

## 1. Inherited exact identity and notation

Use the physical setup and exact identity in `../session189-odd-transitive/PROOF.md`, Sections 1--3. Let the reversed module A have physical transitive order `0 -> 1 -> ... -> m-1`, where `m=2k+1 >= 3`. Let `R=T-A` have `n>=1` physical names. Put

- `O={u in R: A -> u}` and `I={v in R: v -> A}`;
- `alpha=1+sum_(v in I) E_v`, `beta=1+sum_(u in O) S_u`;
- `tau=sum_(u in O,v in I) G_uv`;
- `ell_i=sum_(v in I) G_iv`, `rho_i=sum_(u in O) G_ui`;
- `A_i=alpha S_i-beta E_i`, `B_i=alpha beta(ell_i-rho_i)`.

All products and extractions use the squarefree algebra on the physical names of R. In particular, repeated names annihilate a product. The complete physical named quartic contrast is

    Delta_(01;r,i) = Lambda_i(X_r),
    Lambda_i(f) = [R](A_i f(tau)+B_i f'(tau)).              (1)

Here `X_r=2 E(P_(r-1))`, with

    P_u(s)=(s+1)^u s^(L-u)-s^u(s+1)^(L-u),  L=m-2,
    E(P)(t)=sum_j Delta^j P(0) t^j.                       (2)

Equation (1) compares the source with the tournament obtained by reversing only A. It is not, without a separate restoration theorem, an identity for reversing only a root quotient or an XOR family of original modules.

## 2. Full odd kernel span

Write `w=t(1+t)`.

**Theorem 1 (full span).** The k kernels indexed by

    r in {2,...,k} union {m-1}

form a basis over Q of the vector space

    { f(t(1+t)): f in Q[w], deg f <= k-1 }.              (3)

When k=1, the first set is empty. Consequently, each moment `Lambda_i(w^h)`, `0<=h<=k-1`, is a rational linear combination of those physical quartic contrasts. Clearing denominators gives an integer query certificate.

**Proof.** The transform E is a degree-preserving invertible linear map: the Newton basis `binom(s,j)` maps to `t^j`. For reflection `J_s P(s)=P(-1-s)`, the identity

    binom(-1-s,j)=(-1)^j binom(s+j,j)

and Vandermonde's formula give

    E(J_s P)(t)=E(P)(-1-t).                              (4)

Because L is odd, `P_u(-1-s)=P_u(s)`. Each P_u has degree at most `L-1=2k-2`. The invariant polynomials of that degree are exactly the polynomials in `s(s+1)` of degree at most k-1, a k-dimensional space. (Translate by 1/2 to reduce to even polynomials.)

The selected polynomials are `P_1,...,P_(k-1),P_L`. The polynomial P_L has nonzero constant term 1. For `1<=u<=k-1`, P_u has order of vanishing exactly u at s=0, with coefficient -1 at s^u, since `L-u>u`. They are linearly independent by their distinct lowest nonzero powers. Thus they span the full invariant space. Equation (4) and invertibility of E identify its image with (3). The factor 2 does not affect the rational span. QED.

## 3. Exact truncation and saturation

Define raw moments

    c_(i,0)=[R] A_i,
    c_(i,j)=[R](A_i tau^j+j B_i tau^(j-1)),  j>=1.       (5)

Let `d=floor((n+1)/2)`.

**Lemma 2 (finite jet).** For every i, `Lambda_i(t^j)=c_(i,j)=0` whenever j>d. If `k-1>=d`, then all selected quartics for i vanish if and only if `c_(i,j)=0` for every `0<=j<=d`.

**Proof.** Every monomial of tau uses at least two names: its paths begin in O and end in the disjoint set I. The singleton terms of S_i and E_i cancel in A_i, so A_i has support degree at least 2. B_i has support degree at least 1. The two terms in (5) therefore require at least `2j+2` and `2j-1` names, respectively. This proves the cutoff.

The polynomials `1,w,...,w^d` form a basis of `Q[t]/(t^(d+1))`: `w^h=t^h(1+t)^h` has coefficient 1 at t^h and no lower terms. Thus their moment values determine all raw moments through d. Theorem 1 supplies all these w moments if `k-1>=d`. QED.

More explicitly,

    Lambda_i(w^h)=sum_(j=h)^(2h) binom(h,j-h)c_(i,j),    (6)

with moments beyond d set to zero. Descending substitution in this unitriangular system recovers the raw moments. Equivalently, for 1<=j<=d use the explicit polynomial

    f_j(w)=sum_(h=j)^d (-1)^(h-j) (j/h)
                    binom(2h-j-1,h-j) w^h.             (6a)

The formal inverse of w=t(1+t), or Lagrange inversion for t=w/(1+t), gives `f_j(t(1+t))=t^j mod t^(d+1)`. Hence `Lambda_i(f_j(w))=c_(i,j)`. Use f_0=1. This is an exact structural argument, not extrapolation from finite ranks.

## 4. A matching invariant of the exterior

Let B be the bipartite graph with vertex classes O and I and an edge `{u,v}` precisely when the tournament has `u -> v`, for u in O and v in I. Let q be its maximum matching size.

**Lemma 3.** `tau^(q+1)=0` in the squarefree algebra.

**Proof.** Each directed O-to-I path has a consecutive edge directed from O to I: follow the path until it first enters I. From each of q+1 pairwise disjoint O-to-I paths, select such an edge. They would form a matching of size q+1. Every nonzero product contributing to `tau^(q+1)` would therefore contradict maximality. QED.

This uses the physical exterior graph, not a graph in which multi-name child tournaments have been replaced by unweighted vertices.

## 5. Nonperfect matching: a signed highest moment

Assume q has no perfect matching on all n names, equivalently `2q<n`.

**Lemma 4.** At moment `j=q+1`,

    c_(i,q+1)=-(q+1)[R]alpha beta rho_i tau^q <=0  (i in O),
    c_(i,q+1)= +(q+1)[R]alpha beta ell_i tau^q >=0  (i in I).  (7)

The inequality is strict whenever some maximum matching leaves i unmatched. In particular, at least one moment in (7) is nonzero.

**Proof.** Lemma 3 kills the A_i term. If i is in O, ell_i is a sum of O-to-I paths, and `ell_i tau^q=0` by the same matching argument. If i is in I, `rho_i tau^q=0`. This gives (7), whose remaining coefficients count paths with nonnegative weights.

Choose a maximum matching leaving i unmatched. Use its q oriented edges as the q factors of tau. Use the singleton i as rho_i if i is in O, or as ell_i if i is in I. On the remaining unmatched I names, choose a Hamilton path and put it into alpha; on the remaining unmatched O names, choose a Hamilton path and put it into beta. If either set is empty, use the corresponding constant term 1. Every tournament has a Hamilton path, by elementary vertex insertion into a path. These factors have disjoint support and cover R. They contribute positively to the remaining unsigned extraction in (7). In fact its magnitude is at least `(q+1)!`. QED.

## 6. Perfect matching: the next signed moment

Now `n=2q`, `|O|=|I|=q`, and there is a perfect O-to-I matching. Suppose first q>=2. Let

    J=sum_(u in O,v in I,u->v) z_u z_v,
    M_(i,j)=[R-{i,j}] J^(q-1),
    T_i=[R-{i}] tau^(q-1).                              (8)

The symbol T_i here is a nonnegative integer extraction, not a tournament. For i in O and j in I, M_(i,j) counts ordered cross-edge matchings of the complement.

**Lemma 5.** At moment q,

    c_(i,q)=-q( sum_(j in I,j->i) M_(i,j)+T_i ) <=0  (i in O),
    c_(i,q)=+q( sum_(u in O,i->u) M_(u,i)+T_i ) >=0  (i in I).  (9)

At least one O mark is strictly negative and at least one I mark is strictly positive.

**Proof.** The A_i term has support degree at least `2q+2>n`. Fix i in O. In the positive term `alpha beta ell_i tau^(q-1)`, all q paths must be cross edges and alpha beta must contribute its constant. Its extraction is

    sum_(j in I,i->j) M_(i,j).

In the negative term `alpha beta rho_i tau^(q-1)`, rho_i must be the singleton i: a nonsingleton O-to-O path, together with q-1 disjoint cross paths, would require q+1 distinct O names. The q-1 tau paths either are all edges, leaving one additional I name for the degree-one part of alpha beta, or have total support `2q-1`, in which case alpha beta is constant. The negative extraction is therefore

    sum_(j in I) M_(i,j)+T_i.

Subtracting proves the first formula in (9). The I formula is obtained identically, exchanging the positive singleton ell_i and negative cross-path rho_i channels.

To prove strictness on O, fix a perfect matching `u_a -> v_a`, `a=1,...,q`. Choose any directed edge `v_a -> v_b` between two distinct I names. Leave `i=u_b` out. The path `u_a -> v_a -> v_b`, together with all matching edges except those indexed a,b, is a family of q-1 O-to-I paths covering `R-{i}`. Thus `T_i>0`, indeed `T_i>=(q-1)!`, and `c_(i,q)<=-q!`. For strictness on I, choose an O edge `u_a -> u_b`, omit i=v_a, and use `u_a -> u_b -> v_b` with the other matching edges. The same argument applied to the partners of any two specified names on one side shows that at least one of those two marks is strict. Thus every same-side group with at least two physical names is strict in (9). QED.

There is one small exception to the signed-high-moment construction, not to noncancellation. If q=1 and the matching is perfect, R is the single edge `o -> v`, with o in O and v in I. Direct substitution gives

    c_(o,0)=2,  c_(v,0)=-2.                             (10)

Thus the constant moment already separates this directed three-cycle after contracting A to a root.

## 7. Universal saturated local-converse separation

**Theorem 2.** Let A be an odd transitive module of size m>=3 with a nonempty arbitrary exterior R of n physical names. If

    m >= 2 floor((n+1)/2)+3,                            (11)

then at least one selected named mixed quartic `Delta_(01;r,i)` is nonzero for reversal of A alone.

More generally, partition R into any nonempty groups, each contained entirely in O or entirely in I. At least one selected grouped quartic is nonzero. The groups need not be modules; uniform incidence to A is sufficient. In particular, grouping actual exterior modules that see A uniformly cannot hide the separation.

**Proof.** If the cross-edge maximum matching is not perfect, Lemma 4 supplies a nonzero raw moment with fixed sign on each incidence side. Its index `q+1<=floor((n+1)/2)=d`. If matching is perfect and q>=2, Lemma 5 supplies a nonzero raw moment of index q=d, again with fixed sign on each incidence side. The remaining case is (10). No partition into same-side groups can cancel every such nonzero moment. Lemma 2 and Theorem 1 then supply a nonzero physical or grouped quartic among the selected kernels.

For a group G, the grouped quartic is exactly `sum_(i in G) Delta_(01;r,i)`. The queried coefficient has already used two a variables and two c variables, so the inherited bidegree bound leaves no higher alias route. QED.

The sufficient uniform threshold (11) is `m>=n+3` for even n and `m>=n+4` for odd n. The proof gives sharper graph-dependent thresholds:

- nonperfect cross matching: `m>=2q+5`;
- perfect matching q>=2: `m>=2q+3=n+3`;
- perfect matching q=1: every odd m>=3 already works.

At most k selected kernels per exterior group are needed, each recovered by at most sixteen Boolean finite-difference queries. Nonvanishing of the linear certificate implies at least one actual two-port query differs.

## 8. Scope and limitations

1. Theorem 2 concerns source versus **local converse of A**, with all exterior names and arcs fixed. It is stronger than constant-kernel noncancellation, but is not a claim that the constant contracted skew alone is always nonzero.
2. All matching sizes, degree cutoffs, and the domination condition count **physical names**. Replacing a multi-name module by one quotient vertex would not preserve this argument without a separate profile-weighted theorem.
3. The grouped statement requires each group to have uniform incidence to A. Allowing a group to mix O and I can cancel the signed witnesses. For example, in the three-cycle contraction of (10), grouping both exterior names together kills every raw moment: the two constant moments cancel and all positive moments vanish. This is a genuine counterexample to unrestricted grouping, not a full named-response collision.
4. The nonzero group may lie anywhere in the full exterior R. If the problem only permits marking names inside a designated core, while a fixed outside ambient must be imitated, the present theorem does not ensure an admissible mark. No arbitrary-ambient root-family conclusion is asserted on that basis.
5. A quotient-only source-relative root converse requires its own all-cover restoration and source-family reduction. Applying this local theorem to such a target without those steps would be invalid.
6. The proof requires only the previously established exact identity, standard polynomial algebra, elementary matching reasoning, and Hamilton paths in tournaments. Finite controls verify the implementation, not the all-order theorem.
