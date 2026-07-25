---
id: THM-2195
title: "Transitive quotients exactly control universal tournament substitution products"
status: >
  PROVED + VERIFIED-EXACT. Fix a quotient tournament Q. For all factor tuples of equal
  corresponding positive orders, unlabeled arc-reversal distance between
  Q-substitutions is the sum of corresponding factor distances if and only
  if Q is transitive. The positive direction iterates THM-2183's order-join
  isometry. For every nontransitive Q, a directed triangle permits a cyclic
  three-block transport whose external cost is only linear in the block
  order, while two internal factor distances can be chosen quadratic. For
  Q=C3 there is an explicit nine-vertex witness with product sum two and
  substitution distance zero. This classifies the fixed-correspondence
  product law. Minimizing the factor sum over size-compatible quotient
  automorphisms still fails universally: the cyclic five-tournament gives
  an asymptotic counterexample. The canonical rotation of an equal-size
  directed-triangle block has exact external cost twice its block order times
  the total order of exterior quotient vertices whose triangle-incidence word
  is nonconstant. Thus the cyclic-five witness costs 4N, not merely at most
  6N, and wins once d_iso(A,B)>2N. More generally, every partial automorphism
  of an induced quotient has prescribed-block cost equal to the weighted
  Hamming derivative of its exterior incidence words. Its zero kernel is
  exactly the subgroup extending to the whole quotient, while each cyclic
  orbit records twice its number of nonconstant binary runs. In signed
  coordinates the derivative is a constant minus the permutation trace of
  a positive-semidefinite Gram matrix. Averaging over any automorphism subgroup is exactly the
  exterior-word variance lost under projection to the subgroup-orbit
  constants; this quantitatively detects failure of the whole subgroup to
  extend. These statements are exact only with block markers retained;
  unrestricted vertex transport can be cheaper. The remaining object is an
  optimal unmarked block-transport cost, not a product formula.
source: codex-2026-07-24-tournament-substitution-product
depends_on:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
related:
  - THM-1960
  - THM-2176
  - THM-2216-residual-capacity-hinge-gram-law
script: 04-computation/tournament_partial_automorphism_derivative_thm2195.py
output: 05-knowledge/results/tournament_partial_automorphism_derivative_thm2195.out
script_sha256: 596740b0b85cba16d5a84c82b90dbdebdbc8a7b67ded24aa303cd89e1ed6485b
output_sha256: b4985d152d600cb014e8e79bb2dcefd908a80269144bbbe0f5d32ffd7f30e61e
hash_basis: working-tree bytes (LF)
---

# THM-2195 -- transitive quotients exactly control universal substitution products

For tournaments `T,S` of equal order, let

```text
d_iso(T,S)
 =min_(bijections pi:V(T)->V(S))
   #{unordered {x,y}: pi reverses the orientation of {x,y}}. (1)
```

This is the unmerged unlabeled arc-reversal distance used in THM-2183.

Let `Q` be a tournament on the fixed vertex set `{1,...,q}`. Given
tournaments `T_1,...,T_q`, its substitution

```text
Q[T_1,...,T_q]                                      (2)
```

has the arcs of `T_i` inside block `i`, while every vertex of block `i`
beats every vertex of block `j` exactly when `i->_Q j`.

## 1. Exact classification

> **Theorem.** The identity
>
> ```text
> d_iso(Q[T_1,...,T_q],Q[S_1,...,S_q])
>   =sum_(i=1)^q d_iso(T_i,S_i)                       (3)
> ```
>
> holds for every two factor tuples satisfying
>
> ```text
> |T_i|=|S_i|>=1                 for every i          (4)
> ```
>
> if and only if `Q` is transitive.

The quantifier "every" in this statement is load-bearing. A particular
nontransitive quotient and a restricted factor class can still satisfy (3).

## 2. Transitive quotients

If `Q` is transitive, relabel its vertices in their unique linear order.
Then

```text
Q[T_1,...,T_q]=T_1 join T_2 join ... join T_q,       (5)
```

where every earlier block beats every later block. THM-2183 proves that
order-join is an exact `l1` product for `d_iso`, with no marked-cut
hypothesis. Iterating that theorem gives

```text
d_iso(T_1 join ... join T_q,S_1 join ... join S_q)
 =sum_i d_iso(T_i,S_i).                              (6)
```

Condition (4) makes the two accumulated prefix orders equal at every
iteration. This proves the forward direction of (3).

## 3. Every nontransitive quotient fails

A tournament is transitive exactly when it contains no directed triangle.
Suppose therefore that `Q` is nontransitive, and choose vertices

```text
C={i,j,k},              i->_Q j->_Q k->_Q i.         (7)
```

Let `rho` cyclically permute this triangle:

```text
rho(i)=j,       rho(j)=k,       rho(k)=i.             (8)
```

No claim is made that `rho` extends to an automorphism of `Q`.

Fix a positive integer `N` and two `N`-vertex tournaments `A,B`. Put

```text
(T_i,T_j,T_k)=(A,B,A),
S_(rho(l))=T_l                 for l in C.            (9)
```

At every vertex outside `C`, take

```text
T_l=S_l=the one-vertex tournament.                   (10)
```

Thus the corresponding factor orders agree. On the triangle, (9) reads

```text
(S_i,S_j,S_k)=(A,A,B).                              (11)
```

Consequently the fixed-correspondence right side of (3) is

```text
2d_iso(A,B).                                        (12)
```

Now map every source block `l in C` bijectively to target block `rho(l)`
using an internal isomorphism supplied by (9), and fix the outside
singletons. The internal cost is zero. The cross-block cost inside `C` is
also zero because cyclic rotation preserves the directed triangle.
Outside `C`, all singleton-to-singleton arcs agree.

Only pairs formed by a triangle vertex and an outside singleton can
disagree. There are exactly

```text
3N(q-3)                                             (13)
```

such unordered pairs. Hence

```text
d_iso(Q[T_1,...,T_q],Q[S_1,...,S_q])
 <=3N(q-3).                                         (14)
```

It remains to choose `A,B` so that

```text
2d_iso(A,B)>3N(q-3).                                 (15)
```

This is always possible for sufficiently large `N`. Write

```text
M=binom(N,2),             R=floor(3N(q-3)/2).        (16)
```

Fix any labelled `N`-vertex tournament `A`. For one bijection of its
vertices, at most

```text
sum_(s=0)^R binom(M,s)                               (17)
```

labelled tournaments lie within reversal radius `R`. Taking the union over
all `N!` bijections shows that the number of labelled `B` with
`d_iso(A,B)<=R` is at most

```text
N! sum_(s=0)^R binom(M,s).                           (18)
```

For fixed `q`, the base-two logarithm of (18) is `O_q(N log N)`, whereas
there are `2^M` labelled tournaments and `M` is quadratic in `N`.
Indeed, when `R>=1`, for all large `N` one has `R<=M/2` and

```text
sum_(s=0)^R binom(M,s)<=(R+1)(eM/R)^R,
log_2(N!)<=N log_2 N.                                (18a)
```

Both logarithms on the right are `O_q(N log N)=o(M)`. When `R=0`
(equivalently `q=3`), the sum in (18) is one and the same comparison is
immediate.
Therefore, for all sufficiently large `N`,

```text
N! sum_(s=0)^R binom(M,s)<2^M.                       (19)
```

Choose `B` outside that union. Then `d_iso(A,B)>R`, which implies (15).
Equations (12), (14), and (15) strictly contradict (3). Thus no
nontransitive quotient has the universal fixed-correspondence product law.

## 4. Sharp explicit triangle witness

For `Q=C_3`, no asymptotic choice is needed. Let

```text
A=the transitive three-vertex tournament,
B=the directed three-cycle.                         (20)
```

The two isomorphism classes differ by one arc reversal, so

```text
d_iso(A,B)=1.                                       (21)
```

Use the tuples

```text
T=(A,B,A),               S=(A,A,B).                 (22)
```

Their substitutions have nine vertices. The corresponding-factor sum is
two. Cyclically rotating the three blocks is an isomorphism from `C_3[T]`
to `C_3[S]`, so

```text
d_iso(C_3[T],C_3[S])=0<2=sum_i d_iso(T_i,S_i).       (23)
```

This is minimal in mechanism: block order two has only one tournament
isomorphism class, so order three is the first factor order that can carry
the mismatch.

## 5. The automorphism gauge also fails

Put `n_i=|T_i|=|S_i|` and define the size-compatible automorphisms

```text
Aut_n(Q)={sigma in Aut(Q):n_i=n_(sigma(i)) for every i}.
```

The next natural candidate would be the well-typed automorphism-gauged sum

```text
D_aut(T,S)
 =min_(sigma in Aut_n(Q)) sum_i d_iso(T_i,S_(sigma(i))). (24)
```

The witness (23) is absorbed by (24), because cyclic rotation is an
automorphism of `C_3`. Nevertheless, (24) is not a universal formula.

Take the cyclic five-tournament on `Z/5Z`,

```text
x->_Q y        iff        y-x in {1,2}.              (25)
```

Its automorphism group consists exactly of the five translations. Indeed,
after composing with a translation, an automorphism may be assumed to fix
zero. It must preserve the outneighborhood `{1,2}` and its oriented edge
`1->2`, hence fixes both vertices; the same argument fixes the
inneighborhood `{3,4}`.

The vertices

```text
C={0,1,3}                                             (26)
```

form the directed triangle

```text
0->1->3->0.
```

Give these three blocks order `N` and the other two blocks order one. No
nonidentity translation preserves this block-order vector: a nonempty
proper subset of `Z/5Z` is not invariant under a nonzero translation.
Therefore

```text
Aut_n(Q)={identity}.                                  (27)
```

Apply the construction (9)--(11) on `C`. The automorphism-gauged cost is
then still

```text
D_aut(T,S)=2d_iso(A,B),                              (28)
```

whereas the partial triangle rotation pays at most

```text
3N(5-3)=6N.                                          (29)
```

The counting argument (16)--(19), now with `q=5`, supplies `N`-vertex
tournaments with `d_iso(A,B)>3N`. For those factors,

```text
d_iso(Q[T_i],Q[S_i])<=6N<D_aut(T,S).                 (30)
```

Thus quotient automorphisms do not account for every substitution bypass.

### Exact cost of the partial triangle rotation

The coarse count in (13) can be replaced by an exact quotient derivative.
Allow the outside blocks to have arbitrary positive orders `m_t`, equal in
the source and target, and keep identical factors there. For `t notin C`, put

```text
epsilon_C(t)=
  0  if t beats all of C or loses to all of C,
  1  otherwise.                                      (31)
```

Use the canonical map which rotates the three `N`-vertex blocks by `rho`,
uses the internal isomorphisms from (9), and fixes every outside block. Its
cost is exactly

```text
E_Q(C,rho;m)=2N sum_(t notin C) m_t epsilon_C(t).     (32)
```

Indeed, fix `t notin C` and encode its incidences by the cyclic binary word

```text
b_t(l)=1_(l->_Q t),                  l in C.          (33)
```

The number of quotient incidences changed by the rotation is

```text
#{l in C:b_t(l)!=b_t(rho(l))}.                        (34)
```

A cyclic binary word of length three has zero transitions when constant and
exactly two transitions otherwise. Each changed incidence contributes all
`N m_t` cross-block pairs. Internal pairs, pairs within `C`, pairs wholly
outside `C`, and constant-incidence exterior pairs contribute zero. Summing
over `t` proves (32).

For the cyclic five-tournament in (25)--(26), both exterior vertices `2,4`
have nonconstant incidence to `C`, and both exterior blocks are singletons.
Therefore the displayed block map has exact cost

```text
E_Q(C,rho;m)=4N.                                      (35)
```

Consequently (30) sharpens to

```text
d_iso(Q[T_i],Q[S_i])<=4N<D_aut(T,S)
              whenever d_iso(A,B)>2N.                (36)
```

The construction in Section 3 already supplies factors with
`d_iso(A,B)>3N`, so the witness remains unconditional. Equation (32) is an
exact cost for this integral block map, not a proof that it is globally
optimal among vertex-level bijections. Its preserved sidecar is the exterior
incidence word, compressed losslessly here to the nonconstant mass in (31).

## 6. Partial-automorphism derivative

The triangle formula is the first case of a general exact marked-block law.
Let `I=V(Q)`, let `C subset I`, and let

```text
rho in Aut(Q[C]).
```

Extend it to a permutation `sigma` of `I` by fixing `I\C`. Suppose the
source and target blocks satisfy

```text
|T_i|=|S_(sigma(i))|=:n_i.                           (37)
```

Let `d_sigma` be the minimum reversal cost among bijections which send each
whole source block `T_i` onto the target block `S_(sigma(i))`. For an
unordered quotient pair define the well-defined bit

```text
Delta_Q(i,j;sigma)
 =1_[(i->_Q j) xor (sigma(i)->_Q sigma(j))].          (38)
```

Then

```text
d_sigma(Q[T_i],Q[S_i])
 =sum_(i in I)d_iso(T_i,S_(sigma(i)))
  +sum_({i,j} subset I)n_i n_j Delta_Q(i,j;sigma).   (39)
```

Indeed, unordered vertex pairs partition into internal block pairs and
cross-block pairs. Internal minimizations are independent. On a fixed
quotient pair either all `n_i n_j` cross pairs agree or all reverse. This
also proves the corresponding formula for an arbitrary size-compatible
block permutation; minimizing it over such permutations is the exact
whole-block comparator.

Now assume the blocks over `C` all have order `N`, and put `n_t=m_t` for
`t notin C`. Encode the incidence of an exterior vertex by

```text
b_t(c)=1_[c->_Q t],                    c in C.        (40)
```

The `C`--`C` contribution in (39) vanishes because `rho` is an automorphism
of `Q[C]`; the exterior--exterior contribution vanishes because `sigma`
fixes the exterior. Hence the exact external derivative is

```text
E_Q(C,rho;m)
 =N sum_(t notin C)m_t d_H(b_t,b_t composed rho).    (41)
```

There is no factor `1/2`: every changed coordinate represents exactly
`N m_t` unordered cross-block pairs.

Writing

```text
A_t={c in C:c->_Q t},
```

gives the equivalent set formula

```text
d_H(b_t,b_t composed rho)
 =|A_t symmetric_difference rho^(-1)(A_t)|
 =2(|A_t|-|A_t intersection rho^(-1)(A_t)|).         (42)
```

Thus every row cost is even and `E_Q` is a multiple of `2N`. More
structurally,

```text
E_Q(C,rho;m)=0
 iff every b_t is constant on every rho-orbit
 iff rho extends by the identity to an automorphism of Q.       (43)
```

If `rho` is transitive on `C`, the zero condition says exactly that `C` is
a tournament module.

For a cyclic orbit

```text
O=(c_0,...,c_(ell-1)),       rho(c_j)=c_(j+1),
```

its contribution to the Hamming distance is the cyclic transition count

```text
tau_t(O)=#{j:b_t(c_j)!=b_t(c_(j+1))}.                (44)
```

It is zero for a constant word and otherwise twice the number of cyclic
`1`-runs. The length-three case has only the values zero and two, recovering
(32). Longer cycles retain genuinely more data: a binary word on a
five-cycle can have transition cost two or four. The triangle's
nonconstant-exterior mass is therefore a low-order collapse of the cyclic
run ledger, not the general invariant.

Equivalently, with the weighted exterior signature

```text
kappa(c)=(b_t(c))_(t notin C),
```

equation (41) is

```text
E_Q(C,rho;m)
 =N sum_(c in C)d_m(kappa(c),kappa(rho(c))).          (45)
```

More generally,

```text
D_m(alpha,beta)
 =N sum_(t notin C)m_t
      d_H(b_t composed alpha,b_t composed beta)      (46)
```

is a right-invariant pseudometric on `Aut(Q[C])`. Its zero kernel consists
exactly of the automorphisms extending to `Q` while fixing the complement
pointwise. The derivative therefore measures the precise obstruction to
extending a local symmetry.

The scope marker remains essential. Equations (39)--(46) are exact among
unsplit block bijections inducing the prescribed permutation. They give a
constructive upper bound for unrestricted `d_iso`, not a positive lower
bound. If source and target substitutions are identical, the identity has
cost zero even when a chosen nonextendable partial rotation has positive
cost. Any global optimality claim still needs a block-preservation or
uncrossing theorem.

## 7. Signed Gram trace and subgroup-orbit variance

The Hamming derivative has an exact positive-semidefinite form.  Put

```text
n=|C|,                         M=sum_(t notin C)m_t,
s_t(c)=2b_t(c)-1 in {-1,+1},
G_C=sum_(t notin C)m_t s_t s_t^T.                  (47)
```

Thus `G_C` is positive semidefinite.  For `rho in Aut(Q[C])`, define its
permutation trace against `G_C` by

```text
Tr_rho(G_C)=sum_(c in C)G_C(c,rho(c)).              (48)
```

The binary sign identity

```text
d_H(b_t,b_t composed rho)
 =1/2[n-sum_c s_t(c)s_t(rho(c))]                   (49)
```

and (41) give

```text
E_Q(C,rho;m)
 =N/2[nM-Tr_rho(G_C)].                              (50)
```

This is a genuine Gram representation, not a spectral lower bound by
itself: different permutation traces of the same PSD matrix can still be
close to the identity trace.

There is, however, an exact averaged law.  Let

```text
Gamma <= Aut(Q[C])
```

be any subgroup, let `O` range over its vertex orbits in `C`, and put

```text
r_t(O)=#{c in O:b_t(c)=1}.                          (51)
```

Uniformly averaging (41) over `Gamma` gives

```text
1/|Gamma| sum_(rho in Gamma) E_Q(C,rho;m)
 =2N sum_(t notin C)m_t
      sum_O r_t(O)(|O|-r_t(O))/|O|.                 (52)
```

Indeed, for a fixed `c in O`, the image `rho(c)` is uniform on `O`.
Among all ordered pairs in `O`, exactly
`2r_t(O)(|O|-r_t(O))` have opposite bits.  Summing the resulting
probability over the `|O|` possible starting vertices proves (52).

Equivalently, let

```text
(P_rho v)(c)=v(rho(c)),
Pi_Gamma=1/|Gamma| sum_(rho in Gamma)P_rho          (53)
```

be the orthogonal projection onto vectors constant on every `Gamma`-orbit.
Then

```text
1/|Gamma| sum_rho E_Q(C,rho;m)
 =N/2 sum_(t notin C)m_t
       (||s_t||_2^2-||Pi_Gamma s_t||_2^2).          (54)
```

On an orbit `O`, the projected sign has value
`(2r_t(O)-|O|)/|O|`; expanding (54) therefore recovers (52).

The zero and positive cases are now exact:

```text
average derivative=0
 iff every exterior word b_t is constant on every Gamma-orbit
 iff every rho in Gamma extends by the identity on I\C
     to an automorphism of Q.                       (55)
```

If the subgroup fails to extend, one nonconstant orbit contributes at least
`N m_t`, because

```text
2r(|O|-r)/|O|>=1             for 1<=r<|O|.
```

Hence the average is at least `N min_(t notin C)m_t`, and some nonidentity
element of `Gamma` has derivative at least that average.  The exact formula
(52), rather than this coarse corollary, retains which exterior words and
orbits break the symmetry.

This is the tournament-side transfer of the Gram viewpoint used for
residual capacities in THM-2216.  There the Gram inner product upper-bounds
a selected top-`p` cover; here the Gram trace measures the cost of a
prescribed local symmetry.  In both cases positive semidefiniteness becomes
quantitative only after the relevant sidecar is retained: label alignment
there, subgroup orbits here.

## 8. Boundary and the surviving problem

THM-2183's two-image swap works because a transitive quotient presents
constant cuts in one consistent order. A directed triangle permits a
cyclic block transport. In a larger quotient, its cost against the exterior
is linear in `N`; internal tournament separation is quadratic, so the
triangle transport eventually wins.

The theorem now rules out both the fixed-correspondence product and its
quotient-automorphism repair outside their stated regimes. It does not
classify a formula allowing fractional or vertex-level transport between
quotient blocks. The correct next object is a block transport matrix carrying
internal edit costs, block-size conservation, and quotient-arc disagreement
costs.

This boundary is the tournament analogue of THM-2176's distinction between
decomposition-respecting cost and bypass after forgetting a product marker.
Here transitivity is exactly the condition which makes every such bypass
uncrossable.

## 9. Exact indexing referee

Run

```bash
python3 04-computation/tournament_partial_automorphism_derivative_thm2195.py
python3 -O 04-computation/tournament_partial_automorphism_derivative_thm2195.py
```

The companion exhausts all `1,098` labelled tournaments of orders two
through five, all `32,766` nonempty induced-subtournament automorphism
groups, and every automorphism in them: `41,026` partial-automorphism cases.
It separately checks literal expanded block pairs, quotient weights,
exterior Hamming words, symmetric differences, orbit transitions, parity,
the extension kernel, the right-invariant pseudometric law, the signed Gram
trace, and the full-automorphism-group orbit/projection averages. Ordinary
and optimized runs are byte-identical to the frozen output. The script is an
exact indexing referee; the subgroup averaging proof above applies to every
subgroup and in arbitrary order.

QED.
