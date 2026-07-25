---
id: THM-2195
title: "Transitive quotients exactly control universal tournament substitution products"
status: >
  PROVED. Fix a quotient tournament Q. For all factor tuples of equal
  corresponding positive orders, unlabeled arc-reversal distance between
  Q-substitutions is the sum of corresponding factor distances if and only
  if Q is transitive. The positive direction iterates THM-2183's order-join
  isometry. For every nontransitive Q, a directed triangle permits a cyclic
  three-block transport whose external cost is only linear in the block
  order, while two internal factor distances can be chosen quadratic. For
  Q=C3 there is an explicit nine-vertex witness with product sum two and
  substitution distance zero. This classifies the fixed-correspondence
  product law; it does not classify the automorphism-gauged or optimal
  block-transport distance.
source: codex-2026-07-24-tournament-substitution-product
depends_on:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
related:
  - THM-1960
  - THM-2176
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

## 5. Boundary and the surviving problem

THM-2183's two-image swap works because a transitive quotient presents
constant cuts in one consistent order. A directed triangle permits a
cyclic block transport. In a larger quotient, its cost against the exterior
is linear in `N`; internal tournament separation is quadratic, so the
triangle transport eventually wins.

The theorem classifies exactly the displayed, fixed-correspondence identity
(3). Put `n_i=|T_i|=|S_i|` and define the size-compatible automorphisms

```text
Aut_n(Q)={sigma in Aut(Q):n_i=n_(sigma(i)) for every i}.
```

It does **not** prove or refute the well-typed automorphism-gauged formula

```text
min_(sigma in Aut_n(Q)) sum_i d_iso(T_i,S_(sigma(i))), (24)
```

nor a formula allowing fractional or vertex-level transport between
quotient blocks. The witness (23) is absorbed by (24), because cyclic
rotation is an automorphism of `C_3`. For a general `Q`, the construction
uses a triangle rotation which need not be a global automorphism and pays
the exterior invoice (13). The correct next object is therefore a block
transport matrix carrying both internal edit costs and quotient-arc
disagreement costs.

This boundary is the tournament analogue of THM-2176's distinction between
decomposition-respecting cost and bypass after forgetting a product marker.
Here transitivity is exactly the condition which makes every such bypass
uncrossable.

QED.
