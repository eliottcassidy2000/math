---
id: THM-2185
title: "Rank-two safe-cube floor and height-500 rank-three harvest"
status: >
  PROVED + VERIFIED-EXACT. For every nonzero thirteen-coordinate integer row
  v and every saturated rank-two lattice L of integer relations on v, the
  Haar mass of the 1/14-safe cube inside the relation torus K_L is at least
  15/343. The proof classifies the exact support-at-most-three relation
  geometry and uses a cubic danger-count minorant. THM-2164 and a positive
  Jackson comparison at degree 500 then imply that every zero-safe distinct
  positive thirteen-speed row has dim_Q W_500 at least three. Consequently
  every twelve-speed deletion of such a row has relation rank at least two
  by coefficient height 105000. This removes THM-2178's sparse mod-14 digit
  terminal, but it does not prove LRC(14) or preserve zero-safety under a
  relation-state pump.
source: codex-2026-07-24-rank-two-safe-cube
depends_on:
  - THM-2145-two-block-spectral-crossing-and-the-exact-low-core-relation-carry
  - THM-2164-relative-packet-rank-harvesting
related:
  - THM-2054-relative-fejer-whole-product-decorrelation
  - THM-2169-bounded-relation-on-every-lrc-deletion
  - THM-2171-radix-quotient-order-pumping
  - THM-2178-mod14-transverse-code-rank-harvest
script: 04-computation/lrc14_rank_two_safe_cube_rank_three_harvest_thm2185.py
output: 05-knowledge/results/lrc14_rank_two_safe_cube_rank_three_harvest_thm2185.out
script_sha256: f06880a669c9af0b7edb8b305b5948c929efec59822e41b091bef56a02290891
output_sha256: e03a265487b15d72f754e28915034fa9ae5a253778f271923512fcbd349d005c
hash_basis: working-tree bytes (LF)
---

# THM-2185 -- rank-two safe-cube floor and height-500 rank-three harvest

Put

```text
J=[1/14,13/14],
S(v)={t in R/Z:v_i t in J for every i},
Lambda(v)={m in Z^13:m.v=0}.                          (1)
```

For a saturated lattice `L subset Z^13`, write

```text
K_L={x in (R/Z)^13:l.x=0 mod 1 for every l in L}.     (2)
```

The central point is that a rank-two relation torus cannot avoid the safe
cube. The only possible failures of three-coordinate Haar independence are
so sparse that they can be classified exactly.

## 1. Projection and annihilator bookkeeping

Let `v in Z^13` have

```text
v_i!=0                         for every i,            (3)
```

and let

```text
L subset Lambda(v)
```

be saturated of rank two. Saturation gives

```text
ann(K_L)=L.                                           (4)
```

For a coordinate set `T subset {1,...,13}`, let `pi_T` be projection onto
those coordinates and regard `Z^T` as the integer vectors supported on `T`.
Character duality gives the exact identity

```text
ann(pi_T(K_L))=L intersection Z^T.                   (5)
```

Consequently:

1. if `L intersection Z^T={0}`, then `pi_T(K_L)` is the whole torus
   `(R/Z)^T`, and its Haar pushforward is Haar;
2. if `L intersection Z^T=Z a` for a primitive vector `a` supported on
   `T`, then `pi_T(K_L)` is exactly the subtorus `K_a`;
3. no nonzero member of `L` has support one, because
   `c e_i in Lambda(v)` would imply `c v_i=0`, contrary to (3).

The use of saturation is load-bearing in (4)--(5). It also makes `K_L`
connected, although the projection conclusions already follow directly
from the displayed annihilators.

## 2. A cubic danger-count minorant

On `K_L`, define

```text
D_i={x:||x_i||<1/14},
p=measure(D_i)=1/7,
N(x)=#{i:x in D_i}.                                  (6)
```

Every coordinate character is nontrivial by Section 1, so every `D_i` has
the displayed Haar mass. Consider

```text
P(n)=-(n-1)(n-3)(n-4)/12

    =1-C(n,1)+(5/6)C(n,2)-(1/2)C(n,3).               (7)
```

For every integer `0<=n<=13`,

```text
P(n)<=1_(n=0).                                       (8)
```

Indeed, `P(0)=1`, it vanishes at `1,3,4`, and it is negative at `2` and at
every `n>=5`. Therefore

```text
measure_(K_L)(K_L intersection J^13)
 =Pr(N=0)
 >=E[P(N)].                                          (9)
```

If every set of at most three coordinate characters is Haar independent,
then

```text
E C(N,k)=C(13,k)p^k,                 0<=k<=3.         (10)
```

Substitution in (7) gives the exact baseline

```text
1-13/7+(5/6)(78/49)-(1/2)(286/343)
 =18/343.                                            (11)
```

We next show that every failure of the hypothesis in (10) retains a uniform
positive floor.

## 3. Complete sparse-relation classification

Call a nonzero member of `L` **sparse** if its support has size at most
three. Section 1 rules out support one. Since `L` has rank two, exactly one
of the following cases occurs.

### Case A: no sparse relation

Equation (5) makes every projection onto at most three coordinates Haar.
Thus (10)--(11) apply and

```text
measure_(K_L)(K_L intersection J^13)>=18/343.        (12)
```

### Case B: exactly one rational sparse line, of support three

Let its primitive generator be `a`, with support `A` of size three. Every
one- and two-coordinate projection is still Haar. All three-coordinate
projections are Haar except the projection onto `A`. Put

```text
q=measure_(K_a)(intersection_(i in A) D_i).           (13)
```

Every pair projection of `K_a` is the full two-torus, because a pair
relation would be a second sparse line. Hence

```text
0<=q<=p^2.                                            (14)
```

Only the third factorial moment changes, from its independent value by
`q-p^3`. Equations (7), (11), and (14) give

```text
E[P(N)]
 =18/343-(q-p^3)/2
 >=18/343-(p^2-p^3)/2
 =15/343.                                            (15)
```

### Case C: exactly one rational sparse line, of support two

Let its primitive generator have support `A` of size two, and again let `q`
be the joint danger mass on that pair. The pair factorial moment changes by

```text
delta=q-p^2.                                         (16)
```

For each of the other eleven coordinates, (5) identifies the corresponding
three-coordinate image with

```text
K_a times R/Z.
```

Thus each of those eleven triple masses changes by `p delta`. The coefficient
of `delta` in (7) is

```text
5/6-(11/2)p=1/21.                                    (17)
```

Since `q>=0`,

```text
E[P(N)]
 =18/343+delta/21
 >=18/343-p^2/21
 =53/1029
 >15/343.                                            (18)
```

### Case D: two independent sparse relations

Choose independent sparse relations `a,b` and put

```text
U=supp(a) union supp(b),
m=|U|.                                               (19)
```

They span `L` over `Q`, so every relation in `L` is supported on `U`. We
have

```text
3<=m<=6.                                             (20)
```

The upper bound is immediate. If `m<=2`, the two-dimensional relation plane
on those coordinates would contain a support-one vector; equivalently, it
could not be orthogonal to the nonzero row `v_U`. This proves the lower
bound.

The torus and its Haar measure factor as

```text
K_L=K_(L|U) times (R/Z)^(13-m).                       (21)
```

Every core coordinate is Haar, so the union bound on the core and exact
independence on the free coordinates give

```text
measure_(K_L)(K_L intersection J^13)
 >=(1-m/7)(6/7)^(13-m).                              (22)
```

The right side decreases for `m=3,4,5,6`. At `m=6` it is

```text
(1/7)(6/7)^7
 =279936/5764801
 =15/343+27831/5764801
 >15/343.                                            (23)
```

Cases A--D are exhaustive. We have proved the central theorem:

> **Rank-two safe-cube floor.** If `v_i!=0` for all thirteen coordinates and
> `L subset Lambda(v)` is saturated of rank two, then
>
> ```text
> measure_(K_L)(K_L intersection J^13)>=15/343.       (24)
> ```

This theorem does not use positivity or distinctness. Those hypotheses enter
only through the LRC relation supplier below.

## 4. Positive Jackson comparison and the third relation

Use THM-2145's normalized squared-Fejer kernel

```text
J_N=F_N^2/integral F_N^2
```

and put

```text
q_N=J_N*1_J.                                         (25)
```

Then

```text
0<=q_N<=1,
degree(q_N)<=2N-2,
||q_N-1_J||_1<=eta_N,                                (26)
```

where

```text
eta_N
 =1/2-[4/(pi^2 C_0)]
       sum_(1<=k<=2N-3, k odd) C_k/k^2,              (27)

C_0=N(2N^2+1)/3,

C_k=(4N^3-6Nk^2+2N+3k^3-3k)/6,       0<=k<=N,
C_k=((2N-k)^3-(2N-k))/6,              N<k<=2N-2.     (28)
```

All `C_k` in (28) are positive. The companion independently reconstructs
them by convolving the two triangular Fejer coefficient lists.

Take

```text
N=251,                       H=2N-2=500.              (29)
```

Using `pi<355/113`, the finite exact odd-mode sum in (27) gives

```text
eta_251<21/12500.                                    (30)
```

The adjacent `N=250` rational-pi upper-bound ledger does not clear the
required target `15/(343*26)`. This is a hostile boundary for this exact
certificate, not an optimality theorem for Jackson smoothing or relation
height.

Now let `v` be a distinct positive thirteen-speed row satisfying

```text
measure(S(v))=0.                                     (31)
```

THM-2164 gives

```text
dim_Q W_105(v)>=2.                                   (32)
```

Choose two independent height-`105` relations and let `L` be the saturation
of their integer span. Then `L` is a saturated rank-two sublattice of
`Lambda(v)`.

Suppose for contradiction that

```text
dim_Q W_500(v)=2.                                    (33)
```

Every relation in `Lambda(v)` of coefficient height at most `500` then lies
in the rational plane of `L`; saturation makes it a member of `L`. Hence the
finite Fourier expansions of

```text
I_line=integral_(R/Z) product_i q_251(v_i t) dt,
I_K=integral_(K_L) product_i q_251(x_i) dx            (34)
```

have exactly the same surviving index set, and

```text
I_line=I_K.                                          (35)
```

Every coordinate character on both averaging spaces pushes Haar to Haar.
The product telescope in (26) therefore gives

```text
|I_line-measure(S(v))|<=13 eta_251,                  (36)

|I_K-measure_(K_L)(K_L intersection J^13)|
 <=13 eta_251.                                       (37)
```

Equations (24), (31), and (35)--(37) would imply

```text
15/343<=26 eta_251
       <26(21/12500)
       =546/12500.                                   (38)
```

But exact subtraction gives

```text
15/343-546/12500
 =111/2143750
 >0.                                                  (39)
```

This contradiction proves:

> **Height-500 rank-three harvest.** Every distinct positive
> thirteen-speed row with zero safe measure satisfies
>
> ```text
> dim_Q W_500(v)>=3.                                  (40)
> ```

In particular, the support-at-most-three mod-`14` codeword in THM-2178 is no
longer a terminal branch. The ambient torus has the uniform floor (24)
whether or not that sparse digit occurs.

## 5. Two bounded relations on every deletion

Choose independent relations

```text
a,b in Lambda(v),             ||a||_infinity,||b||_infinity<=105,
c in Lambda(v)\span_Q{a,b},   ||c||_infinity<=500.    (41)
```

Fix a coordinate `i`. If

```text
a_i=b_i=0,
```

then `a,b` already give two independent height-`105` relations on the
deletion of coordinate `i`.

Otherwise choose

```text
r in {a,b} with r_i!=0
```

and let `s` be the other height-`105` relation. Define

```text
u=s_i r-r_i s,
w=c_i r-r_i c.                                       (42)
```

Then

```text
u_i=w_i=0.                                           (43)
```

The basis `r,s,c` is independent, and `r_i!=0`; comparing its `s` and `c`
coefficients shows that `u,w` are independent. Their heights obey

```text
||u||_infinity<=2*105^2=22050,
||w||_infinity<=2*105*500=105000.                    (44)
```

Thus:

> **Deletion rank-two corollary.** For every coordinate `i`, the twelve-speed
> deletion `v\{v_i}` has relation rank at least two inside coefficient height
> `105000`.                                               (45)

This strengthens THM-2169's one-relation-on-every-deletion theorem along the
rank axis, although its first relation has a much smaller height `1247`.

## 6. Scope and exact referee

The conclusions are structural reductions, not LRC(14):

1. rank three does not bound the speed magnitudes;
2. THM-2171's algebraic relation pump still mixes the safe phase predicate;
3. the deletion relations in (45) do not supply a uniform twelve-speed
   stability theorem; and
4. the safe-cube floor belongs to the ambient relation torus, not necessarily
   to the original one-parameter orbit.

The companion uses only exact integer and `Fraction` arithmetic. It:

1. checks (7)--(8) at every danger count `0,...,13`;
2. verifies the exact floors (11), (15), (18), and (22)--(23);
3. independently compares the closed Jackson formula with coefficient
   convolution at every one of the `501` modes for `N=251`;
4. checks (30), (38)--(39), and the hostile adjacent `N=250` ledger; and
5. retains the elementary Fejer/log route as an independent weaker control.

Normal and optimized Python executions are required to agree with the frozen
transcript. QED.
