---
id: THM-2242
title: "Tournament complement transport and knot-kernel Green rigidity"
status: >
  PROVED + HOSTILE-AUDITED. In THM-2221's pinned equal-block tournament
  model, the maximum whole-block exterior response is exactly total possible
  reversal minus an integer Wasserstein distance between the core-pattern
  histogram and its complement. Exterior antipodes exist exactly when every
  pattern and its complement have equal multiplicity; Frobenius--Koenig,
  Hall, the permanent, optimal transport, and the Kantorovich dual give
  equivalent exact certificates and a quantitative imbalance deficit. The
  quotient loses the labelled core kernel G and therefore does not classify
  the full tournament distance. In every jointly nonexpansive commutative
  metric monoid, the continuation-kernel map is an isometric embedding. For knots, Schubert
  prime decomposition then makes every Green class in the kernel image a
  singleton, the unknot the only idempotent and regular element, and prime
  knots exactly the semigroup atoms. The l-infinity kernel embedding alone
  does not imply partial-cube structure.
source: codex-2026-07-25-complement-transport-green-rigidity
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2235-response-antipode-barriers-for-lrc-sheets-tournament-cycles-and-knot-kernels
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
external:
  - "Horst Schubert, Die eindeutige Zerlegbarkeit eines Knotens in Primknoten, Sitzungsberichte der Heidelberger Akademie der Wissenschaften (1949)."
---

# THM-2242 -- complement transport and kernel rigidity

The two parts share one useful correction. A continuation or permutation
response is faithful only when its residual cost is retained. In the
tournament problem that residual is an integer transport deficit plus the
core kernel. In the knot problem the full continuation kernel is already an
isometric copy of the original metric monoid.

## 1. Tournament pattern histograms

Use exactly the pinned equal-order setting of THM-2221. Thus `C` is the core
label set, every core block has order `N>=1`, and the nonempty exterior
ledger is

```text
mu:{0,1}^C -> N,              M=sum_w mu(w)>0.        (1)
```

Let

```text
W=supp(mu).
```

For each core coordinate `c`, form its column pattern across the positive
ledger:

```text
v_c=(w(c))_(w in W) in {0,1}^W.                     (2)
```

Put

```text
n_a=#{c in C:v_c=a},
n^c_a=n_(bar a),                                    (3)
```

where `bar a` is the bitwise complement. Define the weighted Hamming metric

```text
d_mu(a,b)=sum_(w in W) mu(w)|a_w-b_w|.              (4)
```

The total masses of `n` and `n^c` are both `|C|`.

For two nonnegative integer histograms `r,s` of equal mass, let

```text
W_(d_mu)(r,s)
 =min_T sum_(a,b) T_(a,b)d_mu(a,b),                 (5)
```

where `T` ranges over nonnegative integer matrices with row margins `r` and
column margins `s`. This is the unnormalized integer Wasserstein cost.

## 2. Exact exterior maximum

For a whole-block permutation `sigma in Sym(C)`, let

```text
E_mu(sigma)
 =N sum_(c in C) D_mu(c,sigma(c))                   (6)
```

be THM-2221's exterior reversal term. Then

```text
max_(sigma in Sym(C)) E_mu(sigma)
 =N[M|C|-W_(d_mu)(n^c,n)].                          (7)
```

### Proof

A permutation induces the integer coupling

```text
T_(a,b)
 =#{c:v_c=a, v_(sigma(c))=b}                        (8)
```

with both margins `n`. Conversely, every integer self-coupling of `n` is
realized by a permutation: split each pattern class according to the row of
`T` and biject the pieces to the corresponding target classes.

By the definitions of `D_mu` and `v_c`,

```text
E_mu(sigma)/N=sum_(a,b)T_(a,b)d_mu(a,b).            (9)
```

For every pair of patterns,

```text
d_mu(a,b)=M-d_mu(bar a,b).                          (10)
```

The total mass of `T` is `|C|`. Reindexing the rows of (9) by
`a -> bar a` changes the row margin from `n` to `n^c`, so maximizing (9) is
exactly `M|C|` minus the minimum (5). This proves (7). QED.

The exact Kantorovich dual is

```text
W_(d_mu)(n^c,n)
 =max_phi sum_a phi(a)(n^c_a-n_a),                  (11)
```

over real potentials satisfying

```text
|phi(a)-phi(b)|<=d_mu(a,b).                         (12)
```

Thus a feasible potential is a rigorous exterior-response deficit
certificate, while min-cost flow supplies the sharp value.

## 3. Antipodes, the permanent, and the quantitative deficit

An exact exterior antipode has

```text
E_mu(sigma)=NM|C|.                                  (13)
```

Equations (7) and (5) give the equivalent conditions

```text
an exterior antipode exists
 iff W_(d_mu)(n^c,n)=0
 iff n_a=n_(bar a) for every a.                     (14)
```

The last equivalence uses that `d_mu` is a genuine metric on the pattern
set `{0,1}^W`: all coordinates in `W` have positive weight.

There is a literal Frobenius--Koenig/Hall certificate. Define the
`|C| x |C|` zero-one matrix

```text
B_(c,d)=1_[v_d=bar(v_c)].                            (15)
```

Then

```text
per(B)!=0
 iff B has a perfect matching
 iff n_a=n_(bar a) for every a.                     (16)
```

Because (15) is a disjoint union of complete bipartite blocks between
complementary pattern classes, the general Frobenius--Koenig obstruction
collapses exactly to the histogram equalities in (14).

When (14) holds, the number of exterior antipode permutations is

```text
product_a n_a!
 =product_({a,bar a}) (n_a!)^2.                     (17)
```

Let

```text
mu_min=min_(w in W)mu(w).                            (18)
```

Every transport unit sent between distinct patterns costs at least
`mu_min`, and the minimum amount that must move is the total-variation
imbalance. Hence

```text
W_(d_mu)(n^c,n)
 >=mu_min/2 sum_a |n_a-n_(bar a)|.                  (19)
```

Combining (7) and (19) gives

```text
max_sigma E_mu(sigma)
 <=N[M|C|-mu_min/2 sum_a|n_a-n_(bar a)|].           (20)
```

### Information lost by the histogram

The map

```text
(labelled columns v_c) -> histogram n              (21)
```

preserves the unconstrained maximum exterior response and its antipode
count. It forgets which pattern is attached to which core vertex, and
therefore forgets THM-2221's core cost

```text
G(NP_sigma).                                        (22)
```

This loss is real. Put the singleton ledger word `0011` on the four
vertices of a transitive tournament. Its pattern histogram is complement
balanced, so exterior antipodes exist. The transitive four-tournament has
trivial automorphism group, however, so none of those nontrivial antipodes
has zero core cost. The full pinned response still requires the joint
optimization

```text
min_sigma [G(NP_sigma)+E_mu(sigma)]                 (23)
```

or its split-transport extension. Histogram balance is not a classifier of
the whole tournament.

## 4. Isometry of continuation kernels

Let `(M,+,0,d)` be any commutative monoid with a metric satisfying joint
nonexpansivity:

```text
d(a+c,b+d)<=d(a,b)+d(c,d).                          (24)
```

For `x in M`, define its two-ended continuation kernel

```text
P_x(a,b)=d(x+a,b).                                   (25)
```

Then

```text
||P_x-P_y||_infinity=d(x,y).                        (26)
```

Indeed, reverse triangle inequality and simultaneous translation give

```text
|d(x+a,b)-d(y+a,b)|
 <=d(x+a,y+a)
 <=d(x,y).                                          (27)
```

Equality is attained at `(a,b)=(0,x)`, because

```text
P_x(0,x)=0,              P_y(0,x)=d(y,x).           (28)
```

Thus the continuation representation is injective and isometric, not only
faithful. With THM-2176's min-plus convolution convention it is also a
monoid embedding:

```text
P_x tensor P_y=P_(x+y).                             (29)
```

## 5. Gordian catalysts as translated kernel distance

Apply Section 4 to the connected-sum monoid of knots and Gordian distance.
For knots `K,L`, define

```text
d_cat(K,L)=inf_J d_G(K#J,L#J).                      (30)
```

Equations (26) and (29) give the exact kernel formula

```text
d_cat(K,L)
 =inf_J ||P_K tensor P_J-P_L tensor P_J||_infinity. (31)
```

Thus catalytic contraction is distance contraction along common principal
right translates inside one isometric kernel monoid. It is not an inverse:
THM-2235 already proves that the only min-plus unit is the unknot.

## 6. Green rigidity of the knot-kernel image

Let

```text
I_K={P_(K#J):J a knot}                              (32)
```

be the principal ideal generated by `P_K`. Schubert prime decomposition
makes the connected-sum monoid a free commutative monoid on prime knots.
Using injectivity from (26),

```text
I_K subset I_L
 iff K=L#C for some knot C.                          (33)
```

The forward direction follows because `P_K in I_L`; injectivity gives
`K=L#C`. The reverse direction is immediate by adding an arbitrary `J`.

Consequently distinct knots generate distinct principal ideals. Since the
kernel monoid is commutative,

```text
every Green L,R,H,D,J class is a singleton.          (34)
```

The same free-monoid argument gives the complete elementary structure:

```text
only idempotent:       P_U,
only regular element: P_U,
semigroup atoms:       P_K for K prime.              (35)
```

For example, idempotence would give `K#K=K`, and regularity would give
`K=K#J#K`; prime multiplicities force `K=U`. The atom statement is exactly
irreducibility under connected sum.

These classifications do not make unknotting number additive. They locate
nonadditivity in metric contraction along the rigid free-monoid order, not
in hidden algebraic identifications of kernel elements.

## 7. Partial-cube boundary

Equation (26) is a Kuratowski-type `l_infinity` embedding and supplies no
partial-cube conclusion by itself. As a minimal hostile control, take
`M=Z/3Z` with its discrete translation-invariant metric. It satisfies
(24), so (26) holds, but its distance-one graph is `K_3`, which is not
bipartite and therefore not a partial cube.

Any genuine partial-cube claim for a knot or tournament response graph must
construct convex cuts or Djokic--Winkler theta classes. The continuation
kernel alone supplies neither. QED.
