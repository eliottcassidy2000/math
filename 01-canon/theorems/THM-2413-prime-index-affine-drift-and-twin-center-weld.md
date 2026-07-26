---
id: THM-2413
title: "Prime-index affine drift and twin-center summand--multiplicand weld"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE; AWAITING INDEPENDENT HOSTILE
  AUDIT. Operation cospan transitivity is witness composition. After
  deleting equal-parent fibers, the additive and multiplicative local
  transitivity-defect counts at x are respectively x-1 and tau(x)-2;
  hence primes are exactly the positive x>=2 with no multiplicative
  defect. Prime valuations identify multiplication with
  finite-support coordinatewise addition, while ordinary addition has
  exact equal-valuation carry walls. A014574 is equivalently the set
  of centers m whose neighbors form a Boolean divisor diamond for
  m^2-1, or whose sum/product discriminant is four. It is also the set
  of adjacent plateaux of the fixed slope-two prime drift p_k-2k.
  This local plateau is not a plateau of A373813 and says nothing
  about slopes in an optimal line cover. Brun convergence and the
  2026 prime-line asymptotics are CITED, not reproved.
source: codex-2026-07-26-prime-index-drift
depends_on:
  - THM-362-natural-operation-graph-shadows
related:
  - THM-361-product-sum-defect-normal-form
  - THM-446-sidon-additive-relation-ladder-and-the-dyadic-hard-core-of-erdos-64
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - HYP-3003-summand-multiplicand-farey-basis-merge
external: >
  OEIS A014574 and A373813; Kominers--Mrazovic--Pomerance--Sole,
  "Lines in the Prime Number Graph" (2026); Rybin--Zhang--Luo,
  "XX^t Can Be Faster", arXiv:2505.09814v2.
script: 04-computation/prime_drift_twin_center_weld_thm2413.py
output: 05-knowledge/results/prime_drift_twin_center_weld_thm2413.out
script_sha256: 13330b8db5c13b3a789278997f3edfe64a00d60cccae68d9055fdf538482f2ab
output_sha256: 674c3e67a9578b2353a4bd4d2e953a50804cdeac0328d146fe0cd2eb85590254
hash_basis: working-tree bytes (LF)
cite_by_filename: true
---

# THM-2413 -- prime drift and the twin-center operation weld

**PROVED + VERIFIED-EXACT CANDIDATE; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The arrows

```text
x -> z and y -> z  iff  z=x+y,
x -> z and y -> z  iff  z=xy
```

are two-input operation cospans before they are graphs. Their simple
shadows are the total order and divisibility, as THM-362 proves.
Associativity explains why both shadows are transitive: two hidden
witnesses compose by addition or multiplication.

Keeping one extra bit of the forgotten fiber -- whether the two parents
were equal -- exposes a sharper difference. The local failures of
transitivity then number `x-1` additively but `tau(x)-2`
multiplicatively. Primes are exactly the zero-defect multiplicative
vertices. This is a precise sense in which a diagonal irregularity in
the dense summand cospan becomes arithmetic structure in the sparse
multiplicand cospan.

The twin-prime midpoint sequence A014574 supplies an unusually tight
weld between the two operations. If `m-1` and `m+1` are prime, then

```text
(m-1)+(m+1)=2m,
(m-1)(m+1)=m^2-1,
```

and the product has a four-element Boolean divisor interval. In
prime-index coordinates the same event is an adjacent plateau of
`p_k-2k`. This is local fixed-slope geometry, not a characterization
of the global minimum line-cover sequence A373813.

## 1. Operation cospans and transitivity as witness composition

On the positive integers define the labeled fibers

```text
C_+(z)     ={(x,y): x+y=z},
C_times(z) ={(x,y): xy=z}.                                         (1)
```

Forgetting the second parent gives the loopless shadows

```text
x R_+ z      iff  exists c>=1, z=x+c,
x R_times z  iff  exists c>=1, z=xc and z>x.                       (2)
```

Thus

```text
x R_+ z      iff x<z,
x R_times z  iff x|z and x<z.                                      (3)
```

If `x R_+ z` has witness `c` and `z R_+ w` has witness `d`,
then

```text
w=(x+c)+d=x+(c+d),                                                  (4)
```

so `c+d` witnesses `x R_+ w`. Similarly,

```text
w=(xc)d=x(cd),                                                      (5)
```

so `cd` witnesses `x R_times w`. Transitivity is therefore not an
accident of the simple shadows. It is associativity seen after one
input of each operation has been hidden.

The information loss is radically different. The additive shadow
forgets to total order and has successor Hasse edges. The
multiplicative shadow forgets to the divisor poset and has Hasse edges

```text
x -> xp,                         p prime.                            (6)
```

On `{1,...,80}` the additive shadow has `3160` edges. The
multiplicative shadow has `288` when the unit source `1` is retained
and `209` when that source row is deleted. Both conventions occur in
the historical repo; they must not be conflated.

## 2. Equal-parent deletion and the prime defect theorem

Retain from each cospan the information that its two parents are
distinct. Define

```text
x R_+^neq z
  iff exists c>=1 with c!=x and z=x+c
  iff x<z and z!=2x,                                                (7)

x R_times^neq z
  iff exists c>=2 with c!=x and z=xc
  iff x|z, x<z, and z!=x^2.                                        (8)
```

Fix `x>=2`. A two-step additive path rooted at `x` has the form

```text
x -> x+c -> x+c+d.                                                  (9)
```

If its two displayed edges survive (7), its composite direct edge
fails exactly when

```text
x+c+d=2x,                         equivalently c+d=x.                (10)
```

Every ordered positive solution of (10) automatically makes both
steps in (9) legal: `c!=x` and `d!=x+c`. Hence the exact number of
rooted additive transitivity defects is

```text
#{(c,d)>=1: c+d=x}=x-1.                                             (11)
```

The multiplicative calculation is parallel:

```text
x -> xc -> xcd                                                      (12)
```

has no composite direct edge in (8) exactly when

```text
xcd=x^2,                           equivalently cd=x.                (13)
```

The legal witnesses have `c,d>=2`. They are all ordered divisor
pairs except `(1,x)` and `(x,1)`, so the exact number is

```text
#{(c,d)>=2: cd=x}=tau(x)-2.                                         (14)
```

Consequently, for every `x>=2`,

```text
x is prime
  iff tau(x)-2=0
  iff R_times^neq has no two-step transitivity defect rooted at x.  (15)
```

This does not say that addition causes primes. It says that the same
diagonal-deletion operation turns the uniform additive composition
law into the divisor-sensitive multiplicative one; primality is the
zero set of the latter defect statistic.

## 3. Multiplication is addition on a nonuniform prime basis

Unique factorization gives the monoid isomorphism

```text
v:N_(>=1) -> direct-sum_p N,
v(n)=(v_p(n))_p,
v(ab)=v(a)+v(b).                                                     (16)
```

Thus the positive integers under multiplication form the free
commutative monoid on the primes. Divisibility is coordinatewise
order, (6) adds one prime basis vector, and

```text
Omega(n)=sum_p v_p(n)                                               (17)
```

is the unweighted Hasse rank. The ordinary size embedding is instead

```text
log n=sum_p v_p(n) log p.                                           (18)
```

So multiplication has uniform combinatorics in infinitely many prime
coordinates but no uniform metric base: different axes carry weights
`log p`.

Ordinary addition is nonlinear in these coordinates. If

```text
a=p^r u,  b=p^s v,                 p does not divide uv,             (19)
```

then

```text
v_p(a+b)=min(r,s),                          r!=s,
v_p(a+b)=r+v_p(u+v),                        r=s.                     (20)
```

The equal-valuation wall is exactly where carries and cancellation
occur. For `p=2`, `u` and `v` are odd on that wall, so the valuation
strictly rises. Formula (20), rather than a causal slogan, is the
precise way additive diagonal irregularities propagate through the
prime-coordinate model of multiplication.

At the level of functions, the same distinction is

```text
(f *_C g)(n)=sum_(a+b=n) f(a)g(b),          Cauchy convolution,
(f *_D g)(n)=sum_(ab=n) f(a)g(b),           Dirichlet convolution.  (21)
```

Under (16), Dirichlet convolution becomes multivariate Cauchy
convolution on finite-support exponent vectors. The summand and
multiplicand graphs are therefore two convolution geometries, not two
unlabeled pictures of comparable density.

## 4. A014574 as a sum--product discriminant-four weld

Let

```text
T={m>=4: m-1 and m+1 are prime}.                                   (22)
```

This is A014574. For `m>=4`, the following are equivalent:

1. `m in T`;
2. `m^2-1` has the exact divisor set
   `{1,m-1,m+1,m^2-1}`;
3. the divisor interval of `m^2-1` is a Boolean diamond whose two
   middle elements differ by two and sum to `2m`;
4. the unordered pair recovered from

   ```text
   S=2m,  P=m^2-1                                                   (23)
   ```

   has discriminant

   ```text
   S^2-4P=4.                                                        (24)
   ```

Indeed, (1) factors `m^2-1` as the product of the two distinct primes
`m-1,m+1`, proving (2)--(4). Conversely, (2) directly makes those two
neighbors prime. For (3), a positive integer with four divisors is
either a product of two distinct primes or a prime cube. The
prime-cube middle gap can equal two only for
`1,2,4,8`, whose midpoint is `m=3`; the condition `m>=4` excludes this
unique hostile.

The same pair lands in both operation cospans:

```text
m-1 ----+
         +-- sum ----------> 2m
m+1 ----+

m-1 ----+
         +-- product ------> m^2-1 ----(+1)----> m^2.               (25)
m+1 ----+
```

Also `P+1=m^2`, so the product target is one below the square of the
midpoint. Apart from the initial center `m=4`, every center is
divisible by six, because twin primes above three are `-1,+1`
modulo six.

## 5. Reciprocal support

Every `m in T` satisfies the exact identity

```text
1/m
 =1/2(1/(m-1)+1/(m+1))
  -1/[m(m^2-1)].                                                    (26)
```

This follows by putting the first term over denominator `m^2-1`.
Summing (26) over twin centers and writing `B_2` for Brun's constant
in the convention

```text
B_2=sum_(m in T) [1/(m-1)+1/(m+1)]                                 (27)
```

gives

```text
sum_(m in T) 1/m
 =B_2/2-sum_(m in T) 1/[m(m^2-1)].                                 (28)
```

Brun's theorem is an external cited input, so (28) proves that the
reciprocal support of A014574 converges. It does not prove that `T` is
infinite; a finite set would also satisfy the statement.

## 6. Prime-index affine drift and A373813

Let `p_k` be the `k`th prime and

```text
P_k=(k,p_k).                                                        (29)
```

For a fixed slope `a`, define the affine drift

```text
d_k^(a)=p_k-ak.                                                     (30)
```

A set of prime points lies on one line of slope `a` exactly when
their drifts (30) are equal. Equivalently, three indices `i,j,k` are
collinear exactly when their divided differences agree:

```text
(p_j-p_i)/(j-i)=(p_k-p_i)/(k-i).                                   (31)
```

Thus the prime-line problem is a hypergraph cover problem: vertices
are prime indices and every maximal equal-drift fiber, with its slope
retained, is a hyperedge. A373813 is

```text
L(n)=minimum number of affine-line hyperedges covering
     P_1,...,P_n.                                                   (32)
```

Restriction and adjoining a singleton line give

```text
L(n)<=L(n+1)<=L(n)+1,                                               (33)
```

so every increment is zero or one. Calling `p_n` awkward when
`L(n)>L(n-1)`, telescoping gives

```text
L(n)=#{k<=n: p_k is awkward}.                                      (34)
```

The 2026 paper *Lines in the Prime Number Graph* proves the external
asymptotic facts

```text
L(n)=o(n)
```

and convergence of the reciprocal sum over awkward primes, with
quantitative refinements from prime-number-theorem remainder bounds.
Those analytic results are **CITED**, not reproved here.

For the fixed slope `a=2`,

```text
d_(k+1)^(2)-d_k^(2)
 =p_(k+1)-p_k-2.                                                    (35)
```

Therefore

```text
d_(k+1)^(2)=d_k^(2)
  iff p_k,p_(k+1) are twin primes,                                  (36)

m=p_k+1=2k+d_k^(2)+1.                                               (37)
```

So A014574 records the centers of adjacent edges inside the
slope-two drift fibers. The exceptional triple `3,5,7` simply gives
two adjacent edges in one three-point fiber.

Equations (36)--(37) do **not** imply any of the following:

- a plateau of `L(n)`;
- an awkward or nonawkward classification of either endpoint;
- that an optimal A373813 cover uses the slope-two line;
- that adjacent points are the only useful collinearities; or
- infinitude of twin primes.

The exact line-cover prefix through `n=16` is

```text
1,1,2,2,2,3,3,3,3,4,4,4,4,4,4,4,                                 (38)
```

and its increment indices are `1,3,6,10`. This finite check is a
control for (32)--(34), not evidence for an asymptotic pattern.

## 7. Product fibers followed by sum aggregation

Every matrix coefficient has the cospan form

```text
(XX^t)_(ij)=sum_k X_(ik) X_(jk):                                   (39)
```

first create multiplicative atoms indexed by `k`, then aggregate them
in a summand fiber. The `4 x 4` block construction of
Rybin--Zhang--Luo uses eight recursive symmetric calls and 26 general
products; its search first generates product candidates and exact
linear relations, then solves a minimum-cover problem for the target
entries.

This supplies a useful algorithmic pattern for (21), (25), and (32):

```text
labeled product fibers -> exact additive span -> minimum cover.     (40)
```

It proves no prime, twin-prime, tournament, or LRC result. In
particular, collapsing a Farey packet `p/q` to the two scalars `p+q`
and `pq` loses the root packet, order, and phase-owner data required
by HYP-3003.

## 8. Exact companion

Run:

```text
python3 04-computation/prime_drift_twin_center_weld_thm2413.py
python3 -O 04-computation/prime_drift_twin_center_weld_thm2413.py
```

The standard-library integer/Fraction companion checks:

- both unit conventions for the operation shadows through `80`;
- all diagonal-defect counts through `80`;
- `2500` valuation-vector products;
- the first twelve twin centers and every weld equivalence;
- the fixed-slope drift identity on the first forty primes;
- the reciprocal identity exactly;
- the discriminant/diamond converse through midpoint `499`; and
- an exact bitmask dynamic-programming reconstruction of (38).

Both modes reproduce:

```text
05-knowledge/results/prime_drift_twin_center_weld_thm2413.out
```
