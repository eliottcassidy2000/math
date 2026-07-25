---
id: THM-2283
title: "Mixed rank-two safe-torus floors and adaptive scalar rank-three harvest"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every saturated
  rank-two relation torus on one guard and eight ordinary coordinates has
  mixed safe-box Haar mass at least 40/273. When the guard speed is odd,
  as it is on the live scalar LRC(14) section, the floor is 152/1029.
  If the relation plane contains a guard pair or an ordinary pair, the
  respective odd-guard floors are 19/126 and 2060/13377. The proof uses a
  minimal-support-union classification, the exact ordinary and mixed comb
  overlap laws, Hunter trees, and a disjoint-triple factorization at union
  size six. Exact squared-Fejer comparison first crosses the uniform live
  floor at bandwidth 52, with true adjacent failure at 51. Consequently
  every live scalar cover has relation rank at least three by coefficient
  height 594, adaptively by heights 102, 198, or 594 according to the
  THM-2274 branch. The fixed original-row section has rank at least three
  by heights 204, 396, or 1188. No scalar profile is excluded, no owner or
  exact Fourier atom is selected, and LRC(14) remains open. THM-2295
  later supersedes the uniform rank count by proving rank five at height
  196; the mixed-torus floors and support-conditioned certificates here
  remain separate structural inputs.
source: codex-2026-07-25-mixed-safe-torus-rank-three
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2185-rank-two-safe-cube-floor-and-height-500-rank-three-harvest
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2274-mixed-scalar-relative-rank-harvest-and-adaptive-pair-crossing
related:
  - THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
  - THM-2295-weighted-basis-safe-floor-and-scalar-rank-five-harvest
script:
  - 04-computation/lrc14_mixed_rank_two_safe_torus_rank_three_thm2283.py
  - 04-computation/lrc14_mixed_rank_two_hunter_referee_thm2283.py
output:
  - 05-knowledge/results/lrc14_mixed_rank_two_safe_torus_rank_three_thm2283.out
  - 05-knowledge/results/lrc14_mixed_rank_two_hunter_referee_thm2283.out
script_sha256:
  - f2045b18dad1b1106db112ccbc036dfd54f580cb8b5e8dd3596e60805b729586
  - d30b12f3ba6de9e0e1f6f535f624cd324e99240c6ba9f2490b8547ddfc6719a4
output_sha256:
  - 3a378e972fa78b3c06151a705dde1d9840d75678003736432ca03c1e387c1557
  - 6d32e6d361474b5dae66ba063c11a24374b83fb3e64b0b9de07fb306f13d42f2
hash_basis: working-tree bytes (LF)
---

# THM-2283 -- mixed rank-two safe-torus floors and adaptive rank three

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2274 now gives two independent relations on the live nine-coordinate
scalar section by height `594`, with sharper branch heights `35`, `198`,
and `594`. The missing implication was:

```text
rank two on the null scalar orbit
  -> the corresponding rank-two torus has positive mixed safe mass
  -> a bounded positive Fourier comparison sees the discrepancy
  -> a third scalar relation.                                  (1)
```

The old version of this theorem obtained only `72/16807` by pricing one
Haar pair in a six-coordinate core. The load-bearing refinement is to
choose two sparse generators with **minimal support union**. That choice
almost determines the graph of non-Haar pair projections and makes the
six-coordinate branch split into two independent triple subtori.

## 1. Relation tori and pair projections

Let

```text
w=(w_0,w_1,...,w_8) in Z^9,          w_i!=0,          (2)

Lambda(w)={m in Z^9:m.w=0}.
```

Let `L subset Lambda(w)` be saturated of rank two and put

```text
K_L={x in (R/Z)^9:l.x=0 mod 1 for every l in L}.      (3)
```

Saturation gives

```text
ann(K_L)=L.                                           (4)
```

For a coordinate set `T`, character duality gives

```text
ann(pi_T(K_L))=L intersection Z^T.                    (5)
```

No member of `L` has support one: `c e_i in Lambda(w)` would imply
`c w_i=0`.

On `K_L`, define one guard danger and eight ordinary dangers:

```text
D_0={x:||x_0||<=1/7},                 p_0=2/7,

D_i={x:||x_i||<1/14},                 p_i=1/7,
                                      1<=i<=8.         (6)
```

Every coordinate marginal is Haar, so endpoints are null. The mixed safe
box is the atom on which none of the nine dangers occurs.

### Pair-projection lemma

For every ordinary pair `i,j`,

```text
mu_(K_L)(D_i intersection D_j)>=1/91.                 (7)
```

For every guard/ordinary pair,

```text
mu_(K_L)(D_0 intersection D_j)>=1/91.                 (8)
```

If `w_0` is odd, the sharper live bound is

```text
mu_(K_L)(D_0 intersection D_j)>=1/42.                 (9)
```

To prove this, put

```text
M_ij=L intersection Z^{i,j}.
```

The lattice `M_ij` is saturated in `Z^{i,j}` and has rank zero or one.
Rank two would contain a nonzero support-one relation. If its rank is
zero, (5) makes the pair projection the full two-torus, and the
intersection mass is the product `1/49` or `2/49`.

If its rank is one, its primitive generator has two nonzero coefficients.
Because it lies in `Lambda(w)`, it is

```text
+/-(w_j/g,-w_i/g),           g=gcd(|w_i|,|w_j|).      (10)
```

Thus the pair projection, with Haar measure, is the connected orbit

```text
t |-> ((w_i/g)t,(w_j/g)t),                            (11)
```

up to harmless signs. For an ordinary pair, THM-1166 gives (7). Its
statement assumes distinct positive speeds; equal absolute speeds are a
separate easier case, with overlap `1/7`.

For (8), the guard danger contains the radius-`1/14` danger at its first
frequency, so (7) applies to that contained pair. If `w_0` is odd, then
`w_0/g` is odd, and THM-2080's exact unequal-comb theorem applied to
(11) gives (9). Equal absolute speeds again give the larger mass `1/7`.

This proof is about literal symmetric pair overlaps. It introduces no
tournament orientation.

## 2. The cubic minorant and the one-sparse-line cases

Let `N(x)` be the number of active dangers. Use the THM-2185 polynomial

```text
P(n)=-(n-1)(n-3)(n-4)/12

    =1-n+(5/6)binomial(n,2)-(1/2)binomial(n,3).       (12)
```

For every integer `0<=n<=9`,

```text
P(n)<=1_(n=0).                                        (13)
```

If every projection onto at most three coordinates is Haar, its first
three factorial moments are independent. Their elementary symmetric sums
are

```text
sum_i p_i=10/7,

sum_(i<j)p_i p_j=44/49,

sum_(i<j<k)p_i p_j p_k=16/49.                         (14)
```

Hence the three-wise baseline is

```text
B=1-10/7+(5/6)(44/49)-(1/2)(16/49)
 =23/147.                                             (15)
```

Call a nonzero member of `L` sparse when its support has size at most
three. Support one is impossible. The rational span of the sparse members
has dimension zero, one, or two.

### No sparse relation

Every projection onto at most three coordinates is Haar, so

```text
mu_(K_L)(N=0)>=23/147.                                (16)
```

### One sparse line of support three

Let its primitive generator have support `A`, `|A|=3`, and let

```text
q=mu_(K_L)(intersection_(i in A)D_i).
```

Every pair projection is Haar and every triple except `A` is Haar.
The triple intersection lies inside an ordinary pair, so `q<=1/49`.
Only the cubic factorial moment changes. If `A` is ordinary,

```text
E P(N)
 >=B-(1/2)(1/49-1/343)
 =152/1029.                                           (17)
```

If `A` contains the guard,

```text
E P(N)
 >=B-(1/2)(1/49-2/343)
 =307/2058.                                           (18)
```

### One sparse line of support two

Let `A={i,j}`, put

```text
q=mu(D_i intersection D_j),          delta=q-p_i p_j.
```

For every `k` outside `A`, any relation supported on `A union {k}` is
sparse and hence proportional to the line on `A`. Thus (5) identifies the
projection on `A union {k}` with the pair subtorus times a free circle.
Therefore the coefficient of `delta` in `E P(N)` is

```text
11/42,        A ordinary;

1/3,          A contains the guard.                  (19)
```

The pair-projection lemma now matters quantitatively. For an ordinary
pair, `q>=1/91`, so

```text
E P(N)
 >=B+(11/42)(1/91-1/49)
 =2060/13377.                                         (20)
```

For a guard pair, (8) gives the general floor

```text
E P(N)
 >=B+(1/3)(1/91-2/49)
 =40/273.                                             (21)
```

When `w_0` is odd, (9) sharpens this to

```text
E P(N)
 >=B+(1/3)(1/42-2/49)
 =19/126.                                             (22)
```

## 3. Two sparse directions: the minimal-union graph

Suppose the sparse members span `L` over `Q`. Choose independent sparse
relations `a,b` for which

```text
U=supp(a) union supp(b),             m=|U|            (23)
```

is minimal. Then

```text
3<=m<=6.                                              (24)
```

The lower bound excludes two independent relations on the same two
coordinates, since eliminating one coordinate would create a support-one
relation. Because `a,b` span `L_Q`, every member of `L` is supported on
`U`, and

```text
K_L=K_(L|U) times (R/Z)^(9-m).                        (25)
```

Let `G` be the graph on `U` whose edges are the non-Haar pair
projections. By (5), an edge is exactly a pair-supported relation in `L`.
Minimality in (23) forces:

```text
m=3: G is contained in K_3;
m=4: G is a matching;
m=5: G has at most one edge;
m=6: G is empty.                                      (26)
```

Indeed, two adjacent edges are independent sparse relations with union
three. Two disjoint edges are independent with union four. These exclude
adjacent edges when `m>=4` and two edges when `m>=5`. Finally, if `m=6`
and one pair relation `c` existed, at least one of `a,b` would be
independent of `c`, producing an independent sparse pair with union at
most `2+3=5`.

### Hunter trees

For events indexed by the vertices of any tree `T`, the pointwise
inequality

```text
1_(union_i D_i)
 <=sum_i 1_(D_i)-sum_({i,j} in E(T))1_(D_i intersection D_j)   (27)
```

follows because, for a nonempty active vertex set `S`,

```text
|S|-|E(T[S])|=#components(T[S])>=1.
```

For `m=4,5`, the complement of the graph in (26) contains an all-Haar
spanning tree. If the guard lies in `U`, the tree can use at least `m-2`
guard edges and one ordinary edge. For `m=6`, every pair is Haar.

At `m=3`, use any tree and (7)--(9). Combining its Hunter credit with
the exact free-coordinate factors gives the following complete table.
The last column retains the stronger `1/42` credit whenever a chosen tree
uses a non-Haar guard edge.

| `m` | guard outside `U` | guard in `U` | odd guard in `U` |
|---:|---:|---:|---:|
| 3 | `2099520/10706059` | `1912896/10706059` | `155520/823543` |
| 4 | `155520/823543` | `147744/823543` | `149040/823543` |
| 5 | `19440/117649` | `2592/16807` | `18360/117649` |

For example, when the guard is in `U`, the `m=4,5` core floors are

```text
1-2/7-(m-1)/7+(2m-3)/49,
```

and the remaining `9-m` ordinary coordinates contribute
`(6/7)^(9-m)`.

### Why the six-coordinate core factors

When `m=6`, both `a,b` have support three and their supports are disjoint.
Let `a_0,b_0` be the primitive generators of their rational lines.
Since `L` is saturated,

```text
L=(Q a_0 intersection Z^U)
  direct_sum
  (Q b_0 intersection Z^U)
 =Z a_0 direct_sum Z b_0.                             (28)
```

Thus the core torus is a product of two three-coordinate one-relation
subtori. A triple of ordinary dangers has union mass at most `3/7`; a
triple containing the guard has union mass at most `4/7`. Therefore the
two core safe masses are at least

```text
16/49,        guard outside U;

12/49,        guard in U.                             (29)
```

After multiplying the free-coordinate safe masses, the full floors are

```text
(16/49)(5/7)(6/7)^2=2880/16807,

(12/49)(6/7)^3=2592/16807.                            (30)
```

This factorization is what the former one-Haar-pair proof missed.

## 4. The safe-torus floors

Comparing (16)--(22), the table, and (30) gives three useful statements.

> **General mixed floor.** Under only (2),
>
> ```text
> mu_(K_L)(mixed safe box)>=40/273.                   (31)
> ```

> **Odd-guard mixed floor.** If `w_0` is odd,
>
> ```text
> mu_(K_L)(mixed safe box)>=152/1029.                 (32)
> ```

> **Pair-conditioned odd-guard floors.** If `w_0` is odd and `L`
> contains a pair-supported relation, then
>
> ```text
> guard pair in L     -> safe mass>=19/126;
>
> ordinary pair in L  -> safe mass>=2060/13377.       (33)
> ```

For (31), the minimum displayed bound is the general guard-pair value
`40/273`. For (32), oddness raises that value, and the minimum becomes the
ordinary support-three value `152/1029`. The exact gaps above this live
floor include

```text
307/2058-152/1029=1/686,

2060/13377-152/1029=4/637,

19/126-152/1029=19/6174,

2592/16807-152/1029=328/50421.                        (34)
```

No sharpness claim is made. These are universal certified floors.

## 5. Exact squared-Fejer boundary

Use THM-2185's normalized squared-Fejer kernel

```text
J_N=F_N^2/integral F_N^2
```

and convolve each safe-interval indicator with `J_N`. The resulting
polynomial has degree `2N-2`, takes values in `[0,1]`, and has circle
`L^1` error at most

```text
eta_N
 =1/2-[4/(pi^2 C_0)]
       sum_(1<=k<=2N-3, k odd) C_k/k^2,              (35)

C_0=N(2N^2+1)/3,

C_k=(4N^3-6Nk^2+2N+3k^3-3k)/6,       0<=k<=N,
C_k=((2N-k)^3-(2N-k))/6,              N<k<=2N-2.
```

The same error applies to both interval lengths because translating either
indicator by `y` changes it in `L^1` by at most `2||y||`.

THM-2274's independent referee certifies

```text
103993/33102<pi<355/113.                              (36)
```

Substituting these rational bounds into (35), and reconstructing every
`C_k` independently by convolution, gives the true adjacent boundary

```text
eta_51>76/9261=(152/1029)/18,

eta_52<76/9261,                                       (37)

152/1029-18 eta_52>1/600.
```

Thus bandwidth `52`, degree `102`, is the first accepted squared-Fejer
certificate for the uniform live floor. This is a boundary of this
specific comparison, not an optimality theorem for positive
approximations.

The pair-conditioned floors similarly give

```text
guard pair:
  eta_50>(19/126)/18,       eta_51<(19/126)/18;

ordinary pair:
  eta_49>(2060/13377)/18,   eta_50<(2060/13377)/18.   (38)
```

Their first accepted degrees are `100` and `98`.

## 6. Adaptive scalar rank three

Use THM-2274's live scalar row

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3) in Z_(>0)^9,         (39)
```

where `H` is odd and the mixed safe event is null. Write

```text
W_K^*
 =span_Q{m in Z^9:m.w_*=0, ||m||_infinity<=K}.        (40)
```

THM-2274 starts from THM-2275's primitive relation `p` and proves

```text
dim W_35^*>=2,    if |supp(p)|>=3;

dim W_198^*>=2,   if |supp(p)|=2 and H is in supp(p);

dim W_594^*>=2,   if |supp(p)|=2 and H is not in supp(p).       (41)
```

Set

```text
K=102, 198, or 594
```

in the three branches. Let `r` be the independent relation supplied by
THM-2274 in the relevant branch, and let `L` be the saturation of
`Z p+Z r` in `Z^9`. Thus the pair-supported branches visibly retain the
particular guard or ordinary pair `p` needed for (33). Since
`Lambda(w_*)` is saturated,

```text
L subset Lambda(w_*),                 rank(L)=2.       (42)
```

Suppose `dim W_K^*=2`. Use the degree `102` uniform polynomial in the
first branch, the degree `100` guard-pair polynomial in the second, and
the degree `98` ordinary-pair polynomial in the third.

The orbit integral of the nine-factor polynomial retains exactly the
frequencies in `Lambda(w_*)`. The `K_L` integral retains exactly those in
`L`. Every product frequency has coordinate height at most its displayed
degree, hence at most `K`; the assumption `dim W_K^*=2` puts it in
`L_Q`, and saturation gives

```text
L_Q intersection Z^9=L.                              (43)
```

Conversely, every product frequency in `L` lies in `Lambda(w_*)`.
Therefore the two polynomial integrals are equal.

Each of the two nine-factor telescopes costs at most `9 eta_N`. The orbit
safe event has measure zero, while (32) or (33) gives a torus safe mass
strictly larger than `18 eta_N` by (37)--(38). This contradiction proves

```text
|supp(p)|>=3:
  dim W_102^*>=3;

|supp(p)|=2, guard pair:
  dim W_198^*>=3;

|supp(p)|=2, ordinary pair:
  dim W_594^*>=3.                                    (44)
```

In particular every one of the `165` live scalar profiles has

```text
dim W_594^*>=3.                                      (45)
```

This includes the `45` boundary/repeated profiles outside the interior
THM-2266 atlas.

## 7. Fixed-section lift and scope

THM-2203 lifts scalar relation coefficients to the fixed
nine-coordinate section of the original thirteen-speed row by

```text
(r_H,r_1,...,r_8)
 |->(2r_H,r_1,...,r_8),                              (46)
```

with zeros on the four forgotten coordinates. It is integral and
injective, so (44) becomes fixed-section rank three by respective heights

```text
204,             396,             1188.              (47)
```

The theorem improves the previous uniform scalar/fixed heights
`3540/7080` to `594/1188`. It does not exclude a scalar profile, bound the
speed magnitudes, choose a blocker owner, preserve expiration ancestry, or
select one exact energetic integer Fourier atom. The equality of finite
polynomial integrals is a contradiction device; it does not assert that
the original orbit is dense in `K_L`.

THM-2295 subsequently proves the stronger uniform rank-five conclusion at
height `196` by the weighted coordinate-basis mechanism. Thus (45) is no
longer the best uniform rank count. What remains specific here is the
rank-two torus classification, the general/odd/pair-conditioned safe
floors, the true `51/52` Jackson boundary, and the degree-`102` rank-three
certificate on the support-at-least-three branch. These should not be
reinterpreted as a rank-five proof or an exact-frequency selector.

The connection ledger is

```text
source:
  THM-2274's adaptive rank-two scalar packet, THM-1166's ordinary
  overlap, and THM-2080's odd-guard overlap;

target:
  a third relation on the same quotient-faithful scalar section;

map:
  minimal sparse support union -> pair-projection graph -> Hunter/factor
  safe floor -> positive Fourier comparison;

preserved:
  all nine scalar labels, guard parity, interval types, rational rank,
  pair-support branch, and fixed-section lift;

destroyed:
  the chosen starting basis, exact sparse coefficients, Fourier frequency,
  root sheet, owner, and post-expiration ancestry;

hostile controls:
  equal/opposite pair speeds, a support-three line containing the guard,
  a guard support-two line, disjoint support triples, pi-bracket
  directions, the N=51/52 boundary, and saturation in (43);

needed sidecar:
  a bounded anchored rank-three minor together with a legal route from
  its residue address to one exact energetic integer frequency.        (48)
```

The two exact companions use explicit `require` checks in normal and
optimized Python. The primary script checks every sparse-case fraction,
enumerates all admissible non-Haar graphs and Hunter trees, verifies the
disjoint-triple table, certifies all three adjacent Jackson boundaries,
and reconstructs the coefficients by convolution. The independent
referee rebuilds the same claims using a separate spanning-tree
enumerator and direct convolution. Both modes reproduce their stored
transcripts byte for byte. LRC(14) remains open. QED.
