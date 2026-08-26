---
id: THM-4219
title: "No-sink endpoint-energy floor and near-ordinal sharpness"
status: >
  PROVED / FINITE-EXACT / OPEN FIREWALL. Proves the all-order no-sink
  floor Delta_V>=n(n-1)H^2 and the weaker equality classification
  Delta_V=2nH^2 iff C3; gives an exact one-defect-permutation reformulation;
  proves two all-order near-ordinal strong families satisfy the open
  five-copy inequality and make its coefficient sharp; and records a
  FINITE-EXACT strong-class census through order nine. The five-copy
  inequality for every strong tournament and its conditional three-eighths
  consequence remain OPEN in HYP-9081.
source: codex-endpoint-gap-session-20260826
depends_on:
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
related:
  - HYP-9081-strong-tournament-five-copy-endpoint-energy-inequality
script: 04-computation/tournament_no_sink_endpoint_energy_floor_thm4219.py
output: 05-knowledge/results/tournament_no_sink_endpoint_energy_floor_thm4219.out
finite_census_script: 04-computation/tournament_strong_endpoint_five_copy_census_thm4219.cpp
finite_census_output: 05-knowledge/results/tournament_strong_endpoint_five_copy_census_thm4219.out
script_sha256: 0adbbc261ecae9c3401bb8c773bcf2a75a1004a74365cb6da6e8f40e8f21219c
output_sha256: 6ddebc148eee898f8bac0b1577c6d247293612707c13517e3d154bdc2d78d2c6
finite_census_script_sha256: cbaebfdd5b8327e87e3595744adc957560c8187daf4deede02b88b6a146c662b
finite_census_output_sha256: d90fcf6b4f03940451c606c238b4c757d2ebcc000c9f9bfb2305e8a6b90407f2
hash_basis: raw LF bytes
primary_audit: >
  PASS. Literal subset Hamilton-path DP checks the split injection,
  one-defect identity, canonical slacks, exact defect identity, and no-sink
  floor on all 1,099 labelled tournaments through order five, including all
  738 no-sink rows. It independently checks the two tower profiles against
  literal DP and the exact partition sums through orders twelve. Normal and
  python -O streams byte-match.
finite_census_audit: >
  PASS. Homebrew nauty 2.9.3 gentourng -c supplies one representative of
  every strong tournament class through order nine: 184,537 classes total.
  Exact subset DP and checked exact cross multiplication find no five-copy
  failure and equality only at C3. Clang O0/O3 streams byte-match through
  order nine; ASan/UBSan passes order seven.
---

# THM-4219 -- no-sink endpoint-energy floor and near-ordinal sharpness

**PROVED / FINITE-EXACT / OPEN FIREWALL.**

This packet separates its proved candidate content from the open statement
in HYP-9081. In particular, no finite census or near-ordinal calculation
below is promoted into an all-order proof of the five-copy inequality.

## 1. Conventions and statements

All tournaments are finite and nonempty. Retain THM-4208's quantities

```text
H=Hamilton-path count,                  W=capacity mass,
v_i=V_i^0+V_i^1,                       a=sum_i v_i=W+H,
b=W+2H=a+H,                            m=sum_i v_i^2,
Delta_V=a(a+2H)-3m=b^2-H^2-3m.                       (1)
```

For a tournament of order `n`, write `d_i^+` for the outdegree of `i`.
THM-4208 gives the canonical nonnegative endpoint slack

```text
sigma_i=H+sum_(j->i)v_j-v_i>=0,
Delta_V/2=sum_i v_i sigma_i.                          (2)
```

> **Theorem 1 (no-sink endpoint-energy floor).** Every no-sink tournament
> `X` of order `n` satisfies
>
> ```text
> Delta_V(X)>=n(n-1)H(X)^2.                            (3)
> ```
>
> Moreover,
>
> ```text
> Delta_V(X)>=2nH(X)^2,
> Delta_V(X)=2nH(X)^2  iff  X=C3.                      (4)
> ```
>
> No equality classification for the stronger floor `(3)` is claimed.

Call an adjacency in a permutation `(x_1,...,x_n)` **bad** when
`x_(k+1)->x_k`. Let `D_i` count permutations having exactly one bad
adjacency and having it immediately after vertex `i`; put `D=sum_i D_i`.

> **Theorem 2 (exact one-defect model).** For every tournament,
>
> ```text
> v_i=H+D_i,                 W=(n-1)H+D.               (5)
> ```
>
> Consequently the open five-copy inequality
>
> ```text
> 5m<=(W+H)(W+3H)                                      (6)
> ```
>
> is exactly the one-defect dispersion inequality
>
> ```text
> 5sum_i D_i^2
>  <=D^2+2(n-4)HD+n(n-3)H^2.                          (7)
> ```

Theorem 3 below proves `(6)` on two all-order families and shows that the
coefficient `5` cannot be increased in an all-strong version of `(6)`.

Let `T(n,r)` have transitive base `B={0,...,n-2}`, final vertex `z=n-1`,
and

```text
z->i  iff  i in S_r={0,...,r-1,n-3}.                  (8)
```

Every other base-to-`z` edge points toward `z`.
For `r=1,2` in the ranges below these tournaments are strong: `z` reaches
the whole base through `0`, while every base vertex either points directly
to `z` or reaches a base vertex that does.

> **Theorem 3 (two exact near-ordinal towers).** Put
> `G=(W+H)(W+3H)-5m`.
>
> For `T(n,1)`, `n>=5`, set `x=2^(n-4), y=3^(n-4)`. In natural vertex
> order,
>
> ```text
> H=3x+3,
> v_0=3x+3,
> v_i=3^i 2^(n-i-3)+3*2^i,                 1<=i<=n-5,
> v_(n-4)=2y+3x-1,
> v_(n-3)=4y+6x-1,
> v_(n-2)=4y+6x-3,
> v_z=6x-1,                                           (9)
>
> W=14y+18x-12,              b=14y+24x-6,
> m=(609x^2+570xy-330x+196y^2-180y+45)/5,
> Delta_V=2(504x^2+825xy-270x+196y^2-150y)/5,
> G=6(17xy+2y-7x^2+4x-3)>0.                           (10)
> ```
>
> For `T(n,2)`, `n>=6`, set `x=2^(n-5), y=3^(n-5)`. Then
>
> ```text
> H=9x+1,
> v_0=9x+2,                  v_1=12x+5,
> v_i=15*3^(i-2)2^(n-i-3)+2^(i+1),        2<=i<=n-5,
> v_(n-4)=10y+4x-3,
> v_(n-3)=20y+8x-3,
> v_(n-2)=20y+4x-9,
> v_z=12x-3,                                          (11)
>
> W=70y+14x-20,              b=70y+32x-18,
> m=(871x^2+1800xy-540x+2940y^2-1620y+347)/3,
> Delta_V=2(36x^2+1340xy-315x+980y^2-450y-12),
> G=2(2220xy+270y-763x^2-405x-383)/3>0.               (12)
> ```
>
> Finally,
>
> ```text
> (W+H)(W+3H)/m
>   =5+(255/98)(2/3)^(n-4)+o((2/3)^n)   for T(n,1),
>   =5+(74/49)(2/3)^(n-5)+o((2/3)^n)    for T(n,2).    (13)
> ```
>
> Thus `5` is the sharp possible coefficient in `(6)`.

> **Finite-exact census 4.** Among one representative of every strong
> tournament isomorphism class through order nine, `(6)` has no failure and
> equality occurs only at `C3`. The class counts and minimum gaps by order
> `3,...,9` are
>
> ```text
> classes:  1, 1, 6, 35, 353, 6008, 178133,
> min G:    0, 44, 510, 3186, 19842, 122778, 753810.   (14)
> ```
>
> The minimum exact ratios `(W+H)(W+3H)/m` are approximately
>
> ```text
> 5,
> 5.2303664921,
> 5.3944315545,
> 5.3683662851,
> 5.3334509705,
> 5.2841077119,
> 5.2273897553.                                       (15)
> ```

The finite statement is `FINITE-EXACT`, not an all-order theorem.

## 2. Proof of the no-sink floor

Fix `i` and a Hamilton path `R`. Split `R` immediately after its occurrence
of `i`; the prefix is a nonempty directed path ending at `i`, and the suffix
is a Hamilton path of the complement. Concatenation recovers `R`, so this is
an injection into the objects counted by `v_i`. Hence

```text
v_i>=H.                                                (16)
```

Sum the canonical slacks `(2)`. Every `v_j` occurs once for each outneighbor
of `j`, so

```text
sum_i sigma_i
 =nH+sum_i(d_i^+-1)v_i.                                (17)
```

If `X` has no sink, every `d_i^+>=1`; using `(16)` and
`sum_i d_i^+=n(n-1)/2` gives

```text
sum_i sigma_i
 >=nH+H sum_i(d_i^+-1)
 =n(n-1)H/2.                                           (18)
```

Equations `(2)`, `(16)`, and `(18)` yield `(3)`. More precisely, the proof
has the exact nonnegative certificate

```text
Delta_V-n(n-1)H^2
 =2sum_i(v_i-H)sigma_i
  +2H sum_i(d_i^+-1)(v_i-H).                           (19)
```

Thus equality in `(3)` holds exactly when every vertex with `v_i>H` has
both `sigma_i=0` and `d_i^+=1`; this packet does not classify when that can
happen.

A no-sink tournament has `n>=3`. If `n>3`, then `(3)` is strictly stronger
than the lower bound in `(4)`. If `n=3`, no-sink forces every outdegree to be
one, hence `X=C3`. Directly, `H(C3)=3`, `v=(3,3,3)`, and
`Delta_V(C3)=54=2*3*3^2`. This proves `(4)`.

## 3. Proof of the one-defect model

An object counted by `v_i` is a pair `(P,Q)` in which `P` is nonempty and
ends at `i`, while `Q` is a Hamilton path of the complement and may be empty.
Concatenate `P,Q`. If `Q` is empty or `i->first(Q)`, the result is a Hamilton
path. For fixed `i`, splitting that Hamilton path immediately after `i`
recovers the unique preimage. Otherwise `first(Q)->i`, and the concatenation
has exactly one bad adjacency, immediately after `i`. Cutting there is again
the unique inverse. This proves `v_i=H+D_i`; summing and using
`sum_i v_i=W+H` proves `(5)`.

Substitute `v_i=H+D_i` and `a=nH+D` into `(6)`. Its right side minus its left
side is exactly

```text
D^2+2(n-4)HD+n(n-3)H^2-5sum_iD_i^2,                  (20)
```

which proves the equivalence with `(7)`. The same substitution into the
unconditional defect gives the companion identity

```text
Delta_V
 =n(n-1)H^2+2(n-2)HD+D^2-3sum_iD_i^2.                (20a)
```

Thus `(19)` proves that the final three terms in `(20a)` are nonnegative in
every no-sink tournament, while `(20)` isolates the strictly stronger open
five-copy dispersion demand.

## 4. Direct counting of the towers

The tower formulae are geometric sums, not interpolation. Write
`B={0,...,m-1}`, `z=m`, and `S=S_r`. For `U subset B`, every Hamilton path of
the induced tournament on `U union {z}` has a unique form

```text
L^up, z, R^up,               U=L disjoint_union R,     (21)
```

where `up` means increasing base order. It is directed exactly when

```text
L is empty or max(L) notin S,
R is empty or min(R) in S.                             (22)
```

Let `h_S(U)` count the partitions in `(21)--(22)`, and define

```text
F(C)=1+sum_(q in C\S) 2^|C intersect [0,q)|.           (23)
```

The term `1` chooses an empty left part; if its maximum is `q notin S`, every
smaller element is chosen freely and every larger element is forced into the
other part. Hence `(23)` counts precisely the splits `C=L disjoint_union Q`
with `L` empty or `max(L) notin S`.

The endpoint masses now have the exact finite-sum description

```text
H=h_S(B),                       v_z=F(B),

v_i=sum_( {i+1,...,m-1} subset U subset B\{i}) h_S(U)
    +sum_(q in S, q<=i)
       sum_(R subset [q,i], min R=q, max R=i) F(B\R).  (24)
```

In the first sum `z` lies in the complementary path; in the second it lies
in the rooted path ending at `i`. These two cases are disjoint and exhaustive.

Evaluating `(23)--(24)` at `r=1` gives `(9)`. On the ordinary spine the same
geometric summation can be recorded as

```text
2v_(i+1)=3v_i+3*2^i,                 1<=i<=n-6,        (25)
```

with the four boundary values evaluated directly from `(24)`. At `r=2`, the
ordinary-spine recurrence is

```text
2v_(i+1)=3v_i+2^(i+1),               2<=i<=n-6,        (26)
```

and the boundary evaluation gives `(11)`. Summing the displayed endpoint
masses gives `W`; summing their squares gives `m`; substituting in `(1)` and
the definition of `G` gives `(10)` and `(12)`.

For `T(n,1)`, `y/x=(3/2)^(n-4)>=3/2` and `x>=2`; hence
`17xy-7x^2>0` and the remaining `2y+4x-3` is positive. For `T(n,2)`,

```text
2220xy-763x^2>=2567x^2,
270y-405x>=0,                                           (27)
```

and `x>=2`, proving both strict gap claims. Keeping the leading `y^2` and
`xy` terms gives `(13)`.

## 5. Finite-exact universe and controls

The census universe in `(14)` is generated by

```bash
clang++ -std=c++17 -O3 -Wall -Wextra -pedantic \
  04-computation/tournament_strong_endpoint_five_copy_census_thm4219.cpp \
  -o /tmp/census4219
for n in $(seq 3 9); do
  /opt/homebrew/bin/gentourng -q -c "$n" | /tmp/census4219
done
```

Here nauty `gentourng -c` supplies the strong-only universe; the evaluator
does not independently recheck strongness. It uses the literal recurrence

```text
finish[M,j]=sum_(i in M\{j}, i->j) finish[M\{j},i],
h[M]=sum_j finish[M,j],
v_j=sum_(M contains j) finish[M,j]h[V\M].              (28)
```

All inequality arithmetic is exact. The minimum-ratio representative is
selected by checked unsigned-128 cross multiplication, not floating-point
comparison. Positive controls are `C3` word `101`, giving
`(H,W,v,m,G)=(3,6,(3,3,3),27,0)`, and the unique strong order-four word
`100111`, giving `(H,W,m,G)=(5,22,191,44)`. Filter-hostile controls are the
transitive words `111` and `111111`, whose gaps are `-42` and `-170`; thus
strongness in HYP-9081 is load-bearing.

## 6. Logical firewall

This theorem does **not** prove `(6)` for every strong tournament. It also
does not identify `(6)` with the desired three-eighths defect inequality.
For a strong tournament, `(16)` gives `b>=4H`, and if `G>=0` then

```text
Delta_V>=2(b^2-H^2)/5>=3b^2/8.                         (29)
```

But the converse implication has genuine slack. Exactly,

```text
Delta_V-3b^2/8=(b^2-16H^2+24G)/40,                    (30)
```

or equivalently the three-eighths inequality permits `m` to exceed the
five-copy upper bound by

```text
(b^2-16H^2)/120.                                       (31)
```

The distinction is already visible on the conjectural equality family:

```text
G(P_a triangleright C3)=-6(4^a-1)<0,       a>=1,       (32)
```

although `Delta_V=3b^2/8` there. Thus `(6)` is a strong-terminal sufficient
route only. HYP-9081 states the open quantifiers and the conditional ordinal
reduction without blending them into this theorem packet.
