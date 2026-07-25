---
id: THM-2160
title: "Dyadic shell extraction and the optimal half-tail critical-run rule"
status: >
  PROVED. Let independent flips have probabilities p and 1-p, with 0<p<1,
  and let n be the length of the first constant run. If h is the largest
  power of two not exceeding n, exact fixed-composition bisection gives a
  deterministic fair bit for every p. A cyclic checksum proves the requested
  max(2,2n-1) deadline. A sharper owner-stratified rule has T(1)=2, T(2)=3,
  and T(n)<=2n-2 for n>=3; when n=h>=2 it stops at 3h/2, which is optimal
  among branchwise composition-exact h-shell rules. Such shell bisection is
  possible exactly when h is a power of two, and a global shell coloring can
  ignore at most one fixed tail coordinate.
source: codex-2026-07-24-biased-coin-dyadic-shell
depends_on: []
related:
  - THM-2162
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
script: 04-computation/biased_coin_dyadic_half_tail_referee_codex_20260724.py
output: 05-knowledge/results/biased_coin_dyadic_half_tail_referee_codex_20260724.out
---

# THM-2160 -- dyadic critical-run fair extraction

Let

```text
X_1,X_2,... in {0,1}
```

be independent, with `P(X_i=0)=p` and `P(X_i=1)=1-p`, where `0<p<1`.
Outside the two constant rays, define the critical value

```text
n=min{k>=1:X_(k+1)!=X_1}.
```

Put

```text
h=2^floor(log_2 n),                 h<=n<=2h-1.       (1)
```

We give one deterministic rule which proves both requested bounds.

## 1. The dyadic-shell rule

For `h>=2`, write the first `2h` bits as

```text
b^h y_1...y_h.
```

The first half is constant because `h<=n`, and the second half is not `b^h`
because `n+1<=2h`. Define

```text
S(y)=sum_(i=1)^h i y_i mod h,                         (2)
```

using the representative in `{0,...,h-1}`, and output

```text
heads  iff  0<=S(y)<h/2,
tails  iff  h/2<=S(y)<h.                              (3)
```

For `h=1`, output heads on `10` and tails on `01`.

The rule is causal. Wait until the first change reveals `n` and hence `h`;
then read any still-unseen bits required by (2).

## 2. Every composition class is bisected

Fix `h=2^r>=2`. For a tail of Hamming weight `j`, where `1<=j<h`, write

```text
j=2^a u,                         u odd,
delta=h/2^(a+1).                                     (4)
```

Regard the tail positions as the cyclic group `Z/hZ`, with position `h`
represented by zero. Rotation of a weight-`j` tail by `delta` preserves its
weight and changes its checksum by

```text
j delta=(u h)/2=h/2 mod h.                            (5)
```

It therefore bijects the tails in the lower checksum half with those in the
upper half. Equivalently, rotate by `+delta` on the lower half and `-delta`
on the upper half to obtain a fixed-point-free, output-reversing involution.

Now fix the total number `k` of ones in the whole length-`2h` shell prefix.
If `1<=k<h`, the first half is `0^h` and the tail has weight `k`; if
`h<k<=2h-1`, the first half is `1^h` and the tail has weight `k-h`.
These classes are bisected by (4)--(5).

At the middle weight `k=h`, the only shell words not already covered are

```text
0^h 1^h,                       1^h 0^h.               (6)
```

Their checksums are respectively

```text
h(h+1)/2=h/2 mod h,                 0 mod h,           (7)
```

so (3) assigns opposite outputs. Thus every fixed-composition class in this
dyadic shell contains equally many heads and tails prefixes.

All prefixes in a class with `k` ones have the same probability

```text
p^(2h-k)(1-p)^k.                                      (8)
```

Hence heads and tails have equal probability on each shell. The shells,
indexed by the largest power of two below `n`, partition every nonconstant
sequence. The two constant rays have probabilities

```text
lim_(m->infinity) p^m=lim_(m->infinity)(1-p)^m=0.
```

Summing the shell equalities proves

```text
P(heads)=P(tails)=1/2                                  (9)
```

for every unknown `p in (0,1)`.

## 3. The two stopping bounds

Reading through flip `2h` always suffices, and

```text
2h<=2n.                                                (10)
```

This proves the first requested bound.

For the sharper bound, the coefficient of `y_h` in (2) is

```text
h=0 mod h,
```

so the output ignores flip `2h`. That flip is needed only when it is itself
the first change. A valid stopping time is therefore

```text
T(1)=2,
T(n)=max(n+1,2h-1)                    for n>=2.        (11)
```

For `n>=2`,

```text
n+1<=2n-1,                    2h-1<=2n-1,             (12)
```

which gives

```text
T(n)<=max(2,2n-1).                                    (13)
```

For the sample prefix `00001`, one has `n=h=4`. The decision uses only

```text
S=1+2X_6+3X_7 mod 4,
```

so at most two additional flips are required.

## 4. Hostile controls and sharpness

The smallest nontrivial shell is

```text
0001 -> heads,     0010 -> tails,     0011 -> tails,
1100 -> heads,     1101 -> heads,     1110 -> tails. (14)
```

In particular, the two critical-value-two cylinders receive fixed opposite
answers at time three; their conditional imbalance is compensated by the
critical-value-three leaves in the same composition classes. Fairness is not
being asserted conditional on the exact critical value.

The envelope cannot be uniformly improved to `max(2,2n-2)`. Such a rule would
have to stop after two flips on critical values one and two. Causality would
therefore make all four length-two cylinders terminal. The probability of
either output would be a subset sum of

```text
p^2,                    p(1-p), p(1-p),              (1-p)^2,
```

and no such subset sum is identically `1/2`: evaluation at `p=0` already
gives constant term zero or one. Thus the single envelope
`max(2,2n-1)` cannot be shifted down uniformly by one. The obstruction is
only `n=2`; the next construction improves every `n>=3`.

## 5. A faster half-tail construction

The checksum is especially simple, but it is not deadline-optimal on a
dyadic shell. Assume `h>=2`, write

```text
h=2t,
b=X_1,
z_i=X_(h+i) xor b.                                    (15)
```

The shell condition is `z!=0`. We first define the heads set for `b=0`.
On the stratum `z_1=1`, output heads exactly when

```text
z_2+...+z_t=0 mod 2.                                  (16)
```

Thus every sequence with exact critical value `n=h` is decided after the
first `t=h/2` relative tail bits.

It remains to complete (16) on `z_1=0` so that every relative Hamming layer
is bisected. Let `A_j` be the number of heads of relative weight `j` already
prescribed by (16). Its generating polynomial is

```text
A(u)
 =u ((1+u)^(t-1)+(1-u)^(t-1))/2 (1+u)^t
 ={u(1+u)^(h-1)+u(1+u)(1-u^2)^(t-1)}/2.              (17)
```

Write

```text
e_j=[u^j]u(1+u)(1-u^2)^(t-1).                         (18)
```

Among the still-unlabelled words with `z_1=0` and weight `j`, choose exactly

```text
R_j=(binom(h-1,j)-e_j)/2                              (19)
```

as heads, for each `1<=j<h`; a lexicographic choice makes the rule fully
deterministic. The endpoint `z=1^h` can initially be labelled arbitrarily.

This completion is legal. Explicitly,

```text
e_(2a+1)=e_(2a+2)=(-1)^a binom(t-1,a).                (20)
```

Since `h-1` and `t-1` have all binary digits equal to one, Lucas' theorem
makes both binomial coefficients in (19)--(20) odd, so `R_j` is integral.
Moreover

```text
|e_j|<=binom(h-1,j).                                  (21)
```

For example, inject an `a`-subset of `t-1` positions into a `j`-subset of
`2t-1` positions by adjoining `a+1` fixed positions when `j=2a+1`, or
`a+2` fixed positions when `j=2a+2`; the latter case has `a<=t-2` because
`j<h`. Hence `0<=R_j<=binom(h-1,j)`.

Pascal's identity and (17)--(19) now give

```text
A_j+R_j
 ={binom(h-1,j-1)+e_j+binom(h-1,j)-e_j}/2
 =binom(h,j)/2.                                       (22)
```

Thus every nonendpoint relative-weight layer is bisected. For `b=1`, use
the complementary label on the same relative word `z`. Its layers are also
bisected, and the two middle-composition words, both represented by
`z=1^h`, receive opposite labels. The fixed-composition proof of Section 2
therefore proves fairness for every `p`.

If `n=h`, rule (16) stops at

```text
h+t=3h/2.                                             (23)
```

If `h<n<=2h-1`, then `z_1=0`; reading through flip `2h` determines the
completion (19). Consequently

```text
T(1)=2,          T(2)=3,
T(n)<=2n-2                         for every n>=3.     (24)
```

Indeed, (23) is at most `2h-2` for `h>=4`, while `n>h` gives
`2h<=2n-2`. For the sample `00001`, `h=4,t=2`; only flip six is additionally
needed, and the lower-branch convention outputs heads exactly when `X_6=0`.

## 6. Why dyadic shells, and what is optimal

### 6.1 Dyadicity is necessary and sufficient

Consider the same length-`2h` shell for an arbitrary positive integer `h`
and demand equal heads/tails counts in every total-composition class. At
total weights `j=1,...,h-1`, only the prefix `0^h` occurs and the class has
size `binom(h,j)`. Therefore all internal binomial coefficients must be even.
Over `F_2`,

```text
(1+u)^h=1+u^h
```

holds exactly when `h` is a power of two: if the binary expansion of `h`
has two or more nonzero digits, Lucas' theorem supplies an internal `j` with
`binom(h,j)` odd. Conversely, Sections 1--2 construct the bisection when
`h` is a power of two. Hence

```text
composition-exact h-shell extraction exists
iff h is a power of two.                              (25)
```

Equivalently, cyclic rotation partitions every nonconstant binary necklace
of dyadic length into even orbits; an alternating coloring of each orbit,
plus the middle pair, gives a noncanonical sufficiency proof.

### 6.2 The half-tail deadline is shell-optimal

Fix either initial-bit branch and suppose every word with exact critical
value `n=h` is decided from only the first `d` relative tail bits. In the
relative-weight-`h-1` layer, the `h` words have one zero and must split
`h/2`--`h/2`. The exceptional word with `z_1=0` contributes either zero or
one head, so the `h-1` exact-`n=h` words must contribute `h/2` or `h/2-1`.

Among those words, zeros in positions `2,...,d` yield `d-1` distinguishable
prefixes, while all `h-d` zeros in later positions share one prefix.
Therefore their possible head counts lie in

```text
[0,d-1] union [h-d,h-1].                              (26)
```

Neither required count belongs to (26) when `d<h/2`. Thus

```text
d>=h/2,                                               (27)
```

and (16) attains equality. The time `3h/2` in (23) is optimal within this
branchwise composition-exact shell model; (27) is not asserted as a lower
bound for every conceivable cross-shell stopping rule.

### 6.3 One globally ignored coordinate is maximal

Fix one branch and encode a composition-bisecting coloring by
`f:{0,1}^h->{-1,+1}`, extending it to the omitted zero word if necessary.
If `f` is independent of a fixed set of `r` tail coordinates, its signed
layer enumerator factors as

```text
F(u)=sum_z f(z)u^|z|=(1+u)^r G(u).                    (28)
```

Composition balance kills every coefficient of degree `1,...,h-1`, so

```text
F(u)=epsilon_0+epsilon_h u^h.                         (29)
```

For even `h`, divisibility by `1+u` forces
`epsilon_h=-epsilon_0`. The root `u=-1` of
`epsilon_0(1-u^h)` is simple, proving

```text
r<=1.                                                 (30)
```

The checksum's ignored final coordinate is therefore maximal among global
fixed-coordinate quotients. The faster rule evades (30) by stratified
obliviousness: it ignores half the tail on `z_1=1` but may inspect the whole
tail on `z_1=0`, where (19) pays the exact layer defect.

## 7. Transfer boundary

The mechanism is a composition-preserving cyclic symmetry. It works because
Bernoulli cylinder mass depends only on Hamming weight. A target distribution
which is not exchangeable inside composition classes cannot inherit the
involution without an additional weight or orbit-incidence sidecar. This
distinction is essential in LRC incidence fibers, GMC root packets, and
valuation-sensitive polynomial descent.

There is nevertheless one exact algebraic bridge to THM-2201. For
`h=2^r`,

```text
F_2[C_h]
 =F_2[g]/(g^h-1)
 =F_2[epsilon]/(epsilon^h),          epsilon=g-1.     (31)
```

Thus the same Pascal/Hasse change of basis which triangularizes the
thirteen-sheet LRC fibre also triangularizes a dyadic shell. Dyadicity in
(25) is exactly the characteristic-two identity
`g^h-1=(g-1)^h`; the binomial coefficients off the two endpoints vanish.
The two problems then ask different questions of the local algebra:

```text
biased coin:
  find one balanced composition coloring, so a carefully chosen quotient
  may ignore one fixed coordinate;

LRC root fibre:
  reconstruct every labelled owner incidence under deletion/insertion, so
  the full Hasse jet is required and the socle epsilon^12 records the
  occupancy-thirteen state lost by augmentation.                    (32)
```

This explains both the usefulness and the limit of the puzzle analogy. The
critical run suggests looking for a first invisible filtration coordinate;
it does not license replacing THM-2198's guard, deep-owner bits, or winding
word by an unweighted fair coloring.

QED.
