---
id: THM-2160
title: "A dyadic checksum extracts a fair bit within the sharp critical-run envelope"
status: >
  PROVED. Let independent flips have probabilities p and 1-p, with 0<p<1,
  and let n be the length of the first constant run. If h is the largest
  power of two not exceeding n, a cyclic checksum on the second half of the
  first 2h flips splits every fixed-composition shell exactly in half. The
  resulting deterministic bit is fair for every p. It is decided after at
  most 2n flips, and in fact after at most max(2,2n-1). The stronger uniform
  envelope max(2,2n-2) is impossible.
source: codex-2026-07-24-biased-coin-dyadic-shell
depends_on: []
related: []
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
gives constant term zero or one. Thus (13) is globally one-step sharp.

## 5. Transfer boundary

The mechanism is a composition-preserving cyclic symmetry. It works because
Bernoulli cylinder mass depends only on Hamming weight. A target distribution
which is not exchangeable inside composition classes cannot inherit the
involution without an additional weight or orbit-incidence sidecar. This
distinction is essential in LRC incidence fibers, GMC root packets, and
valuation-sensitive polynomial descent.

QED.
