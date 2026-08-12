---
id: THM-3340
title: "A single delayed donor and cyclic rotation prove every pointwise AMM 12592 floor"
status: >
  PROVED + VERIFIED-EXACT. For every dyadic horizon M there is a deterministic
  exactly fair unknown-bias coin extractor with T_M(1)=M and T_M(n)=n+1 for
  every 2<=n<M, retaining THM-2225 for n>=M. The rule labels every even
  critical value n>=2 heads and every odd n>=3 tails at the first
  disagreement; one lexicographic subset of the delayed n=1 branch repairs
  each Hamming layer. Left cyclic rotation bijects odd branches with the even
  branches ending in their initial bit. The unmatched even branches inject
  into the n=1 layer, proving the exact repair capacity 0<=r_w<=C(M-2,w-1).
  A second polynomial factorization sharpens the unmatched defect to at most
  one donor orientation, hence C(M-2,w-1)/2<=r_w<=C(M-2,w-1).
  Choosing M>n proves T_opt(n)=n+1 for every n>=2; M=2 handles n=1. Thus the
  pointwise optimum is exactly n+1 for every positive integer n. This does
  not construct one extractor attaining all floors simultaneously.
source: codex-kps-2026-08-12-single-donor-rotation
depends_on:
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
related:
  - THM-2966-spine-normal-form-for-critical-run-fair-extractors
  - THM-3032-sharpened-half-tail-extractor-and-shell-four-pareto-frontier
  - THM-3337-cross-shell-compression-attains-the-T4-floor
  - THM-3338-horizon16-cross-shell-surgery-closes-the-first-fifteen-pointwise-AMM-values
script: 04-computation/amm12592_single_donor_all_pointwise_floors_thm3340.py
output: 05-knowledge/results/amm12592_single_donor_all_pointwise_floors_thm3340.out
script_sha256: 961b0cedfc2aa8f18247c8f58a5629882f0c27809cfec4fe40ad206287b4d1cc
output_sha256: 67381c3442c2bb8ca87e93a0e86f3793c5824855cbe396d3ae2060edd1822841
hash_basis: working-tree bytes (LF)
---

# THM-3340 -- one donor proves all pointwise floors

Let independent bits satisfy

```text
P(X_i=0)=p,       P(X_i=1)=q=1-p,       0<p<1,
```

and let

```text
n=min{k>=1:X_(k+1)!=X_1}
```

be the critical value of a nonconstant stream. For an extractor, `T(n)` is
the worst stopping time over streams of critical value `n`.

## 1. The finite dyadic rule

Fix `M=2^r`, `r>=1`. On a branch with `2<=n<M`, stop at flip `n+1` and use

```text
n even:     heads on both initial-bit orientations;
n odd:      tails on both initial-bit orientations.                 (1)
```

On `n=1`, read through flip `M`; its finite rule is specified in section 3.
If the first `M` bits are constant, so `n>=M`, continue the THM-2225
cyclic-checksum extractor without modification. These events are disjoint.

We prove that the `n=1` branch has enough capacity to make every nonconstant
length-`M` Hamming layer exactly half heads.

## 2. The rotation-capacity lemma

Fix a weight `1<=w<M`. In the length-`M`, weight-`w` layer write

```text
A_w = binom(M-2,w-1),
E_w = number of words with even critical value n>=2,
O_w = number of words with odd critical value n>=3.                 (2)
```

There are `2A_w` words of critical value one: `A_w` beginning `01` and
`A_w` beginning `10`. Thus

```text
binom(M,w)=2A_w+E_w+O_w.                                           (3)
```

**Rotation lemma.**

```text
0 <= E_w-O_w <= 2A_w.                                              (4)
```

*Proof.* Let `L` move the first bit of a word to the end. If a word has odd
critical value `n>=3`, then `L` has critical value `n-1`, which is even, and
its last bit equals its initial bit. Conversely, right rotation sends every
even-critical word whose last bit equals its initial bit back to a unique
odd-critical word. Therefore `L` is a weight-preserving bijection

```text
{odd n>=3}  <-->  {even n>=2 and X_M=X_1}.                         (5)
```

The defect `E_w-O_w` consequently counts the even-critical words with
`X_M!=X_1`, so it is nonnegative.

For such an unmatched word, let its even critical value be `n` and rotate
left by `n-1` places. The image begins with two different bits, hence has
critical value one, and rotation preserves weight. This map is injective:
the image ends in a run of exactly `n-1` copies of its initial bit, preceded
by the opposite bit `X_M`; that terminal run recovers `n`, and rotating back
recovers the word. The target contains only `2A_w` critical-value-one words,
proving the upper bound in (4). QED

This is the exact capacity sidecar. A first-moment comparison alone would
show neither injection nor the one-sided repair below.

## 3. Repairing the donor layer

Because `M` is a power of two, Lucas's theorem gives

```text
binom(M,w) even,       1<=w<M.                                    (6)
```

Equations (3) and (6) imply that

```text
D_w:=E_w-O_w
```

is even. Define

```text
r_w := A_w-D_w/2 = binom(M,w)/2-E_w.                              (7)
```

The rotation lemma gives

```text
0<=r_w<=A_w.                                                       (8)
```

There is also a sharper algebraic certificate. Put `m=M-2` and let

```text
F_M(x)=sum_(n=2)^(M-1) (-1)^n
       (1+x^(n-1))(1+x)^(M-n-1)=sum_(k=0)^m f_k x^k.              (8a)
```

Thus `f_(w-1)=E_w-O_w`. A geometric-series calculation factors

```text
F_M(x)=A_m(x)+x^m A_m(1/x),
A_m(x)=((1+x)^m-1)/(x+2)
      =x sum_(i=0)^(m/2-1) (1+x)^(2i).                           (8b)
```

If `a_k=[x^k]A_m(x)`, then

```text
a_k=sum_(i=0)^(m/2-1) binom(2i,k-1),
f_k=a_k+a_(m-k).                                                   (8c)
```

The hockey-stick identity gives

```text
a_k<=binom(m-1,k),       a_(m-k)<=binom(m-1,k-1),
```

so `0<=f_k<=binom(m,k)=A_(k+1)`. Moreover `f_k` is even: modulo
two, multiplying (8a) by `x` and using
`(1+x)^M=1+x^M` makes every coefficient vanish. Consequently (7) actually
obeys the stronger interval

```text
A_w/2<=r_w<=A_w.                                                   (8d)
```

The rotation proof supplies the conceptual injection; (8a)--(8d) supply an
independent closed coefficient certificate and show that at least half of
the chosen `01` donor orientation is always used.

There are exactly `A_w` words of weight `w` beginning `01`. Order their
last `M-2` bits lexicographically within their weight class and declare the
first `r_w` of them heads. Declare every `n=1` word beginning `10` tails.
This is a literal deterministic rule, since (8) places the requested count
inside the available class.

By (1), the total number of heads in layer `w` is

```text
E_w+r_w=binom(M,w)/2.                                              (9)
```

Thus every nonconstant Hamming layer is bisected.

For reference, the quantities in (2) can be calculated without enumerating
words. The number in layer `w` with critical value `n` is

```text
binom(M-n-1,w-1)+binom(M-n-1,w-n),                                (10)
```

with out-of-range binomial coefficients interpreted as zero. Summing (10)
over even and odd `n` gives (7).

## 4. Exact fairness and deadlines

Equation (9) makes the finite head probability

```text
sum_(w=1)^(M-1) [binom(M,w)/2] p^(M-w)q^w
  = (1-p^M-q^M)/2.                                                (11)
```

THM-2225 has the same aggregate head probability on the event that the first
`M` bits are nonconstant. Replacing that finite prefix and retaining its
continuation from `0^M` and `1^M` therefore preserves total head probability
`1/2` for every `p`. The resulting extractor is exactly fair, deterministic,
and causal. Its finite deadline profile is

```text
T_M(1)=M,
T_M(n)=n+1,             2<=n<M.                                  (12)
```

No verdict is made on a constant prefix; branches in (1) stop exactly when
their first disagreement appears, while the donor verdict uses only the
first `M` bits.

## 5. Every pointwise optimum

No exactly fair rule can stop on a constant prefix. If it did so at `0^k`
or `1^k`, one verdict would have probability at least `p^k` or `q^k`, which
exceeds `1/2` as the corresponding bias tends to one. Hence universally

```text
T(n)>=n+1.                                                         (13)
```

For any fixed `n>=2`, choose a power of two `M>n`. The extractor above has
`T_M(n)=n+1`, meeting (13). For `n=1`, take `M=2`; then the donor rule is the
usual `01/10` split and stops at flip two. Therefore

```text
                         T_opt(n)=n+1

for every positive integer n.                                    (14)
```

The quantifiers matter: (14) is pointwise. The horizon `M`, and hence the
extractor, may depend on the target `n`. This theorem does **not** assert one
extractor with `T(n)=n+1` for every `n`, nor does it settle the optimal
asymptotic envelope `C*` of HYP-9061.

## 6. Exact audit

```bash
python 04-computation/amm12592_single_donor_all_pointwise_floors_thm3340.py
python -O 04-computation/amm12592_single_donor_all_pointwise_floors_thm3340.py
```

Both runs use only standard-library integer arithmetic. They check (2)--(10)
for every dyadic horizon through `M=512`, and independently enumerate every
word for `M=2,4,8,16`, checking the two rotation maps, their inverses and
injectivity, the lexicographic donor rule, and exact layer bisection. The
proof itself is the rotation argument above; the computation is a hostile
referee, not a dependency. QED.
