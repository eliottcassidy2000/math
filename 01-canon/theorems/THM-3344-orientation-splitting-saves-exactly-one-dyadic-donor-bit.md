---
id: THM-3344
title: "Orientation splitting saves exactly one dyadic donor bit"
status: >
  PROVED + VERIFIED-EXACT.  There is one deterministic exactly fair
  unknown-bias coin extractor with T(1)=2, T(n)=n+1 at every nonpower of two,
  and T(n)=2n-1 at every power of two n>=2.  On each dyadic annulus it labels
  all nonboundary critical values at their pointwise floor, splits the two
  top orientations, and repairs the lower donor one bit before the annulus
  endpoint.  Among shell-composition-exact rules that stop every non-donor
  critical value at its pointwise floor, the donor bound 2n-1 is optimal:
  its required enumerator can have (1+x)-adic valuation at most one.
source: codex-kps-2026-08-12-orientation-split-donor
depends_on:
  - THM-3343-shifted-donor-rotation-bisects-exactly-the-dyadic-annuli
related:
  - THM-2160-dyadic-checksum-extracts-a-fair-bit-under-the-critical-run-deadline
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors
script: 04-computation/amm12592_orientation_split_saves_exactly_one_donor_bit_thm3344.py
output: 05-knowledge/results/amm12592_orientation_split_saves_exactly_one_donor_bit_thm3344.out
script_sha256: 14b4e433652f8fb1bdb17697a1cc41129f85e4266e8cf36630742c47234231e4
output_sha256: 0dd7112179e299ea2c059d89361a6a1dede4c8472e713c15597cdf6933af28fd
hash_basis: working-tree bytes (LF)
---

# THM-3344 -- one orientation split saves exactly one donor bit

Fix a power of two `d>=2` and put

```text
M=2d,                H=M-1=2d-1.                                  (1)
```

We refine THM-3343 on the annulus `S_(d,M)` of words whose first `d` bits
are equal and whose first `M` bits are nonconstant.

## 1. The odd-prefix construction

First work with nonconstant prefixes of length `H` and critical value at
least `d`.  At every critical value above the donor prescribe

```text
n-d odd:       heads,
n-d even:      tails,                       d<n<H,                  (2)
```

on both initial-bit orientations.  For prefix Hamming weight `w`, let

```text
A_w = number of critical-d prefixes,
E_w = number of prescribed heads,
O_w = number of prescribed tails,
D_w = E_w-O_w,
B_w = A_w+E_w+O_w.                                               (3)
```

The rotation proof of THM-3343 applies at length `H` without change:
rotation by one bijects the negative branches with the positive branches
ending in their initial bit, and each unmatched positive branch injects into
the donor by rotation through `n-d`.  Hence

```text
0<=D_w<=A_w.                                                       (4)
```

Every interior prefix layer has odd size.  Indeed its enumerator is

```text
(1+x)^(H-d)(1+x^d)-1-x^H,
```

and, modulo two,

```text
(1+x)^(d-1)(1+x^d)=(1+x)^(2d-1)=(1+x)^H.                         (5)
```

Since the binary expansion of `H=2d-1` consists entirely of ones, Lucas's
theorem makes every `binom(H,w)` odd.  Thus

```text
A_w-D_w=B_w-2E_w                                                   (6)
```

is a positive odd integer.

Set `tau_w=(-1)^w` and choose exactly

```text
q_w=(A_w-D_w+tau_w)/2                                             (7)
```

critical-`d` prefixes of weight `w` as heads, lexicographically.  Equation
(4) and the positive oddness in (6) give

```text
0<=q_w<=A_w.                                                       (8)
```

The signed heads-minus-tails count of the completed prefix layer is

```text
2q_w-A_w+D_w=tau_w.                                               (9)
```

## 2. The final orientation split

Every nonconstant `H`-prefix is now terminal; its two length-`M`
continuations receive the same label.  The only remaining annulus words are

```text
0^H 1,                1^H 0.                                     (10)
```

Label the first heads and the second tails.  The signed prefix enumerator is

```text
P(x)=sum_(w=1)^(H-1)(-1)^w x^w=(x^H-x)/(1+x).                    (11)
```

Refining a terminal `H`-prefix by its last bit multiplies its enumerator by
`1+x`.  The two words in (10) contribute `x-x^H`, so

```text
(1+x)P(x)+x-x^H=0.                                                (12)
```

Thus every length-`M` Hamming layer in the annulus has the same number of
heads and tails.

## 3. The uniform profile

Apply the construction independently to every consecutive dyadic annulus
`S_(2^r,2^(r+1))`; handle `d=1` by the ordinary `01/10` split.  The annuli
partition all nonconstant streams, and composition bisection proves exact
fairness for every `0<p<1`.  The stopping profile is

```text
T(1)=2,
T(n)=2n-1,          n>=2 a power of two,
T(n)=n+1,           otherwise.                                    (13)
```

For `n=d`, the donor verdict uses the first `H=2d-1` bits.  If
`d<n<H`, (2) decides at the first disagreement, time `n+1`.  At `n=H`,
the first disagreement occurs at `M=n+1` and (10) decides it.

This strictly Pareto-improves the cyclic-checksum profile of THM-2225 while
retaining its worst-case envelope `max(2,2n-1)`: all non-dyadic critical
values now meet their individual floor.

## 4. Why a second donor bit cannot be saved

The following sharpness statement is restricted but architecture-complete.
Consider any composition-exact coloring of `S_(d,M)` in which every
non-donor critical value `n>d` stops at its pointwise floor `n+1`.  Its label
may depend on the initial orientation.  Write these signs as
`s_(n,0),s_(n,1) in {+-1}` and define

```text
C_(n,0)(x)=x(1+x)^(M-n-1),
C_(n,1)(x)=x^n(1+x)^(M-n-1),
A(x)=C_(d,0)(x)+C_(d,1)(x),
S(x)=sum_(n>d,o) s_(n,o) C_(n,o)(x).                              (14)
```

Layer balance forces the donor-head enumerator

```text
R(x)=(A(x)-S(x))/2.                                                (15)
```

If all donor words stopped by `M-2`, refining their terminal prefixes to
length `M` would give each contribution at least two free bits.  Therefore

```text
(1+x)^2 divides R(x).                                              (16)
```

We show this is impossible for every orientation-sign word.  Both `A(-1)`
and `A'(-1)` vanish.  At `x=-1`, only the row `n=M-1` survives in `S`.
Thus `R(-1)=0` first forces

```text
s_(M-1,1)=-s_(M-1,0)=:-a.                                        (17)
```

The derivative of that top row at `-1` is `-(M-2)a`.  Rows at most `M-3`
have a square factor `1+x`.  The row `M-2`, when present, can alter the
derivative only by one of `-2,0,2`.  For `d>=4`, `M-2>=6`, so cancellation
is impossible.  For `d=2`, the `M-2` row is the donor rather than a fixed
row and the top derivative has magnitude two.  Hence

```text
R(-1)=0  implies  R'(-1)!=0,                                     (18)
```

so `v_(1+x)(R)<=1`.  This contradicts (16).  Some donor branch must reach
time `M-1`, and (13) is sharp in the stated shell-local floor-interior class.

## 5. Exact audit

```bash
python 04-computation/amm12592_orientation_split_saves_exactly_one_donor_bit_thm3344.py
python -O 04-computation/amm12592_orientation_split_saves_exactly_one_donor_bit_thm3344.py
```

Both modes pass.  The companion checks (4)--(9) through `d=256`, constructs
the literal extractor on every word through `M=16`, checks (12) symbolically
through `d=256`, and exhausts all 2,956 feasible orientation-sign
configurations for `d=2,4,8`, always finding maximal `(1+x)`-adic valuation
one.  The derivative proof is uniform; the census is hostile corroboration.
QED.
