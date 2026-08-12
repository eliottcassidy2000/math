---
id: THM-3337
title: "Cross-shell compression attains the AMM 12592 floor T(4)=5"
status: >
  PROVED + VERIFIED-EXACT. There is one deterministic exactly fair extractor
  for an unknown-bias Bernoulli coin whose deadline vector on critical values
  1 through 7 is (2,3,5,5,7,7,8), and which agrees with the THM-2225 cyclic-
  checksum extractor for every critical value n>=8. Consequently the
  pointwise optimum at n=4 is exactly 5. The construction deliberately
  violates separate dyadic-shell balance: the doubled Bernstein deficits of
  shells {2,3} and {4,5,6,7} are nonzero exact opposites. Thus the first open
  value left by THM-3032 closes by cross-shell cancellation, not by a better
  coloring of the m=4 shell.
source: codex-kps-2026-08-12-anthropic-inertia-transfer
depends_on:
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-3032-sharpened-half-tail-extractor-and-shell-four-pareto-frontier
related:
  - THM-2966-spine-normal-form-for-critical-run-fair-extractors
  - "Claude, More than two thirds of the zeros of the Riemann zeta function lie on the critical line (2026), finite-compression inspiration only"
script: 04-computation/amm12592_cross_shell_T4_floor_thm3337.py
output: 05-knowledge/results/amm12592_cross_shell_T4_floor_thm3337.out
script_sha256: af051babc893801ffc133448fbb15ba758f0199a2f27d387a9f3dfedc27660fd
output_sha256: b7f80f0c98a07f08b49ae0940d5836410d9de32af6b2ea8a737c44e02b2c8685
hash_basis: working-tree bytes (LF)
---

# THM-3337 -- cross-shell compression gives `T(4)=5`

Let independent bits satisfy

```text
P(X_i=0)=p,       P(X_i=1)=q=1-p,       0<p<1,
```

and let the critical value of a nonconstant stream be

```text
n=min{k>=1:X_(k+1)!=X_1}.
```

For one extractor write `T(n)` for its worst stopping time over streams of
critical value `n`. The pointwise optimum is the infimum of `T(n)` over all
deterministic exactly fair extractors.

## 1. The finite replacement

On the branches `1<=n<=7`, use the following rule.

| `n` | stop at | heads exactly when |
|---:|---:|---|
| 1 | 2 | `X_1=0` |
| 2 | 3 | `X_1=1` |
| 3 | 5 | `X_5=0` |
| 4 | 5 | `X_1=0` |
| 5 | 7 | `X_7=1` |
| 6 | 7 | `X_1=1` |
| 7 | 8 | `X_1=0` |

If the first eight bits are constant, so `n>=8`, use the THM-2225 cyclic-
checksum extractor without modification. The cases are disjoint and every
displayed read occurs after the first disagreement, so this is one causal,
deterministic rule on every nonconstant stream. On `n>=8` it inherits the
proved finite THM-2225 deadline.

## 2. Exact fairness of the replacement

The head probability contributed by the seven displayed branches is

```text
H_<(p)
 = pq + q^2 p
   + (p^3q+q^3p)p
   + p^4q
   + (p^5q+q^5p)q
   + q^6p + p^7q.                                  (1)
```

The terms follow directly from the table. For example, on `n=3` either
prefix is `0001` or `1110`, followed by the required bit `X_5=0`, giving
`(p^3q+q^3p)p`. Expanding with `q=1-p` gives

```text
H_<(p)=(1-p^8-q^8)/2.                              (2)
```

There is also a composition-level proof: among all nonconstant eight-bit
words, the rule labels exactly

```text
4, 14, 28, 35, 28, 14, 4
```

heads in Hamming weights `1,...,7`, exactly half of
`binom(8,1),...,binom(8,7)`. This is verified exhaustively and independently
of floating point by the companion script.

THM-2225 bisects every composition class inside each dyadic shell. Its
aggregate head probability on the same event `n<=7` is therefore also the
right side of (2), half the probability that the first eight bits are not
constant. Let `R(p)` be the head probability of the unchanged THM-2225 rule
on `n>=8`. Exact fairness of that rule says

```text
(1-p^8-q^8)/2 + R(p)=1/2.                           (3)
```

The replacement changes the first term of (3) to an equal polynomial and
leaves `R(p)` unchanged. Hence the spliced extractor is exactly fair for every
`0<p<1`.

## 3. Why shell-by-shell search missed it

Refine every stopped prefix to length eight and, in each Hamming layer, write
`2*(heads)-(total)` for the doubled deficit. The shell `n=1` remains balanced,
while the other two dyadic shells have

```text
n in {2,3}:       (0,-2,-4, 0, 4, 2,0),
n in {4,5,6,7}:   (0, 2, 4, 0,-4,-2,0).             (4)
```

Thus neither shell is composition-balanced, but their defects cancel in
every global Hamming layer. THM-3032 proved `T(4)>=6` only inside the class of
shell-balanced rules and explicitly left cross-shell compensation open. The
new rule occupies exactly that missing class.

This is the finite analogue of retaining an indefinite block sidecar: the
positive and negative shell defects are useful separately and vanish only
after the correct global compression. The analogy to the 2026 zeta-zero
finite-compression argument is motivational; no analytic claim from that
paper is a dependency.

## 4. Optimality at four

No exactly fair online rule can stop on a constant prefix `b^k`. If it did,
one verdict would have probability at least `p^k` (for `b=0`) or `q^k` (for
`b=1`), exceeding `1/2` as the corresponding bias tends to one. Before flip
`n+1`, every critical-value-`n` stream still has a constant observed prefix.
Therefore every exactly fair extractor has

```text
T(n)>=n+1.                                           (5)
```

The displayed rule has `T(4)=5`, so (5) is sharp and the pointwise optimum is

```text
                         T_opt(4)=5.                 (6)
```

The same single extractor has the simultaneous finite profile
`(2,3,5,5,7,7,8)` on `1<=n<=7`; no claim of global Pareto optimality for that
whole vector is needed.

## 5. Exact audit

```bash
python 04-computation/amm12592_cross_shell_T4_floor_thm3337.py
python -O 04-computation/amm12592_cross_shell_T4_floor_thm3337.py
```

Both runs check causality, all 254 nonconstant length-eight words, every
Hamming layer, the two opposite shell-deficit vectors, the power-basis
identity (2), and an independent THM-2225 checksum-prefix control. QED.
