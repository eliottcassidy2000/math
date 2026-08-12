---
id: THM-3338
title: "Horizon-16 cross-shell surgery closes every pointwise AMM 12592 value through n=15"
status: >
  PROVED + VERIFIED-EXACT. There is one deterministic exactly fair unknown-
  bias coin extractor with deadline vector (2,3,4,5,15,7,8,9,10,11,12,13,
  14,15,16) on critical values 1 through 15, agreeing with THM-2225 for
  n>=16. Thus it stops at the information floor n+1 simultaneously for all
  n<=15 except n=5. THM-2253 supplies a different floor-attaining rule at
  n=5, so the pointwise optimum is T_opt(n)=n+1 for every 1<=n<=15. In
  particular the previously unsettled even values 8,10,12 close. Exact
  length-16 Hamming-layer bisection proves fairness. The lower and upper
  dyadic shells are not separately balanced: their doubled deficit vectors
  are nonzero exact opposites. The companion standard-library audit checks
  all 65536 words, causality, the Bernoulli identity, and an independent
  THM-2225 checksum-prefix control.
source: codex-kps-2026-08-12-anthropic-finite-compression-transfer
depends_on:
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-2253-dyadic-contrast-extractor
  - THM-3337-cross-shell-compression-attains-the-T4-floor
related:
  - THM-3032-sharpened-half-tail-extractor-and-shell-four-pareto-frontier
  - "Claude, More than two thirds of the zeros of the Riemann zeta function lie on the critical line (2026), finite-compression inspiration only"
script: 04-computation/amm12592_horizon16_floor_thm3338.py
output: 05-knowledge/results/amm12592_horizon16_floor_thm3338.out
script_sha256: 01fbcab5f1d3406fbad837a4c22e2e07c7116402dfc154be8b38225191449f9f
output_sha256: efdc09355398216c6001e2ab68f84c8a46ec4cbcd1cb347e976b6cace9de76b9
hash_basis: working-tree bytes (LF)
---

# THM-3338 -- horizon-16 cross-shell surgery

Let independent bits satisfy

```text
P(X_i=0)=p,       P(X_i=1)=q=1-p,       0<p<1,
```

and let the critical value of a nonconstant stream be

```text
n=min{k>=1:X_(k+1)!=X_1}.
```

For one extractor, `T(n)` is its worst stopping time over streams with
critical value `n`.

## 1. A finite-prefix surgery lemma

Fix a horizon `M`. Every nonconstant length-`M` word has a unique critical
value `n<M`. Suppose a deterministic causal rule on these branches labels
exactly half of every nonconstant Hamming layer as heads:

```text
#{x in {0,1}^M: |x|=k and H(x)} = binom(M,k)/2,   1<=k<M.       (1)
```

For `M` a power of two, all the right sides are integers by Lucas. The head
probability of the finite rule is then

```text
sum_(k=1)^(M-1) [binom(M,k)/2] p^(M-k)q^k
  = (1-p^M-q^M)/2.                                           (2)
```

THM-2225's shell-balanced rule has the same aggregate probability on the
event that the first `M` bits are nonconstant. Replace its verdicts on that
event by the finite rule and leave its continuation from `0^M` and `1^M`
unchanged. Equation (2) shows that the total head probability remains
exactly `1/2` for every `p`. Notice that no conditional-fairness assertion on
the two constant prefixes is used.

This is the finite-prefix surgery principle. The present theorem is its
`M=16` instance.

## 2. The rule away from the donor branch

For every `n` in `{1,...,15}\{5}`, stop as soon as the first disagreement is
seen, at flip `n+1`. Declare heads exactly for the displayed initial bit:

| `n` | 1 | 2 | 3 | 4 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heads for `X_1` | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 1 | 1 | both | neither | neither | both | 0 |

The words `both` and `neither` mean exactly what they say; the verdict still
occurs after the disagreement, not on a constant prefix.

## 3. The single donor branch

If `n=5`, the first six bits are `000001` or `111110`. Read the next nine
bits, positions 7 through 15. Let `k` be their number of ones and rank the
nine-bit continuation lexicographically among all continuations with that
weight, starting at rank zero. Declare heads iff its rank is smaller than
the corresponding entry below:

```text
initial bit 0:  (1,9,19,25,10,  0, 0, 9,9,0),
initial bit 1:  (1,0,27,84,126,116,58,15,0,0),                 (3)
```

where entries are indexed by `k=0,...,9`. Every entry lies between zero and
`binom(9,k)`, so (3) defines literal subsets of continuation words. This
branch stops at flip 15. If no disagreement has appeared by flip 16, retain
THM-2225 from the current constant prefix.

All cases are disjoint. Every verdict is a function of the prefix available
at its stated deadline, so the spliced rule is deterministic and causal.

## 4. Exact layer bisection

Refine each stopped prefix by all free completions to length 16. Direct
integer counting gives:

| weight `k` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heads | 8 | 60 | 280 | 910 | 2184 | 4004 | 5720 | 6435 | 5720 | 4004 | 2184 | 910 | 280 | 60 | 8 |

These are exactly `binom(16,k)/2`. Thus (1) holds, the finite-prefix surgery
lemma applies, and the infinite spliced extractor is exactly fair.

The cancellation is genuinely cross-shell. In weights `1,...,15`, the
doubled deficits `2*(heads)-(total)` are zero on shells `{1}` and `{2,3}`.
On the remaining two shells they are

```text
{4,5,6,7}:
  (0, 4, 2,-8,-16,-14,-6,0, 6,14,14, 2,-6,-4,0),

{8,9,10,11,12,13,14,15}:
  (0,-4,-2, 8, 16, 14, 6,0,-6,-14,-14,-2, 6, 4,0).           (4)
```

Neither shell is balanced, but (4) cancels coordinatewise. This is the same
missing degree of freedom exposed at horizon eight by THM-3337, now acting
between two larger dyadic shells.

## 5. Deadlines and pointwise optima

The finite profile is

```text
(T(1),...,T(15))=(2,3,4,5,15,7,8,9,10,11,12,13,14,15,16).   (5)
```

No exactly fair online rule can stop on a constant prefix. Otherwise one
verdict has probability at least `p^r` or `q^r`, exceeding `1/2` as the
corresponding bias tends to one. Hence every rule obeys

```text
T(n)>=n+1.                                                    (6)
```

Equations (5)--(6) prove `T_opt(n)=n+1` simultaneously for every
`1<=n<=15`, `n!=5`. THM-2253 attains the same floor at every odd `n`, in
particular at `n=5`. Therefore

```text
                    T_opt(n)=n+1,       1<=n<=15.             (7)
```

This closes the previously open even values `n=8,10,12`; `n=4` is also an
independent second floor witness after THM-3337. The first pointwise value
not settled by (7) and the earlier odd/near-shell-top constructions is now
`n=16`.

**Later supersession.** The last sentence was the frontier at this theorem's
checkpoint. THM-3340 now proves `T_opt(n)=n+1` for every positive integer;
there is no unresolved pointwise value. The simultaneous horizon-16 profile
and opposite-shell defect cancellation proved here remain strictly stronger
finite Pareto information than that pointwise statement.

## 6. Exact audit and scope

```bash
python 04-computation/amm12592_horizon16_floor_thm3338.py
python -O 04-computation/amm12592_horizon16_floor_thm3338.py
```

Both runs use standard-library integer arithmetic. They check all `65536`
length-16 words, every lexicographic rank and count box in (3), causality at
every deadline, all fifteen layer equalities, the Bernoulli polynomial
identity, (4), and an independently implemented THM-2225 checksum prefix.

The finite integer search used to discover (3) is not a proof dependency.
The theorem is the displayed rule plus the exact elementary audit. No claim
is made that (5) is globally Pareto-optimal or that the donor cost at `n=5`
is minimal. QED.
