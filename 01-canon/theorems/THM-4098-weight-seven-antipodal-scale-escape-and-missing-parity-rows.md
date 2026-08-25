---
id: THM-4098
title: "Weight-seven antipodal scale escape and the missing parity rows"
status: >
  PROVED + PROVED RELATIVE TO THM-2061/2072/4092 + VERIFIED-EXACT +
  INDEPENDENTLY VERIFIED-EXACT. If a closed two-phase safe core contains an
  interval J of length L, every finite outlier bank of antipodal weight seven
  and cardinality below seven escapes its first-moment wall after any common
  dilation q with qL>=1. For odd q this is a weight-seven overlap theorem; for
  even q the half phase collapses to a sub-seven single-phase bank. Applied to
  the exact AP intervals of THM-4092, this supplies all previously missing
  AP7/AP6/AP5 parity rows at q>=55,35,35 and closes their dyadic two-tail seams
  through every further common dilation. This is not an arbitrary-core
  supplier or LRC(14).
source: codex-arithmetic-boundary-breakthrough-20260825
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
  - THM-4092-parity-weighted-antipodal-k-comb-density-compiler
related:
  - THM-1155-threespeed-exhaustive-and-ceiling
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-4081-lrc-antipodal-height12-obstruction-and-six-speed-floor
script: 04-computation/lrc_weight_seven_scale_escape_thm4098.py
output: 05-knowledge/results/lrc_weight_seven_scale_escape_thm4098.out
independent_audit_script: 04-computation/lrc_weight_seven_scale_escape_thm4098_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_weight_seven_scale_escape_thm4098_independent_audit.out
script_sha256: a6f77908408ae854cdc3e24f7dd7de0e330d032a9e9893a5d0cdce3a1daabb46
output_sha256: 00a2a184154ed668d69e18babad51f4a99cae16c9fa69b34cb6b85b5a8a9db47
independent_audit_script_sha256: 929edc49da3b0ac83979faceed33fd2b7c5b48e8b81a5b564158a723ce92deeb
independent_audit_output_sha256: eb33a5c20defc168d6200ddaccd1eaa95eddd031d8d3620f3bd272d649f40f9e
hash_basis: raw LF bytes
---

# THM-4098 -- weight-seven antipodal scale escape

Put `delta=1/14`. For a finite set `A` of positive integer speeds, retain the
closed two-phase safe set

```text
G_A^+-={theta in R/Z:
          ||a theta||>=delta and ||a(theta+1/2)||>=delta
          for every a in A}.                              (1)
```

THM-4092 proves a first-moment absorption criterion only when the total
antipodal weight is strictly below seven. The present theorem treats the
boundary `W=7` by using the scale coordinate discarded there. Its mechanism
is not a better first moment: overlap makes a base survivor, and dilation
turns that survivor into a periodic net.

## 1. The parity-sensitive pullback identity

For one speed put

```text
U_v={theta:min(||v theta||,||v(theta+1/2)||)<delta},
omega(v)=1 if v is even, and 2 if v is odd.                (2)
```

As in THM-4092, `U_v` has measure `omega(v)/7`. Also define the one-phase
closed safe set

```text
G_A^0={theta:||a theta||>=delta for every a in A}.          (3)
```

Let `mu_q(theta)=q theta mod 1`. Directly from the half phase,

```text
boxed:
G_(qA)^+- = mu_q^(-1)(G_A^+-)   if q is odd,
G_(qA)^+- = mu_q^(-1)(G_A^0)    if q is even.              (4)
```

Indeed

```text
q(theta+1/2) = q theta+1/2 mod 1   for odd q,
q(theta+1/2) = q theta       mod 1 for even q.             (5)
```

Thus odd dilation preserves the two base phases, while even dilation merges
them. Equation `(4)` is equality of closed sets, including every wall.

## 2. Every weight-seven bank has a positive base interval

Let `V` be a finite set of distinct positive speeds, put

```text
k=|V|,                 W=sum_(v in V) omega(v),             (6)
```

and assume

```text
boxed: k<7 and W=7.                                        (7)
```

Then both `G_V^+-` and `G_V^0` contain nondegenerate closed intervals.

For the two-phase assertion, subadditivity and `(2)` give

```text
sum_(v in V) mu(U_v)=W/7=1.                                (8)
```

This is the equality wall that defeats the THM-4092 first moment. But the
danger combs cannot be disjoint: if `B=max(V)`, every `U_v` contains the
open circular interval

```text
||theta||<1/(14B).                                         (9)
```

Choose two speeds in `V`; there are at least two because every individual
weight is at most two. Their intersection has positive measure. Hence

```text
mu(union_(v in V) U_v)
 <= sum_v mu(U_v)-mu(U_u intersect U_v)<1.                 (10)
```

The complement `G_V^+-` therefore has positive measure. It is the complement
of a finite union of open intervals, hence is a finite union of closed
intervals and points. Positive measure forces one interval of length
`ell_odd>0`.

For the one-phase assertion, each one-phase danger comb has measure `1/7`,
so

```text
mu((G_V^0)^c)<=k/7<1.                                     (11)
```

The same finite-wall argument supplies an interval of length
`ell_even>0`. Notice the division of labor: overlap repairs the odd branch
at exact weight seven, while ordinary subadditivity repairs the even branch
after its half-phase collapse.

## 3. A dilated interval becomes a fine safe net

Let `K` be any closed circular interval of length `ell>0`. The pullback
`mu_q^(-1)(K)` consists of `q` equally spaced closed intervals of length
`ell/q`. Every complementary gap has length

```text
(1-ell)/q < 1/q.                                           (12)
```

Consequently every closed circular interval `J` of length `L` meets this
pullback whenever

```text
boxed: qL>=1.                                              (13)
```

Indeed `J` is longer than every complementary gap and therefore cannot lie
in one. Equality in `(13)` is allowed because `(12)` is strict.

Combine Sections 1--3. In the odd case use the positive interval in
`G_V^+-`; in the even case use the positive interval in `G_V^0`. We obtain
the general scale-escape statement.

> **Theorem 3.1 (weight-seven scale escape).** Let `D,V` be finite sets of
> distinct positive speeds, with `D` disjoint from `qV`, and suppose
> `J subset G_D^+-` is a closed interval of length `L>0`. If `(7)` holds and
> `q` is any positive integer with `qL>=1`, then
>
> ```text
> boxed: J intersect G_(qV)^+- is nonempty,                (14)
> ```
>
> and hence `G_(D union qV)^+-` is nonempty.

The conclusion is uniform over the actual magnitudes, ratios, and residue
classes in `V`. The sufficient scale depends only on the inherited core
length. It is not asserted optimal for an individual bank.

## 4. The three missing AP parity rows

THM-4092 proves the following exact core intervals:

| core `D_m` | safe interval `J_m` | length `L_m` |
|---|---|---:|
| `{1,...,7}` | `[4/35,13/98]` | `9/490` |
| `{1,...,6}` | `[4/35,1/7]` | `1/35` |
| `{1,...,5}` | `[4/35,1/7]` | `1/35` |

Let `k=11-m` and take any `k` distinct positive base speeds `V`. The rows
missing from the strict first-moment table are exactly

| core | `k` | required odd members of `V` | `W=k+o` | sufficient `q` |
|---|---:|---:|---:|---:|
| AP7 | 4 | 3 | 7 | **55** |
| AP6 | 5 | 2 | 7 | **35** |
| AP5 | 6 | 1 | 7 | **35** |

The last column is the first integer satisfying `(13)`:

```text
55*(9/490)>1>54*(9/490),
35*(1/35)=1>34*(1/35).                                    (15)
```

At those scales every member of `qV` exceeds the AP core, so the eleven
speeds are distinct. Theorem 3.1 gives

```text
boxed:
G_({1,...,m} union qV)^+- is nonempty                     (16)
```

for every bank in the displayed row and every integer `q` at least the
displayed threshold. For odd `q`, parity is preserved, so `(16)` literally
fills each missing weight-seven row. For even `q`, all scaled outliers become
even; `(16)` remains true through the separate one-phase branch of `(4)`.

This is genuinely beyond the selected-interval first moment. At unit scale,
the AP7 interval is literally covered by the THM-4092 hostile bank

```text
V={8,9,11,13},       W=1+2+2+2=7.                         (17)
```

After every `q>=55`, Theorem 3.1 certifies a survivor for the same shape
`qV`. The exact referees also replay `q=55` and `q=56`, separating the odd
and even pullback branches. No claim is made that `55` is the first successful
dilation for `(17)`; it is the uniform all-bank threshold.

## 5. Adaptive endpoint and dyadic-seam consequence

Put

```text
C={1,...,m} union qV                                      (18)
```

in any row of Section 4. The survivor set inside `J_m` is a finite union of
closed intervals and points. Choose a component endpoint. It is either a
core-interval endpoint or a literal tooth endpoint of a scaled outlier.
Exactly as in THM-4092, it has an even presentation

```text
theta=r/N,             N<=14q max(V).                     (19)
```

At the stated thresholds the right side dominates the inherited endpoint
denominators `98` and `70`. Both `theta` and `theta+1/2` lie in the weak-safe
set `G_C` of THM-2061. THM-2072's antipodal-safe-pair certificate therefore
gives, for every pair of distinct positive odd tails `x,y`,

```text
boxed: 2C union {x,y} is 1/14-lonely.                     (20)
```

The certificate transports through every further common dilation `d>=1`.
Use `theta/d`: at phase zero the old products are unchanged; at the half
phase odd `d` preserves the old half phase and even `d` collapses it to the
old zero phase. Thus

```text
boxed:
2d({1,...,m} union qV) union {x,y} is 1/14-lonely          (21)
```

for every `d>=1`. This is a new infinite structured supplier for the dyadic
two-tail seam. It is not a claim that an arbitrary eleven-core has an AP
subcore or a weight-seven dilation representation.

## 6. Exact audits and loss map

Primary reproduction:

```bash
python3 -B 04-computation/lrc_weight_seven_scale_escape_thm4098.py
python3 -B -O 04-computation/lrc_weight_seven_scale_escape_thm4098.py
```

The primary referee merges literal open intervals over the declared universe
of all row-compatible banks drawn from `{1,...,10}`. It checks `155` base
banks, `310` scaled rows at the first sufficient scale and its successor,
and `10,108` owner-complement gates. The smallest direct uncovered margin is
`29/7700`. It also verifies the exact q=1 hostile `(17)`.

Independent reproduction:

```bash
python3 -B 04-computation/lrc_weight_seven_scale_escape_thm4098_independent_audit.py
python3 -B -O 04-computation/lrc_weight_seven_scale_escape_thm4098_independent_audit.py
```

The independent referee does not merge danger intervals. It partitions by
all rational walls and classifies walls and open cells directly. On its
separate `{1,...,9}` bank universe it checks `80` banks, `240` target cells,
and `66,240` exact pointwise instances of both pullback branches in `(4)`.
Normal and optimized transcripts are byte-identical for both referees; every
gate is an explicit exception check and neither script uses floating point.

The typed connection is

```text
source:        one positive two-phase core interval J and a W=7 bank V
operation:     common dilation V |-> qV
map:           multiplication mu_q on the circle
preserved:     target two-phase safety and closed wall conventions
destroyed:     the base half-phase distinction when q is even
needed sidecar:dilation parity plus one positive base safe interval
recovery:      the q-periodic preimage net, then an arrangement endpoint
hostile:       AP7 with {8,9,11,13} at q=1
target:        an antipodal safe pair for an eleven-speed dyadic core.       (22)
```

The theorem crosses the exact `W=7` wall only after adding the scale
coordinate. It does not cover `W>7`, remove the requirement of a positive
core interval, classify low dilations, imply physical entry for an arbitrary
LRC residual, or prove LRC(14). **QED.**
