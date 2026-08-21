---
id: THM-3597
title: "AMM dyadic Rule-A transition and adjoint horizons through R=1024"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
  AUDIT.  At the exact golden-floor epochs R=256,512,1024, Rule A first
  survives among the tested consecutive nonnegative offsets at D0=1,5,15.
  Every smaller tested offset dies at the listed row, and positive integer
  truncated-adjoint Pascal multipliers force any surviving alternative to
  depart from that Rule-A prefix by the listed earlier row.  This is a
  fixed-policy finite atlas, not alternative feasibility, a uniform
  extractor, or an AMM constant.
source: kps-s188 / AMM adjoint scaling continuation, 2026-08-21
audit: >
  Author exact audit only.  The optimization-safe standard-library companion
  reconstructs all golden floors, Rule-A rows, 21 deaths, three surviving
  prefixes, positive dual multipliers, strict sign boundaries and monotonic
  walls.  Ordinary and optimized replay are byte-identical to the stored
  output.  Independent hostile audit remains pending.
depends_on:
  - THM-3588-amm-r512-truncated-adjoint-pascal-repair-horizons
related:
  - THM-3577-amm-r512-offset-transition-and-causal-horizon
  - THM-3330-growth-law-y-box-lucas-diagonal-and-attainment-to-2047
  - THM-3331-elin-conditional-theorem-eps-star-and-the-superblock-law
script: 04-computation/amm12592_dyadic_rule_a_adjoint_atlas_thm3597.py
output: 05-knowledge/results/amm12592_dyadic_rule_a_adjoint_atlas_thm3597.out
script_sha256: bf6a429fb9785a02f3433d1bd2cdca8f3b8f2872bbbdc61ea8f6bcf514d5cc1e
output_sha256: fcd21410dc86a5c974b13ce30f5d889e7c38b2042a056111a1db6f411f274c21
semantic_sha256: 4169c97f724ae4a510dbc2b258f2f0039387d7758a552ee783d30d6b24b29fbe
hash_basis: LF-normalized bytes
---

# THM-3597 -- AMM dyadic Rule-A transition and adjoint horizons through 1024

**PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
AUDIT.**  THM-3588's positive dual certificate is now reconstructed at three
dyadic epochs.  The result separates two coordinates which must not be
conflated: the row where Rule A first dies and the much earlier row by which
any surviving alternative must already have departed from it.

## 1. Exact finite atlas

For `R in {256,512,1024}` and `0<=i<R`, put

```text
d_i(D0)=floor(log_5(phi^(2(R+i))))+D0,                 (1)
```

where the floor is decided exactly by the Fibonacci--Lucas comparison
`5^d<=phi^(2m)`.  Run the inherited integer Rule A on this degree word.  The
first consecutive nonnegative offset at which the policy survives all `R`
rows, together with the endpoint degrees of that surviving word, is

| `R` | failed offsets | first surviving offset | survivor degree endpoints |
|---:|---:|---:|---:|
| 256 | `0` | 1 | `154,306` |
| 512 | `0,...,4` | 5 | `311,616` |
| 1024 | `0,...,14` | 15 | `627,1239` |

For every failed offset, let `j` be its first Rule-A death row.  The exact
death and adjoint-horizon atlas is:

| `R` | `D0` | death `j` | first negative cut `s` | every survivor differs by row `<=s-1` |
|---:|---:|---:|---:|---:|
| 256 | 0 | 61 | 25 | 24 |
| 512 | 0 | 107 | 36 | 35 |
| 512 | 1 | 110 | 39 | 38 |
| 512 | 2 | 113 | 40 | 39 |
| 512 | 3 | 116 | 43 | 42 |
| 512 | 4 | 121 | 45 | 44 |
| 1024 | 0 | 207 | 64 | 63 |
| 1024 | 1 | 209 | 66 | 65 |
| 1024 | 2 | 211 | 67 | 66 |
| 1024 | 3 | 213 | 69 | 68 |
| 1024 | 4 | 216 | 70 | 69 |
| 1024 | 5 | 219 | 72 | 71 |
| 1024 | 6 | 221 | 74 | 73 |
| 1024 | 7 | 223 | 76 | 75 |
| 1024 | 8 | 227 | 78 | 77 |
| 1024 | 9 | 229 | 81 | 80 |
| 1024 | 10 | 233 | 82 | 81 |
| 1024 | 11 | 236 | 84 | 83 |
| 1024 | 12 | 240 | 88 | 87 |
| 1024 | 13 | 244 | 90 | 89 |
| 1024 | 14 | 250 | 95 | 94 |

The `R=512` row agrees exactly with independently audited THM-3588.  The
other two epochs are new finite continuations of the same mechanism.

## 2. Positive truncated-adjoint wall

Use THM-3588's halved junk states `y_(i,t)`, corrections `u_(i,t)`, Lucas
box bounds `L_(i,t)<=u_(i,t)<=U_(i,t)`, and transition kernels

```text
K_i(x)=(1+x)^(1+d_i-d_(i-1)).                          (2)
```

If an admissible alternative agrees with Rule A through row `s-1`, its
difference state `z` satisfies the one-sided cell inequalities

```text
z_(i,t)-(K_i z_(i-1))_t <= u^A_(i,t)-L_(i,t)           (3)
```

from the cut onward, together with the fatal-top lower inequality at row
`j`.  Propagate that top evaluation backwards by the transposed Pascal
kernels, discarding ancestry outside the actual state width.  The resulting
nonnegative integer multipliers cancel every free `z` coefficient and give

```text
0 <= B_s,
B_s=load^A_(j,d_j)+
    sum_(i=s)^(j-1) sum_t lambda_(i,t)(u^A_(i,t)-L_(i,t)). (4)
```

For every row in the table, the companion proves

```text
B_(s-1)>0>B_s.                                         (5)
```

Moreover

```text
B_s-B_(s+1)=sum_t lambda_(s,t)(u^A_(s,t)-L_(s,t))>=0. (6)
```

Thus all later cuts remain negative.  If an alternative agreed through row
`s-1`, `(4)--(5)` would say `0<=B_s<0`; it must differ by row `s-1` or
earlier.  The positive value `B_(s-1)` is only a hostile boundary for this
dual family and is not a feasible-prefix witness.

## 3. What scales and what does not

The three exact departure vectors are

```text
R=256:  (24),
R=512:  (35,38,39,42,44),
R=1024: (63,65,66,68,69,71,73,75,77,80,81,83,87,89,94). (7)
```

They are proof-obligation horizons, not death predictions.  The archived
all-scale transient work in THM-3330/3331 supplies a different coordinate:
front/stall dynamics for this fixed greedy policy.  Neither coordinate alone
proves infeasibility of the full entry polytope.  In particular, an adjoint
wall says where a successful policy must cease agreeing with Rule A; a
superblock certificate says why a specified continuation class then dies.
Their lawful composition requires an overlap theorem on the same state and
inequality class, which is not asserted here.

## 4. Verification and exact pins

Reproduce with

```bash
python3 04-computation/amm12592_dyadic_rule_a_adjoint_atlas_thm3597.py
python3 -O 04-computation/amm12592_dyadic_rule_a_adjoint_atlas_thm3597.py
```

The standard-library companion hash-pins THM-3588's parent implementation,
then independently rebuilds `(1)--(6)`.  For each of the 21 failed offsets it
pins the two boundary signs by their bit lengths and SHA-256 digest, records
the active dual-cell count, multiplier mass and maximum multiplier, and
checks monotonicity before the wall.  It also replays the three surviving
prefixes.  The complete integer ledger has semantic digest

```text
4169c97f724ae4a510dbc2b258f2f0039387d7758a552ee783d30d6b24b29fbe.
```

This theorem proves no feasible alternative, no infeasibility before the
displayed horizon, no rule uniform in `R`, no characteristic-free asymptotic
claim, and no value or bound for the uniform AMM constant `C*`.

**QED.**
