---
id: THM-3577
title: "AMM R=512 offset transition and causal horizon"
status: >
  PROVED + VERIFIED-EXACT.  At the exact golden-floor epoch R=512, Rule A
  fails before the first feed-free F1--F3 entry for offsets D0=0,...,4 and
  succeeds at D0=5.  The death rows are 107,110,113,116,121; the minimal
  coefficient-capacity horizons are 58,58,60,60,62, forcing any surviving
  alternative to differ by rows 49,52,53,56,59 respectively.  At D0=5 the
  row-126 state of degree 386 enters the row-127 degree-387 cone with support
  18 and exact F2/F3 margins 14 and 1229.  This is a fixed-policy finite
  atlas, not an alternative-prefix infeasibility or uniform-extractor result.
source: kps-s188 + subagent-laplace
depends_on:
  - THM-3374-amm-r512-rule-a-causal-repair-horizon
related:
  - THM-3332-s-cone-one-row-certificate-32768-closures-and-the-256-floor-point
companion: 04-computation/amm12592_r512_offset_transition_causal_horizon_kps_s188.py
output: 05-knowledge/results/amm12592_r512_offset_transition_causal_horizon_kps_s188.out
script_sha256: c782f3236b5f3815a98b8b78bf55481c6e26a3d8d431ec8b40fbab42ed7c5210
output_sha256: b34896bb94b83a1f9d14e296bd4b5f7aefb1634dc412d93abf069843b193017d
hash_basis: LF-normalized bytes
---

# THM-3577 -- AMM R=512 offset transition and causal horizon

**PROVED + VERIFIED-EXACT.**  The first six integer translates of the golden
deadline word have a sharp Rule-A transition at `D0=4/5`.  Every quantity in
the result is reconstructed with integer arithmetic; no LP solver status or
floating candidate is used.

## 1. Exact profiles and policy

Fix `R=512`.  For `0<=i<R`, let

```text
d_i(D0)=floor(log_5(phi^(2(R+i))))+D0,                 (1)
```

where `phi=(1+sqrt(5))/2`.  The floor in `(1)` is decided without floating
point by the Fibonacci--Lucas comparison

```text
5^d <= phi^(2m).                                       (2)
```

Run the inherited integer Rule A on this degree word.  At each row it clamps
every reversed-Bernstein coefficient to the exact Lucas box and retains the
halved junk row.  The forced top coefficient has no state variable: it must
lie in its two-point admissible interval.  A nonzero residual there is the
death certificate.

Let `e(D0)` be the first feed-free row.  The state tested for THM-3332 entry
is row `e-1`; when the degree jumps, its degree and the next cone degree are
different.  Keeping these two numbers separate removes the apparent
`386/387` off-by-one at `D0=5`.

## 2. The six-offset atlas

Exact replay gives:

| `D0` | first feed-free row | state/next degree | Rule-A death | capacity horizon | any survivor differs by row |
|---:|---:|---:|---:|---:|---:|
| 0 | 130 | 383/383 | 107 | 58 | `<=49` |
| 1 | 129 | 383/384 | 110 | 58 | `<=52` |
| 2 | 128 | 384/384 | 113 | 60 | `<=53` |
| 3 | 128 | 385/385 | 116 | 60 | `<=56` |
| 4 | 127 | 385/386 | 121 | 62 | `<=59` |
| 5 | 127 | 386/387 | survives | -- | -- |

Thus Rule A fails for every `D0<=4` in the tested cell and reaches the
certified tail cone at `D0=5`.

For the surviving `D0=5` row-126 state, the exact THM-3332 data are

```text
support_max=18,
F2_margin=14,
F3_min_margin=1229.                                   (3)
```

In particular all three entry conditions hold with positive nontrivial
margins.

## 3. Causal repair horizons

Suppose Rule A dies at row `j` with negative forced-top residual `-N`.  If an
alternative admissible prefix agrees with Rule A until only the preceding
`L` rows may change, THM-3374's coefficient-capacity lemma bounds its possible
change at the fatal coordinate by

```text
C_L=sum_(s=1)^L 2^(s+1) binom(d_(j-s),s).              (4)
```

For each failing offset, the horizon in the table is the least `h` with

```text
C_(h-1)<N<=C_h.                                       (5)
```

The exact fatal residuals and both integers in every bracket `(5)` are
printed in the stored output.  Consequently an alternative agreeing through
row `j-h` could change only the last `h-1` rows and cannot survive.  It must
differ at some row at most `j-h`, giving precisely `49,52,53,56,59`.

The upper inequality in `(5)` is only the point where the triangle-inequality
obstruction becomes inconclusive.  It does not construct a repair.  In
particular, the `D0=4` row says that a viable alternative must change the
policy by row 59; no late clamp near the row-121 death can decide feasibility.

## 4. Scope and next exact target

This theorem concerns one explicit policy at one finite epoch and six integer
offsets.  It proves neither:

- infeasibility of another causal prefix at `D0<=4`;
- feasibility of the full entry polytope at `D0=4`;
- a single extractor valid at every scale;
- a deadline slope below two; nor
- the value of the uniform AMM constant `C*`.

Its constructive use is to shrink the honest `D0=4` search: an active-face,
integer-lattice, or Farkas computation must expose changes no later than row
59.  Its hostile use is the exact `D0=5` positive control; any claimed exact
solver upgrade must recover or certify that known feasible entry state.

## 5. Exact verification

Run

```bash
python3 04-computation/amm12592_r512_offset_transition_causal_horizon_kps_s188.py
python3 -O 04-computation/amm12592_r512_offset_transition_causal_horizon_kps_s188.py
```

The standard-library companion reconstructs `(1)--(2)`, all feeds and Lucas
caps, every Rule-A row, the five deaths, all ten capacity brackets in `(5)`,
and the three entry margins `(3)`.  Ordinary and optimized transcripts agree.

**QED.**
