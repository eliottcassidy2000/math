---
id: THM-3601
title: "AMM R=2048 terminal Rule-A failure and adjoint horizon"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
  AUDIT.  At the exact golden-floor epoch R=2048, the last failing Rule-A
  offset D0=37 dies at row 508.  A positive integer truncated-adjoint
  certificate is nonnegative at cut 195 and negative at cut 196, so every
  admissible continuation which survives this fatal inequality must differ
  from Rule A by row 195.  The adjacent offset D0=38 survives the fixed Rule-A
  epoch.  This is not alternative feasibility, a uniform extractor, or an
  AMM constant.
source: kps-s188 / THM-3597 dyadic adjoint continuation, 2026-08-21
audit: >
  Author exact audit only.  The optimization-safe standard-library companion
  hash-pins THM-3597, rebuilds the golden degree words and Rule-A traces,
  performs the complete one-sweep adjoint ledger, proves its single strict
  sign wall and monotonicity, and independently compares the sweep with the
  legacy coefficient-cancellation certificate at cuts 195 and 196.  Ordinary
  and optimized replay are byte-identical to the stored output.  Independent
  hostile audit remains pending.
depends_on:
  - THM-3597-amm-dyadic-rule-a-transition-and-adjoint-horizons-through-R1024
related:
  - THM-3588-amm-r512-truncated-adjoint-pascal-repair-horizons
  - THM-3577-amm-r512-offset-transition-and-causal-horizon
  - THM-3330-growth-law-y-box-lucas-diagonal-and-attainment-to-2047
  - THM-3331-elin-conditional-theorem-eps-star-and-the-superblock-law
script: 04-computation/amm12592_r2048_terminal_failure_adjoint_horizon_thm3601.py
output: 05-knowledge/results/amm12592_r2048_terminal_failure_adjoint_horizon_thm3601.out
script_sha256: 071db60cd42ce1d1fff7c693e2b2baeeb782046cd93169401636022da06b0205
output_sha256: 51caae743b5fa55a4c30840baf83232ddd62626c16607179c2e1b9f129d79dc8
semantic_sha256: aac4e238138a4acde9618854eee894d08131830e4f6745dd9351432c434c1247
hash_basis: LF-normalized bytes
---

# THM-3601 -- AMM R=2048 terminal Rule-A failure and adjoint horizon

**PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
AUDIT.**  The dyadic adjoint atlas reaches its first scale beyond the finite
range of THM-3597.  The new point also exposes a possible golden-ratio scaling
coordinate, recorded only as a hypothesis in Section 4.

## 1. Exact terminal pair

For `0<=i<2048`, put

```text
d_i(D0)=floor(log_5(phi^(2(2048+i))))+D0,              (1)
```

with the floor decided by the exact Fibonacci--Lucas comparison inherited
from THM-3597.  On these degree words, the fixed integer Rule A has the exact
adjacent behavior

```text
D0=37: first death j=508, (d_0,d_j,d_2047)=(1261,1565,2485),
D0=38: survives all 2048 rows, (d_0,d_2047)=(1262,2486).             (2)
```

Thus `37` is the terminal failing offset and `38` the first surviving offset
within the already-established consecutive Rule-A scan.  Neither assertion
says that a different admissible policy fails or succeeds.

## 2. Exact positive adjoint wall

Use the halved junk states, Lucas box bounds and Pascal transitions of
THM-3597.  If an admissible alternative agrees with the failing Rule-A trace
through row `s-1`, transposed propagation of the fatal-top evaluation yields
nonnegative integer multipliers and the necessary inequality

```text
0 <= B_s,
B_s=load^A_(508,d_508)
    +sum_(i=s)^507 sum_t lambda_(i,t)(u^A_(i,t)-L_(i,t)).             (3)
```

The complete backward sweep proves

```text
B_195 > 0 > B_196,
B_s >= B_(s+1) for 1<=s<507,                                      (4)
```

and its negative cuts are exactly `196,...,507`.  Therefore agreement
through row `195` would imply `0<=B_196<0`.  Every surviving admissible
alternative must depart from Rule A at or before row `195`.

At the first negative cut the multiplier metadata is

```text
active cells = 48828,
multiplier mass =
383906924668156805798067812893197399113080992965104527379505084660104758411637611723600286539507103880742384734627445584706149259400431166642522540376,
maximum multiplier =
10248578385150114262857028294077844120966011618655540198805391799146850862676148277125424815376100473925059468454633373051883752561439644748608883448. (5)
```

The two boundary signs have bit lengths `(1211,1212)` and joint digest

```text
f38b3838b74bc6ab557fd6919efcb42901350073745abd16d95d9f4fe7d2ab95. (6)
```

The companion separately assembles the older coefficient-cancellation
certificate at cuts `195` and `196`; both values and all multiplier metadata
agree exactly with the one-sweep ledger.  This is a hostile indexing control,
not a second mathematical proof of the inherited inequalities.

## 3. Why this is a proof obligation, not an AMM bound

The three relevant quantities remain distinct:

```text
D0*=38       first offset where the fixed Rule A survives,
j=508        first death row of the adjacent failing trace,
q=195        latest row by which every surviving repair must differ.       (7)
```

The adjoint certificate controls only continuations agreeing with that trace
up to a cut.  It does not construct a continuation after departure, prove
infeasibility of the full entry polytope, or show that an early departure
enters THM-3331's superblock basin.  The missing theorem is still a typed
overlap statement between the dual departure obligation and a primal
absorbing region.

## 4. A sharpened scaling target, not a theorem

For the hardest failed offsets at the four audited dyadic epochs, `(R,j,q)` is

```text
(256,61,24), (512,121,44), (1024,250,94), (2048,508,195).             (8)
```

At the new scale

```text
j/R = 0.248046875,
q/j = 0.3838582677...,
q/R = 0.09521484375.                                                 (9)
```

These are close respectively to `1/4`, `phi^-2=0.381966...`, and
`phi^-2/4=0.0954915...`.  The factorization is structurally attractive:
`phi^-2` is already the exact tangency coordinate in the AMM archimedean
capacity calculation.  Four scales, however, do not prove a limit, and the
`R=512` ratios are visibly nonmonotone.  The pre-registered hostile test is
the terminal failed offset at `R=4096`; only after that test should one seek a
continuum adjoint equation whose zero lies at `theta=phi^-2/4`.

## 5. Verification and scope

Reproduce with

```bash
python3 04-computation/amm12592_r2048_terminal_failure_adjoint_horizon_thm3601.py
python3 -O 04-computation/amm12592_r2048_terminal_failure_adjoint_horizon_thm3601.py
```

The standard-library companion hash-pins the exact THM-3597 implementation,
then reconstructs `(1)--(6)` rather than trusting archived state.  Its pinned
semantic digest is

```text
aac4e238138a4acde9618854eee894d08131830e4f6745dd9351432c434c1247.
```

This theorem is confined to `R=2048`, offsets `37,38`, and the fixed Rule A.
It proves no alternative feasibility, no asymptotic limit in Section 4, no
uniform extractor, and no value or improved upper bound for the AMM constant
`C*`.

**QED.**
