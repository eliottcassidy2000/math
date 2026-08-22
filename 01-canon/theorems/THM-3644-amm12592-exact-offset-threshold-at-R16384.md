---
id: THM-3644
title: "AMM12592 exact Rule-A offset transition at R=16384"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
  ADJACENT-AUDITED.  For the fixed golden-floor Rule-A policy at R=16384,
  every offset D0=401,...,412 dies and every offset D0=413,...,416 closes.
  Thus D0=413 is the first closure inside this exact consecutive bracket.
  No monotonicity below D0=401, global threshold, alternative feasibility,
  asymptotic law, or AMM bound is claimed.
source: kps-s189 + agent Harvey / exact clamped-Pascal continuation, 2026-08-21
audit: >
  The full 16-row durable ledger was checked against pinned exact degree
  coordinates and status semantics.  An independent FLINT/fmpz
  implementation reproduced the stable fields at offsets 408,410,411,412,
  413, including the load-bearing adjacent death/closure pair.  Live small
  controls independently exercise both terminal branches of the pinned
  engine in ordinary and optimized Python modes.
depends_on:
  - THM-3597-amm-dyadic-rule-a-transition-and-adjoint-horizons
  - THM-3633-amm-r16384-fixed-failed-trace-adjoint-phase-shock
related:
  - THM-3626-amm-r8192-adjoint-horizon-and-phase-rebound
script: 04-computation/amm12592_R16384_offset_transition_thm3644.py
output: 05-knowledge/results/amm12592_R16384_offset_transition_thm3644.out
script_sha256: ddf248e8939bbdefa7b2544bbd3df1c23e47e53e00aae4013d09149eae2ca59c
output_sha256: 28cedef2dc0b176b62f0633d93ea23a7c316993ed197743bd54f4466eb860c21
semantic_sha256: d7d96200feff9a6e453da08f6716fb60d34ab2f8f1cd53998c5bcd0af5a92df4
hash_basis: raw LF bytes for files; canonical JSON for semantic ledger
---

# THM-3644 -- AMM12592 exact Rule-A offset transition at `R=16384`

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
ADJACENT-AUDITED.**  This is a local transition theorem for one fixed policy
and one exact consecutive bracket.  It deliberately does not promote the
observed transition to a global offset threshold.

## 1. Fixed profile and exact transition

For `0<=i<16384`, set

```text
d_i(D0)=floor(log_5(phi^(2(16384+i))))+D0.              (1)
```

The floor in (1) is decided by the inherited exact Fibonacci--Lucas
comparison, never by floating point.  Running the pinned clamped-Pascal
Rule-A engine gives

```text
D0=401,...,412 : DIE,
D0=413,...,416 : CLOSED,
OPEN_RESIDUAL   : none.                                 (2)
```

Consequently `D0=413` is the first closing offset **inside the tested bracket
`401<=D0<=416`**.  THM-3633 separately proves that `D0=400` dies, but neither
the present computation nor its dependencies prove that closure is monotone
as a function of `D0`.  In particular, (2) is not a global least-offset claim.

## 2. Complete finite ledger

The columns are the offset, terminal status, terminal event row, degree at
the start/event/last row, number of `-2` rows, the `T6b` bound, and overflow
bit length when the policy dies.

| `D0` | status | event | `(d_start,d_event,d_last)` | `-2` rows | `T6b` | overflow bits |
|---:|:---|---:|:---|---:|---:|---:|
| 401 | DIE | 4058 | `(10198,12625,19995)` | 2130 | 4014 | 9879 |
| 402 | DIE | 4061 | `(10199,12627,19996)` | 2120 | 4016 | 9881 |
| 403 | DIE | 4066 | `(10200,12631,19997)` | 2138 | 4018 | 9882 |
| 404 | DIE | 4072 | `(10201,12636,19998)` | 2132 | 4020 | 9887 |
| 405 | DIE | 4075 | `(10202,12639,19999)` | 2147 | 4022 | 9893 |
| 406 | DIE | 4077 | `(10203,12641,20000)` | 2138 | 4024 | 9893 |
| 407 | DIE | 4080 | `(10204,12644,20001)` | 2153 | 4026 | 9897 |
| 408 | DIE | 4082 | `(10205,12646,20002)` | 2143 | 4028 | 9894 |
| 409 | DIE | 4085 | `(10206,12649,20003)` | 2159 | 4030 | 9898 |
| 410 | DIE | 4091 | `(10207,12653,20004)` | 2153 | 4032 | 9901 |
| 411 | DIE | 4102 | `(10208,12661,20005)` | 2176 | 4034 | 9913 |
| 412 | DIE | 4116 | `(10209,12670,20006)` | 2179 | 4036 | 9920 |
| 413 | CLOSED | 10116 | `(10210,16259,20007)` | 8191 | 4038 | -- |
| 414 | CLOSED | 10126 | `(10211,16266,20008)` | 8191 | 4040 | -- |
| 415 | CLOSED | 10115 | `(10212,16261,20009)` | 8191 | 4042 | -- |
| 416 | CLOSED | 10126 | `(10213,16268,20010)` | 8191 | 4044 | -- |

Thus the load-bearing adjacent comparison is

```text
D0=412 : DIE at row 4116, overflow bits 9920,
D0=413 : CLOSED at capture row 10116.                  (3)
```

The nonmonotone variation of the event rows within each block is retained in
the table; only the finite status assertion in (2) is proved.

## 3. Engine semantics and independent controls

The archival companion source-pins the exact engine and its independent floor
referee by the respective hashes

```text
8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad,
c679d5c1546160acb9ea5a71c5365178737ef2af5ba36a01eafdb1759c58aa75. (4)
```

It checks the full sixteen-row ledger against digest

```text
56640dafae5f05d891c2d6243c875be618d989006cdc5655df58348638955801. (5)
```

As live positive and hostile controls, the same pinned engine reconstructs

```text
(R,D0)=(512,4)   : DIE at 121,
(R,D0)=(512,5)   : CLOSED at 312,
(R,D0)=(1024,14) : DIE at 250,
(R,D0)=(1024,15) : CLOSED at 639.                      (6)
```

Separately, a FLINT/fmpz implementation reproduced the stable fields for
offsets `408,410,411,412,413`.  Its durable transcript segment has SHA-256

```text
aaec1a93c67fdad94a6a56bd9fa96edabfd1693945b84aec022db53a5d4b7ab5. (7)
```

This gives an independent adjacent audit of (3), not merely a replay through
the same Python implementation.

## 4. Reproduction and strict boundary

The ordinary verification, including all live controls, is

```bash
python3 04-computation/amm12592_R16384_offset_transition_thm3644.py
python3 -O 04-computation/amm12592_R16384_offset_transition_thm3644.py
```

Both streams must match the stored transcript after LF normalization.
Optional `--replay-adjacent` recomputes offsets `412,413` from the pinned
engine (about twenty-two minutes on the audited machine), while
`--replay-full` recomputes all sixteen expensive rows.

The theorem certifies only the fixed Rule-A status pattern in (2).  It proves
no monotonicity outside the bracket, no feasible alternative policy, no
uniform extractor, no asymptotic threshold, and no numerical bound for
AMM12592.  **QED.**
