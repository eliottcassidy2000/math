---
id: THM-3633
title: "AMM R=16384 fixed-failed-trace adjoint phase shock"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + FINITE-EXACT + VERIFIED-EXACT /
  PENDING INDEPENDENT HOSTILE AUDIT.  At the fixed archived failed trace
  (R,D0)=(16384,400), the exact Rule-A profile dies at row 4055 and its
  positive truncated-adjoint wall changes sign between cuts 1516 and 1517.
  Hence every admissible continuation avoiding that fatal inequality must
  depart by row 1516.  The normalized horizon lies below the golden candidate
  and is exactly 15/8192 farther away than the R=8192 control.  A local phase
  invoice isolates the dyadic q-defect -30.  D0=400 is not known terminal;
  no AMM bound or asymptotic law is claimed.
source: kps-s189 + agent Newton / THM-3626 dyadic continuation, 2026-08-21
audit: >
  PENDING -- author/agent exact reconstruction passed the complete profile,
  fatal, checkpointed adjoint, contiguous-wall, boundary, active-ledger,
  phase-invoice, dyadic-defect, parent-pin, no-assert, semantic, and normal/-O
  gates.  Independent hostile reconstruction remains required.
depends_on:
  - THM-3626-amm-r8192-adjoint-horizon-and-phase-rebound
related:
  - THM-3626-amm-r8192-adjoint-horizon-and-phase-rebound
  - THM-3616-amm-R4096-adjoint-horizon-and-golden-finite-scale-obstruction
script: 04-computation/amm_binary_sturmian_R16384_D400_adjoint_phase_shock_thm3633.py
output: 05-knowledge/results/amm_binary_sturmian_R16384_D400_adjoint_phase_shock_thm3633.out
script_sha256: a4888d29f0eb1f59debea6f1545b16928bd9278af77bfa9ae5f5480ee08bdab4
output_sha256: 6754b09009d72ecee1fb5e5e3230adf2a660578cfa34a047dd154ddb28088104
semantic_sha256: 08450a048b8e1379d3945e9464021b53203bb10f82638daf6bbfca076bd71a00
hash_basis: LF-normalized bytes for files; canonical JSON for semantic ledger
---

# THM-3633 -- AMM `R=16384` fixed-failed-trace adjoint phase shock

**RESERVED / PROVISIONAL PROOF CANDIDATE + FINITE-EXACT + VERIFIED-EXACT /
PENDING INDEPENDENT HOSTILE AUDIT.**  This certifies one fixed failed trace at
the next dyadic scale after THM-3626.  It does not claim that the archived
offset is the last failure.

## 1. Exact fatal row and adjoint wall

For `0<=i<16384`, use the inherited Fibonacci--Lucas comparison to set

```text
d_i(D0)=floor(log_5(phi^(2(16384+i))))+D0.              (1)
```

At the fixed archived offset `D0=400`, the exact Rule-A trace has

```text
first death j=4055,
(d_0,d_4055,d_16383)=(10197,12622,19994).               (2)
```

The fatal integer is negative and has `9874` bits.  Using the same positive
truncated adjoint as THM-3626, let `B_s` be the necessary right-hand side for
an admissible continuation agreeing with the failed trace through cut `s-1`.
The complete checkpointed sweep gives the contiguous wall

```text
B_1516>0>B_1517,
B_s<0 exactly for s=1517,1518,...,4054.                 (3)
```

Therefore every admissible continuation that avoids this fatal inequality
must depart at or before

```text
q=1516.                                                 (4)
```

The bit lengths of the two boundary integers are respectively

```text
bits(B_1516)=9870,                 bits(B_1517)=9874,   (5)
```

and their exact ordered digest is

```text
6ecf448b01a6225777b389d9293084e5c706c7cdf80fdfdf87f859b9e05d787f. (6)
```

At the first negative cut the positive multiplier ledger contains exactly
`3,221,991` active cells, with digest

```text
de05e4c1cd3748884a8866af9fb274ae01b3d749c80bc63dcfe63806f73a87e6. (7)
```

The enormous profile, fatal, and full-wall ledgers are independently pinned
by the companion digests rather than printed.

## 2. Golden-coordinate deterioration

Let

```text
theta=(3-sqrt(5))/8.                                    (8)
```

The new normalized departure coordinate obeys

```text
q/R=1516/16384=379/4096,
q/R-theta=(-1157+512sqrt(5))/4096<0.                    (9)
```

The sign is exact because

```text
1157^2-5*512^2=27929>0.                                (10)
```

Putting the absolute discrepancies at the two latest scales over denominator
`8192` gives

```text
E_8192 =(2299-1024sqrt(5))/8192,
E_16384=(2314-1024sqrt(5))/8192.                        (11)
```

Consequently

```text
E_16384-E_8192=15/8192>0.                              (12)
```

Thus the improvement seen from `R=4096` to `R=8192` does not persist at this
fixed `R=16384` failed trace.  This is an exact finite-scale deterioration,
not evidence for or against an eventual limit.

## 3. Exact local phase invoice

For either audited scale put

```text
h=5R/8-d_0,
b=2d_0-R+2,
m=j-b,
ell=j-s,                       s=first negative cut,
beta=(sqrt(5)-1)/8,
e=ell-beta R.                                           (13)
```

Since `theta+beta=1/4` and `q=j-ell-1`, direct algebra gives the exact invoice

```text
q-theta R=m+1-2h-e.                                    (14)
```

The two certified records are

```text
R=8192:  (h,b,m,ell)=(31,1988,57,1271),
          e=2295-1024sqrt(5),

R=16384: (h,b,m,ell)=(43,4012,43,2538),
          e=4586-2048sqrt(5).                           (15)
```

Their raw dyadic defects are

```text
q_16384-2q_8192       =-30,
j_16384-2j_8192       =-35,
ell_16384-2ell_8192   = -4.                             (16)
```

More finely,

```text
m_16384-2m_8192=-71,
h_16384-2h_8192=-19,
e_16384-2e_8192= -4,                                   (17)
```

so subtracting twice `(14)` gives the exact integer decomposition

```text
-30=(-71-1)-2(-19)-(-4).                               (18)
```

This identifies the phase shock: the margin term collapses by `71`, partially
cancelled by the headroom and depth sidecars.  It does not prove that these
three variables form a closed recurrence at later scales.

## 4. Verification and strict boundary

The deterministic companion hash-pins THM-3626 and THM-3616, independently
reconstructs the complete `R=8192` control, then rebuilds the `R=16384` trace.
Its block-checkpointed adjoint uses bounded storage while retaining every exact
charge pairing.  Normal and optimized runs execute the same `27,594`
optimization-safe gates; the source has raw LF bytes and no Python `assert`
node.

Reproduce with

```bash
python3 04-computation/amm_binary_sturmian_R16384_D400_adjoint_phase_shock_thm3633.py
python3 -O 04-computation/amm_binary_sturmian_R16384_D400_adjoint_phase_shock_thm3633.py
```

Both streams must agree with the stored transcript after LF normalization.
The result concerns only the fixed failed prefix `(R,D0)=(16384,400)`.
This package does not identify the terminal failing offset.  The theorem
proves no adjacent survivor, alternative feasibility, asymptotic phase law,
uniform extractor, value of the AMM constant, or improvement to its bound.

Independent hostile audit remains the promotion gate.
