---
id: THM-3029
title: "AMM 12592: the archimedean floor rate log_5(phi^2) is ATTAINED by the construction through R = 32, and profile monotonicity makes beam failures non-evidence"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. Improves the
  standing construction bound from C <= 8/5 = 1.6 to C = log_5(5 phi^2) =
  1.5979874356654402 on the epochs where it is verified, i.e. the construction
  MEETS the proved archimedean floor (THM-3009/THM-3017/THM-3024) exactly.
  (M) PROFILE MONOTONICITY. If an epoch closes with degree profile (d_i), it
  closes with EVERY pointwise-larger profile (d'_i >= d_i): lift each block by
  convolving its deltas with [binom(d'-d,k)]_k, which by x + (1-x) = 1 is the
  admissible block representing the CONSTANT 1, so the epoch identity is
  unchanged and admissibility is preserved (THM-3026 lemma L). CONSEQUENCE: a
  beam search that FAILS at a larger profile while SUCCEEDING at a smaller one
  has produced a pure SEARCH ARTEFACT, and such failures are not evidence.
  (D0) THE SLACK TRAP. With d_i = floor(gamma(R+i)) + D0 the EFFECTIVE rate at
  epoch R is max_i d_i/(R+i), which exceeds gamma by up to D0/R. So "gamma
  closes at R with slack D0" does NOT witness rate gamma at that epoch -- at
  R = 32, gamma = 0.48 with D0 = 4 closes, but its effective rate is ~0.56.
  Quantitatively a candidate gamma needs D0(R) >~ (gamma_c(R) - gamma) R, which
  is BOUNDED iff gamma >= lim gamma_c(R) = gamma*. gamma = 1/2 needs D0 ~ 95 by
  R = 1024 and is therefore not asymptotically viable -- exactly as the floor says.
  (C) THE FINITE-R CAPACITY FLOOR gamma_c(R) -- the least rate satisfying (ARCH)
  at epoch R -- rises to gamma* FROM BELOW: 0.375, 0.4412, 0.5000, 0.5435,
  0.5677, 0.5808, 0.5886, 0.5929 for R = 8..1024. So a construction beating 3/5
  at SMALL R is expected and says nothing asymptotic; only large R discriminates.
  (A) THE RESULT. The gamma* FLOOR PROFILE d_i = floor(gamma*(R+i)), D0 = 0,
  CLOSES at R = 8, 16, 32 -- directly by beam at R = 8, 16, and at R = 32 by
  lifting the (gamma = 1/2, D0 = 3) solution, whose blocks are all admissible at
  the target and whose epoch identity is exact. At these R the profiles for
  0.599, 0.598, 0.59799 and gamma* COINCIDE, so what closes is literally the
  floor profile. Effective rates 0.583333, 0.592593, 0.596774 against the 3/5
  profile's 0.600000.
  Hence C = 1 + gamma* = log_5(5 phi^2) is ATTAINED for n <= 63, matching the
  proved lower bound exactly on that range. UPDATE 2026-08-03: R = 64 is now
  CLOSED by two independent, distinct, hostile-re-verified witnesses
  (THM-3302), extending attainment to n <= 127 and superseding C <= 8/5
  there; R = 128 and all-R remain OPEN. Also note MISTAKE-361: the "proved
  lower bound" matched here is the BALANCED-BLOCK floor (THM-3009); the
  general-class promotion (THM-3024) is demoted. Toolchain note: the lost
  liftrate/gammac imports were reconstructed 2026-08-03 with byte-identical
  referee output (05-knowledge/results/amm12592-thm3029-toolchain-repair-boxeph.md).
source: death-star-2026-08-01-coinC2
depends_on:
  - THM-3026
  - THM-3002
related:
  - THM-3009
  - THM-3017
  - THM-3024
  - THM-3028
  - HYP-9061
script: 04-computation/amm12592_floor_rate_attained_thm3029.py
output: 05-knowledge/results/amm12592_floor_rate_attained_thm3029.out
---

# THM-3029 -- the construction meets the floor

## 1. Profile monotonicity (M)

A block is **admissible at `d`** when its coefficients in the basis
`B_{d,k}(x) = x^{d-k}(1-x)^k` satisfy `|delta_k| <= binom(d,k)` and
`delta_k = binom(d,k) mod 2`. THM-3026 (L) observed that, because
`x + (1-x) = 1`,

```text
sum_k binom(r,k) B_{r,k}(x) = 1      identically,
```

so `[binom(r,k)]_k` is an admissible block representing the **constant 1**.
Convolving a block's deltas with it re-expresses the *same polynomial* at
degree `d + r`, still admissible (Vandermonde). Therefore

```text
if an epoch closes with profile (d_i), it closes with every profile
(d'_i) satisfying d'_i >= d_i pointwise.                                  (M)
```

**This is the methodological point of the file.** A beam search that fails at a
*larger* profile while succeeding at a *smaller* one has proved nothing about
the larger profile -- it has exhibited a search artefact. Several such failures
were recorded in this lane and must not be read as evidence. (M) also converts
the search problem into: *find any sufficiently small profile that closes*, and
lift.

**A demonstrated instance, not merely a caution.** Running the direct beam at
`R = 32` over sub-`3/5` rates produces a uniform negative:

```text
rate                          R = 8        R = 16       R = 32          R = 64
3/5 (control)                 CLOSES       CLOSES       CLOSES          CLOSES
0.599                         CLOSES       CLOSES       "did not close" "did not close"
0.598                         CLOSES       CLOSES       "did not close" "did not close"
100/167 (simplest in gap)     CLOSES       CLOSES       "did not close" "did not close"
0.59799 (just above gamma*)   CLOSES       CLOSES       "did not close" "did not close"
```

Every `R = 32` entry above is **false**. Section 4 proves those very profiles
close, by lifting the `(gamma = 1/2, D0 = 3)` solution. So the beam's negatives
are not weak evidence, they are *wrong*, and no negative from this solver may be
cited as evidence about a profile without first checking whether some smaller
profile closes.

## 2. The slack trap (D0)

With `d_i = floor(gamma(R+i)) + D0`, the **effective rate** of the epoch is

```text
gamma_eff(R) = max_i  d_i / (R+i)   <=   gamma + D0/R.
```

At finite `R` a constant `D0` is *not* negligible. Concretely `gamma = 0.48`
with `D0 = 4` closes at `R = 32`, but its effective rate is about `0.56`, not
`0.48`; the apparent violation of the capacity floor `gamma_c(32) = 0.5`
disappears once the slack is accounted for.

Quantitatively, a candidate `gamma` needs

```text
D0(R)  >~  (gamma_c(R) - gamma) * R,
```

which stays **bounded** exactly when `gamma >= lim gamma_c(R) = gamma*`:

```text
   R    gamma_c(R)   D0 needed at 1/2   at 0.55   at 0.58   at 3/5   at gamma*
  32     0.500000          0.0             -         -         -         -
  64     0.543478          2.8             -         -         -         -
 128     0.567669          8.7            2.3        -         -         -
 512     0.588619         45.4           19.8       4.4        -         -
1024     0.592891         95.1           43.9      13.2        -         -
```

`gamma = 1/2` needs unbounded slack and is therefore **not** asymptotically
viable -- an independent confirmation of the proved floor from the construction
side. `gamma = 3/5` and `gamma*` need none at any `R`.

## 3. The finite-R capacity floor (C)

`gamma_c(R)` is the least rate satisfying (ARCH) at epoch `R`. It rises to
`gamma*` **from below**:

```text
R          8        16        32        64       128       256       512      1024
gamma_c  0.3750   0.4412   0.5000   0.5435   0.5677   0.5808   0.5886   0.5929   -> 0.5979874
```

So the necessary condition is *weak* at small `R`, and a construction beating
`3/5` at `R = 8` or `16` is expected and carries no asymptotic information.
Only large `R`, where `gamma_c(R)` approaches `gamma*`, discriminates.

## 4. The result (A)

Let `gamma* = log_5(phi^2) = 0.5979874356654402` and take the **floor profile**
`d_i = floor(gamma*(R+i))`, `D0 = 0`. At `R = 8, 16, 32` the profiles for
`0.599`, `0.598`, `0.59799` and `gamma*` are **identical**, so this is exactly
what a sub-`3/5` construction must produce.

```text
R = 8   direct beam solve, verify() True,             effective rate 0.583333
R = 16  direct beam solve, verify() True,             effective rate 0.592593
R = 32  by (M): lift the (gamma=1/2, D0=3) solution   effective rate 0.596774
        [source verified; target pointwise larger; all lifted blocks admissible;
         epoch identity sum_i x^i Delta_i = q^{R-1} exact]
                                        (the 3/5 profile has rate 0.600000)
```

Hence

```text
C = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654402   is ATTAINED for n <= 63,
```

matching the proved archimedean lower bound (THM-3009, THM-3017, and THM-3024
for the general class) **exactly** on that range. The previous best
construction constant was `C <= 8/5 = 1.6`.

## 4b. The trend: the construction-capacity gap is shrinking

Section 4 shows the `gamma*` profile closes; a separate question is how close the
construction gets to the *finite-R* capacity floor. Minimising the **effective
rate** `max_i d_i/(R+i)` over all `(gamma, D0)` pairs that close and verify:

```text
   R    gamma_c(R)    best effective rate achieved   witness        gap
   8     0.375000            0.444444              107/240, D0=0   0.0694
  16     0.441176            0.500000              120/240, D0=0   0.0588
  32     0.500000            0.545455              110/240, D0=3   0.0455
```

The gap is **monotonically shrinking**: `0.069, 0.059, 0.045`. Since both
`gamma_c(R)` and any achievable rate must converge to `gamma*` if the floor is
tight, a shrinking gap is precisely the behaviour required for
`C* = log_5(5 phi^2)`; a gap bounded away from zero would instead point to
`C* > gamma* + 1`.

This is **EVIDENCE, not proof** -- three values of `R`, and the "best achieved"
figures are upper bounds from a heuristic search (by (M) they can only improve,
never worsen, as the search gets better). It is nevertheless the cleanest
asymptotic signal available from the construction side, and it points the same
way as section 4.

## 5. Scope

(M) and the `D0` analysis are proofs. (C) is a numerical evaluation of (ARCH)
at the stated `R`. (A) is **VERIFIED-EXACT but finite**: it establishes the
floor rate on the epochs `R = 8, 16, 32`, i.e. for `n <= 63`, and **proves
nothing for larger `n`**. In particular this does **not** show
`C* = log_5(5 phi^2)`; that needs the floor profile to close at every `R`,
which is open. The standing global statement is still `C <= 8/5` for
`n <= 127` (THM-3002 5b), and `R = 128` at `3/5` remains open (THM-3028).
Extending (A) to `R = 64` is in progress and not claimed here.

Referee: `amm12592_profile_monotonicity_thm3029.py` -- the lifting construction
and its admissibility/identity checks; the `gamma_c` table is
`amm12592_finite_R_capacity_threshold.py`.
