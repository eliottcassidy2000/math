---
id: THM-3015
title: "The fifth, sixth and seventh resultant log-jets, and the corner-slice laws"
status: VERIFIED-EXACT (interpolated under a probed degree ansatz; L_5 identity PROVED) / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428 (two-agent build + adversarial rebuild; shape, predictions and new law re-checked here)
depends_on:
  - THM-3013-symbolic-fourth-resultant-log-jet-P4
  - THM-3014-jet-D-adic-structure-discriminant-value-and-factorial-readout
related:
  - THM-3011-fourth-resultant-jet-and-the-third-edge-invoice
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
table_p5: 05-knowledge/results/gmc_first_gap_resultant_jet_P5_table_thm3015.json
table_p6: 05-knowledge/results/gmc_first_gap_resultant_jet_P6_table_thm3015.json
table_p5_sha256: 8cf6ed5cfca3b92a9229f79aae26f87ab1e65db1cf288b4299fe37b4a47b1de9
table_p6_sha256: 7004c579194e13f10aa03ceb26707adbaeae01e64b5be85d76792987f20c150e
table_p7: 05-knowledge/results/gmc_first_gap_resultant_jet_P7_table_thm3015.json
table_p7_sha256: bc53e1a23a9694c277de3aa1e9f6c4401be585452b7a20e35abc0a7a050fb287
script: 04-computation/gmc_fifth_and_sixth_resultant_jets_verify_thm3015.py
output: 05-knowledge/results/gmc_fifth_and_sixth_resultant_jets_verify_thm3015.out
script_sha256: e394453cbbd23295c3860b38e8fa6325680c56a63bbb8ba6829ef73edc528218
output_sha256: 5407a624633dd4262edd739e366f313bd14559da12a42c434cde28e9e7ae9533
hash_basis: LF-normalized bytes
---

# THM-3015 -- the fifth and sixth resultant log-jets

**VERIFIED-EXACT.  The `L_5` trace identity is PROVED; the tables are exact
interpolations inside a degree ansatz that was probed with oversampling.**

## 1. The fifth line (PROVED)

With `X_k = m_0^(-1)m_k`, `Y = sum_(k>=1) X_k t^k`, and
`log det(I+Y) = sum_(p>=1) (-1)^(p+1) tr(Y^p)/p`, extracting `[t^5]` over
compositions of `5` with cyclic-invariance multiplicities gives

    L_5 = tr X_5 - tr(X_1X_4) - tr(X_2X_3) + tr(X_1^2X_3) + tr(X_1X_2^2)
          - tr(X_1^3X_2) + tr(X_1^5)/5.                                (1)

The same extraction at `[t^4]` returns THM-3011 (1) verbatim, which is the
internal check.  The general line is
`L_n = sum_p ((-1)^(p+1)/p) sum_(k_1+..+k_p=n) tr(X_(k_1)...X_(k_p))`.

## 2. The tables

| | terms | degrees `(M,U,V)` | support | absences |
|---|---|---|---|---|
| `P_5` | **1313** | `(10,20,10)` | `b+2c <= 20` | 18, all at the three corners |
| `P_6` | **2176** | `(12,24,12)` | `b+2c <= 24` | 21, all at the three corners |

`P_5`'s absences sit at `M`-powers `{1,3,7,8,9,10}` on each of `(20,0)`, `(0,10)`,
`(0,0)`.  Shape ledger `j=1..6`: terms `27/122/333/717/1313/2176`, absences
`0/3/10/12/18/21`, and

    |P_j| = (2j+1)^3 - 3(2j - 2 - floor(j/2)) - [j=3]                  (2)

reproduces all six.  The `[j=3]` correction is THM-3013's `(5,0,3)`, still the
**only** non-corner absence anywhere; a `j>=7` sporadic is not excluded.

`c_j` remains **convention-dependent** (THM-3013 section 2).  The contents
`1/content(D^(2j)L_j)` for `j=1..6` are `1, 2^4*3, 2^7*3^2, 2^12*3^4*5,
2^15*3^5*5, 2^19*3^8*5*7`: the `2`-exponents `0,4,7,12,15,19` and `3`-exponents
`0,1,2,4,5,8` show **no** pattern.  The one regularity is that the prime set is
exactly `{p <= j+1}` -- `5` entering at `j=4`, `7` at `j=6` -- consistent with the
`1/(j(j+1))` equidistribution denominators.


## 2a. `P_7`, and the off-corner test

**`P_7` has exactly `3348` terms**, degrees `(M,U,V) = (14,28,14)`, support
`b+2c <= 28`, `27` absences, and **zero off-corner absences**.  Equation (2)
predicts `3375 - 3*(14-2-3) = 3348` -- confirmed.  So the `j=3` sporadic
`(5,0,3)`, the only non-corner absence anywhere, does **not** recur at `j=7`;
(2) now holds on `j=1..7`.

`c_7`: `1/content = 2889476997120 = 2^22 * 3^9 * 5 * 7`.  The `2`- and
`3`-exponents `0,4,7,12,15,19,22` and `0,1,2,4,5,8,9` still show no pattern, but
the prime set is `{2,3,5,7} = {p <= j+1}`, continuing the one regularity there is.

**All 21 pre-registered `j=7` coefficients confirmed**, seven on each of the three
corner slices:

    A_7 = (23/28)M^8 - (12/7)M^7 + (47/24)M^6 - (707/768)M^4 + (11267/43008)M^2 + 6657676188189/14336
    E_7 = (23/28)M^8 - (11/7)M^7 + (37/18)M^6 - (301/324)M^4 +  (2677/10206)M^2 +   72615186151/54432
    C_7 = (23/28)M^8 -  (3/7)M^7 +  (23/6)M^6 -   (23/12)M^4 +       (23/42)M^2 + 102839873747325/224

so the six laws of section 4 now hold on `j=1..7`, including the `[M^(j-4)]=0`
law at a third point and the `[M^(j-5)]` falling-factorial law at a second.

**Provenance note.**  The `j=7` build agent stalled before returning, but it had
already written the tables; they were recovered from disk and verified here
independently against the pre-registered predictions and equation (2).  The
digest is recorded above.  This is a recovered artifact checked after the fact,
not a completed audited run -- treat it as one notch weaker than `P_5`/`P_6`
until someone reproduces the build.

## 3. Seven pre-registered predictions, all CONFIRMED

THM-3011 section 2a and THM-3014 section 2 recorded these **before** `P_5`
existed.  Re-checked here against the table in all three normalisations:

    A_5: [M^6] = 23/15,  [M^5] = -12/5,  [M^4] = +47/24     all CONFIRMED
    E_5: [M^6] = 23/15,  [M^5] = -11/5,  [M^4] = +37/18,  [M^3] = 0   all CONFIRMED

Full slices:

    A_5 = (23/15)M^6 - (12/5)M^5 + (47/24)M^4 - (707/1920)M^2 + 1184426373/640
    E_5 = (23/15)M^6 - (11/5)M^5 + (37/18)M^4 - (301/810)M^2 + 16866997/648
    C_5 = (23/15)M^6 -  (3/5)M^5 +  (23/6)M^4 -   (23/30)M^2 + 8248311801/40

with `C_j := Q_j(M,0,0)` the third corner.  Hence
`p_5(R_M) -> (23/3)M^6 - 12M^5 + (235/24)M^4 - (707/384)M^2 + 1184426373/128`,
confirming equidistribution `46/(j+1)` at `j=5`.

## 4. Two new laws, and what they explain

Extending THM-3011/THM-3014's four laws to six, verified on all three corner
slices at `j=1..6` (and, by the auditor, at `j=7`):

    [M^(j+1)] = (-1)^(j+1) * 46/(j(j+1))                    all three slices
    [M^j]     = (-1)^j * kappa/j,     kappa = 12 (A), 11 (E), 3 (C)
    [M^(j-1)] = (-1)^(j+1) * lambda,  lambda = 47/24 (A), 37/18 (E), 23/6 (C)   j>=2
    [M^(j-2)] = 0                                                              j>=3
    [M^(j-3)] = (-1)^j * mu * (j-1)(j-2),                                      j>=4
                mu = 707/23040 (A), 301/9720 (E), 23/360 (C)          ** NEW **
    [M^(j-4)] = 0                                                     j>=5  ** NEW **

> **The `[M^(j-3)]` law retro-explains two constants the canon had flagged as
> unexplained.**  THM-3011's `707/3840` and THM-3014's `301/1620` are exactly the
> `j=4` members: `(707/23040)*3*2 = 707/3840` and `(301/9720)*3*2 = 301/1620`.
> Both re-derived here.

`A_j` and `E_j` are now pinned to two unknowns each: `[M^(j-5)]` (one data point)
and the constant term.

## 5. Structural findings

- **`E_j` needs no D-adic re-expansion.**  `D` is linear in `V` with slope `-3`
  and `deg_V Q_j = 2j`, so `E_j = [V^(2j)]Q_j / 9^j` exactly (verified `j=1..6`).
  This simplifies THM-3014 section 2: the discriminant regular value is just the
  `V`-corner rescaled.
- **THM-3014's "arithmetic step 18" carries no information.**  It is the `j=4`
  instance of step `2(2j+1)`, a pure degree artifact.  The real D-adic laws,
  measured at `j=1..6`, are `|q_k| = (4j-2k+1)(2j+1)` for `1<=k<=2j-1` (fully
  dense), `|q_0| = (4j+1)(2j+1)-j`, `|q_(2j)| = floor(j/2)+3`.
- **No `P_j -> P_(j+1)` recursion was found, and the absences are a basis
  artifact.**  Newton's identities give a triangular re-encoding into the
  `D`-scaled elementary symmetric functions, and that re-encoding fills the box
  `(2j+1)^3` with **zero** absences at `j=1..6`.  So there is no sparsity to
  exploit: `E_(j+1)` is `(2j+3)^3` genuinely new rationals that can only come
  from the `(j+1)`-st Macaulay slice.  **This is FINITE-EXACT on `j<=6`, not a
  proof that no transfer exists** -- the build agent's word "provably" was
  removed by the auditor and is not used here.
- **But `P_6` is already cheap**, which is the operative point.  Replacing the
  hand-written trace formulas by `d/dt log det A(t) = tr(A(t)^(-1)A'(t))` makes
  the jet order a runtime parameter.  Measured on 8 cores: `P_5` grid 21 s + fit
  4 s; `P_6` grid 100 s + fit 6 s; extrapolating, `P_7 ~ 5-8 min`, `P_8 ~ 20-30
  min`.  There is no wall.  THM-3014 section 3's "a table is a cache, not the
  object" is now implemented rather than asserted.

## 6. Scope

- The degree ansatz was **probed** (last nonvanishing divided difference exactly
  `20/10/10` for `P_5` with 9/13/9 spare vanishing orders), not proved.  Same
  status as THM-3013.
- All `j`-pattern statements are measured on `j<=6` (the auditor adds `j=7` for
  the corner laws only).  Equation (2) and the two new laws are extrapolations
  beyond that.
- Everything inherits THM-2997's encoded-wall continuation for `M>=27`.
- Nothing here proves no-return, ULC, or GMC(2).
