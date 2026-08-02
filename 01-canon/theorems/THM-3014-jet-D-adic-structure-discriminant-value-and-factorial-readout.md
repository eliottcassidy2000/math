---
id: THM-3014
title: "D-adic structure of the resultant log-jets, the discriminant regular value, and the factorial readout"
status: PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428
depends_on:
  - THM-3013-symbolic-fourth-resultant-log-jet-P4
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
related:
  - THM-3011-fourth-resultant-jet-and-the-third-edge-invoice
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-2812-consecutive-three-slot-factorial-moment-six-detection
  - THM-2836-sfc3-arbitrary-support-shifted-window-census
  - THM-2173-sparse-projective-factorial-moment-floor
script: 04-computation/gmc_jet_D_adic_structure_and_factorial_readout_thm3014.py
output: 05-knowledge/results/gmc_jet_D_adic_structure_and_factorial_readout_thm3014.out
script_sha256: 1a94956f603f4c5b7d9b7d7c17903d51b64222893670acb43cdc1185c2b35133
output_sha256: ea8e58163523fe39a35db0bd08e9fb64613fa850f4c7d9c17b171ab878dd1bab
hash_basis: LF-normalized bytes
---

# THM-3014 -- D-adic jet structure, the discriminant value, and the factorial readout

**PROVED + VERIFIED-EXACT.**

## 1. The D-adic re-expansion is exact but NOT a compression

`D=U^2+3U-3V-1` is **linear in `V`**, so `(U,V) -> (U,D)` is an invertible
polynomial change of variables, `V=(U^2+3U-1-D)/3`.  Re-expanding THM-3013's
`P_4 = sum_k q_k(M,U) D^k`:

| `k` | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0 |
|---|---|---|---|---|---|---|---|---|---|
| terms in `q_k` | 5 | 27 | 45 | 63 | 81 | 99 | 117 | 135 | 149 |
| `deg_U q_k` | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 |

- **`deg_U q_k = 2(2j-k)` exactly** -- forced, since `deg_U D = 2`.
- Interior term counts are an arithmetic progression of step `18`.
- **Total `721` against `717`**: D-coordinates are *not* sparser.
- `P_4` is **not** divisible by `D` (`q_0 != 0`), as it must not be: the `D^{2j}`
  factor exists precisely to clear the pole.

So the answer to "factor the coefficients D-adically" is: the expansion is clean
and graded, but it buys no compression.  Its value is section 2.

## 2. The discriminant regular value `E_j(M)`, and its laws

Since `deg_U P_j = 4j = 2j * deg_U D`, the top coefficient `[D^{2j}]P_j` is
**forced `U`-free**.  Hence `L_j` has a *regular value at the discriminant locus*
`D=0` that is a pure polynomial in `M`:

    E_j(M) := [D^(2j)] P_j / c_j.

Computed exactly (`U`-freeness confirmed at every `j`):

    E_1 = 23M^2 - 11M + 10
    E_2 = -(23/3)M^3 + (11/2)M^2 - (37/18)M - 553/6
    E_3 = (23/6)M^4 - (11/3)M^3 + (37/18)M^2 + 31709/54
    E_4 = -(23/10)M^5 + (11/4)M^4 - (37/18)M^3 + (301/1620)M - 1651481/432

**Laws** (all verified on `j=1..4`):

    [M^(j+1)] E_j = (-1)^(j+1) * 46/(j(j+1))     (equidistribution)
    [M^j]     E_j = (-1)^j * 11/j
    [M^(j-1)] E_j = (-1)^(j+1) * 37/18            (j >= 2)
    [M^(j-2)] E_j = 0                             (j >= 3)

`E_j` and THM-3011's top-`U` slice `A_j` share **only** the leading
equidistribution coefficient; the second differs (`11/j` vs `12/j`) and the third
differs in value *and* sign (`37/18` vs `47/24`).  They are two different limits
of the same jet, agreeing exactly where the limit-independent law lives.

## 3. `P_5`: forced shape, predicted slices, and a better representation

Shape is forced by the pattern established at `j=1..4`: max terms
`(2j+1)^3 = 1331`, degrees `(M,U,V)=(10,20,10)`, support `b+2c <= 20`.

Predicted extreme slices, from section 2 and THM-3011 section 2a:

    A_5:  [M^6]=23/15,  [M^5]=-12/5,  [M^4]=+47/24
    E_5:  [M^6]=23/15,  [M^5]=-11/5,  [M^4]=+37/18,  [M^3]=0

**A better representation than a fifth table.**  All `L_j` are coefficients of
ONE object, `log(R_M(n)/(LC n^d)) = sum_j L_j n^(-j)`.  Tabulating `P_5` as `1331`
integers repeats work: the right artifact is the generating series together with
the `c_j` normalisation, from which every `P_j` follows.  Concretely, the
THM-3013 pipeline already computes `L_1..L_4` from a single `36x36` slice tower;
extending the log-determinant recursion by one line yields `L_5` at no structural
cost, and the interpolation cost grows only as the ansatz `(2j+1)^3`.  A table is
a cache, not the object.

## 4. The factorial readout, and the SFC lane

THM-3000's first-occurrence law gives jet `J_j` the coefficient
`-(j-1)! * C(k-2,j-3)` at edge `k` -- a **factorial** times a **binomial**.
Exactly:

    Lambda_k(J) := sum_(j>=3) [J_j] W_(j-1)^(k) * J_j
                 = -sum_(i>=0) (i+2)! C(k-2,i) J_(i+3)
                 = -L( z^2 * Phi_k(z) ),   Phi_k(z) = sum_i C(k-2,i) J_(i+3) z^i,

where `L(z^n)=n!` is the **Edo--van den Essen factorial functional** that defines
the Strong Factorial Conjecture lane (THM-2173, THM-2812, THM-2836).  Verified
for `k=2..8`.

> The Newton edge reads the log jets **through `L`**, twisted by a binomial
> transform whose index is the edge.  The edge index sets the binomial weight;
> the jet index sets the factorial.

**Loss ledger -- this is a shared functional, NOT a reduction.**

| | SFC lane | circuit lane |
|---|---|---|
| functional | `L(z^n)=n!` | same |
| argument | `f^m` for a `t`-term `f` | `z^2 Phi_k(z)`, `Phi_k` a binomial transform of the jets |
| dependence on data | **nonlinear** (an `m`-th power) | **linear** in the jets |
| depth parameter | power `m`, up to `t` | edge index `k` |
| atom parameter | term count `t` | cluster count `K` (THM-3004) |

The decisive difference is the middle row: `SlotSFC_1(t)` asks whether
`L(f), ..., L(f^t)` can all vanish for a `t`-term `f`, which is a question about
`L` on **powers**; `Lambda_k` is `L` on an object **linear** in the jets.  Nothing
here proves, reduces to, or is implied by SlotSFC_1(2) (proved, Edo--van den Essen 2013)
or SlotSFC_1(3) (THM-2836's certified census box, THM-2812's consecutive case).
MISTAKE-237 retracted one NC2 -> JC(2) "bridge"; this file deliberately claims a
common evaluation functional and an index dictionary, nothing more.

**What the dictionary does suggest, as questions rather than claims.**
The parallel `t <-> K` and `m <-> k` is suggestive because both lanes already have
a depth law: SFC needs depth `t` for a `t`-atom object, and THM-3004 gives
`2K-3` sign changes for a `K`-cluster object.  Two cheap tests follow:
**Test (i) is answered, NEGATIVELY, in this file.**  Does `L(f^m)` obey a
sign-change bound of the THM-3004 shape `2t-3`?  **No.**  Randomised exact search
over `t=2,3,4` with supports in `[0,8)` and integer coefficients finds sequences
`L(f), ..., L(f^8)` that alternate at **every** step -- `7` of `7` possible sign
changes -- at every `t`, including `t=2` where `2t-3=1`.  Example: `t=2`,
support `(1,5)`, coefficients `(2,-9)`.  The reason is structural and kills the
analogy cleanly: the sign of `L(f^m)` is driven by which term dominates the
`m`-th power, so a single negative coefficient of large modulus alternates the
whole sequence irrespective of `t`.  Sign-change count is therefore **not** the
right SFC-side invariant, and the `t <-> K` half of the dictionary does not
transport THM-3004's law.

That leaves (ii): does the binomial transform `Phi_k` have an SFC-side meaning --
is `L(z^2 Phi_k)` ever an `L(f^m)` for an explicit `f`?  Given (i), a negative
answer to (ii) would pin the dictionary as cosmetic, and that is now the honest
default expectation rather than an open hope.

## 5. Boundaries

- Section 1's counts are for `j=4`; the `deg_U q_k = 2(2j-k)` law is degree-forced
  and therefore general, the term counts are not.
- Section 2's laws are verified at `j=1..4` only.  The `[M^(j-2)]=0` law rests on
  two data points (`j=3,4`) and is the weakest of the four.
- Section 3's `P_5` predictions are **extrapolations** of section 2 and THM-3011
  section 2a, not computations.  Nothing about `P_5` is verified here.
- Section 4 is an identity plus a dictionary.  It is not a bridge and must not be
  cited as one.
- Everything inherits THM-3013's status (exact interpolation under a measured
  ansatz) and THM-2997's continuation wall for `M>=27`.
