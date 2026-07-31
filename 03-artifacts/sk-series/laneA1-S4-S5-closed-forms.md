# Lane A1 — S(4), S(5) for S(k) = Σ_{n≥0} C(2n,n)C(4n,2n) / ((kn+1) 64^n)

Session: coinC2 / lane A1.
Canonical script: **`/tmp/math-wt-coinC2/04-computation/sk_series_pslq_laneA1.py`**
```
python3 sk_series_pslq_laneA1.py table 520     # rebuild the 515-digit S(k) table
python3 sk_series_pslq_laneA1.py verify        # check every representation below
python3 sk_series_pslq_laneA1.py schwarz       # the k<=3 / k>=4 regime diagnostic
python3 sk_series_pslq_laneA1.py pslq          # the exclusion runs + k=3 sanity control
```
Supporting scripts (same directory): `sk_pslq_engine.py`, `sk_probe.py`,
`sk_final_pslq.py`, `sk_kitchen_sink.py`, `sk_verify_forms.py`.
All numerics mpmath, working precision ≥ 350 dps for PSLQ unless stated.

---

## 0. Headline

**S(4) and S(5) were NOT reduced to a classical closed form, and I now believe none
exists of the type S(1), S(2), S(3) have.**  The reason is structural and sharp:

> The hypergeometric that controls S(k) has Schwarz angle sum **> 1 exactly for
> k = 1, 2, 3** and **= 1 for every k ≥ 4**.  k = 4 is precisely the wall.

What *was* obtained: a uniform new integral representation valid for all k that
reproduces every known closed form mechanically, several new exact reductions of
S(4) (including a **rational** two-dimensional integral), 510-digit values, and a
battery of PSLQ exclusions that pin down what S(4) and S(5) are *not*.

---

## 1. High-precision values  (VERIFIED-NUMERIC, ≈ 514 digits)

Computed as `S(k) = ∫_0^1 2F1(1/4,3/4;1;t^k) dt` with the AGM evaluation
`2F1(1/4,3/4;1;z) = sqrt(2/(1+w))·(2/π)·K(m)`, `w = sqrt(1−z)`, `m=(1−w)/(1+w)`
(mpmath parameter convention), tanh–sinh quadrature.
Cross-checked at dps 520 vs dps 560: agreement to **5·10⁻⁵¹⁵** for every k.
Reproduces S(1), S(2), S(3) closed forms to 500 digits (residuals 0, 2.4e-501, 2.2e-500).

```
S( 4) = 1.06935255444126805858296193953427806132593201450630803756288610958327684527563737450723312250703705075321091427672643087533912066206933253937630168946...
S( 5) = 1.05708716209477162686360610231100103612637197979788061199083857250115970978291491185127374646059607375264843457690897777701370523254730527086260687956...
S( 6) = 1.04852004921802894368141152435909714731278245795491273017979672433994043531184186774349725838658678299586863921052748465530034477216288324578724024534...
S( 8) = 1.03733013650043673193740244495426303855919241047157137205810623740479266004296547004778133125115625864838780918167542057364216956636180660043447850234...
S(12) = 1.02555584240313734228603437132357545378261992725825348560111919978824587075317200154367997247536375697492962098438232655285773037684242919806933184921...

π·S(4) = 3.3594701291301671589900224730224240656919222546691...
π·S(5) = 3.3209372626410174933367229658774746263264441740161...
```
Full 515-digit values: `/tmp/math-wt-coinC2/04-computation/Sk_520.json`, `Sk_560.json`.

---

## 2. New exact results

### 2.1  Simplification of the known S(3)  — **PROVED**

The stated form `π S(3) = −√3 log(5−2√6) − 2 arctan(√2/5)` hides an exact identity:

> **arctan(√2/5) = π − 3 arctan(√2)**,
> because `(1+i√2)³ = −5 + i√2`, and both sides lie in (π/2, π).

Hence
```
π S(3) = 2√3 · log(√3+√2)  −  2π  +  6 arctan(√2)
       = 2√3 · arccosh(√3)  −  2π  +  6 arccos(1/√3).
```
Confirmed by PSLQ at 400 dps, residual 1.7e-400, coefficients (1, +2, −6, −2√3).
This is the right normal form: the alphabet is `{π, arctan √2, √3·log(√3+√2)}`,
i.e. `arg` and `log` of units of Q(√−2) and Q(√6) — **weight 1**.
(For comparison: `π S(2) = 4 log(1+√2)`, `π S(1) = 8√2/3`.)

### 2.2  Uniform representation for all k  — **PROVED**, verified 50 digits for k = 1…12

Start from the session's (v),
`S(k) = (2/π) ∫_0^{π/2} 2F1(1/2, 2/k; 1+2/k; cos 2θ) dθ`.
Its parameters satisfy `c − a − b = 1 + 2/k − 1/2 − 2/k = 1/2` **for every k**, so the
quadratic transformation `2F1(a,b;a+b+½;z) = 2F1(2a,2b;a+b+½;(1−√(1−z))/2)` applies.
With `z = cos 2θ`, `√(1−z) = √2 sin θ`:

> **S(k) = (2/π) ∫_0^{π/2} 2F1( 1, 4/k ; 1 + 2/k ; (1 − √2 sin θ)/2 ) dθ**

Numerical check (dps 50): exact agreement for k = 1,2,3,4,5,6,8; 2.7e-51 for k = 12.

This single formula *mechanically* reproduces the known cases:
* k = 2: `2F1(1,2;2;y) = 1/(1−y)` ⟹ `S(2) = (4/π)∫_0^{π/2} dθ/(1+√2 sin θ) = (4/π)log(1+√2)` — exactly the session's route (vi), now derived rather than assumed.
* k = 1: `2F1(1,4;3;y)` is rational + log ⟹ elementary.
* k = 4: `2F1(1,1;3/2;y) = arcsin√y /(√y√(1−y))`; with `y = (1−cos ν)/2` this is `ν / sin ν`.

### 2.3  Thomae normal form  — **PROVED**, verified 50 digits

From `S(k) = 3F2(1/4, 3/4, 1/k; 1, 1+1/k; 1)` (session (iii)), Thomae's relation
singling out the parameter `c = 1/k` gives

> **S(k) = (16 / (3√2 · k · π)) · 3F2( 1, 1, 1 − 1/k ; 5/4, 7/4 ; 1 )**

and the companion (singling out `a = 1/4`)

> **S(k) = (4 / (3π√2)) · 3F2( 1, 3/4, 3/4 + 1/k ; 7/4, 1 + 1/k ; 1 )**.

Writing `T(a) = 3F2(1,1,a;5/4,7/4;1)`, so `S(k) = (16/(3√2 kπ)) T(1−1/k)`, one has the
clean symmetric double integral (verified numerically for k = 1,2,3):

> **T(a) = 3 ∫_0^1 ∫_0^1 ρ² [ 1 − (1−ρ⁴)(1−v⁴) ]^{−a} dρ dv = ∫_0^1 2F1(a,1;7/4;1−v⁴) dv.**

`T(0) = 1` gives `S(1) = 8√2/(3π)` in one line.
Numeric: `T(1/2) = 1.86967572042069…`, `T(2/3) = 2.72005441344966…`, `T(3/4) = 3.5631…`.

### 2.4  Four new exact forms for S(4)  — **PROVED**, verified to 40 digits (diff 0.0)

Taking k = 4 in 2.2, `2F1(1,1;3/2;(1−√2 sinθ)/2) = ω/sin ω` where `cos ω = √2 sin θ`:

| # | form | check |
|---|------|-------|
| F1 | `S(4) = (2/π) ∫_0^{π/2} (ω/sin ω) dχ`,  `cos ω = √2 cos χ` | diff 0.0 @40d |
| F2 | `S(4) = (4/π) ∫_0^1 ∫_0^{π/2} dσ dχ / (1 + 2√2 σ cos χ + σ²)` | diff 0.0 @40d |
| F3 | `S(4) = (2/π) ∫_0^{π/2} [ ω + arcsinh(sin ω) ] dω / √(1+sin²ω)` | 2.3e-41 |
| F4 | `S(4) = (√2/π) ∫_0^1 [artanh b − arctan b] (1−b⁴)^{−3/4} db` | 2.4e-10 (quad-limited) |
| F5 | `S(4) = (2√2/(3π)) · 3F2(1,1,3/4; 5/4,7/4; 1)` | exact (2.3) |

**F2 is the sharpest new object: a completely rational 2-form.**  Its inner σ-integral is
`∫_0^1 dσ/(σ²+2σ cos ω+1) = ω/(2 sin ω)` exactly, which is how F1 arises.

F3 is the session's (vii) in the `ω` chart: it says `S(4) = (2/π)(A+B)` with
```
A = ∫_0^{π/2} ω dω / √(1+sin²ω) = 0.95766408075068216581951644768633921482910423797955893192957579894541160142...
B = ∫_0^{π/2} arcsinh(sin ω) dω / √(1+sin²ω) = 0.72207098381440141367549478882487281801685688935499576958870909530230239891...
```
(190-dps values; `(2/π)(A+B) − S(4) = −1.8e-152`.)  Integration by parts turns A into
**`A = (1/√2) ∫_0^{π/2} F(u | 1/2) du`** — the mean of the *incomplete* elliptic integral
of the first kind at the lemniscatic modulus.  That is the obstruction: `dω/√(1+sin²ω)`
is the lemniscatic elliptic measure and `ω` is the *angle*, so A, B are quasi-period /
"amplitude-moment" objects, not periods.

### 2.5  Why k ≤ 3 closes and k ≥ 4 does not  — **CONJECTURAL (strong)**

For the integrand of 2.2, `2F1(1, 4/k; 1+2/k; ·)`, the Schwarz triple
`(|1−c|, |a−b|, |c−a−b|) = (2/k, |1−4/k|, 2/k)` has angle sum

```
 k :  triple                    sum     regime
  1 : (2,      3,      2   )    7       spherical
  2 : (1,      1,      1   )    3       spherical
  3 : (2/3,    1/3,    2/3 )    5/3     spherical
  4 : (1/2,    0,      1/2 )    1       EUCLIDEAN
  5 : (2/5,    1/5,    2/5 )    1       EUCLIDEAN
  6 : (1/3,    1/3,    1/3 )    1       EUCLIDEAN
  8 : (1/4,    1/2,    1/4 )    1       EUCLIDEAN
 12 : (1/6,    2/3,    1/6 )    1       EUCLIDEAN
```
(identically 1 for every k ≥ 4, since `4/k + (1 − 4/k) = 1`).  The same computation on
the original `2F1(1/2, 2/k; 1+2/k; ·)` gives sum `4/k` for k ≤ 4 and `1` for k ≥ 4 — the
boundary sits at k = 4 either way.

Spherical (sum > 1) = finite/small projective monodromy = the regime where these
θ-integrals collapse to logs and arctans.  Sum = 1 is the parabolic regime attached to
elliptic curves; the k = 4 period is literally an integral against the lemniscatic
measure (2.4).  **This is exactly the observed boundary: closed forms are known for
k = 1, 2, 3 and none exists for k ≥ 4.**  It is evidence, not a proof — a proof would
need a differential-Galois / Nesterenko-type argument.

---

## 3. PSLQ exclusions  (all at 380–500 dps, ≈ 510 correct digits available)

Every run reports `none` unless noted.  Runs that *did* return a relation returned one
with coefficient **0 on the target** — i.e. an internal degeneracy of my own basis, not
an evaluation.  Those degeneracies (all correct, all known) were:

* `π² − 4 log²(1+√2) − 8 Li₂(√2−1) + 8 Li₂(1−√2) = 0`  (Landen at the fixed point of x↦(1−x)/(1+x))
* `Li₂(1/φ) = π²/10 − log²φ`
* `2 E(1/2) K(1/2) − K(1/2)² = π/2`  (Legendre)
* `π/2 = arctan√2 + arctan(1/√2)`, `π = 2arctan√2 + arctan(2√2)`
* `8G + 3 Cl₂(π/3) − 12 Cl₂(π/6) = 0`
* `Ti₂(√2−1) = Cl₂(π/4) − G/4 − (π/8) log(1+√2)`  (rediscovered at residual 0.0; also derivable by hand from `Ti₂(tanθ) = θ log tanθ + Σ_{n odd} sin 2nθ/n²` at θ = π/8)
* `π² − 6 log²(1+√2) − 18 Li₂((√2−1)²) + 12 Li₂(−(√2−1)²) = 0`  (Abel at the same fixed point)

After removing these, **no relation involving π·S(4) or π·S(5) survives** in:

| target | basis | n | maxcoeff | verdict |
|--------|-------|---|----------|---------|
| π S(4) | weight 1: {π, log2, log3, log(1+√2), log(2+√3), log(√3+√2), arctan√2} × {1,√2,√3,√6} | 32 | 1e5 | none |
| π S(4) | weight 1 wide (+ arctan√6, arctan½, arctan⅓, arctan2√2, arctan(√2/3)) × {1,√2,√3,√6} | 52 | 1e5 | none |
| π S(4) | weight 2 at 8th roots: {π², G, Cl₂(π/4), log²2, πlog2, log²(1+√2), πlog(1+√2), log2·log(1+√2), arctan²√2, π arctan√2, arctan√2·log2, arctan√2·log(1+√2), Li₂(√2−1)} | 14 | 1e5 | none |
| π S(4) | same × {1,√2} | 28 | 1e5 | none |
| π S(4) | weight 2 at 12th roots (+ Cl₂(π/3), Cl₂(π/6), log²(2+√3), π log(2+√3), log²3, π log3) | 21 | 1e5 | none |
| π S(4) | elliptic {K, E, πK, πE, K·log(1+√2), K·arctan√2, K·log2, E·log(1+√2)} × {1,√2} at m=1/2 | 26 | 1e5 | none (only Legendre) |
| π S(4) | Γ(1/4): {1,√2,π,√2π,log2,log(1+√2), Γ(1/4)²/√π, √π³/Γ(1/4)², Γ(1/4)⁴/π³, π³/Γ(1/4)⁴} × {1,√2} | 14 | 1e12 | none (only 1e12-size junk) |
| π S(5) | weight 1 Q(√2,√5): {π,log2,log5,log(1+√2),log φ,log(√5+2),log(√5+√2),arctan√2,arctan√5} × {1,√2,√5,√10} | 39 | 1e5 | none |
| π S(5) | weight 2 at 20th roots (+ Cl₂(π/5), Cl₂(2π/5), log²φ, π log φ, log2·logφ, log²5, π log5) | 21 | 1e5 | none |
| π S(5) | same × {1,√5} | 42 | 1e5 | none |
| π S(6), π S(8), π S(12) | weight 1 Q(√2,√3); weight 2 at 8th and 12th roots | ≤42 | 1e5 | none |
| π S(4) | kitchen sink: weights 0/1/2 + K,E,πK,πE,K·log(1+√2),K·arctan√2,K·log2,E·log(1+√2) all in one basis | 26 | 1e5 | none |
| A alone, B alone (halves of S(4)) | {K, E, πK, πE, K·log(1+√2), K·arctan√2, K·log2, E·log(1+√2), E·arctan√2, π, π², G} | 13 | 1e5 | none |
| π S(4), π S(5), π S(6) | fully-clean dilog basis {π², G, Cl₂(π/4), log²2, πlog2, log²(1+√2), πlog(1+√2), log2·log(1+√2), Li₂((√2−1)²)} × {1,√2} | 20 | 1e5 | none |

**Sanity control** (proving the search engine works): with basis `{π, log(√3+√2), arctan√2} × {1,√3}`
PSLQ instantly recovers
`π S(3) + 2π − 6 arctan√2 − 2√3 log(√3+√2) = 0`, residual **1.7e-400**.

**False-positive scale.**  With D ≈ 400 significant digits and n ≤ 52 basis vectors, a
spurious relation needs coefficients of order 10^{D/(n+1)} ≥ 10^{7.5}; I capped
`maxcoeff` at 10⁵, so the "none" verdicts are safe.  The two survival tests actually
run (the k=3 sanity relation at 400 and 500 dps; the S(1)/S(2)/S(3) closed forms at
60 and 500 dps) both survived precision increases of >400 digits.

---

## 4. What I did *not* do / open ends

* No proof that S(4), S(5) are non-elementary.  2.5 is a heuristic (Schwarz/monodromy).
  A real proof route: show `A = (1/√2)∫_0^{π/2} F(u|1/2) du` is not in the ring generated
  by π, K(1/2), E(1/2), logs of algebraic numbers — plausible via Nesterenko/quasi-period
  arguments but not attempted here.
* F2 (the rational 2-form) has not been attacked with a genuine dilogarithm/Abel-map
  calculation.  If any lane wants one more shot at S(4), that is the object to hit:
  `S(4) = (4/π) ∫_0^1 ∫_0^{π/2} dσ dχ / (1 + 2√2 σ cos χ + σ²)`.
  Doing the χ-integral first produces `log(algebraic)/√(σ⁴−6σ²+1)`, and
  `σ⁴−6σ²+1 = (σ²−(√2−1)²)(σ²−(√2+1)²)` — a genuine elliptic curve, which is the
  same wall from a second direction.
* Only PSLQ over the listed alphabets.  Not tried: weight-3 constants, Ti₂/Cl₂ at
  arguments that are *not* rational multiples of π (e.g. `Ti₂(3−2√2)` at
  `δ = 2 arctan(3−2√2)`, which arises naturally from F1 but is not a torsion angle),
  or bases with degree-4 algebraic coefficients (2^{1/4}, 3^{1/4}).

---

## 5. One-line summary for the orchestrator

`S(1), S(2), S(3)` close because their controlling `2F1` is **spherical**; from `k = 4`
on it is **Euclidean**, and `S(4) = (2/π)∫_0^{π/2}(ω/sin ω)dχ` with `cos ω = √2 cos χ`
is an integral against the lemniscatic elliptic measure.  510-digit values are in hand;
extensive PSLQ finds nothing in weight-1 (Q(√2,√3), Q(√2,√5)), weight-2 (8th, 12th,
20th roots of unity), elliptic (K, E at m = 1/2) or Γ(1/4) bases.
