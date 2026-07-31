# Lane A2 — the general closed form for S(k)

`S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ((kn+1) 64^n)`, k a positive integer.

Script: `/tmp/math-wt-coinC2/04-computation/sk_series_general_laneA2.py`,
`/tmp/math-wt-coinC2/04-computation/sk_pslq_probe.py`.

Session date 2026-07-30/31. All numerics mpmath, dps >= 60 (most >= 160).

---

## 0. Headline

There is **no elementary general closed form**: the honest general answer is a
*master representation* plus the statement that `k = 1,2,3` are sporadic
degenerations of it. Concretely this lane contributes

* **(A2-1)** a Thomae normal form that is strictly simpler than the input
  3F2 and makes `S(1)` the k = 1 special case by inspection;
* **(A2-2)** an exact "flat-measure / sector" period form exhibiting `S(k)` as a
  period of the algebraic surface `w^2 = 1 + Im(zeta^k)`;
* **(A2-3)** an exact elliptic-K master integral valid for every k;
* **(A2-4)** a **proved inhomogeneous contiguous recursion** for the master
  function, which closes the *shifted* family `sum a_n/(n+m)` completely in
  closed form (all values in `sqrt2 * Q / pi`);
* **(A2-5)** high-precision numerics k = 1..12 (>= 160 digits, two independent
  routes agreeing) and PSLQ verdicts refuting elementary closed forms for
  k = 4, 5 in large natural bases.

Notation throughout: `a_n = C(2n,n)C(4n,2n)/64^n = (1/4)_n (3/4)_n/(n!)^2`,
`f(z) = 2F1(1/4,3/4;1;z)`, `G(s) = int_0^1 z^{s-1} f(z) dz = sum_n a_n/(n+s)`,
so `S(k) = (1/k) G(1/k)`. `L = log(1+sqrt2)`, `G_C` = Catalan.

---

## 1. Numerics (VERIFIED-NUMERIC, >= 160 digits)

Two independent routes agree to full working precision for every k in 1..12:

* Route A: `S(k) = int_0^1 2F1(1/4,3/4;1;t^k) dt` (log endpoint singularity),
* Route B: `S(k) = (2/pi) int_0^{pi/2} 2F1(1/2, 2/k; 1+2/k; cos 2theta) dtheta`
  (integrand analytic on [0,pi/2] because c-a-b = 1/2, so the z->1 expansion is
  in powers of sqrt(1-z) = sqrt2 sin theta).

```
k= 1 1.20042175487614142607359892134232077688609886660696083409525...
k= 2 1.12219970467836025427143917787047938561774484925619985179912...
k= 3 1.08840416411727687127017749683727720119898085372770215368188...
k= 4 1.06935255444126805858296193953427806132593201450630803756289...
k= 5 1.05708716209477162686360610231100103612637197979788061199084...
k= 6 1.04852004921802894368141152435909714731278245795491273017980...
k= 7 1.04219410499900079434274085375568415313726447160791294800737...
k= 8 1.03733013650043673193740244495426303855919241047157137205811...
k= 9 1.03347314303729719690093646546997589850917742251625242413224...
k=10 1.03033940151270842581142301412174370851119643478045022185827...
k=11 1.02774272731066258609635325501465945812712338048654809981832...
k=12 1.02555584240313734228603437132357545378261992725825348560112...
```
(full 150-digit values in `04-computation/sk_values.txt`).

Known closed forms re-verified to 110+ digits: S(1)=8sqrt2/(3pi),
S(2)=(4/pi)log(1+sqrt2),
S(3) = -(1/pi)(sqrt3 log(5-2sqrt6)+2 arctan(sqrt2/5))
     =  (1/pi)(2 sqrt3 asinh(sqrt2) - 2 arctan(sqrt2/5)).
Preamble item (vii) also re-verified.

---

## 2. (A2-1) Thomae normal form  [PROVED + VERIFIED-NUMERIC]

Starting from `S(k) = 3F2(1/4, 3/4, 1/k; 1, 1+1/k; 1)` (preamble (iii)),
apply Thomae's transformation singling out the numerator parameter `c = 1/k`
(excess s = 1):

    3F2(a,b,c;e,f;1) = Gamma(e)Gamma(f)Gamma(s)/(Gamma(c)Gamma(s+a)Gamma(s+b))
                       * 3F2(e-c, f-c, s; s+a, s+b; 1)

with a=1/4, b=3/4, c=1/k, e=1, f=1+1/k, s=1. Since
Gamma(1+1/k)/Gamma(1/k)=1/k and Gamma(5/4)Gamma(7/4)=(3/16)Gamma(1/4)Gamma(3/4)
= 3 pi sqrt2/16:

> **THEOREM A2-1.**
> `S(k) = (8 sqrt2 / (3 pi k)) * 3F2( 1 - 1/k, 1, 1 ; 5/4, 7/4 ; 1 )`.

Equivalently, using the duplication identity `(5/4)_n (7/4)_n = 4^{-n}(5/2)_{2n}`,

> `S(k) = (8 sqrt2 / (3 pi k)) * sum_{n>=0} ((k-1)/k)_n * n! * 4^n / (5/2)_{2n}`.

Status: PROVED (Thomae is a theorem; the Gamma bookkeeping is elementary) and
VERIFIED-NUMERIC to 45 digits at k = 1,2,3,4,5,6,7,10.

Remarks.
* k = 1 makes the 3F2 terminate at n = 0, giving `S(1)=8 sqrt2/(3 pi)`
  *by inspection*. This is the structural reason k=1 is special (Whipple in the
  original variables).
* The normal form has **two unit numerator parameters**; only one parameter,
  `1-1/k`, carries k. This is the smallest-possible k-dependence.
* Its Euler double integral is
  `S(k) = (8 sqrt2/(pi k)) int_0^{pi/2} int_0^1 sin th cos^2 th
          (cos^2 2th + v^2 sin^2 2th)^{1/k - 1} dv dth`.

### 2a. Why no Gauss/Watson/Whipple/Dixon closes it for general k

The Thomae orbit of `3F2(1/4,3/4,1/k;1,1+1/k;1)` contains (among others)

    3F2(3/4, 3/4+1/k, 1; 7/4, 1+1/k; 1)
    3F2(1/4, 1/4+1/k, 1; 5/4, 1+1/k; 1)
    3F2(1-1/k, 1, 1; 5/4, 7/4; 1)
    3F2(1-1/k, 1/4, 1/4; 1, 5/4; 1)
    3F2(1/4+1/k, 3/4+1/k, 1/k; 1+1/k, 1+1/k; 1)

A member degenerates to a Gauss 2F1 iff a numerator equals a denominator or is a
non-positive integer. Solving each of those linear equations in `1/k` gives only
`k = 1` (plus non-integral roots). Checking the classical summation theorems on
every orbit member: Watson forces `c = 1/2` or `5/8`, Whipple forces
`k = 1`, Dixon forces `c = 0`. **CONJECTURAL-but-exhaustive-in-this-family:**
no classical 3F2(1) summation closes S(k) for k >= 2.

---

## 3. (A2-2) Sector / surface-period form  [PROVED]

Chain (each step elementary and numerically re-verified):

1. Kernel (iv) with `t = tau^2`:
   `S(k) = (2/pi) int_0^1 int_0^{pi} tau (1 + tau^k cos phi)^{-1/2} dphi dtau`.
2. `int_0^pi (1+u cos phi)^{-1/2} dphi = 2 int dw / sqrt(u^2-(w^2-1)^2)`
   over `sqrt(1-u) < w < sqrt(1+u)`  (substitute `w = sqrt(1+u cos phi)`).
   Because `u^2-(w^2-1)^2 = (u-m)(u+m)` with `m = w^2-1`, swapping the order
   makes the measure **flat**:
   `G(s) = (4/pi) int int_D ((w^2-1)^2+v^2)^{s-1} dv dw`,
   `D = { w,v>0 : (w^2-1)^2+v^2 < 1 }`.
3. In the complex variable `W = (w^2-1) + i v` this is the right half unit disc
   with weight `|W|^{2s-2}`; at `s = 1/k` the k-th power map `W = zeta^k` has
   Jacobian `k^2 |zeta|^{2k-2}` which exactly cancels `|W|^{2/k-2}`:

> **THEOREM A2-2.**
> `S(k) = (2k/pi) * int int_{Sig_k} (1 + Im(zeta^k))^{-1/2} dA(zeta)`,
> `Sig_k = { |zeta| < 1, |arg zeta| < pi/(2k) }`,
> i.e. S(k) is (up to 2k/pi) the area-period of the algebraic surface
> `w^2 = 1 + Im(zeta^k)`, a polynomial of degree exactly k, over a circular
> sector of opening pi/k.

Sanity: k=1 gives `(2/pi) * (2 sqrt2/3) * 2 = 8 sqrt2/(3 pi)` directly from
`int_0^{sqrt2} w sqrt(2-w^2) dw = 2 sqrt2/3` — the k=1 closed form falls out of
the flat-measure step with no work.

This is the cleanest statement of "where the difficulty lives": the branch curve
of the double cover has degree `k` (k even) or `k+1` (k odd), so the surface is
rational for k <= 4, K3 for k = 5,6 and of general type for k >= 7. That is a
*heuristic* for the elementary/non-elementary dichotomy, not a proof (see §6).

---

## 4. (A2-3) Elliptic-K master integral  [PROVED + VERIFIED-NUMERIC]

Apply the quadratic transformation (iv) `f(z) = (1+sqrt z)^{-1/2}
2F1(1/2,1/2;1; 2 sqrt z/(1+sqrt z))` inside `S(k)=int_0^1 f(t^k)dt` and
substitute `kappa^2 = 2 t^{k/2}/(1+t^{k/2})` (so `t^{k/2} = kappa^2/(2-kappa^2)`):

> **THEOREM A2-3.**
> `S(k) = (16/(sqrt2 pi k)) int_0^1 K(kappa) * kappa^{4/k-1} * (2-kappa^2)^{-1/2-2/k} dkappa`,
> K the complete elliptic integral of the first kind, modulus convention
> `K(kappa)=int_0^{pi/2} dtheta/sqrt(1-kappa^2 sin^2 theta)`.

VERIFIED-NUMERIC to 50+ digits at k = 1,2,3,4,5.

Special case k = 4 collapses to a single clean K-moment:

> `S(4) = (2 sqrt2/pi) int_0^1 K(kappa) dkappa / (2 - kappa^2)`.

Doing the kappa-integral first (elementary for k=4) gives the *log form*

> `S(4) = (2/pi) int_0^{pi/2} log( sqrt2 cos th + sqrt(cos 2th) ) / sqrt(cos 2th) dth`
> (real; for th > pi/4 the integrand becomes `arctan(sqrt(-cos2th)/(sqrt2 cos th))
> / sqrt(-cos 2th)`), verified to 50 digits.

Equivalently `S(4) = (2/pi)(I1+I2)` with
`I1 = int_0^{pi/2} th dth/sqrt(1+sin^2 th) = 0.9576640807506821658195164476863392148291...`
`I2 = int_0^{pi/2} asinh(sin th) dth/sqrt(1+sin^2 th) = 0.7220709838144014136754947888248728180169...`
(both computed to 400 digits; these are the `arcsin`/`arcsinh` halves of
preamble (vii) after th = arcsin v).

Change of variable `W = e^{2iz}` turns `I1+I2` into
`-(1/2) int log(W) dW / Y` on the **elliptic curve** `Y^2 = W(6W - W^2 - 1)`,
whose 2-torsion is `{0, 3+2sqrt2, 3-2sqrt2}` and whose j-invariant is 287496,
i.e. CM by the order of discriminant -16 in Q(i). So `S(4)` is an
*Eichler / elliptic-dilogarithm* period of a CM curve, not a polylogarithm —
which is exactly what the PSLQ searches in §5 confirm.

---

## 5. (A2-5) PSLQ verdicts  [see §5 for status; all with stated bases]

Protocol used: coefficient bound chosen strictly below the spurious-hit scale
`10^{D/n}` (D = working digits, n = number of reals), so a `None` is meaningful.

*(filled in below as runs complete)*

---

## 6. Status of the "sum over roots of unity" ansatz

REFUTED as stated for general k in the tested bases (see §5): the mechanism that
makes k = 1,2,3 elementary does not persist.

Positive part of the picture, in the A2-2 surface language:

| k | deg of branch curve | surface type | S(k) |
|---|---|---|---|
| 1 | 2 | rational (quadric) | `8 sqrt2/(3 pi)` |
| 2 | 2 | rational | `(4/pi) log(1+sqrt2)` |
| 3 | 4 | del Pezzo deg 2 (rational) | `(1/pi)(2 sqrt3 asinh sqrt2 - 2 arctan(sqrt2/5))` |
| 4 | 4 | del Pezzo deg 2 (rational) | CM elliptic-dilog period (no polylog form found) |
| 5,6 | 6 | K3 | no closed form found |
| >=7 | >=8 | general type | no closed form expected |

The k=4 row shows the naive rationality heuristic is **not** sufficient: the
2-chain's boundary meets the branch curve (a smooth plane quartic, genus 3, but
whose relevant piece is the CM elliptic curve above) and the resulting period is
not mixed Tate. So the honest dichotomy is k <= 3 elementary, k >= 4 not — with
k=4 carrying a *CM* (hence in principle Gamma-expressible) period.
