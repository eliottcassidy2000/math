# Lane A2 — the general closed form for S(k)

`S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ((kn+1) 64^n)`, k a positive integer.

Script (all claims reproducible): `/tmp/math-wt-coinC2/04-computation/sk_series_general_laneA2.py`
(helpers: `sk_pslq_probe.py`; data: `i12_400.txt`, `s35_250.txt`, `s5_pslq.log`).

Notation. `a_n = C(2n,n)C(4n,2n)/64^n = (1/4)_n (3/4)_n/(n!)^2`,
`f(z) = 2F1(1/4,3/4;1;z)`, `G(s) = int_0^1 z^{s-1} f(z) dz = sum_n a_n/(n+s)`,
so `S(k) = (1/k) G(1/k)`.  `L = log(1+sqrt2)`, `G_C` = Catalan,
`varpi = int_0^1 dv/sqrt(1-v^4) = Gamma(1/4)^2/(4 sqrt(2 pi)) = 1.31102877714605990523...`.

---

## 0. Headline

**There is no elementary general closed form.**  `k = 1, 2, 3` are three separate
degenerations, and the mechanism dies at `k = 4`: `S(4)` is provably an
*Eichler / elliptic-dilogarithm* period of the CM elliptic curve
`Y^2 = W(6W - W^2 - 1)` (j = 66^3 = 287496, CM by the order of discriminant −16
in Q(i)), and it survives an 8-digit-coefficient PSLQ search at 380 digits
against a 25-element polylogarithmic + lemniscatic basis.

What this lane *does* deliver as the honest answer to "Question 2" is a
**master normal form with the minimum possible k-dependence** plus three exact
representations and one exact recursion:

| # | statement | status |
|---|---|---|
| A2-1 | `S(k) = (8 sqrt2/(3 pi k)) 3F2(1-1/k, 1, 1; 5/4, 7/4; 1)` | PROVED + numeric to 251 digits |
| A2-2 | `S(k) = (2k/pi) ∫∫_{|zeta|<1, |arg zeta|<pi/2k} (1+Im zeta^k)^{-1/2} dA` | PROVED + numeric |
| A2-3 | `S(k) = (16/(sqrt2 pi k)) ∫_0^1 K(kappa) kappa^{4/k-1}(2-kappa^2)^{-1/2-2/k} dkappa` | PROVED + numeric |
| A2-4 | `(s-1)^2 G(s-1) = (s-1/4)(s-3/4) G(s) - 1/(pi sqrt2)` (Re s > 1) | PROVED + numeric |
| A2-5 | S(4) is a CM elliptic Eichler period; no polylog form | CONJECTURAL (strong) |

---

## 1. Numerics  [VERIFIED-NUMERIC, >= 160 digits]

Two independent routes agree to full working precision, k = 1..12:
Route A `S(k) = int_0^1 2F1(1/4,3/4;1;t^k)dt`;
Route B `S(k) = (2/pi) int_0^{pi/2} 2F1(1/2,2/k;1+2/k;cos 2th)dth` (integrand
analytic on [0,pi/2] since c-a-b = 1/2 makes the z->1 expansion run in powers of
`sqrt(1-z) = sqrt2 sin th`).

```
k= 1 1.20042175487614142607359892134232077688609886660696083409525
k= 2 1.12219970467836025427143917787047938561774484925619985179912
k= 3 1.08840416411727687127017749683727720119898085372770215368188
k= 4 1.06935255444126805858296193953427806132593201450630803756289
k= 5 1.05708716209477162686360610231100103612637197979788061199084
k= 6 1.04852004921802894368141152435909714731278245795491273017980
k= 7 1.04219410499900079434274085375568415313726447160791294800737
k= 8 1.03733013650043673193740244495426303855919241047157137205811
k= 9 1.03347314303729719690093646546997589850917742251625242413224
k=10 1.03033940151270842581142301412174370851119643478045022185827
k=11 1.02774272731066258609635325501465945812712338048654809981832
k=12 1.02555584240313734228603437132357545378261992725825348560112
```
Given closed forms re-verified to 110+ digits, including
`S(3) = (1/pi)(2 sqrt3 asinh(sqrt2) - 2 arctan(sqrt2/5))`
(`asinh sqrt2 = log(sqrt2+sqrt3)`, `5-2sqrt6 = (sqrt3-sqrt2)^2`), and preamble (vii).

---

## 2. THEOREM A2-1 — Thomae normal form  [PROVED]

From `S(k) = 3F2(1/4, 3/4, 1/k; 1, 1+1/k; 1)` (excess s = 1), apply Thomae
singling out the numerator `c = 1/k`:

    3F2(a,b,c;e,f;1) = Γ(e)Γ(f)Γ(s)/(Γ(c)Γ(s+a)Γ(s+b)) · 3F2(e-c, f-c, s; s+a, s+b; 1)

with `a=1/4, b=3/4, c=1/k, e=1, f=1+1/k, s=1`.  Then `Γ(1+1/k)/Γ(1/k)=1/k`
and `Γ(5/4)Γ(7/4) = (3/16)Γ(1/4)Γ(3/4) = 3 pi sqrt2/16`, giving

> ### `S(k) = (8 sqrt2 / (3 pi k)) · 3F2( 1 - 1/k, 1, 1 ; 5/4, 7/4 ; 1 )`

Equivalently, using `(5/4)_n (7/4)_n = 4^{-n} (5/2)_{2n}` (Legendre duplication),

> `S(k) = (8 sqrt2 / (3 pi k)) · sum_{n>=0} ((k-1)/k)_n · n! · 4^n / (5/2)_{2n}`.

VERIFIED-NUMERIC: residual 0 to 45 digits at k = 1,2,3,4,5,6,7,10, and to
**251 digits** at k = 3 (against the known S(3) closed form) and k = 5.

Why this is the right normal form.
* Only **one** parameter, `1-1/k`, carries k; the other four are absolute.  This
  is the minimum k-dependence achievable in the 3F2(1) Thomae orbit.
* `k = 1` makes the series terminate at n = 0, so `S(1) = 8 sqrt2/(3 pi)` *by
  inspection*.  That is the structural reason k = 1 is special.
* Euler double integral of the same object (VERIFIED-NUMERIC to quadrature
  accuracy at k = 1..5):
  `S(k) = (8 sqrt2/(pi k)) ∫_0^1 ∫_0^{pi/2} sin th cos^2 th
          (cos^2 2th + v^2 sin^2 2th)^{1/k-1} dth dv`.
  (Route: `n!/(7/4)_n = (3/4)B(n+1,3/4)`, `(A)_n/(5/4)_n` a Beta, then
   `x = sin^2 th`, and `2F1(A,1;3/2;z) = ∫_0^1 dv/(1-z+z v^2)^A`.)

### 2a. No classical summation closes S(k) for k >= 2  [PROVED for the orbit]

Thomae orbit members include

    3F2(3/4, 3/4+1/k, 1; 7/4, 1+1/k; 1)
    3F2(1/4, 1/4+1/k, 1; 5/4, 1+1/k; 1)
    3F2(1-1/k, 1, 1; 5/4, 7/4; 1)
    3F2(1-1/k, 1/4, 1/4; 1, 5/4; 1)
    3F2(1/4+1/k, 3/4+1/k, 1/k; 1+1/k, 1+1/k; 1)

Degeneration to a Gauss 2F1 needs a numerator = a denominator or a non-positive
integer; every such linear equation in 1/k has only the root `k = 1` (plus
non-integral roots).  Watson forces `c ∈ {1/2, 5/8, 7/8}` with an incompatible
`(a+b+1)/2`; Whipple forces `k = 1`; Dixon forces `c = 0`.  So Watson/Whipple
(and Dixon) really do apply only at k = 1, on the *whole* orbit — not just on
the original form.

---

## 3. THEOREM A2-2 — flat measure and the sector/surface period  [PROVED]

Two exact steps, both of independent use.

**(i) Flattening.**  In
`G(s) = (2/pi) ∫∫_{X^2+Y^2<1, Y>0} (X^2+Y^2)^{s-1}(1+X)^{-1/2} dX dY`
(this is preamble (iv) with `z = r^2`, `phi` traded for `w = sqrt(1+r cos phi)`),
substitute `X = w^2 - 1`.  Then `(1+X)^{-1/2} dX = 2 dw`, i.e. **the measure
becomes flat**:

    G(s) = (4/pi) ∫∫_D ((w^2-1)^2 + v^2)^{s-1} dv dw,
    D = { w,v > 0 : (w^2-1)^2 + v^2 < 1 }.

At `s = 1` the integrand is 1 and `int_0^{sqrt2} w sqrt(2-w^2) dw = 2 sqrt2/3`
returns `S(1) = 8 sqrt2/(3 pi)` with no work — a one-line re-derivation of k = 1.

**(ii) k-th power map.**  Write `W = (w^2-1) + i v`, so the region is the right
half unit disc with weight `|W|^{2s-2}`.  At `s = 1/k` put `W = zeta^k`; the
Jacobian `k^2|zeta|^{2k-2}` cancels `|W|^{2/k-2}` exactly, and the k-to-1 cover
contributes `1/k`:

> ### `S(k) = (2k/pi) ∫∫_{Sigma_k} (1 + Im(zeta^k))^{-1/2} dA(zeta)`,  `Sigma_k = { |zeta| < 1, |arg zeta| < pi/(2k) }`

VERIFIED-NUMERIC (30 digits, k = 1..4).

So **S(k) is the area-period of the algebraic surface `w^2 = 1 + Im(zeta^k)`**
— a *polynomial of degree exactly k* — over a circular sector of opening pi/k.
The singular locus `1+Im zeta^k = 0` sits exactly at the two corners
`zeta = e^{±i pi/(2k)}` of the sector.  For k = 3 the polynomial factors as
`1 + y(sqrt3 x - y)(sqrt3 x + y)`, which is where the `sqrt3` in S(3) comes from;
for k = 2 it is `1 + 2xy` (rotate 45°: `1 + u^2 - v^2`), the `sqrt2`.

Branch-curve degrees: `k` (k even) or `k+1` (k odd), so the double cover is
rational for k <= 4, K3 for k = 5,6, of general type for k >= 7.  This is a
*heuristic*, and §5 shows it is **not sufficient**: k = 4 is already
non-polylogarithmic.

---

## 4. THEOREM A2-3 — elliptic-K master integral  [PROVED]

Insert the quadratic transformation (iv) into `S(k) = ∫_0^1 f(t^k) dt` and
substitute `kappa^2 = 2 t^{k/2}/(1+t^{k/2})`, i.e. `t^{k/2} = kappa^2/(2-kappa^2)`:

> ### `S(k) = (16/(sqrt2 pi k)) ∫_0^1 K(kappa) · kappa^{4/k-1} · (2-kappa^2)^{-1/2-2/k} dkappa`

with `K(kappa) = ∫_0^{pi/2} dtheta/sqrt(1-kappa^2 sin^2 theta)`.
VERIFIED-NUMERIC to 50 digits at k = 1,2,3,4,5.
(Note `2 - kappa^2 = 1 + kappa'^2`, so the weight is `kappa^{4/k-1}(1+kappa'^2)^{-1/2-2/k}`.)

The inner `kappa`-integral (after swapping in K's integral) is elementary
exactly when `4/k - 1` and `-1/2-2/k` cooperate: k = 1, 2, 4 only.  In
particular the beautifully short

> `S(4) = (2 sqrt2 / pi) ∫_0^1 K(kappa) dkappa / (2 - kappa^2)`,

and doing that integral first gives the **log form** (verified to 50 digits)

> `S(4) = (2/pi) ∫_0^{pi/2} log( sqrt2 cos th + sqrt(cos 2th) ) / sqrt(cos 2th) dth`

(real throughout: for th > pi/4 the integrand equals
`arctan(sqrt(-cos 2th)/(sqrt2 cos th))/sqrt(-cos 2th)`, because
`|sqrt2 cos th + i sqrt(-cos 2th)| = 1`).

---

## 5. THEOREM A2-4 — contiguous recursion for the master function  [PROVED]

Integrate `z^{s-1} · [ z(1-z) f'' + (1-2z) f' - (3/16) f ] = 0` (the
hypergeometric ODE for `f = 2F1(1/4,3/4;1;·)`) over (0,1).  Every boundary term
collects into `z^{s-1}(1-z)[ z f' - (s-1) f ]`, which -> `1/(pi sqrt2)` at z = 1
(since `f ~ -(1/(pi sqrt2)) log(1-z)`) and -> 0 at z = 0 for Re s > 1.  Hence

> ### `(s-1)^2 G(s-1) = (s - 1/4)(s - 3/4) G(s) - 1/(pi sqrt2)`,  Re s > 1.

VERIFIED-NUMERIC: residual at quadrature accuracy for
s = 1.5, 1.7, 2, 2.2, 2.25, 7/3, 9/4, 11/5.  (My first attempt dropped the
boundary term and produced a *constant* defect −0.22508…, which is exactly
−1/(pi sqrt2) — that is how the term was found.)

Corollaries.
* **Shifted family closed in full**: `sum_n a_n/(n+m) = sqrt2 · r_m / pi` with
  `r_1 = 8/3`, `r_{m+1} = (m^2 r_m + 1/2)/((m+3/4)(m+1/4))`.  First values
  `8/3, 152/105, 10568/10395, 178328/225225, 47453768/72747675,
   782539544/1405485081, 51331850888/105411381075, 838519635608/1933976154825`
  (checked to 60 digits).  Only at m = 1 does this equal
  `Γ(m)^2/(Γ(m+1/4)Γ(m+3/4))`; the inhomogeneity breaks the Gamma pattern for m >= 2.
* Expanding the recursion at s = 1 (order eps^1) gives
  `sum_n a_n/(n+1)^2 = (16/3)(8 sqrt2/(3 pi) - 1) = 1.06891602600608760572586...`
  (VERIFIED-NUMERIC to 11 digits by direct summation; the series converges like
  `1/(pi n^3)`).
* Structurally: writing `G(s) = h(s) P(s)` with
  `h(s) = Γ(s)^2/(Γ(s+1/4)Γ(s+3/4))`, the recursion becomes
  `P(s) - P(s-1) = (1/(pi sqrt2)) Γ(s-3/4)Γ(s-1/4)/Γ(s)^2`.
  The summand is a Gamma ratio, not a rational function of s, so `P` is a
  *generalised* digamma — direct evidence that `G(1/k)` has no elementary form.

---

## 6. THEOREM A2-5(conj) — what S(4) actually is

`S(4) = (2/pi)(I1 + I2)` where (preamble (vii) after `th = arcsin v`; both
integrands analytic, so 400 digits are cheap)

```
I1 = ∫_0^{pi/2} th dth / sqrt(1+sin^2 th)          = 0.95766408075068216581951644768633921482910423797955893192957579894541...
I2 = ∫_0^{pi/2} asinh(sin th) dth / sqrt(1+sin^2 th)= 0.72207098381440141367549478882487281801685688935499576958870909530230...
I1+I2 = pi S(4)/2                                   = 1.67973506456508357949501123651121203284596112733...
```

Both equal `∫ z dz/sqrt(2-cos^2 z)`: `I1` along `[0, pi/2]`, `I2` along
`[0, iL]` (`L = log(1+sqrt2)`, verified numerically), so

> `I1 + I2 = ∫_{iL}^{pi/2} z dz / sqrt(2 - cos^2 z)` ,

whose endpoint `z = iL` is a branch point (`cos(iL) = cosh L = sqrt2`).
Substituting `W = e^{2iz}` (so `z = log(W)/(2i)`, `2-cos^2 z = (3-cos 2z)/2`):

> `I1 + I2 = -(1/2) ∫ log(W) dW / Y`  on the elliptic curve  `Y^2 = W(6W - W^2 - 1)`,

whose 2-torsion is `{0, 3+2 sqrt2, 3-2 sqrt2}`, `lambda = 17-12 sqrt2`, and

> **j = 256(1-λ+λ²)³/(λ²(1-λ)²) = 287496 = 66³ exactly** (verified to 44 digits)

i.e. **CM by the order of discriminant −16 in Q(i)** (the lemniscatic family, as
one expects from `y^2 = 1-v^4`).  So `S(4)` is an Eichler integral / elliptic
dilogarithm of a CM curve — a weight-3 CM L-value type constant, *not* a
polylogarithm.  Consistent with §7: no polylogarithmic relation exists.

---

## 7. PSLQ verdicts  [protocol stated; all negative results are bounded]

Protocol.  With n reals at D digits, PSLQ manufactures spurious relations of
size ~`10^{D/n}`.  I therefore always (a) run PSLQ on the basis alone first — a
hit with zero leading coefficient is an internal basis identity, not a hit — and
(b) cap `maxcoeff` far below `10^{D/n}`, so a `None` is meaningful.

**Internal identity caught this way** (it poisoned an early run):
`-pi^2 + 4 L^2 + 8 Li2(sqrt2-1) - 8 Li2(1-sqrt2) = 0`,
i.e. `chi_2(sqrt2-1) = pi^2/16 - (1/4) log^2(1+sqrt2)` — the classical Landen
value.  Removing `Li2(1-sqrt2)` made the 24-element basis independent (PSLQ on
the basis alone returns `None` at 380 digits, maxcoeff `10^8`).

### 7a. k = 4 — REFUTED within a 24-element basis
Target: each of `I1`, `I2`, `I1+I2 = pi S(4)/2`.
Basis (24 = 12 x {1, sqrt2}):
`pi^2, pi log2, pi L, log^2 2, log2·L, L^2, G_C, varpi·pi, varpi·log2, varpi·L,
 varpi^2, Li2(sqrt2-1)`, each also times `sqrt2`.
Working precision **380 digits** (I1, I2 known to 400+; stability 1e-422).
Spurious scale `10^{380/25} = 10^{15.2}`.

| maxcoeff | I1 | I2 | I1+I2 |
|---|---|---|---|
| 1e5 | None | None | None |
| 1e8 | None | None | None |

Also **None** for the smaller bases `{pi^2, pi L, L^2, G_C}×{1,sqrt2}` and
`{varpi pi, varpi L, varpi log2, varpi^2}×{1,sqrt2}` at 280 digits with
maxcoeff `1e12`, and for the 14-element weight-2 basis at 150/180 digits (whose
apparent hits grew with precision — textbook false positives, reported as such).

Multiplicative test (is `pi S(4)/2` a Gamma/pi monomial?): PSLQ on
`log|X|` against `{log pi, log 2, log Γ(1/8), log Γ(3/8)}` at 280 digits,
maxcoeff `1e8` → **None** for I1, I2, I1±I2.  (`log Γ(1/4)` and `log(1+sqrt2)`
had to be dropped first: `pi Γ(1/8)^2/(Γ(3/8)^2 Γ(1/4)^2 (1+sqrt2)) = 1`.)

`findpoly(pi S(4), deg<=8, coeff<=1e8)` → None; same for `S(4)`.

**Verdict: `pi S(4)` is CONJECTURALLY not a polylogarithmic (weight <= 2) or
lemniscatic-product constant.**  Caveat stated honestly: PSLQ refutations are
basis-limited, and the S(3) answer shows the required algebraic numbers
(`sqrt3`, `sqrt2+sqrt3`, `5 + i sqrt2` of norm 27) are *not* guessable a priori.
The independent structural evidence in §6 is what makes this verdict strong.

### 7b. The "2 c_k log(sqrt2 + c_k) − 2 arctan(t_k)" ansatz — REFUTED
The S(3) form fits `pi S(k) = 2 c_k log(sqrt2 + c_k) - 2 arctan(t_k)` with
`c_k = 2 cos(pi/2k)` (`c_3 = sqrt3`).  Solving for `t_k` at 100 digits:

| k | c_k | implied t_k | minimal polynomial (deg<=8, coeff<=1e6) |
|---|---|---|---|
| 2 | sqrt2 | −0.30098471622657998077… | **none** |
| 3 | sqrt3 | 0.28284271247461900976… | `25 t^2 - 2` → `t = sqrt2/5` ✓ |
| 4 | sqrt(2+sqrt2) | 0.55272092804101081737… | **none** |
| 5 | sqrt((5+sqrt5)/2) | 0.71374939860309665581… | **none** |
| 6 | sqrt(2+sqrt3) | 0.81905141422841897727… | **none** |

The ansatz reproduces `sqrt2/5` exactly at k = 3 and fails **even at k = 2**
where the answer *is* elementary.  So the S(3) shape is a coincidence of k = 3,
not a template.  **REFUTED.**

### 7c. k = 4 — the "right" cyclotomic weight-1 basis also fails
The sector for k = 4 has corners at `zeta = e^{±i pi/8}`, so the natural field is
`Q(zeta_16)^+ = Q(sqrt(2+sqrt2))`.  Put `c_4 = 2cos(pi/8) = sqrt(2+sqrt2)` and

* coefficients `{1, sqrt2, c_4, sqrt2 c_4}` (note `sqrt(2-sqrt2) = (sqrt2-1)c_4`,
  which PSLQ caught as an internal identity on the first attempt — removed),
* logs `{pi, log 2, log(1+sqrt2), log(sqrt2 + c_4), log(1 + c_4)}`.

20-element independent basis, **360 digits**, spurious scale `10^{17}`:

| maxcoeff | `pi S(4)` | `pi S(4)/2` |
|---|---|---|
| 1e5 | None | None |
| 1e7 | None | None |

So the S(3)-style "algebraic coefficient × log(algebraic)" shape is absent at
k = 4 even in the cyclotomic field the geometry hands you.

### 7d. k = 5
`S(5)` computed to 250 digits from A2-1.  All **None**:

| basis | size | digits | spurious scale | maxcoeff | result |
|---|---|---|---|---|---|
| `{1,sqrt2,sqrt5,sqrt10,c_5,sqrt2 c_5}×{pi,log2,L,log5,log phi,log(sqrt2+c_5)}` | 36 | 240 | 1e6.5 | 1e5 | None |
| `{pi^2,pi log2,pi L,pi log5,pi log phi,L^2,log^2 2,log^2 phi,log2 L,L log phi,log2 log phi,G_C}×{1,sqrt5}` | 24 | 240 | 1e9.6 | 1e5 | None |
| `{1,sqrt2,sqrt5,c_5}×{pi,log2,L,log(sqrt2+c_5)}` | 16 | 235 | 1e13 | 1e6 | None |
| `{1,sqrt2,sqrt5,sqrt10}×{pi,log2,L,log phi,log5}` | 20 | 235 | 1e11 | 1e6 | None |
| `{pi^2,pi L,pi log2,pi log phi,L^2,log^2 phi,G_C}×{1,sqrt5}` | 14 | 235 | 1e15 | 1e6 | None |

(`c_5 = 2cos(pi/10) = sqrt((5+sqrt5)/2)`, `phi = (1+sqrt5)/2`; internal-relation
pre-check `None` for every basis.)  Consistent with the K3 prediction of §3.

---

## 8. Answer to Question 2, as I would state it

> **S(k) has no closed form in elementary functions for general k.**  The
> sharpest general statements are
>
>   `S(k) = (8 sqrt2/(3 pi k)) · 3F2(1 - 1/k, 1, 1; 5/4, 7/4; 1)`
>        `= (16/(sqrt2 pi k)) ∫_0^1 K(kappa) kappa^{4/k-1}(2-kappa^2)^{-1/2-2/k} dkappa`
>        `= (2k/pi) ∫∫_{|zeta|<1, |arg zeta|<pi/(2k)} (1 + Im zeta^k)^{-1/2} dA(zeta)`,
>
> which realise `S(k)` as an area-period of the degree-k surface
> `w^2 = 1 + Im(zeta^k)`.  The k = 1, 2, 3 evaluations are the three
> degenerations of that surface that stay rational *with a rational boundary
> configuration*: k = 1 because the Thomae normal form terminates, k = 2 because
> `1 + Im zeta^2 = 1 + 2xy` is a quadric, k = 3 because
> `1 + Im zeta^3 = 1 + y(sqrt3 x - y)(sqrt3 x + y)` splits over `Q(sqrt3)`
> (whence the `sqrt3` in the answer).  At k = 4 the period is already an
> Eichler integral of the CM elliptic curve `Y^2 = W(6W-W^2-1)`, `j = 66^3`;
> for k >= 5 the surface is K3 or of general type.

Companion exact results proved along the way (these *are* general and closed):

* `sum_n a_n/(n+m) = sqrt2 r_m/pi`, `r_m ∈ Q` explicit, for every integer m >= 1;
* `sum_n a_n/(n+1)^2 = (16/3)(8 sqrt2/(3 pi) - 1)`;
* `(s-1)^2 G(s-1) = (s-1/4)(s-3/4) G(s) - 1/(pi sqrt2)` for Re s > 1.

## 9. Open / next

1. Identify `pi S(4)/2` as an explicit weight-3 CM L-value (level 16 or 64,
   `eta(4 tau)^6`-type); if that works, `S(4)` gets a genuine closed form of a
   *different* kind, and the k=4 row of the table is filled.
2. `S(6)`: `1 + Im zeta^6` gives a K3; check whether it is *singular* (rank 20),
   in which case `S(6)` should also be a CM period.
3. The A2-1 normal form `3F2(1-1/k,1,1;5/4,7/4;1)` is the right object to feed to
   any 3F2(1) contiguous-relation / creative-telescoping machinery: it has one
   free parameter, so Zeilberger in that parameter is cheap.
