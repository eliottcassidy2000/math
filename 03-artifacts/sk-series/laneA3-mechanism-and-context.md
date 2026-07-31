# Lane A3 — the k = 3 mechanism, and the context/provenance cross-checks

`S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ((kn+1) 64^n)`, k a positive integer.

Script (everything below is reproducible from it):
**`/tmp/math-wt-coinC2/04-computation/sk_series_mechanism_laneA3.py`**

```
python3 sk_series_mechanism_laneA3.py verify    # the reduction chain, dps 80 (ids) / 30 (2-D quad)
python3 sk_series_mechanism_laneA3.py xdeg      # THE MECHANISM (one table)
python3 sk_series_mechanism_laneA3.py k3        # the explicit k=3 evaluation, 8 steps
python3 sk_series_mechanism_laneA3.py pencil    # the one-parameter fibre pencil
python3 sk_series_mechanism_laneA3.py sep       # separation normal form + the red herring
python3 sk_series_mechanism_laneA3.py schwarz   # Schwarz triples WITH the reduction
python3 sk_series_mechanism_laneA3.py known     # which neighbours are classically summable
python3 sk_series_mechanism_laneA3.py thomae    # the S_5 orbit and its summability scan
```

---

## 0. Headline

**The mechanism is one-variable elementary calculus, and the wall at k = 4 is Liouville's
theorem, not a Schwarz-list, monodromy, or CM phenomenon.**

Take lane A2's sector form (A2-2), which I re-verified independently:

> `S(k) = (2k/pi) ∫∫_{Sigma_k} dx dy / sqrt(f_k)`,  `f_k = 1 + Im (x+iy)^k`,
> `Sigma_k = { x^2+y^2 < 1, |y| < x tan(pi/2k) }`.

`Sigma_k` is **x-simple**: for each `y`, `x` runs from `x0 = |y| cot(pi/2k)` (a ray) to
`x1 = sqrt(1-y^2)` (the arc). So **integrate in x first**. Then everything follows from

> ### `deg_x ( 1 + Im (x+iy)^k ) = k - 1` exactly

| k | `f_k = 1 + Im (x+iy)^k` | `deg_x` | `∫ dx / sqrt(f_k)` |
|---|---|---|---|
| 1 | `1 + y` | 0 | algebraic |
| 2 | `1 + 2xy` | 1 | algebraic, `2 sqrt(1+2xy)/(2y)` |
| **3** | `1 + 3x^2 y - y^3` | **2** | **arcsinh / arcsin — LAST elementary case** |
| 4 | `1 + 4x^3 y - 4x y^3` | 3 | incomplete **elliptic F** — the wall |
| 5 | `1 + 5x^4y - 10x^2y^3 + y^5` | 4 | elliptic (genus 1) |
| 6 | `1 + 6x^5y - 20x^3y^3 + 6xy^5` | 5 | hyperelliptic (genus 2) |
| 7 | … | 6 | hyperelliptic (genus 2) |

`∫ dx/sqrt(P(x))` is elementary iff `deg P <= 2`. Hence the inner integral is elementary
**iff k <= 3**, and at k = 4 it is an incomplete elliptic integral of the **first** kind —
which is *precisely* lane A1's independently-found `A = (1/sqrt2) ∫_0^{pi/2} F(u | 1/2) du`.
That agreement between two completely different routes is the strongest confirmation
available that this is the right mechanism.

This also predicts the observed **weight-1 alphabet** `{pi, log(alg), arctan(alg)}`: the
inner integral produces exactly an `arcsinh`/`arcsin`, and a single integration by parts in
`y` turns the outer integral into a purely **algebraic** one.

---

## 1. The explicit k = 3 evaluation (every step numerically verified)

Target `pi S(3)/6 = 0.5698870876878961109766952540969…`.

**[1]** `pi S(3)/6 = ∫_{-1/2}^{1/2} ∫_{sqrt3|y|}^{sqrt(1-y^2)} dx dy / sqrt(1 + 3x^2y - y^3)`.
The y-range is `[-1/2, 1/2]` because the ray and arc meet at the sector corners
`zeta = e^{±i pi/6}` — which are exactly the zeros of `1 + Im zeta^3` on `∂Sigma_3`.
**VERIFIED-NUMERIC** (resid 1.7e-39).

**[2]** `f_3 = (1-y^3) + 3y x^2` is quadratic in x, so

```
 y > 0 :  inner = (3y)^{-1/2} [ asinh sqrt(3y(1+y)/(1+y+y^2)) − asinh sqrt(9y^3/(1-y^3)) ]
 y =-z<0: inner = (3z)^{-1/2} [ asin  sqrt(3z(1-z)/(1-z+z^2)) − asin  sqrt(9z^3/(1+z^3)) ]
```
**VERIFIED-NUMERIC** (resid 0.0 at dps 40). The first argument is the **arc** endpoint,
the second the **ray** endpoint.

**[3–5]** IBP in y against `∫(3y)^{-1/2}dy = (2/sqrt3) sqrt y`. The boundary terms at
`y = ±1/2` **cancel** (both arguments are equal at the corner) and vanish at 0. The
remainder is *purely algebraic*, because of four exact collapses (all verified):

```
 1 + 3y(1+y)/(1+y+y^2) = (1+2y)^2/(1+y+y^2)        1 + 9y^3/(1-y^3) = (1+2y)(1-2y+4y^2)/(1-y^3)
 1 - 3z(1-z)/(1-z+z^2) = (1-2z)^2/(1-z+z^2)        1 - 9z^3/(1+z^3) = (1-2z)(1+2z+4z^2)/(1+z^3)
```

> `pi S(3)/6 = −(2/sqrt3) ∫_0^{1/2} [ (A1−A2) + (B1−B2) ] dy` ,
> `A1 = sqrt3/(2(1+y+y^2) sqrt(1+y))`   `A2 = 9y/(2(1-y^3) sqrt(1+8y^3))`
> `B1 = sqrt3/(2(1-y+y^2) sqrt(1-y))`   `B2 = 9y/(2(1+y^3) sqrt(1-8y^3))`

**VERIFIED-NUMERIC to ~32 digits** (2.4e-32 at dps 60; quadrature-limited by the integrable
endpoint singularity of `B2` at `y = 1/2`). **PROVED** step by step (each collapse is an
algebraic identity).

**[6–7] Where `2 sqrt3 log(sqrt3+sqrt2)` is born — the ARC half.** `A1`, `B1` carry only
`sqrt(linear)`, so `1 ± y = tau^2` rationalises them:

> `∫_0^{1/2}(A1+B1) dy = sqrt3 ∫_{1/sqrt2}^{sqrt(3/2)} d tau/(tau^4 − tau^2 + 1)`,
> `tau^4 − tau^2 + 1 = (tau^2 + sqrt3 tau + 1)(tau^2 − sqrt3 tau + 1)`  ← **the sqrt3**
> upper endpoint `tau = sqrt(3/2) = sqrt3/sqrt2`                       ← **the sqrt2**

Partial fractions over `Q(sqrt3)` then give, **PROVED, resid 0.0 at dps 60**:

> `∫_0^{1/2}(A1+B1)dy = (1/2)log(5+3sqrt2) − (1/2)log(sqrt3+sqrt2) − (1/4)log 7`
> `  + (sqrt3/2)[ arctan(sqrt6−sqrt3) + arctan(sqrt6+sqrt3) − arctan(sqrt2−sqrt3) − arctan(sqrt2+sqrt3) ]`

The coefficient of `log(sqrt3+sqrt2)` is exactly `−1/2`, and the prefactor `−(2/sqrt3)·6`
turns it into `+6/sqrt3 = 2 sqrt3`. **That is literally the `2 sqrt3 log(sqrt3+sqrt2)` of
the answer**, and it comes from the ARC `|zeta| = 1`, not from the rays.

So the provenance is now pinned:
* `sqrt3` = the discriminant of `1 + y + y^2`, the cyclotomic factor of `1 − y^3`
  (equivalently `cot(pi/6)`, the slope of the sector rays);
* `sqrt2` = the corner `y = 1/2`, i.e. `zeta = e^{i pi/6}`, entering as `tau = sqrt3/sqrt2`.

**[8] The RAY half.** `A2 = (9/2) y /((1-y^3) sqrt(1+8y^3))`; with `u = 2y` it is
`(9/2) u du /((8-u^3) v)` on the **equianharmonic** curve `v^2 = u^3 + 1` — a differential
of the **third kind** with poles over `u^3 = 8`. Over `u = 2` the two points are `(2, ±3)`,
and `(2,3)` is the classical **6-torsion** point of `y^2 = x^3+1`:

```
 1P=(2,3)  2P=(0,1)  3P=(-1,0)  4P=(0,-1)  5P=(2,-3)  6P=O      [verified by group law]
```

so the polar divisor `(2,3) − (2,−3) ~ 2[(2,3)]` is **3-torsion** — exactly the condition
for a third-kind abelian integral to collapse to a **logarithm of an algebraic function**.
STATUS: the torsion is **PROVED**; that the residual first-kind (`du/v`) component of this
particular differential vanishes is **inferred** from the verified total, not proved
independently, and I did not exhibit the explicit logarithm. PSLQ at 200 dps against
`{pi, arctan sqrt2, log(1+sqrt2), log(sqrt3+sqrt2), log2, log3, ∫_0^1 du/sqrt(1-u^3)} ×
{1, sqrt3}` (14 elements, non-degenerate control, maxcoeff 1e8, spurious scale 1e13) returns
**None** for the ray half alone — so its closed form needs algebraic arguments outside that
basis, and the "no equianharmonic period" conclusion is basis-limited, not settled.

---

## 2. Two red herrings, both REFUTED rigorously

### 2a. "The inner integrand is elementary" is neither necessary nor sufficient  [PROVED]

I first derived a **separation-of-variables normal form** valid for every k. From the
quadratic transformation (`c = a+b+1/2` holds for all k) plus Euler,
`S(k) = (2^{1+4/k}/pi) ∫_0^{sqrt2} ∫_0^1 D^{-4/k} d rho ds / sqrt(2-s^2)`,
`D = (1+s) + (1-s) rho^{k/2}`; the substitution `rho = P ((1+s)/(1-s))^{2/k}` (so
`D = (1+s)(1+P^{k/2})`) **separates the integrand for every k**:

> **`S(k) = (2^{1+4/k}/pi) ∫_0^{sqrt2} (1-s^2)^{-2/k}(2-s^2)^{-1/2}
>   [ ∫_0^{Lam_k(s)} (1+P^{k/2})^{-4/k} dP ] ds`,  `Lam_k(s) = ((1-s)/(1+s))^{2/k}`**

(PROVED; VERIFIED-NUMERIC for k = 2,4,6, resid ~2e-16 at dps 30 — quadrature-limited by the fractional-power endpoint singularities, not an identity defect). All k-dependence of the *region* is the
single curve `P = Lam_k(s)`.

Now apply **Chebyshev's binomial-differential theorem** to the separated inner integral
`∫ (1+P^{k/2})^{-4/k} dP` (`m=0, n=k/2, p=-4/k`): `(m+1)/n = 2/k ∈ Z` iff `k ∈ {1,2}`;
`(m+1)/n + p = −2/k ∈ Z` iff `k ∈ {1,2}`; `p ∈ Z` iff `k ∈ {1,2,4}`. Identically, in the
polar chart, `∫_0^1 r dr/sqrt(1+sigma r^k)` (`m=1,n=k,p=-1/2`) is elementary iff
`k ∈ {1,2,4}`. Therefore:

> * **k = 3: BOTH iterated integrals are non-elementary, yet S(3) IS elementary.**
> * **k = 4: the inner integral IS elementary, yet S(4) is NOT.**

This is a theorem (Chebyshev), not a heuristic. It **refutes** any mechanism story of the
form "S(k) is elementary because the inner 2F1 degenerates", in *both* directions. In
particular it refutes the general claim of the prior repo reflection
`07-reflections/the-binomial-product-series-Sk-generating-function-and-general-closed-form-opus-S4.md`
(see §4) and shows lane A1's §2.2 "this formula mechanically reproduces the known cases"
narrative is a coincidence of k ∈ {1,2,4}, not the mechanism.

### 2b. The Schwarz angle-sum criterion, as used in lane A1 §2.5, is not valid  [PROVED]

The Schwarz triple is only defined **up to `lam -> ±lam + integer` with an EVEN total
shift** (these moves do not change the projective monodromy); membership in Schwarz's list
must be tested on the *reduced* triple. Doing that:

| k | raw triple of `2F1(1,4/k;1+2/k;·)` | raw sum | reduced triple | reduced sum |
|---|---|---|---|---|
| 1 | (2, 3, 2) | 7 | (0,0,1) | 1 (degenerate) |
| 2 | (1, 1, 1) | 3 | (0,0,1) | 1 (degenerate) |
| 3 | (2/3, 1/3, 2/3) | 5/3 | (1/3,1/3,1/3) | **1 — EUCLIDEAN** |
| 4 | (1/2, 0, 1/2) | 1 | (0,1/2,1/2) | 1 |
| 5 | (2/5,1/5,2/5) | 1 | (1/5,2/5,2/5) | 1 |
| … | … | | | 1 for every k |

The reduced sum is **exactly 1 for every k**, including k = 1,2,3. The same holds for the
`2F1(1/2,2/k;1+2/k;·)` of preamble (v). So the "spherical for k ≤ 3 / Euclidean for k ≥ 4"
separation of lane A1 §2.5 is an artefact of not reducing the triple. **The verdict
(k ≤ 3 vs k ≥ 4) happens to be right, but the stated criterion does not produce it.** The
criterion that does is §0.

---

## 3. New exact statements collected on the way (all PROVED + verified)

**A3-1 (the exact-derivative property, and why k = 1 is free for the whole family).**
For `y = 2F1(a,b;c;z)`, the ODE `z(1-z)y'' + (c-(a+b+1)z)y' - ab y = 0` is an *exact*
derivative `d/dz[z(1-z)y'] = ab y` **iff `a+b = c = 1`** — exactly the signature family
`2F1(lam, 1-lam; 1; z)`. Hence, with `f = 2F1(1/4,3/4;1;z)`,
`d/dz[z(1-z) f'] = (3/16) f` (verified to 5e-51), and

> `S(1) = ∫_0^1 f = (16/3) lim_{z->1}(1-z)f'(z) = 16/(3 pi sqrt2) = 8 sqrt2/(3 pi)` — one line;
> more generally `S_lam(1) = sin(pi lam)/(pi lam (1-lam))` (verified for lam = 1/2,1/3,1/4,1/6,2/5).

**A3-2 (one clean IBP).** `S(k) = 8 sqrt2/(3 pi k) + (1 − 1/k) T(k)`,
`T(k) = ∫_0^1 2F1(1/4,3/4;2;t^k)dt`, using `(1-z)f'(z) = (3/16) 2F1(1/4,3/4;2;z)`.
Verified for k = 1,2,3,4,5,7 to 1e-41. (This is lane A2's recursion A2-4 in the form that
makes the `k=1` degeneration transparent.)

**A3-3 (odd-k isotrivial rescaling).** Because the fibre curves `w^2 = 1 + sigma r^k` are
all isomorphic under `r -> sigma^{-1/k} r`,

> **`S(k) = (2k/pi) ∫_{-1}^{1} m^{k-3} J_k(m) dm / sqrt(1 - m^{2k})`, `J_k(m) = ∫_0^m u du/sqrt(1+u^k)`**, k odd.

Verified for k = 1,3,5,7. **k = 3 is the unique k for which the weight `m^{k-3}` is
trivial** — a striking marker for the special case, though §0 shows it is not the cause.

**A3-4 (separation of variables for all k)** — §2a above.

---

## 4. (b) Prior work in the repo — what is genuinely reusable

There *is* substantial prior work on this exact problem. In order of relevance:

| file | verdict |
|---|---|
| `07-reflections/central-binomial-3F2-quarter-series-S-of-k-kps-S146.md` | **The origin of the problem.** Establishes (i)–(iii) of the preamble, gives `S(k) = (2/pi)∫_0^{pi/2} g_k(cos^4 phi) dphi`, `g_k(y)=∫_0^1 ds/sqrt(1-y s^k)`, decodes `arctan(sqrt2/5) = pi − 3 arctan sqrt2` via `(1+i sqrt2)^3 = −5+i sqrt2`, and proves S(4) contains `varpi`. **Reusable and correct.** Its k=3 remark ("connects to `w=1+i sqrt2`, special to k=3") is a numerological observation, not the mechanism. |
| `07-reflections/the-binomial-product-series-Sk-generating-function-and-general-closed-form-opus-S4.md` | **Contains an over-claim that should be flagged.** Its Wallis route `S(k) = (2/pi)^2 ∫∫ G_k(sin^2 t sin^4 u) dt du`, `G_k(c)=∫_0^1 dx/(1-cx^k)`, is correct and useful (I re-verified the equivalent `S(k)=(2/pi)∫∫ dx du/sqrt(1-x^k sin^4 u)`). But its conclusion — "*each S(k) has an elementary `(1/pi)*(logs+arctans of algebraics)` form by the `G_k` partial-fraction structure*" — is **REFUTED**: `G_k` elementary does not make the 2-D integral elementary, and §2a shows the analogous inference fails in both directions. This directly contradicts kps-S146 (same repo) and lanes A1/A2. Recommend a correction note. |
| `07-reflections/signature-series-family-CM-and-the-lemniscate-bridge-to-LRC-kps-S147.md` | Generalises to `S_lam(k)`, `lam ∈ {1/2,1/3,1/4,1/6}`, and proves `S_lam(1) = sin(pi lam)/(pi lam(1-lam))` (which I re-derive in one line as A3-1). **Reusable.** Its "elementary reduction is a CM/singular-value phenomenon" explanation is a heuristic that §0 replaces; and my k=4 inner curve `w^2 = 4y x^3 − 4y^3 x + 1` has `j(y) = 110592 y^8/(64y^8 − 27)`, i.e. is **non-isotrivial**, so the k=4 obstruction is not a single CM curve at that stage (the CM curves that A1/A2 find appear later in their routes — not a contradiction, but the "CM causes it" story is not supported at the level where the failure actually happens). |
| `01-canon/theorems/THM-2000` (support-harmonic Abel–Dini / figurate mass surface) and `THM-2005` | **Not reusable.** They concern reciprocal sums of *integer* sequences (`∑ 1/a_n`), Abel–Dini/Kakeya tail criteria, and digamma evaluation of `∑ 1/(kn+r)`. The shared idiom is only "arithmetic-progression denominator ⇒ k-th-root partial fractions", which is precisely the step §2a shows to be insufficient here. No theorem transfers. |
| `04-computation/lean/TournamentH7/TournamentH7/ArtanhSandwich.lean` | **Reusable, narrowly but genuinely.** It kernel-proves `2(t+t^3/3+t^5/5) ≤ log((1+t)/(1-t)) ≤ 2(t+t^3/3+t^5/(5(1-t^2)))` for `0 ≤ t < 1` in Lean/Mathlib. Both logs in this problem are `artanh` values at algebraic points — `log(1+sqrt2) = artanh(1/sqrt2)` and `log(sqrt3+sqrt2) = artanh(sqrt(2/3))` — so the sandwich yields *certified* rational enclosures of `pi S(2)` and of the log part of `pi S(3)`, if a formal verification of those closed forms is ever wanted. That is the only genuinely transferable asset I found. |
| the "Abel–Dini", "artanh", "3F2" grep hits elsewhere (`HYP-9023`, `THM-2143`, `THM-1985/1990`, the AMM-12592 body) | **Not reusable.** They are certified-inequality / reciprocal-series objects; the `3F2` hits are incidental string matches in unrelated LRC/tournament theorems. |

**Honest summary of (b):** the repo already contains this problem twice (S146, opus-S4) and a
generalisation (S147); it does **not** contain the mechanism, and one of the three
reflections states a general conclusion that is now refuted.

---

## 5. (c) Is `3F2(1/4,3/4,c; 1, c+1; 1)` a known evaluable family?

### 5a. The Thomae orbit, computed independently  [PROVED computationally]

Generating the full `S_5` orbit (basic Thomae step + all 12 trivial permutations, BFS to
closure) of `3F2(1/4,3/4,c;1,1+c;1)` gives **7 distinct members** (the orbit collapses from
the generic 10 because of the repeated `1, 1+c` pattern):

```
 1. 3F2(1-c, 1/4, 1/4;      1,     5/4;   1)
 2. 3F2(1-c, 3/4, 3/4;      1,     7/4;   1)
 3. 3F2(1-c, 1,   1;        5/4,   7/4;   1)      <-- lane A2's normal form A2-1
 4. 3F2(1/4, c+1/4, 1;      5/4,   c+1;   1)
 5. 3F2(3/4, c+3/4, 1;      7/4,   c+1;   1)
 6. 3F2(1/4, 3/4, c;        1,     c+1;   1)      <-- the original
 7. 3F2(c+1/4, c+3/4, c;    c+1,   c+1;   1)
```
(the basic Thomae step itself checked numerically at c = 1/3 and 1/5, resid 4.6e-41).

Solving, over the **whole** orbit, the parameter equations for every classical summation:

| summation | values of c that make it apply, anywhere on the orbit |
|---|---|
| Gauss degeneration (a numerator = a denominator) | `c ∈ {1, 0, −1/4, −3/4}` |
| terminating (a numerator a non-positive integer) | `c ∈ {1,2,3,4,5} ∪ {0,−1,−2,…} ∪ {−1/4,−3/4,−5/4,…}` |
| Watson | `c = 1` only |
| Whipple | `c = 1` only |
| Dixon | `c = 1` only |

`S(k)` sits at `c = 1/k`. **The only positive integer k reachable by any classical
summation anywhere on the 120-element Thomae group orbit is k = 1.** This is an
independent confirmation of lane A2 §2a (which reported 5 of these 7 members).

### 5b. The nearest family that IS evaluable  [PROVED]

The **reflected** Mellin transform is a Gamma quotient, and the proof is one line of Gauss:

```
 ∫_0^1 x^{c-1} 2F1(1/4,3/4;1; 1-x) dx = (1/c) 3F2(1/4,3/4,1; 1, 1+c; 1)
                                      = (1/c) 2F1(1/4,3/4; 1+c; 1)      [the 1 cancels]
                                      = Gamma(c)^2 / (Gamma(c+1/4) Gamma(c+3/4)).
```
Verified numerically (exact at c = 1, 7/5; the c = 1/2, 1/3 checks are quadrature-limited
by the `x^{c-1} log x` singularity, resid 3e-21 / 1e-13, and the identity is proved anyway).

In Legendre language: `2F1(1/4,3/4;1;z) = P_{-1/4}(1-2z)` (verified), and the classical
tabulated formula is `∫_0^1 x^{s-1} P_nu(2x−1) dx = Gamma(s)^2/(Gamma(s-nu)Gamma(s+nu+1))` —
**argument `2x−1`, i.e. the reflected one.** `S(k)` is the same Mellin transform with
argument `1−2x`, for which no classical evaluation exists.

> **Answer to (c): no.** `3F2(1/4,3/4,c;1,c+1;1)` is *not* a known evaluable family; the
> free parameter sits in the pair `(c ; 1+c)` instead of the pair `(1 ; 1+c)`, and it is
> exactly that one displacement that destroys Gauss. The Thomae orbit does not repair it —
> no image is a summable Gauss/Watson/Whipple/Dixon/Saalschütz case for general c, and
> **Question 2 is not settled by Thomae in one stroke.**

---

## 6. Claims, with status

| # | claim | status |
|---|---|---|
| A3-M1 | `deg_x(1 + Im(x+iy)^k) = k−1`, so the x-first inner integral of the sector period is elementary iff k ≤ 3, and is an incomplete elliptic F of the first kind at k = 4 | **PROVED** |
| A3-M2 | full explicit k=3 chain: sector → arcsinh/arcsin → IBP → purely algebraic `−(2/sqrt3)∫_0^{1/2}[(A1−A2)+(B1−B2)]dy` | **PROVED** (steps) + **VERIFIED-NUMERIC (~32 digits)** |
| A3-M3 | arc half `= (1/2)log(5+3√2) − (1/2)log(√3+√2) − (1/4)log7 + (√3/2)[atan(√6−√3)+atan(√6+√3)−atan(√2−√3)−atan(√2+√3)]`; its `−1/2` is the source of `2√3 log(√3+√2)` | **PROVED**, resid 0.0 at dps 60 |
| A3-M4 | ray half = third-kind differential on `v^2=u^3+1` with 3-torsion polar divisor (`(2,3)` is 6-torsion) | torsion **PROVED**; "hence elementary" **CONJECTURAL** (first-kind component not independently checked; explicit log not exhibited) |
| A3-R1 | "inner integral elementary" is neither necessary (k=3) nor sufficient (k=4): the separated / polar inner integral is elementary exactly for k ∈ {1,2,4} | **PROVED** (Chebyshev) |
| A3-R2 | the raw Schwarz angle sum is not the Schwarz criterion; the reduced triple has sum 1 for every k, so lane A1 §2.5's criterion is invalid (its verdict is nonetheless right) | **PROVED** |
| A3-N1 | separation normal form `S(k) = (2^{1+4/k}/pi)∫(1-s^2)^{-2/k}(2-s^2)^{-1/2}[∫_0^{Lam_k}(1+P^{k/2})^{-4/k}dP]ds` | **PROVED** + VERIFIED-NUMERIC |
| A3-N2 | `d/dz[z(1-z)y'] = ab y` iff `a+b=c=1`; gives `S_lam(1)=sin(pi lam)/(pi lam(1-lam))` and `S(k)=8√2/(3πk)+(1−1/k)T(k)` | **PROVED** + VERIFIED-NUMERIC |
| A3-N4 | `x = yX` gives `f_k = 1 + y^k Q_k(X)`, `Q_k = Im((X+i)^k)` of degree `k−1` with roots `cot(j pi/k)`: a one-parameter hyperelliptic pencil of genus `floor((k−2)/2)` | **PROVED** + VERIFIED-NUMERIC (k=3..7) |
| A3-N3 | odd-k isotrivial rescaling `S(k)=(2k/pi)∫_{-1}^1 m^{k-3}J_k(m)dm/sqrt(1-m^{2k})`; k=3 unique with trivial weight | **PROVED** + VERIFIED-NUMERIC |
| A3-C1 | Thomae orbit has 7 distinct members; only c = 1 is classically summable anywhere on it | **PROVED** (computational, exact in sympy) |
| A3-C2 | the reflected family `3F2(1/4,3/4,1;1,1+c;1) = Gamma(c)^2/(Gamma(c+1/4)Gamma(c+3/4))` by Gauss; S(k) is its unrepairable neighbour | **PROVED** |
| A3-P1 | prior repo reflection opus-S4's general claim ("every S(k) is elementary by the `G_k` partial fraction") is refuted | **REFUTED** (by A3-R1, and by kps-S146 / lanes A1,A2) |

---

## 7. What this hands the other lanes

1. **The criterion to quote is A3-M1**, not Schwarz and not CM. It is elementary, sharp, and
   independently corroborated by lane A1's incomplete-`F` finding at k = 4.
2. **A one-parameter normal form for the fibre curve (A3-N4, PROVED, trivial but useful).**
   Because `y` is real, `x = yX` gives

   > `f_k(x,y) = 1 + y^k · Im((X+i)^k)`,  so  `∫ dx/sqrt(f_k) = y ∫ dX / sqrt(1 + c Q_k(X))`,
   > `c = y^k`,  `Q_k(X) = Im((X+i)^k)` — a degree `k−1` polynomial with the `k−1` distinct
   > real roots `cot(j pi/k)`, `j = 1..k−1`.

   So the whole family collapses to the **one-parameter** hyperelliptic pencil
   `w^2 = 1 + c Q_k(X)`, genus `floor((k−2)/2)`:
   `Q_3 = 3X^2−1` (conic — rational, hence k = 3 elementary);
   `Q_4 = 4X(X^2−1)` (elliptic pencil, non-isotrivial, `j(y)=110592y^8/(64y^8−27)`);
   `Q_5 = 5X^4−10X^2+1` (genus 1); `Q_6 = 2X(3X^2−1)(X^2−3)` (genus 2);
   `Q_7` (genus 2).
   This turns lane A2's open question 2 into a concrete one: **is the genus-2 pencil
   `w^2 = 1 + 2c X(3X^2−1)(X^2−3)` of CM type?** Its root set `{0, ±1/sqrt3, ±sqrt3, ∞}`
   is stable under `X -> −1/X` (checked), which is exactly the kind of extra involution that
   splits a genus-2 Jacobian — a promising, cheap computation for whoever picks S(6) up.
3. **A3-M4 is the one soft joint** in my chain. Making it hard means exhibiting the explicit
   algebraic `g` with `div(g) = 3(2,3) − 3(2,−3)` on `v^2 = u^3+1` and checking
   `(9/2) u du/((8-u^3)v) = c · dlog g` with no `du/v` residue. That is a finite,
   mechanical computation (Hermite reduction on a genus-1 curve) and would upgrade the
   whole k=3 evaluation to a from-scratch proof independent of the quoted closed form.
4. **A3-R1 is a warning to keep**: any future "S(k) is elementary because the inner function
   degenerates" argument must be checked against k = 4, where the inner function *does*
   degenerate and the answer is still transcendental.
