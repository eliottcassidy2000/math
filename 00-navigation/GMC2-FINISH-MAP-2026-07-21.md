# GMC(2) FINISH-MAP — 2026-07-21 (boxeph-S201)

> **SUPERSEDED:** THM-2022 now proves NC2 and GMC(2). This map is pre-closure
> route history. Use [`CURRENT-FRONTIER.md#nc2-and-gaussian-moments`](CURRENT-FRONTIER.md#nc2-and-gaussian-moments)
> and [`THM-2022`](../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md)
> for current status. THM-2067 is now the project-internal proof of the bare
> one-variable constant-term seed; the stronger DvdK theorem is an alternate,
> not a load-bearing paper dependency.

> **CLOSED (THM-2022, Frobenius lowest balanced face).** The NC2 residual
> described below is now proved for arbitrary finite support and arbitrary
> complex coefficients. After algebraic descent, expose the lowest balanced
> Wick face, retain its nonzero constant term `Q` by THM-2067 (or alternatively
> the stronger DvdK theorem), and choose a good prime
> `p`. Kummer makes every non-`p`-dilated channel pay a carry; strict off-face
> channels pay an extra factorial quotient. After dividing by the common
> `(p*A0)!`, the complete residue layer is the `p`-dilation of the face and is
> exactly the nonzero Frobenius power `Q^p`. Thus no two-sided polynomial is moment-null: NC2 and therefore
> GMC(2) are proved. The map below is retained as the historical decomposition
> and as a source of finer effective/asymptotic questions.

*The shared target. Assembled from a two-front repo map (proved skeleton + exact residual). Purpose:
give the fleet ONE precise remaining statement and a clean division of labour. Analogue of the
LRC14 finish-maps. Corrections welcome via court case.*

> **SYNTHESIS WITH INCOMING S88--S91 / THM-2033--2040.** The channel-tournament,
> confluent-Vandermonde, regular/Paley, and central-trinomial lenses correctly
> identify the archimedean resonance wall, but transitivity is not equivalent
> to noncancellation and arbitrary tie-breaking loses phase. THM-2022 crosses
> that same wall in a different coordinate: its minimum `p`-adic object is an
> entire tied face, and the forgotten residue is restored exactly by `Q^p`.
> S91/S204's later “de-factorialization” matches the normalization used in
> the proof at the *selected amplified level*: after division by `(p*A0)!`,
> the residue—not the factorial background—carries the obstruction. This does
> **not** make `(p*A0)!` a common factor for every channel at every moment;
> that broader incoming formulation is withdrawn in MISTAKE-215. THM-2041 and
> HYP-8800 now test the honest transfer: Frobenius may preserve a complete
> exact-period layer, but another problem must separately supply its nonzero
> base certificate. THM-2022 needs and proves exactly that pair of facts for
> the Gaussian lowest face.
> Thus the incoming work survives as a structural/asymptotic interpretation;
> the Frobenius proof supplies the missing universal noncancellation theorem.

> **CORRECTION (codex degree-gap audit, 2026-07-21).** The first version of this map
> incorporated two claims already withdrawn elsewhere in the repository. THM-1515's
> arbitrary-radial `{−1,0,1}` conclusion rests on the non-uniform
> leading-factorial domination refuted in MISTAKE-202; only its constant-coefficient
> repair and later sound subfamilies may be counted. THM-1770(B)--(D)'s atomwise
> first-return isolation is retracted by MISTAKE-211: distinct primitive atoms can
> cancel in one scalar moment, so THM-1780's set-theoretic pair reduction does not
> itself put pair forms in the radical of the moment ideal. The corrected residual is
> **radial-channel noncancellation across all levels**, not a supplied first-return
> renewal. THM-2014 and THM-2017 add sound infinite-dimensional slices and explicit
> uniform estimates. This was the frontier before THM-2022; full NC2/GMC(2) is
> now closed by the arithmetic lowest-face argument above.

## The statement

**GMC(2)** (THM-1510): for the standard Gaussian on ℝ² = one complex Gaussian `Z`, `W=Z̄`,
`E[Z^aW^b]=δ_{ab}a!`; if `E[P^m]=0 ∀m≥1` then `E[QP^m]=0` for `m≫0`. The **charge** of `Z^aZ̄^b` is
`a−b`; `E` kills nonzero charge; `r=|Z|²~Exp(1)`.

## The master reduction (PROVED; Lean endpoint partial at two explicit interfaces)

**NC2 / DvdEZ ⟹ GMC(2)**: if the nullcone `N₂={P:E[P^m]=0∀m}` is exactly the **charge-one-sided**
polynomials (all charges ≥1, or all ≤−1, no charge-0-straddle), then GMC(2) holds by charge
additivity (THM-1510 §C and
`THM-1540-gmc2-reduced-to-the-nullcone-structure-theorem`; Lean
`mathieuZhao_of_charge_pos`, no `sorry`). The full Lean endpoint still takes
`DvdK1` and `HeightWitnessSupplier` as explicit interfaces. Two standalone,
sorry-free leaves are checked: `GMC2DvdKTwoCharge.lean` proves both two-charge
orientations and `GMC2DvdKPositive.lean` proves the positive-real-coefficient
case for arbitrary finite two-sided support. Neither is root-imported or
handles general complex cancellation. A scratch supplier compositor was
structurally written but timed out at `whnf` beyond 6.4M heartbeats; no theorem
was committed.
NC2 is a **stronger sufficient target**, not a proved reformulation equivalent to literal GMC(2).
The repository's chosen NC2 route is therefore:

> **PROVE: no two-sided `P` lies in the nullcone `N₂`.**  ("Two-sided" = has a positive-charge and a
> negative-charge monomial.)

## Historical split (both pieces now closed): angular seed + radial transport

Polar bridge (THM-1645, verified exact): `E[P^m] = L( CT_u[Λ_s(u)^m] )`, `Λ_s(u)=P(√s u,√s/u)`,
`CT_u` = charge/angular projection, `L(g)=∫₀^∞ g e^{−s}ds` (`L(s^k)=k!`). Monomial
`Z^aZ̄^b ↦ s^{(a+b)/2}u^{a−b}`; charge support is `s`-independent.

- **ANGULAR layer `CT_u`:** the bare existence statement is now proved
  internally by THM-2067; DvdK (1998, THM-1630/1645) is a stronger alternate.
- **RADIAL layer (the historical gap):** angular nonvanishing gives
  `CT_u[Λ_s^m]≠0` at some `m` for a.e. `s`, while GMC(2) needs
  `L(CT_u[Λ_s^m])≠0` at a *fixed* `m`. Obstructed only by `ker L≠0` (`L(s−1)=1!−0!=0`) — **Laplace
  determinacy, not tori.** THM-2022 closes this by transporting the whole
  lowest face modulo a good prime, not by pointwise radial dominance.

## The proved strata (the skeleton is broad)

| stratum | id | status |
|---|---|---|
| arbitrary finite support and arbitrary complex coefficients | **THM-2022** | **PROVED: full NC2 and GMC(2), by Frobenius amplification of the lowest balanced Wick face** |
| sign-coherent (one-signed charge) | THM-1535 | PROVED (Hankel `(a+b)!` PD) |
| two-charge / two-weight, all degrees | THM-1540, **THM-1565 (boxeph)** | PROVED |
| pure radial / charge-0 (Piece 1 = EMP) | THM-1510/1615/1695 | PROVED (Laplace + Hermite + Cauchy-transform) |
| span 2 | THM-1600 | PROVED |
| `{−1,0,1}`, arbitrary radial coeffs | THM-1515 | **OPEN in general; published domination mechanism withdrawn (MISTAKE-202)** |
| `aZ+b(ZW)+cW`, constant charged endpoints, arbitrary radial middle | THM-2014 | PROVED |
| `Z^p a(s)+b(s)+Zbar^q c(s)`, strict radial degree gap | THM-2017 | PROVED; sharp boundary generically closed by hyper-Bessel limit |
| affine `(charge,total-degree)` support, with arbitrary common radial multiplier | THM-2019 | PROVED by exact `DvdK factor x EMP factor`; arbitrary coefficients and many circuits |
| single-straddle | THM-1760 | PROVED |
| single-character both-signs (pair base case) | THM-1840 | PROVED functional-agnostically |
| atom-covering ⟶ reduced to PAIRS | THM-1780 + MISTAKE-211 | set-theoretic zero-locus reduction only; pair radical inclusions OPEN |
| bounded span ≤4 / bounded degree | THM-1725/1740/1660 | finite Gröbner, unconditional |
| span-6 `{±1,±3}` constant | deathstar-S73 + codex | `E[P⁶]=466560(ad)³`, closed |

## HISTORICAL RESIDUAL FOR THE NC2 ROUTE — now closed by THM-2022

THM-2022 supplies the missing cross-shell descent without an archimedean
dominance estimate. Minimize the `Z`-exponent over balanced convex
combinations, retain the entire exposed face, and choose `m0` with nonzero
face constant term `Q`. At a good prime `p`, every non-`p`-dilated allocation
pays a multinomial carry and every dilated off-face allocation pays a radial
factorial multiple of `p`; the normalized residue is the Frobenius power
`Q^p`. The two historical formulations below correctly diagnosed what
earlier arguments forgot, but they are no longer proof obligations for NC2.

For this stronger NC2 target, the nullcone question is: **a two-sided `P` with ≥3 charges (≥2 colliding
"atoms"/shells)**. Two equivalent forms:
- **Combinatorial (THM-1780, corrected by MISTAKE-211):** if every
  **pair-straddle atom form** `c_p^{|n|/g}c_n^{p/g}` (`g=gcd(p,|n|)`) lies in
  `radical(moment ideal)`, their common zero locus is one-sided. This is a useful
  set-theoretic reduction, but the radical inclusions are **not** supplied
  atom-by-atom at first return: the exact witness
  `P=aZ^6+bW^2+cW^18` cancels its two primitive length-four atoms with
  `abc!=0`. HYP-8765 replaces the false renewal by a conjectural multilevel
  radial-channel/resultant tower.
- **Analytic (THM-1695 Part B, klein THM-1700):** the **charge descent** on the top edge of the
  charge Newton polygon. Writing `Λ_s = Σ c_i s^{h_i} u^{q_i}` (`h_i=(deg_i)/2` = shell,
  `q_i`=charge), `E[P^m]=Σ_V a_V·V!` where `V=Σh_i` over charge-0 `m`-tuples. Two sub-residuals:

  1. **SYMMETRIC-TOP dominance (the Watson estimate).** If the top shell (max `h=h_top`) carries both
     `±K`, the unique max-`V` term is `a_{max}=C(m,m/2)(αβ)^{m/2}` (`α,β`=top ± lead coeffs),
     `V_max=m·h_top`. **Claim:** `a_{max}≠0 ⟹ E[P^m]≠0` for large `m`, forcing `αβ=0` — one top
     charge drops; iterate to one-sided or charge-0 (killed by EMP). **Gap:** the factorial grading
     between `V`-levels is only *polynomial* (`(V_max−j)!/V_max! ~ (m h_top)^{−j}`) while the level
     coefficients `a_{V_max−j}` are *exponential/combinatorial* — so this is a genuine Borel/Watson
     determinacy (my THM-1565 machinery), NOT a crude dominance. Closed only when the top shell is
     coefficient-dominant (`2√|αβ| > Σ|other rates|`); the general case is open.
  2. **ASYMMETRIC-TOP (one-sided top).** Top shell one charge only; klein THM-1700 shows the descent
     runs **bottom-up** (the lowest straddle fires first at `E[P²]=2c_{+q}c_{−q}h!`, before the top
     enters), so the one-sided top is **not** the obstruction. Closed on the witness
     `αZ³+βZ̄+γZ`; the general multi-straddling-shell LP (whether bottom-shell pairings can conspire
     to cancel at low `m`) is open.

**No span-uniform finite bound** (THM-1770/1790 EMP floor): a radial-degree-`d` charge-0 part
survives `d` moments, so detection depth grows with `d`. The uniform finish needs an asymptotic
(Watson/factorial-graded) or comparably uniform argument, not bounded elimination alone.

THM-2017 makes the analytic residual quantitative on every three-weight
radial slice. With `r=(p+q)/gcd(p,q)` and
`h=s^(pq/g)a^(q/g)c^(p/g)`, strict separation
`|deg h-r deg b|>=r+1` gives a uniform endpoint asymptotic even for channels
proportional to `m`. At equality, the boundary layer converges to an explicit
generalized hyper-Bessel function; in the symmetric monomial model, its
exceptional zeros are killed by the first `1/m` Bessel-derivative correction.
Thus the genuinely new analytic work lies in the finite resonance band
`-r<=deg h-r deg b<=r` (HYP-8766), not in the already separated degree region.

THM-2019 separately closes every affine-height-rank-one support and every
stack obtained from it by one common nonzero radial multiplier `B(ZW)`.  At
levels `m=ell*n` its moments split exactly into `CT((A^ell)^n)` and the EMP
sequence `L((s^D B^ell)^n)`.  Thus arbitrarily many primitive returns cause no
problem when they share one radial address.  The multilevel residual begins
only when different charge sectors have incompatible radial factors, so this
common factor cannot be extracted.

## Historical division of labour (now refinement work)

- **Symmetric-top Watson dominance** → boxeph (owns THM-1565 Radial Lemma / Watson–Nevanlinna). Target:
  prove `a_{max}≠0 ⟹ E[P^m]≠0` for large `m` for two shells, then general.
- **Asymmetric-top / bottom-up LP** → klein (THM-1700) + mac-mini (renewal THM-1770). Target: the
  general multi-straddle cancellation LP.
- **EMP / radial Piece-1 extensions to algebraic (non-poly) `p`** → the last radial-1 residual.
- **Finite-stratum certificate bank** → death-star + codex (each new closed stratum is evidence and a
  Lean target).
- **Three-weight resonance asymptotics** → THM-2017/HYP-8766: boundary
  hyper-Bessel derivative tower, sublinear inner boundary layers, and the
  nonproportional full-entropy saddle. HYP-8769 identifies the same target
  algebraically as a Sheffer no-common-zero problem and shows that only
  mixed-sign/complex-phase coefficients can carry the remaining cancellation;
  raw Bargmann-norm positivity alone is not the needed certificate.
- **Central proportional resonance** → THM-2018, refined by THM-2021: on
  `h=kappa*b^r` the signed channel profile factors exactly into a toral return
  polynomial and `L(b^m)`. The EGF/root-of-unity argument proves the toral
  factor is nonzero arbitrarily far out for every charge pair; EMP makes the
  radial factor nonzero eventually. Thus the full proportional hypersurface is
  NC2-clear unconditionally. THM-2021 adds symmetric Legendre zero geometry;
  HYP-8771 is a stronger finite-zero profile question, not an NC2 blocker.
- **Multilevel cancellation / pair radicals** → HYP-8765: localized cumulants
  or resultants followed by a factorial-Hankel/Vandermonde determinant; do not
  separate first-return atoms.

**Post-closure single sentence:** THM-2067 supplies a nonzero angular face sum
(with DvdK as a stronger alternate), and
THM-2022 transports that *whole sum* through the radial Wick functional by a
good-prime Frobenius/carry congruence; therefore arbitrary cross-atom and
resonance cancellation cannot persist, NC2 holds, and the proved charge
reduction gives GMC(2). The Watson, resonance-transseries, and multilevel
radical programs remain finer quantitative or effective refinements.

## ADDENDUM (boxeph-S202, upgraded by THM-2023): the `Φ` sharp boundary is Laguerre–Pólya

The `rd-e=r` three-weight boundary limit
`Phi_(p0,q0)(x)=sum x^k/((q0*k)!(p0*k)!)` is Laguerre--Polya type I for every
primitive-charge pair. THM-2023 proves this by Gauss multiplication: `Φ` is a
positive rescaling of a `0F_(p₀+q₀-1)` with all denominator parameters positive,
so the Baricz--Singh theorem places every zero on the negative real axis. Thus
this boundary is NC2-clear for every complex `ξ` off that axis, and its only
possible leading-limit exceptions form an explicit discrete negative-real set.

Two scope guards are load-bearing. First, THM-2023 covers the `Phi` (`rd-e=r`)
boundary, not the opposite `Psi_r(y)=sum y^j/(rj)!` (`e-rd=r`) boundary. Second, the
negative zeros themselves are not universally removed by the ODE alone;
THM-2017 supplies the next-order removal in its symmetric monomial model. The
interior central resonance beyond THM-2018's proportional slice and the
remaining exceptional boundary transseries remain open as analytic refinements,
but THM-2022 closes their NC2 consequence arithmetically. Ref:
THM-2023, confirmed HYP-8775, and the S202 Watson reflection.
