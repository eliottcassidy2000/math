# GMC(2) FINISH-MAP — 2026-07-21 (boxeph-S201)

*The shared target. Assembled from a two-front repo map (proved skeleton + exact residual). Purpose:
give the fleet ONE precise remaining statement and a clean division of labour. Analogue of the
LRC14 finish-maps. Corrections welcome via court case.*

## The statement

**GMC(2)** (THM-1510): for the standard Gaussian on ℝ² = one complex Gaussian `Z`, `W=Z̄`,
`E[Z^aW^b]=δ_{ab}a!`; if `E[P^m]=0 ∀m≥1` then `E[QP^m]=0` for `m≫0`. The **charge** of `Z^aZ̄^b` is
`a−b`; `E` kills nonzero charge; `r=|Z|²~Exp(1)`.

## The master reduction (PROVED, Lean-formalized)

**NC2 / DvdEZ ⟹ GMC(2)**: if the nullcone `N₂={P:E[P^m]=0∀m}` is exactly the **charge-one-sided**
polynomials (all charges ≥1, or all ≤−1, no charge-0-straddle), then GMC(2) holds by charge
additivity (THM-1510 §C, THM-1535, THM-1830; Lean `mathieuZhao_of_charge_pos`, no `sorry`). So the
entire problem is:

> **PROVE: no two-sided `P` lies in the nullcone `N₂`.**  ("Two-sided" = has a positive-charge and a
> negative-charge monomial.)

## The split (PROVED): angular = DvdK, gap = radial

Polar bridge (THM-1645, verified exact): `E[P^m] = L( CT_u[Λ_s(u)^m] )`, `Λ_s(u)=P(√s u,√s/u)`,
`CT_u` = charge/angular projection, `L(g)=∫₀^∞ g e^{−s}ds` (`L(s^k)=k!`). Monomial
`Z^aZ̄^b ↦ s^{(a+b)/2}u^{a−b}`; charge support is `s`-independent.

- **ANGULAR layer `CT_u` = the Duistermaat–van der Kallen theorem** (1998), applied uniformly in `s`
  (THM-1630/1645). **CLOSED.**
- **RADIAL layer** = the whole gap: DvdK gives `CT_u[Λ_s^m]≠0` at some `m` for a.e. `s`; GMC(2) needs
  `L(CT_u[Λ_s^m])≠0` at a *fixed* `m`. Obstructed only by `ker L≠0` (`L(s−1)=1!−0!=0`) — **Laplace
  determinacy, not tori.**

## The proved strata (the skeleton is broad)

| stratum | id | status |
|---|---|---|
| sign-coherent (one-signed charge) | THM-1535 | PROVED (Hankel `(a+b)!` PD) |
| two-charge / two-weight, all degrees | THM-1540, **THM-1565 (boxeph)** | PROVED |
| pure radial / charge-0 (Piece 1 = EMP) | THM-1510/1615/1695 | PROVED (Laplace + Hermite + Cauchy-transform) |
| span 2 | THM-1600 | PROVED |
| `{−1,0,1}`, arbitrary radial coeffs | THM-1515 | PROVED |
| single-straddle | THM-1760 | PROVED |
| single-character both-signs (pair base case) | THM-1840 | PROVED functional-agnostically |
| atom-covering ⟶ reduced to PAIRS | THM-1780 | multi-charge atoms redundant |
| bounded span ≤4 / bounded degree | THM-1725/1740/1660 | finite Gröbner, unconditional |
| span-6 `{±1,±3}` constant | deathstar-S73 + codex | `E[P⁶]=466560(ad)³`, closed |

## THE ONE RESIDUAL — the cross-shell descent

After the above, the nullcone question is: **a two-sided `P` with ≥3 charges (≥2 colliding
"atoms"/shells)**. Two equivalent forms:
- **Combinatorial (THM-1780):** every **pair-straddle atom form** `c_p^{|n|/g}c_n^{p/g}` (`g=gcd(p,|n|)`)
  lies in `radical(moment ideal)` — supplied level-by-level by the first-return renewal (THM-1770),
  well-founded but not uniformly closed.
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
survives `d` moments, so detection depth grows with `d`. The finish must be an asymptotic (Watson /
factorial-graded) argument, not elimination.

## Division of labour (proposed)

- **Symmetric-top Watson dominance** → boxeph (owns THM-1565 Radial Lemma / Watson–Nevanlinna). Target:
  prove `a_{max}≠0 ⟹ E[P^m]≠0` for large `m` for two shells, then general.
- **Asymmetric-top / bottom-up LP** → klein (THM-1700) + mac-mini (renewal THM-1770). Target: the
  general multi-straddle cancellation LP.
- **EMP / radial Piece-1 extensions to algebraic (non-poly) `p`** → the last radial-1 residual.
- **Finite-stratum certificate bank** → death-star + codex (each new closed stratum is evidence and a
  Lean target).

**The single sentence:** GMC(2) is one theorem away — *the symmetric-top factorial-graded Watson
dominance* `a_{max}≠0 ⟹ E[P^m]≠0` (`m≫0`) — above a fully proved DvdK-angular + broad-stratum base
and a Lean-checked reduction; the asymmetric case reduces bottom-up.
