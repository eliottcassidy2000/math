---
id: THM-1715
title: "THE STRUCTURE OF TRINOMIAL TNC: positivity + recurrence + roots of unity. (A) POSITIVITY, PROVED: in the gauge r_0=r_d=1 (middle coefficient = a), CT(Lambda^m) = sum_y (POSITIVE multinomial) a^y, so real a>0 gives CT(m)>0 -- a nullcone point needs NONREAL phase. This strengthens the common-ray cone (THM-1705) in the middle-coefficient direction and localises any violator to complex a. (B) RECURRENCE: CT(Lambda^m) is P-recursive (a diagonal of an algebraic GF); at nondegenerate saddles CT(m) = sum_j c_j w_j^m with w_j = R(u_j)/u_j^N the saddle values, a linear recurrence whose characteristic roots are the w_j. DISTINCT values => Vandermonde => TNC (THM-1625 recast). CT(3)=0 at the witness 1+u^3-u^6 is ONE linear condition on the c_j and cannot force CT(6)=0 -- indeed CT(6)=-30. (C) ROOTS OF UNITY live in the BRANCH SYMMETRY, not the tuned points (THM-1710 refuted the latter): the N small branches form a mu_N orbit to leading order, and a symmetric R (R(u)=S(u^g)) makes saddle VALUES collide in mu_g orbits. When g | N these collisions DESCEND to a smaller instance; the generic trinomial collision is mu_g-symmetric (18/19 in the sweep). The rare ASYMMETRIC collision (e.g. 1+3u+u^3, g=1) is closed by the resultant Res_a(CT(m0),CT(2m0)) != 0 (THM-1710). So trinomial TNC (already proved per-pattern by THM-1680) decomposes structurally as: positivity forces complex phase, the recurrence gives distinct-value Vandermonde, mu_g roots of unity handle symmetric collisions, the resultant handles the rest"
status: PROVED (A); recast/VERIFIED (B, C). Structural consolidation of the trinomial closure (THM-1680); the k-nomial uniform version is the open target.
author: opus-2026-07-20-S424
depends_on: [THM-1680 (trinomial gcd -- the closure), THM-1625 (Vandermonde/collision), THM-1710 (resultant; cyclotomic-tuned-point refutation), THM-1705 (common-ray cone), THM-415 (vanishing sums)]
---

# THM-1715 — Trinomial TNC: positivity, recurrence, roots of unity

The owner's steer was to "keep roots of unity and recurrence." This note places both
correctly within the trinomial closure (THM-1680), after THM-1710 refuted the naive
cyclotomic-tuned-point route.

## A. Positivity (PROVED)

In the gauge `r_0 = r_d = 1` a trinomial `R = 1 + a u^j + u^d` has, for the charge
`Λ = u^{−N}R`,

```
CT(Λ^m) = Σ_{(x,y,z): x+y+z=m,\ −Nx+(j−N)y+(d−N)z=0}  \binom{m}{x,y,z}\, a^{y}
        = Σ_y c_{m,y}\, a^y ,   c_{m,y} = \binom{m}{x,y,z} > 0 .
```

**All coefficients `c_{m,y}` are positive multinomials** (verified: `{−2,1,4}` gives
`3+3a²`, `15+60a²+15a⁴`; `{−2,−1,2}` gives `3`, `15+20a⁴`). Hence:

> **For real `a > 0`, `CT(Λ^m) > 0` at every reachable level.** A nullcone point `a*` must
> be **non-real** (genuinely complex phase). This strengthens the common-ray cone (THM-1705
> §1) along the middle-coefficient axis and confines any violator to `arg(a) ∉ πℤ`.

## B. The recurrence, and the Vandermonde step

`CT(Λ^m) = [u^{Nm}]R^m` is a **diagonal of an algebraic generating function**, hence
**P-recursive** (satisfies a linear recurrence with polynomial-in-`m` coefficients). At
nondegenerate saddles the closed form is

```
CT(Λ^m) = Σ_j c_j\, w_j^m ,     w_j = R(u_j)/u_j^N  (saddle values, THM-1615) ,
```

a linear recurrence whose **characteristic roots are the saddle values** `w_j`. If the `w_j`
are **distinct**, the `c_j` are constants and `Σ c_j w_j^m = 0 ∀m` forces `c_j = 0`
(Vandermonde) — impossible at a genuine saddle. This is THM-1625 §1 in recurrence form.

**Witness.** `R = 1 + u^3 − u^6`, `N = 2`: `CT(m) = 0,0,0,0,0,−30,0,0,126,0,0,1386,…`
(supported on `m ≡ 0 \bmod 3`, P-recursive). `CT(3) = 0` is **one** linear condition on the
`c_j`; it cannot force `CT(6) = 0`, and `CT(6) = −30`.

## C. Where the roots of unity actually live

THM-1710 refuted the idea that *tuned points* are roots of unity. The correct home is the
**branch symmetry**:

- The `N` small branches of `u^N = tR(u)` are a **`μ_N` orbit** to leading order
  (`u_i ≈ ω^i (r_0 t)^{1/N}`, `ω = e^{2πi/N}`). Roots of unity are intrinsic to the covering.
- A **symmetric** `R` (`R(u) = S(u^g)`, `g = gcd` of exponents `≥ 2`) makes the **saddle
  values collide in `μ_g` orbits**. When `g | N` these collisions **descend** to a smaller
  TNC instance (THM-1625 §3).
- **The generic trinomial collision is `μ_g`-symmetric** — 18 of 19 collision cases in the
  coefficient sweep. Roots of unity handle the bulk.
- The **rare asymmetric collision** (`g = 1`, e.g. `1 + 3u + u^3`) has no `μ_g` symmetry and
  is closed instead by the **resultant** `Res_a(CT(m_0), CT(2m_0)) ≠ 0` (THM-1710) — the
  angular separation of the two levels' root sets.

## D. The trinomial closure, decomposed

`R` a trinomial ⟹ not a nullcone element, by:

| step | mechanism | tool |
|---|---|---|
| real-positive `a` | positivity `CT > 0` | §A |
| complex `a`, distinct saddle values | Vandermonde on `w_j^m` | §B / THM-1625 |
| symmetric collision (`μ_g`, `g\|N`) | descent to smaller instance | §C / THM-1625 §3 |
| asymmetric collision (`g=1`) | angular level-separation | §C / THM-1710 |

All four together re-prove THM-1680's per-pattern closure, now with the *reason* exposed:
positivity kills the reals, Vandermonde kills the distinct-value bulk, `μ_g` roots of unity
kill the symmetric collisions, and the resultant kills the asymmetric remainder.

## E. Toward the uniform k-nomial completion

The decomposition generalises, and each step has a `k`-nomial analogue:

1. **Positivity** holds verbatim: `CT(Λ^m)` has positive coefficients in the `k−2` gauge
   parameters, so a violator lies in `(ℂ ∖ ℝ_{≥0})^{k−2}` — the phase must be tuned in
   *every* coordinate. This is a genuine constraint worth exploiting.
2. **Vandermonde** on the `w_j` closes the distinct-value locus in every dimension.
3. **`μ_g` descent** closes symmetric `k`-nomials.
4. The residual is asymmetric-collision `k`-nomials, where the **amoeba/multinomial-radius
   separation** (HYP-8520) is the uniform statement: the root-amoebae of `CT(Λ^{ℓ m_0})` sit
   at radii growing with `ℓ`, so no common zero. Positivity (§1) plus this radius growth is
   the concrete finish.

## Verification

`04-computation/tnc_roots_of_unity_recurrence_opus_S424.py` (positivity; the witness
recurrence sequence; the `μ_3` symmetry of `1+u^3−u^6`), `tnc_collision_symmetry_opus_S424.py`
(18/19 collisions `μ_g`-symmetric; the `1+3u+u^3` asymmetric exception),
`tnc_annulus_test_opus_S424.py` (Eneström–Kakeya annuli nest — the separation is angular, not
radial). Outputs in `05-knowledge/results/`.
