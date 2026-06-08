# THM-441 — Cross-correlation is the adjoint of convolution; this duality unifies the repo's clock/shell faces, the tournament converse, and additive energy

**Status:** PROVED (the operator identities — classical, recomputed exactly) + VERIFIED (the repo
identifications, exact/numeric) + SYNTHESIS (the unifying map). A unifying lens, not a new bound.
**Source:** opus-2026-06-07-S706, from the user's seed "for complex-valued functions the
cross-correlation operator is the adjoint of the convolution operator." Ties together THM-430/S702
(antipodal σ), THM-420/S700 (shell/sums), THM-421/S701 (clock/differences), HYP-2283 (converse),
S599 (unit-distance = additive energy on a Cayley graph), S599g (spectral unification), THM-402/403.

## The operator duality (PROVED, classical)

On a finite abelian group `G` with the Hermitian inner product `⟨a,b⟩=Σ conj(a)·b`:
```
   convolution        (h*g)(x) = Σ_y h(y) g(x−y)        C_h : g ↦ h*g
   cross-correlation  (h⋆g)(x) = Σ_y conj(h(y)) g(x+y)  R_h : g ↦ h⋆g
```
> **(D1) Adjoint.** `⟨h*a, b⟩ = ⟨a, h⋆b⟩`, i.e. `C_h^* = R_h`. (Verified `|LHS−RHS|~10⁻¹⁴`, N=7,12,27.)
> **(D2) Reflection link.** `h⋆g = (σh̄) * g`, where `(σf)(x)=f(−x)` is the **antipodal involution**
> (S702) and `h̄=conj h`. (Verified exactly.) So **correlation = convolution precomposed with σ (and
> conjugation)** — the adjoint is the σ-reflected, conjugated convolution.
> **(D3) Fourier / "adjoint = conjugate the symbol".** `\hat{h*g}=\hat h·\hat g`,
> `\hat{h⋆g}=conj(\hat h)·\hat g`. Convolution operators are exactly the Fourier-diagonal (circulant)
> operators; the adjoint conjugates the symbol per frequency. (Verified `symbol(A^T)=conj(symbol(A))`.)

## The repo map (the synthesis)

> **(M1) LRC clock ⟂ shell are an adjoint pair.** For a speed set `S` (indicator `1_S`):
> - `1_S * 1_S =` **sumset** multiplicities (`v_i+v_j`) `=` the **SHELL face** (sums mod `2n−1`,
>   shell-partners `a+b≡0`, THM-420/S700) — the **convolution**.
> - `1_S ⋆ 1_S =` **difference-set** multiplicities (`v_i−v_j`) `=` the **CLOCK face** (differences
>   mod `n`, torsion-leak, THM-421/S701) — the **cross-correlation / autocorrelation**.
> - By (D2), `1_S⋆1_S = 1_{σS}*1_S`: the clock face is the shell face precomposed with `σ`. So
>   **clock and shell are adjoint operators related by the antipodal `σ` of THM-430** — the S700/S701/
>   S702 trilogy is one statement: *shell = convolution, clock = correlation, σ = the adjoint.*
>   (Verified: `conv==sumset-mult`, `corr==diff-mult`, exact.)
>
> **(M2) The tournament converse is the adjoint.** For a circulant tournament `A=Cay(ℤ/m,H)` (its
> adjacency = convolution by `1_H`), the **converse** `T↦T*` reverses arcs, giving `A^T=` convolution
> by `1_{−H}` `= C_h^* = R_h` — the adjoint. **Self-converse** (the LRC worry-set, THM-402) `=`
> self-adjoint-up-to-multiplier (`∃λ: λH=−H`); the `{±1}` skew-adjacency `S` is **skew-Hermitian**
> `S^*=−S`. (Verified n=5,7: `symbol(A^T)=conj(symbol(A))`, `S^*=−S` exact, self-converse via
> multiplier.) This is the operator content of HYP-2283 (converse swaps `H±Pf`).
>
> **(M3) Additive energy = autocorrelation norm = 4th Fourier moment (Wiener–Khinchin).**
> `E(S) = #{a+b=c+d} = ‖1_S⋆1_S‖² = Σ_ξ |\hat{1_S}(ξ)|⁴`. (Verified three ways equal: 53,85,53.)
> The unit-distance count `U(P)` is the **autocorrelation `1_P⋆1_P` summed over the unit sphere**;
> maximising unit distances = maximising autocorrelation mass on the sphere. This is exactly the
> repo's *unit-distance = additive energy on a cyclotomic Cayley graph* (S599) and the *spectral
> unification* (S599g) — convolution operators diagonalised by one Fourier transform, adjoint =
> conjugate symbol.

## Why it matters

- **One duality under the trilogy.** S700 (shell/sums), S701 (clock/differences), S702 (antipodal σ)
  are now a single fact: the clock and shell faces of LRC are the **cross-correlation and
  convolution** of the speed set, which are **mutually adjoint**, and `σ` is precisely the reflection
  that realises the adjoint (D2). The "two-tower" coprimality (THM-428, `n` vs `2n−1`) is the
  difference-modulus vs sum-modulus of the same set.
- **Positivity engine.** `1_P⋆1_P` has Fourier transform `|\hat{1_P}|² ≥ 0` (positive-definite,
  Wiener–Khinchin). This non-negativity is the analytic spine of the additive-energy/unit-distance
  bounds (S599g Bessel/Dirichlet) and of the LRC moment method (THM-406): the covering-depth factorial
  moments are autocorrelation integrals, and the Vitali wall is their failure to truncate.
- **Self-adjoint = worry-set = tight.** Across all three problems the **self-adjoint** locus is the
  rigid/extremal one: self-converse round tournaments (THM-402 worry-set), real (`|\hat{1_S}|²`)
  symmetric autocorrelations (tight additive energy), the cyclotomic (`σ`-symmetric) configs. Adjoint
  symmetry = the extremal rigidity (the THM-430 "tight = cyclotomic" inversion, now operator-theoretic).

## Scope / honesty

- (D1)–(D3) are classical Fourier-analysis identities (verified exactly here). The **content is the
  repo map** (M1)–(M3): identifying the repo's clock/shell/converse/additive-energy objects as the
  convolution–correlation adjoint pair. This is a unifying reframing; it introduces no new bound and
  resolves no open case.
- (M1)'s `1_S⋆1_S=1_{σS}*1_S` is exact; the "adjoint pair" reading of clock⟂shell is the new lens.

**Tournament reading:** the adjacency of any Cayley/circulant object is a convolution; its converse/
transpose is the cross-correlation adjoint; Hermitian (self-converse) ⟺ real spectrum ⟺ the worry-set.
The whole repo lives among **convolution operators on cyclotomic groups**, with `σ`/conjugation as the
single adjoint.

**Artifacts:** `04-computation/convolution_correlation_adjoint_s706.py` (+`.out`). Reflection
`07-reflections/convolution-correlation-adjoint-unifies-clock-shell-converse-s706.md`. New:
**HYP-2311**. Builds on THM-420/421/428/430, THM-402/403, HYP-2283, S599/S599g, S700/S701/S702.
