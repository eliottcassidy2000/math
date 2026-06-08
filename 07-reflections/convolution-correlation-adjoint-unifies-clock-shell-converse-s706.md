---
source: opus-2026-06-07-S706 (user seed: cross-correlation is the adjoint of convolution)
status: SYNTHESIS — the convolution/cross-correlation ADJOINT duality is the single operator-theoretic
  root of my recent trilogy and several older repo threads. (D1) ⟨h*a,b⟩=⟨a,h⋆b⟩ so correlation =
  adjoint of convolution; (D2) h⋆g=(σh̄)*g, so the antipodal σ (S702) IS the adjoint reflection.
  REPO MAP: shell face (sums v_i+v_j mod 2n−1, S700) = convolution 1_S*1_S; clock face (differences
  v_i−v_j mod n, S701) = cross-correlation 1_S⋆1_S; they are ADJOINT, related by σ (S702). Tournament
  converse (HYP-2283) = adjoint of the circulant adjacency-convolution; self-converse worry-set
  (THM-402) = self-adjoint; skew-adjacency S*=−S. Additive energy = ‖autocorrelation‖² = Σ|F|^4
  (Wiener–Khinchin) = the unit-distance/Cayley object (S599). All Fourier-diagonal; adjoint = conjugate
  the symbol (S599g). THM-441, HYP-2311. All identities verified exactly.
tags: [convolution, cross-correlation, autocorrelation, adjoint, fourier, wiener-khinchin, additive-energy,
  antipodal, sigma, clock-shell, converse, self-adjoint, worry-set, unit-distance, cayley, spectral,
  positive-definite, synthesis, lonely-runner, tournaments]
---

# Cross-correlation is the adjoint of convolution — and that is the whole engine room

**Prompt (user):** for complex-valued functions the cross-correlation operator is the adjoint of the
convolution operator; integrate both, their relation, and how they fit in and extend the repo.

The seed turned out to be the **operator-theoretic root** of my last three sessions and of several
older repo threads at once. Convolution and correlation are not two tools; they are an operator and
its adjoint, and the repo's two LRC faces, its tournament converse, and its additive-energy/
unit-distance counting are all instances of that one duality.

## 1. The duality, exactly

With `⟨a,b⟩=Σ conj(a)b` on a finite abelian group:
- **adjoint:** `⟨h*a,b⟩=⟨a,h⋆b⟩`, so `C_h^* = R_h` — correlation is the adjoint of convolution
  (verified `~10⁻¹⁴`).
- **reflection link:** `h⋆g=(σh̄)*g`, `(σf)(x)=f(−x)` — the **antipodal involution σ** (S702) is
  precisely the reflection that turns convolution into its adjoint (verified exactly).
- **Fourier:** `\hat{h*g}=\hat h\,\hat g`, `\hat{h⋆g}=\overline{\hat h}\,\hat g` — convolution
  operators are exactly the Fourier-diagonal ones; **the adjoint conjugates the symbol** (verified
  `symbol(A^T)=conj(symbol(A))`).

## 2. The trilogy collapses to one statement

For a speed set `S` with indicator `1_S`:
- `1_S * 1_S` = **sumset** multiplicities (`v_i+v_j`) = the **SHELL face** (sums mod `2n−1`,
  shell-partners `a+b≡0`, S700). *Convolution = sums = shell.*
- `1_S ⋆ 1_S` = **difference-set** multiplicities (`v_i−v_j`) = the **CLOCK face** (differences mod
  `n`, torsion-leak, S701). *Correlation = differences = clock.*
- `1_S⋆1_S = 1_{σS}*1_S`: the clock face **is** the shell face precomposed with `σ`.

> So S700 (shell), S701 (clock), S702 (the antipodal σ) are **one fact**: the clock and shell faces of
> LRC are the **correlation and convolution** of the speed set, which are **mutually adjoint**, and σ
> is the adjoint's reflection. THM-428's "two coprime towers `n` vs `2n−1`" is just *difference-modulus
> vs sum-modulus of the same set.* I had been describing the two faces and their σ-link three times; the
> user's seed names the single object: **adjoint duality.**

## 3. The tournament converse is literally the adjoint

A circulant tournament's adjacency is convolution by its connection set `1_H`. The **converse**
(reverse all arcs) is `A^T = ` convolution by `1_{−H}` `= C_h^*` — the adjoint (verified). Hence:
- **self-converse round tournament (the LRC worry-set, THM-402) = self-adjoint** (up to a multiplier
  `λH=−H`);
- the `{±1}` skew-adjacency `S` is **skew-Hermitian** `S^*=−S` (verified exact);
- HYP-2283's "converse swaps `H±Pf`" is the spectral action of the adjoint (conjugate the symbol;
  `H`= the swap-even part, `Pf`= the swap-odd √).

So the repo's central tournament symmetry (converse) and the LRC worry-set (self-converse) are the
adjoint and the self-adjoint locus of the same convolution algebra.

## 4. Additive energy, unit distance, and positivity

`E(S)=#{a+b=c+d}=‖1_S⋆1_S‖²=Σ_ξ|\hat{1_S}(ξ)|⁴` (verified three ways). The **autocorrelation**
`1_S⋆1_S` has Fourier transform `|\hat{1_S}|²≥0` (Wiener–Khinchin, positive-definite). The
**unit-distance count** `U(P)` is the autocorrelation summed over the unit sphere; maximising unit
distances = maximising autocorrelation mass on the sphere — exactly the repo's *unit-distance =
additive energy on a cyclotomic Cayley graph* (S599) and the *spectral unification* (S599g, Bessel/
Dirichlet). The non-negativity of `|\hat{1_S}|²` is the analytic spine under those bounds and under
the LRC moment method (THM-406: the covering-depth factorial moments are autocorrelation integrals; the
Vitali wall is their non-truncation). My S705 `U_count` (unit-distance census for u(22)) is literally
an autocorrelation-on-the-sphere evaluation.

## 5. The recurring punchline: self-adjoint = extremal/rigid

Across all three problems the **self-adjoint** locus is the rigid/tight one:
- self-converse round tournaments = the LRC worry-set (THM-402, `M=1/n` tight);
- real symmetric autocorrelations / `|\hat{1_S}|²` = tight additive energy (the AP/lattice extremals);
- `σ`-symmetric (cyclotomic) configs = the tight/hard LRC (THM-430's "solvable=cyclotomic=tight"
  inversion), now seen as **self-adjointness**.

> Adjoint symmetry is extremal rigidity. The hard, tight, worry-set cases are exactly where the
> convolution operator is self-adjoint (real spectrum, `σ`-invariant) — the operator-theoretic form of
> the repo's recurring "too symmetric to be loose."

## 6. Honest status

- **Proved (classical, verified exactly):** the adjoint identity, the σ-reflection link, the Fourier
  conjugate-symbol law, additive energy = ‖autocorrelation‖² = 4th Fourier moment, converse = adjoint,
  skew-Hermitian skew-adjacency.
- **New = the synthesis (the map):** identifying the repo's clock/shell faces, tournament converse, and
  additive-energy/unit-distance counting as the *single* convolution–correlation adjoint pair, with σ
  as the adjoint reflection and self-adjoint = worry-set/extremal. A unifying lens; **no new bound, no
  open case resolved.**
- **What it buys:** a clean place to stand — every repo object is a convolution operator on a
  cyclotomic group; "what is its adjoint?" (the converse / the σ-reflected face) and "is it
  self-adjoint?" (worry-set?) are now first questions. Concrete next probe: read the LRC covering-depth
  moments (THM-406) explicitly as autocorrelation integrals and ask whether the Vitali-wall
  non-truncation is the failure of the autocorrelation to be a *finite* positive combination of
  characters (a finite-rank/positivity statement) — the analytic twin of the S704 depth wall.

**Artifacts:** `04-computation/convolution_correlation_adjoint_s706.py` (+`.out`). Theorem **THM-441**.
New **HYP-2311**. Builds on THM-420/421/428/430 (clock/shell/two-tower/antipodal), THM-402/403 &
HYP-2283 (converse/worry-set), S599/S599g (additive energy / spectral), THM-406 (Vitali wall), S705
(u(22) autocorrelation census).
