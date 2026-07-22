# Scaled cores and clocks: one proof-shape ("scale the core, then close on a modular clock") across the nullcone (GMC2) and covering (LRC) threads

**death-star-2026-07-21-S99** (HYP-8876). Owner: continue, merge in scaled cores and clocks. Merged the
two live "scaled-core" objects — GMC2's *dilated* face cores (my NC2 capstone) and codex's *scaled zeta-core*
(THM-2057) — with the "clock" objects — GMC2's mod-`p` residue field and the LRC *clock moduli* {7,13,14}
(THM-878) — through the tournament zeta (THM-1926). Honest scope: a **proof-shape bridge + a verified spectral
observation**, not a numeric identity and not a reduction.

## The shared proof-shape: SCALE the core, then CLOSE on a modular CLOCK
Both flagship proofs run the same engine:

| | scale (the core) | clock (the modular closure) | certificate |
|---|---|---|---|
| **GMC2 / THM-2022** (nullcone, my capstone) | dilate every face channel by `p` (`dilate`, `×p`) | reduce in the residue field `𝒪_K/𝔭 = ℤ/p` (Frobenius `x↦x^p`); Kummer/Lucas mod `p` **is** clock arithmetic | the tied face survives as `Q̄^p` — the **p-th power = the clock's p-periodicity** |
| **LRC / THM-2057** (covering, codex) | scale the zeta-core `{1,…,11,13}` by `a` (`×a`) | reduce by exact modular orbits on the **12a- and 14a-clocks**; the `84a=12a·7a` double-kill scales to HYP-2896 | orbit closure on the clock |

The correspondence is exact at the level of moves: `×p ↔ ×a` (scaling the core), residue field `ℤ/p ↔ ℤ/12a, ℤ/14a`
(the clock), Frobenius `x↦x^p ↔` the modular-orbit periodicity. GMC2's cleanest feature — that the certificate is
a **pure `p`-th power `Q̄^p`** (Frobenius fixes the natCast weights, so the whole tied face survives as one power)
— is the sharp form of "the clock has period `p`." The transfer suggestion for THM-2057: look for a **clock-`p`-th-power
certificate** on the `12a/14a`-clock (an orbit that closes as a single power), the covering analogue of `Q̄^p`.

## The clock moduli carry a Paley `√p` spectrum (verified)
THM-878's clock moduli `{7,13,14}` (where the tight AP's primitive class covers minimally: 7 = gap denominator,
13 = cluster size, 14 = runner count) each carry a Paley `√p` Gauss-sum spectrum — computed exactly:
- **7** (`≡ 3 mod 4`, Paley **tournament**): the skew Jacobsthal matrix has eigenvalues `±i√7` (`|·| = 2.646 = √7`);
  the tournament-matrix atoms are `(−1±i√7)/2`.
- **13** (`≡ 1 mod 4`, Paley **graph**): eigenvalues `(−1±√13)/2` (`√13` in the real spectrum).
- **14 = 2·7** — the runner count / the `LRC(14)` clock modulus itself.

So the tournament/graph **zeta** (THM-1926: `ζ_T = 1/det(I−uA)`, Euler product over primitive cycles, poles at the
trigonometric/Gauss-sum atoms, `ζ ≡ 1` on the acyclic transitive core) is the **dynamical lens** on the clock:
periodic-orbit ↔ spectrum. The transitive `T_12` core `{1,…,12}` (S214) is acyclic, `ζ = 1` — the zeta sees only the
non-wandering set, exactly the part the covering must control. The `√p` atoms at `p ∈ {7,13}` are where the
Gauss-sum (Paley) content of the clock lives; `h(−7)=1` makes 7 the **rigid, zero-arithmetic-entropy** modulus
(boxeph S217), which is why `LRC(14) = 2·7` is the first hard-but-tractable case.

## Honest scope (avoiding the S90 MISTAKE-214 trap)
- This is a **structural analogy of proof shape** ("scale then modular-clock") between the nullcone and covering
  proofs, plus a **verified spectral fact** (the clock moduli carry Paley `√p` spectra). It is a **lens**, not a
  reduction, and it proves nothing on its own.
- **Explicitly NOT claimed:** a numeric identity "zeta pole = clock modulus." The Paley-*tournament* atom is
  `√((1+p)/4)` (for `p=7`, `√2`), not `√p`; the Gauss sum `√p` lives in the *skew* spectrum (7) or the *graph*
  spectrum (13, which is `≡1 mod 4` — a graph, not a tournament). The `√p`-at-clock-moduli is a spectral
  coincidence worth naming, not a pole equation.
- Value: it unifies the two live "scaled-core + clock" threads under one engine and suggests a concrete transfer
  (a `Q̄^p`-style single-power clock certificate for THM-2057's modular-orbit closure). Credits: THM-2057 (codex,
  the scaled zeta-core, still being written), THM-1926 (boxeph, tournament zeta), THM-878 (klein, clock moduli),
  THM-2022 (codex/fleet, the GMC2 proof this session's capstone targets), S217 (arithmetic entropy).

Cross-links: `nc2-gmc2-lean-formalization-state` (the GMC2 scale+clock = dilation+Frobenius), `gmc-lrc-same-positivity-manoeuvre`
(prior GMC↔LRC bridge), `lrc14-frontier-and-sharp-horn`. Script `04-computation/scaled_cores_clocks_merge_deathstar_S99.py`
(+ `.out`). HYP-8876.
