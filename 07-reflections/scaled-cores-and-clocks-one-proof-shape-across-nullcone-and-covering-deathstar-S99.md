# Scaled cores and clocks: one proof-shape ("scale the core, then close on a modular clock") across the nullcone (GMC2) and covering (LRC) threads

> **REFUTED / HISTORICAL (MISTAKE-232).** The shared wording “scale” and
> “clock” does not define a map from GMC(2) Frobenius to the additive LRC
> clocks, and the original text incorrectly assigned a Paley spectrum to the
> composite modulus `14`. Only the exact, independent Paley identities at the
> primes `7` and `13` and THM-2057's separate clock proof survive.

**death-star-2026-07-21-S99** (HYP-8876). Owner: continue, merge in scaled cores and clocks. Merged the
two live "scaled-core" objects — GMC2's *dilated* face cores (my NC2 capstone) and codex's *scaled zeta-core*
(THM-2057) — with two prime Paley examples at `7` and `13`. The proposed
bridge through tournament zeta (THM-1926) is retained only as historical
motivation, not a reduction or verified proof-shape correspondence.

## The shared proof-shape: SCALE the core, then CLOSE on a modular CLOCK
Both flagship proofs run the same engine:

| | scale (the core) | clock (the modular closure) | certificate |
|---|---|---|---|
| **GMC2 / THM-2022** (nullcone, my capstone) | dilate every face channel by `p` (`dilate`, `×p`) | reduce in the residue field `𝒪_K/𝔭 = ℤ/p` (Frobenius `x↦x^p`); Kummer/Lucas mod `p` **is** clock arithmetic | the tied face survives as `Q̄^p` — the **p-th power = the clock's p-periodicity** |
| **LRC / THM-2057** (covering, codex) | scale the zeta-core `{1,…,11,13}` by `a` (`×a`) | reduce by exact modular orbits on the **12a- and 14a-clocks**; `lcm(12a,14a)=84a` forces the binding branch | orbit closure on the clock |

The two rows share vocabulary but not a typed correspondence. In particular,
`Z/pZ` is a field of characteristic `p`, whereas `Z/12aZ` and `Z/14aZ` are
additive composite clocks, and Frobenius has no supplied analogue preserving
the LRC witness predicate. GMC2's cleanest feature — that the certificate is
a **pure `p`-th power `Q̄^p`** (Frobenius fixes the natCast weights, so the whole tied face survives as one power)
— has no known LRC transfer. A single-power clock certificate may be explored
as a heuristic prompt only after specifying a source, target, map, preserved
predicate, lost data, and decisive test.

## The prime examples `7` and `13` carry Paley spectra (verified)
The two prime moduli have exact Paley spectra:
- **7** (`≡ 3 mod 4`, Paley **tournament**): the skew Jacobsthal matrix has eigenvalues `±i√7` (`|·| = 2.646 = √7`);
  the tournament-matrix atoms are `(−1±i√7)/2`.
- **13** (`≡ 1 mod 4`, Paley **graph**): eigenvalues `(−1±√13)/2` (`√13` in the real spectrum).
- **14 = 2·7** is only the runner count. It is not a prime power and no
  Paley graph or tournament at modulus `14` exists in this argument.

So the tournament/graph **zeta** (THM-1926: `ζ_T = 1/det(I−uA)`, Euler product over primitive cycles, poles at the
trigonometric/Gauss-sum atoms, `ζ ≡ 1` on the acyclic transitive core) is the **dynamical lens** on the clock:
periodic-orbit ↔ spectrum. That classical lens does not identify a tournament
with an LRC core or prove which part of a covering must be controlled.
MISTAKES-227, -228, -230, and -231 block the former transitive-core,
Paley-LRC, class-number, and entropy explanations.

## Honest scope (avoiding the S90 MISTAKE-214 trap)
- This is only a terminological analogy ("scale then modular-clock") between
  the nullcone and covering proofs, plus verified Paley facts at `7` and `13`.
  It proves nothing and currently supplies no transfer.
- **Explicitly NOT claimed:** a numeric identity "zeta pole = clock modulus." The Paley-*tournament* atom is
  `√((1+p)/4)` (for `p=7`, `√2`), not `√p`; the Gauss sum `√p` lives in the *skew* spectrum (7) or the *graph*
  spectrum (13, which is `≡1 mod 4` — a graph, not a tournament). The `√p`-at-clock-moduli is a spectral
  coincidence worth naming, not a pole equation.
- Surviving value: the exact matrix identities are clean hostile controls for
  Paley computations, while THM-2057 independently exemplifies a successful
  modular-orbit proof. Credits: THM-2057 (codex), THM-1926 (boxeph), THM-878
  (klein), and THM-2022 (codex/fleet).

Cross-links: `nc2-gmc2-lean-formalization-state` (the GMC2 scale+clock = dilation+Frobenius), `gmc-lrc-same-positivity-manoeuvre`
(prior GMC↔LRC bridge), `lrc14-frontier-and-sharp-horn`. Script `04-computation/scaled_cores_clocks_merge_deathstar_S99.py`
(+ `.out`). HYP-8876.
