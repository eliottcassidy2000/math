---
id: HYP-2140
status: EXPLORATION + FRAMEWORK — "rigidity = orbit-type": a taxonomy of all the repo's rigidities
  as orbit types of group actions, unified at the 2-adic seam (even n) where they degenerate
  together. Verified orbit facts; a creative conjecture set (C1–C6) for the LRC frontier.
source: claudebox-2026-06-03-S596
related:
  - HYP-2124
  - HYP-2135
  - HYP-2130
  - HYP-2117
  - HYP-2063
  - HYP-2120
  - HYP-2085
---

# HYP-2140: rigidity is orbit-type — everything as orbits, degenerating at the 2-adic seam

A unifying frame for the many "rigidities" in the repo (local/global from HYP-2135; pinning,
automorphism, spectral, cascade, dynamical). **Each rigidity is a group `G` acting on a space, and
the rigidity TYPE is the ORBIT TYPE.** LRC hardness is *orbit degeneracy*; the cure's arity is the
*orbit size*; and the frontier sits exactly where the orbits collapse — the **2-adic seam**.

## The phase ring and its group actions

Runner phases live in `Z/n`. Three actions (all verified, `lrc_rigidity_as_orbit_type_s596.py`):

- **`⟨−1⟩` reflection** (the witness symmetry, HYP-2135): orbits are `±`-pairs `{x, n−x}`; fixed
  points are `{x : 2x ≡ 0}` — `{0}` for odd `n`, and `{0, n/2}` for even `n`. **The apex `n/2` is
  exactly the second `⟨−1⟩`-fixed point, and it exists iff `n` is even.**
- **`⟨2⟩` doubling** (the binary phase dynamics, HYP-2117): a permutation iff `n` odd; **2-adic
  collapse iff `n` even**. On odd `n` its orbits are the cyclotomic cosets — a single maximal-mixing
  cycle iff 2 is a primitive root (`n=5,11,13,19`).
- **`(Z/n)*` dilation** (`G` is dilation-invariant): orbits are the gcd/divisor classes;
  `#orbits = d(n)`.

## Rigidity = orbit-type (the dictionary, verified facts)

| rigidity (repo) | group `G` | rigid configuration | LRC reading |
|---|---|---|---|
| **local / pinning** (THM-398) | `⟨−1⟩` | a **fixed point** (size-1 orbit) — the apex | single corrector handles it |
| **global / cascade** (HYP-2135) | `⟨−1⟩` | a **free orbit** (`±`-pair) | cure arity = orbit size (2 ⇒ pair-sum) |
| **stabilizer** (HYP-2130) | `Aut ⊂ S_n` | **large stabilizer** (symmetric/tight) ↔ small orbit | `|orbit|·|stab|=|G|`: conservation |
| **dynamical / orbit-closure** | rotation flow on `𝕋^k` | **confined** (rational) vs dense (Q-indep) | Lemma B vs Lemma A; dim = `k−rankΛ` |
| **spectral / determination** | `S_n` | invariant **fiber** is a single orbit (score-rigid) | score- vs cycle-determined boundary |
| **phase / doubling** (HYP-2117) | `⟨2⟩` | **2-adic collapse** (even `n`) vs mixing cycle | the apex / hard frontier vs easy |

## The 2-adic seam: where all rigidities degenerate together (verified)

`n` even `⟺` `⟨−1⟩` gains the apex fixed point `n/2` `⟺` `⟨2⟩` collapses (2-adically) `⟺` the hard
LRC frontier. The **apex is one point wearing every hat**: the `⟨−1⟩`-fixed point (HYP-2135), the
`⟨2⟩`-collapse target (HYP-2117), the field zero-divisor (HYP-2063), and the local pin (THM-398).
Classification of the frontier as orbit-degeneracy:

- **odd `n`, 2 a primitive root** (5,11,13,19): one maximal-mixing `⟨2⟩`-orbit on the units — the
  most ergodic, "fractal-regular," easiest case.
- **odd `n`, `⟨2⟩` splits** (7,17,23): several equal cyclotomic cosets — more orbits to resolve.
- **`n = 2·prime`** (10,14,22): apex + a single 2-collapse — the named hard frontier (n=14).
- **`n = 2^a·m`** (8,16,18,...): deeper 2-adic seam.

## Creative conjecture set

- **C1 (conservation of rigidity = orbit–stabilizer).** Witness-rigidity and structure-rigidity are
  the two factors of `|orbit|·|stab| = |G|`: hard (tight) configs maximize the stabilizer
  (symmetry); easy configs maximize the orbit (genericity). LRC difficulty `∝ |stabilizer of the
  witness configuration|`.
- **C2 (doubling-orbit frontier).** LRC difficulty at `n` is a monotone function of `(v₂(n),
  #cyclotomic cosets of ⟨2⟩ on the odd part)`: pure `2^a` and `2·prime` are the extreme-hard cases;
  `2`-primitive-root odd primes the easiest. (Matches: n=13 easy, n=14 hard.)
- **C3 (sieve arity = orbit size; refines HYP-2135).** The multi-sieve arity needed = the size of
  the largest orbit of the witness symmetry group on the binding runners. `⟨−1⟩` (order 2) ⇒
  pair-sum sieve clears `2·prime`. Richer symmetry (units mod `n`) ⇒ higher arity required — test
  against n=18 (HYP-1992 gate ladder).
- **C4 (orbit-closure controls the margin; refines HYP-2120).** The LRC margin `G−δ` is governed by
  the dynamical orbit-closure dimension `k − rank(Λ)` (Λ = relation lattice): full-dimensional
  (circuit-free) ⇒ the `(1−2δ)^k ≈ e^{−2}` plateau; confined (relations) ⇒ margin `→ 0` (tight).
- **C5 (spectral-rigidity threshold).** The score-determined → cycle-determined boundary (first at
  n=5, where the per-path identity and the perspective-coincidence both break, S590) is where the
  invariant orbit-fibers first carry more than one orbit — a spectral rigidity onset aligned with
  the small-n symmetry collapse.
- **C6 (the apex is the unique multi-hat fixed point).** For every even `n`, `n/2` is simultaneously
  fixed/degenerate under `⟨−1⟩`, `⟨2⟩`, and the field structure; it is the only such point. Hence
  any method keyed to a single fixed point (single corrector) necessarily breaks there, and the
  resolution must change to an orbit-sized (pair, then higher) tool — the structural reason for the
  single→multi corrector transition.

## Convergence with HYP-2124 (opus-S584, independent)

opus independently arrived at the same core from the witness side: "symmetry ⟺ witness-rigidity
duality; AP lonely at `t=j/n` iff `gcd(j,n)=1` (the witness is the `φ(n)` unit clock points = a
rigid `(Z/n)*`-orbit); local pinch rigidity cascades simply-transitively via the unit group and up
the doubling tower; breakage only at a non-unit = a multiple of `n` = the `2·7` seam = n=14's
hardness." That is exactly C1 (conservation) + the `(Z/n)*`/doubling cascade + the 2-adic seam,
reached separately. This document adds the full **taxonomy** (six rigidities as orbit types,
incl. spectral and dynamical/orbit-closure), the explicit **apex-wears-every-hat** unification
(`⟨−1⟩`-fixed ∧ `⟨2⟩`-collapse ∧ zero-divisor ∧ pin, all at even `n`), and the **sieve-arity =
orbit-size** conjecture (C3). Two independent agents + the human's intuition converging on
rigidity = orbit-type is itself evidence the frame is the right one.

## Open / next

- Make C2 quantitative: a difficulty index `R(n) = f(v₂(n), Σ coset sizes)` and regress it against
  the historical/empirical LRC effort per `n`.
- C3 theorem: "minimal sieve arity = max orbit size of the witness symmetry."
- Tie the dynamical orbit-closure (C4) to the relation-lattice theta of HYP-2120 rigorously.

**Artifacts:** `04-computation/lrc_rigidity_as_orbit_type_s596.py` (+`.out`),
`07-reflections/lrc-everything-as-orbits-s596.md`. Builds on HYP-2135 (local/global rigidity),
HYP-2130 (perspectives=rigidity), HYP-2117 (doubling/IFS), HYP-2063 (apex), HYP-2120
(circuit-free/3-term = orbit-closure), HYP-2085 (time-Burnside).
