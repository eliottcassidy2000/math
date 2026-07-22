# The cusp frame is a difficulty-locator — applied across the repo

> **SCOPE CORRECTION 2026-07-21 (codex MISTAKE-226).** “Computable main term
> plus hidden residual” is a useful diagnostic template, but it does not make
> each residual a cusp form or give it a cusp dimension. In particular the
> HYP-8880 identification of LRC clocks with modular cusps and of the LRC
> obstruction with `14.2.a.a` is retracted. The finite tournament counts,
> explicit radial-kernel elements, and figurate identities below survive as
> separate computations; their Eisenstein/cusp names are metaphors unless an
> explicit modular or cohomological carrier is provided.

*boxeph-2026-07-21-S221. Owner: go back through the repo and apply the cusp frame to under-attended problems;
show its power. Builds on S220 (LRC obstruction = f₁₄, the first cusp form), S218 (arithmetic entropy =
genus = cusp dim), the-modular-tournament (Eisenstein 97%/cusp 3% at n=5), THM-1830 (the 3-cycle atom),
THM-1645/S211 (the GMC polar bridge), S207 (figurate/Fibonacci). Verified in
`04-computation/the_cusp_frame_as_a_diagnostic_across_the_repo_boxeph_S221.py`.*

## The frame as a diagnostic

The S220 modular decomposition generalizes into a **difficulty-locator** applicable to any problem with a
computable "main term" and a hidden obstruction:

> **object = EISENSTEIN (the computable floor / main term / local, spectral-or-genus data)
>          ⊕ CUSP (the hidden obstruction = the genus = the deep arithmetic entropy, S218).**
> The **difficulty is always the cusp.** Applying the frame LOCALIZES it to a small/finite object and
> predicts the **"first hard case"** as the first positive cusp dimension.

For LRC(14) the cusp is `f₁₄=14a` (genus 1, the apex 7). The point of this session: **the same frame,
applied to under-attended problems, immediately locates their difficulty** — and they turn out to share the
structure.

## Sweep 1 (under-attended): tournament cospectrality is a reconstruction cusp

Take the char_A spectrum as the *Eisenstein* (local) invariant of a tournament. The **cospectral fiber** —
non-isomorphic tournaments with the same spectrum — is the *cusp*: the hidden data that spectral (local)
information cannot resolve. Verified:

| `n` | #iso classes | #distinct spectra | **cospectral cusp dim** | cospectral family sizes |
|---|---|---|---|---|
| 4 | 4 | 3 | **1** | one pair |
| 5 | 12 | 9 | **3** | a pair + a triple |
| 6 | 56 | 28 | **28** | many pairs, a 3, two 4s, a 5 |

Most tournaments are **spectrally determined** (the Eisenstein/rigid bulk — the transitive one is uniquely
`char = xⁿ`); the cospectral classes are the **reconstruction cusp**. The cusp dimension is the "genus" of
tournament reconstruction, it **first turns positive at `n=4`** (the first cospectral pair), and it grows —
exactly kps's reconstruction wall in cusp language, and the fiber whose `log` is my S218 reconstruction
entropy. This problem had received almost no cusp-frame attention; the frame instantly gives it a genus.

## Sweep 2: intransitivity `c₃` is the tournament's cusp form

The 3-cycle count `c₃(T) = C(n,3) − Σ C(sᵢ,2)` is the tournament's own cusp form (verified): it **vanishes
on the transitive/AP vertex** (`c₃=0` — pure Eisenstein, a gradient flow) and **peaks on the regular/Paley
pole** (`c₃ = 5,14,30` for `n=5,7,9` — the intransitive obstruction). The **3-cycle atom** (THM-1830) is the
minimal cusp — the elementary unit of the intransitive obstruction, the tournament analogue of "the first
newform." So intransitivity is cuspidal: the Eisenstein floor is the score-ordered transitive part; the cusp
is the cyclic deviation.

## Sweep 3: GMC(2) splits Eisenstein (angular, DvdK-closed) ⊕ cusp (radial, Laplace)

The polar bridge `E = L ∘ CT` (THM-1645/S211) *is* an Eisenstein/cusp split: the **angular** constant-term
`CT` (the DvdK integral) is **closed** — the computable, proved *Eisenstein floor* — while the **radial**
Laplace `L` has a **nonzero kernel** (verified: `L(t−1)=0`, `L(t²−3t+1)=0`, `L(t²−2t)=0`) — the *cusp*
obstruction (Laplace determinacy). GMC(2) is rigid because Frobenius (THM-2022) defeats the radial cusp;
**GMC(n≥3) is false because the cusp grows** (the extra counterexamples are "extra cusp forms"). (Honest: `L`
is a linear functional so its kernel is large — the "cusp" here is the *radial-determinacy residual*, an
analogy to the finite genus of the LRC case, not a literal finite cusp space.)

## Sweep 4: figurate cutting = Eisenstein polynomial ⊕ Fibonacci cusp

The cake/bagel cutting sequences (S207) are **smooth Eisenstein polynomials** (`cake(n)=ΣC(n,k)_{k≤3}`, a
degree-3 "ball"), and the **Fibonacci shallow-diagonal** of the same Pascal triangle is the *oscillating,
number-theoretic cusp* reading. The `2ⁿ−cake` deficit (the dropped `C(n,≥4)`) and the Fibonacci wobble are
the cuspidal/arithmetic part of an otherwise smooth polynomial count — S207 recast in the frame.

## The power: one difficulty-locator

```
  problem          EISENSTEIN (computable floor)    CUSP (hidden obstruction = genus = deep entropy)
  ───────          ────────────────────────────     ────────────────────────────────────────────────
  LRC(14)          floor 3/pi^2                      f14=14a (genus 1, first cusp, apex 7)    [S220]
  tournaments      char_A spectrum                  the COSPECTRAL fiber (dim 1,3,28,…)       [Sweep 1]
  intransitivity   transitive (c3=0)                c3 = the 3-cycle count (atom = 1st cusp)  [Sweep 2]
  GMC(2)           angular CT (DvdK-closed)         radial ker L (Laplace determinacy)        [Sweep 3]
  figurate         smooth polynomial                Fibonacci / deviation                     [Sweep 4]
```

The frame is powerful for three reasons, all demonstrated:
1. **It localizes difficulty.** Each hard problem's obstruction collapses to a *small, nameable* cusp object
   — a genus-1 newform, a cospectral fiber, a 3-cycle count, a radial kernel — rather than the whole object.
2. **It predicts the first hard case.** The difficulty appears exactly when the cusp dimension first turns
   positive: LRC at `p=7` (genus 0→1), tournament reconstruction at `n=4` (first cospectral pair),
   intransitivity at `n=3` (the first 3-cycle). Below that, everything is Eisenstein/rigid.
3. **It unifies them.** Every cusp is a **dim = the deep arithmetic entropy (S218)** — the genus, the fiber
   `log`, the kernel: the *hidden* information beyond the local floor. The Eisenstein floor is always the
   easy, computable, local half; **the cusp is always where the proof must go.**

## Scope

Two sweeps are literally modular (LRC `f₁₄`, and the-modular-tournament's `H`); the others (cospectral fiber,
`c₃`, radial `ker L`, Fibonacci) are the *analogous* "computable main term ⊕ hidden obstruction" structure =
the S218 arithmetic entropy. So the contribution is a **diagnostic lens** — find the Eisenstein floor
(easy), localize the cusp (hard), measure its dimension (the genus/entropy) — shown to work and to *reveal
the shared difficulty structure* across LRC, tournaments, GMC, and figurate. Under-attended problems
(tournament cospectrality especially) gain a genus and a "first hard case" for free. Not a new proof step;
a map of where every problem's difficulty lives.

Links: HYP-8885, HYP-8880, HYP-8875, THM-1830, THM-1645,
[[the-lrc14-obstruction-is-the-first-cusp-form-scaled-cores-clocks-are-cusps-boxeph-S220]],
[[arithmetic-entropy-is-a-repo-wide-invariant-the-rigid-extremum-is-the-zero-entropy-point-boxeph-S218]].
