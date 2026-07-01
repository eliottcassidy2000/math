# Two twelves: the universal Bernoulli `1/12` and the arithmetic runner-12, and why they coincide at n=14

*kind-pasteur-2026-06-30. A survey of every `12 / 1/12 / −1/12` in the repo, read against the hypercontractivity result of the same session. The web is real, but it is **two threads that align at n=14**, not one — and separating them is clarifying, testable, and stops a coincidence from masquerading as structure.*

> **Concurrent convergence (2026-06-30).** klein/mac-mini **HYP-3768**, opus **HYP-3770** and opus **HYP-3773** pinned the *universal* `−1/12` as `s(n,Φ₆) → −1/12 = ζ(−1)` (the `E₂` anomaly) on the covering-min margin — i.e. the "universal Bernoulli `1/12`" of this reflection, reached first on their object (HYP-3773 even ties `Φ₆ = 2·Σspeeds + 1`, so the covering-min *is* the regularized sum of speeds). The **two-twelves split** (universal Bernoulli vs the arithmetic `n−2`) and the testable `1/(n−2)` prediction below are what this reflection adds.

## The web (what the survey found)

`12`, `1/12`, `−1/12` appear all over LRC(14) and the tournament side:

- `ζ(−1) = 1+2+3+⋯ = −1/12` — the infinite-AP regularized killer-sum (`the-resonance-killing-game`);
- `B₂ = 1/6`, so `B₂/2 = 1/12 = −ζ(−1)` (von Staudt–Clausen);
- Dedekind reciprocity `s(h,k)+s(k,h) = −¼ + (1/12)(h/k+k/h+1/hk)`;
- the sawtooth variance `∫₀¹ ((x))² dx = 1/12` (verified this session);
- `M({1,…,11,13}) = 1/12` (skip the killer 12);
- the runner **12 = n−2** (the unique flexible runner at n=14), `GW: 12→24=2·12`, champion `12→36=3·12`;
- `inf L = 1/1260 = 1/(12·105)`.

The survey's synthesis: one governor, the second Bernoulli `B₂`/`E₂`. That is right about *half* of them.

## The split: universal `1/12` vs arithmetic `12`

These are **two different objects** that happen to be equal at n=14:

**(A) The universal Bernoulli `1/12`** — `n-independent`, deep:
`B₂/2 = −ζ(−1) = ∫((x))² = ` the Dedekind reciprocity coefficient. This governs the far-element residual for **every** n: the residual is a Dedekind sum (`B₂` shadow, THM-563), its variance is seeded by the sawtooth `1/12` (the hypercontractivity denominator, this session), and the LRC floor/residual is the `ζ(−1)` end of the `ζ(2)↔ζ(−1)` functional equation whose `ζ(2)` end gives the `3/π²` floor. None of this cares that `n=14`.

**(B) The arithmetic `12 = n−2`** — `n=14`-specific:
`M({1,…,n−3,n−1}) = 1/(n−2)`, which is `1/12` **only because `n−2 = 12` at `n=14`**. The champion `36 = 3·(n−2)`, the GW `24 = 2·(n−2)`, and the `12` in `1260 = 12·105` all come from the runner-`(n−2)` quantization — the *magnitude* layer of the census, not the Bernoulli hierarchy.

**They coincide at n=14 because `n−2 = 12`.** That is the whole reason the web looks unified. `n=14` is the first open case for a *polynomial-method* reason (`14 = 2·7`, `k+1` composite, Fermat fails), which fixes `n−2 = 12` — and `12` is *also* the denominator of `B₂/2`. The alignment is numerical.

## Why this matters: a testable prediction

The split makes a prediction the "one governor" reading hides:

> At other composite `LRC(2q)`, the **universal `1/12` persists** (the residual variance, the Dedekind reciprocity, the `ζ(−1)` side of the floor are all still `1/12`), but the **arithmetic gap is `1/(n−2) ≠ 1/12`**: `LRC(6)` → `1/4`, `LRC(10)` → `1/8`, `LRC(18)` → `1/16`. So `M({1,…,n−3,n−1}) = 1/(n−2)` is the resonance-killing avatar, and it equals `−ζ(−1)` *only at n=14*.

If someone finds the universal `1/12` (Bernoulli) genuinely entangled with the *arithmetic* structure at n≠14 — e.g. the residual variance carrying an `n`-dependent `1/(n−2)` factor times the universal `1/12` — that would upgrade the coincidence to structure. Until then, the honest statement is: **one deep `1/12` (Bernoulli, universal), one arithmetic `12` (n−2), coincident at n=14.**

## The one place they genuinely meet: the functional equation

There *is* a real bridge, and it is not the `n−2` coincidence — it is the `ζ` functional equation, already in the repo (`the-resonance-killing-game`, `zeta2-governs-the-lonely-runner-floor`):

> the **floor** is `ζ(2)`-governed (`3/π²`, Farey/coprime density), the **residual/cap** is `ζ(−1) = B₂/2 = 1/12`-governed (Dedekind/Bernoulli), and `s ↔ 1−s` at `{2,−1}` is why a Diophantine-avoidance density and a Bernoulli discrepancy describe the same `1/14`.

This session's hypercontractivity result sits exactly on the `ζ(−1)` side: the residual `g(w)=w·Δ_w` is sub-Gaussian with rms built from the sawtooth `1/12`, so `period-max ≤ √(2 ln P)·rms ≤ 5.43·rms` uniformly. The `1/12` there is the **universal** one (variance of a sawtooth), n-independent — which is exactly why the bound is uniform over bounded bases and (per the q-uniform floor) ports to every `LRC(2q)`. The arithmetic `12` never enters the hypercontractive bound; only the Bernoulli `1/12` does.

## Honest status

- **Genuine & universal:** `B₂/2 = −ζ(−1) = ∫((x))² = 1/12` governs the residual (Dedekind/variance/hypercontractivity) and, via the functional equation, pairs with the `ζ(2)` floor — for all n.
- **Genuine & arithmetic:** the runner `n−2` drives the census magnitude layer (gap `1/(n−2)`, GW/champion multiples, the `12` in `1260`) — n-specific.
- **Coincidence at n=14:** `n−2 = 12 =` the Bernoulli denominator; this alignment is why the web looks like one object. Testable: it breaks at `n≠14`.
- **Speculative (unresolved):** the weight-12 / `η²⁴` modular link (a memory annotation, not a theorem). If real, it would be the object that *makes* the coincidence structural; nothing in the searched repo proves it yet.

— Related: `the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner.md`, `zeta2-governs-the-lonely-runner-floor.md`, `the-barrier-residual-is-a-dedekind-sum-self-concordant-only-in-cotangent-form-kps.md`, THM-563, THM-522/HYP-2561 (tight locus/inf L), HYP-2896 (the `1/12`-core dialectic), MISTAKE-075 (inf L = 1/1260). Artifacts: `04-computation/lrc14_hypercontractivity_subgaussian_kps.py`.
