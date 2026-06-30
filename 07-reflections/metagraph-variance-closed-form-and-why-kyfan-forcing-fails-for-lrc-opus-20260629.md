# Two threads that explain each other: the metagraph variance has a succession-count closed form (proved), and Ky-Fan topological forcing FAILS for the LRC (the even/Borsuk–Ulam regime) — which is exactly why the moment machinery is the only route to the floor

*opus-2026-06-29. Owner: pursue both the closed form for the metagraph variance AND the Ky-Fan
forced-existence count, working them back and forth. Thread A solved cleanly; Thread B hit the
fundamental even/odd wall — and B's failure is the proof that A's machinery is necessary.*

## Thread A (PROVED): the metagraph variance is a succession-weighted permutation count
`H = Σ_π 1[π HP]`, and the pair-covariance (previous reflection) reduces — via the left-regular
symmetry `Cov(π,π')` depends only on `π^{-1}π'` — to a single-permutation sum. With `g(σ)=2^{\#\{k:
σ_{k+1}=σ_k+1\}}` if `σ` has **no descending value-adjacency** (`σ_{k+1}=σ_k−1` never), else `0`:
> **`Var(H) = \dfrac{n!}{4^{n-1}}\big(G(n)-n!\big)`,  equivalently  `\dfrac{E[H^2]}{E[H]^2} = \dfrac{G(n)}{n!}`,
>  where `G(n)=Σ_{σ∈S_n} g(σ)`.**

Inclusion–exclusion on the `n−1` value-bonds (`w_i = 1+[asc_i]−[desc_i]`) collapses this to a
**composition closed form** (each bond-run must be even ⟺ each consecutive-value block is odd):
> **`G(n) = Σ_{\text{compositions of }n\text{ into ODD parts}} k!\cdot 2^{\#\{\text{parts}\ge 3\}}`**,
> `k=` #parts; OGF `Σ_n G(n)x^n = Σ_{k\ge0} k!\,W(x)^k`, `W(x)=x(1+x^2)/(1-x^2)`.

Verified (`Var(H)=3/4,3,285/16,585/4,206325/128`; `G=1,2,8,32,158,928,6350,49752,439670,…`). `G(n)` is a
**succession-weighted count** — Riordan's classical permutation-succession theory (Hertzsprung,
[A002464](https://oeis.org/A002464); cf. [A001100](https://oeis.org/A001100)). And `G(n)/n! =
E[H^2]/E[H]^2 → 1^+` (decreasing `4/3, 4/3, 79/60, 58/45, …`), so **`H` concentrates** — the regular-end
upper tail (idea 6) shrinks relative to the mean. *A clean, classical home for the metagraph second
moment.*

## Thread B (the obstruction, made precise): Ky-Fan forcing FAILS for the LRC
The danger labeling `d(t)=(1[‖s_i t‖<1/14])_i` walks on the cube `Q_{n-1}`, starting and ending at the
all-danger vertex `1…1` (everyone at the observer at `t=0`); **lonely = visits to `0…0`.** Since
`‖s_i(−t)‖=‖s_i t‖`, the walk is **even** (`d(−t)=d(t)`), so:
- **The lonely count is EVEN** — mirror pairs `{t,−t}` (fixed points `t=0,1/2` are non-lonely for
  covering sets, which carry even speeds). Borsuk–Ulam regime, *not* Rédei.
- **The signed/degree count is 0** — entering and leaving `0…0` cancel along a closed walk.
- **The "factor out the 2 (14=2·7)" hope fails:** the HALF-circle `(0,1/2)` lonely-interval count is
  **NOT always odd** — computed parity is *mixed* (5/12 covering sets odd: `2,8,4,24,16,28,30,30,…`).

> **So there is NO parity/degree forcing of a lonely point.** The metagraph's `P_n(−1)=SC>0` forces
> self-converse classes (Ky-Fan, an ODD alternating count = Rédei's `H` odd); the LRC's analogous count
> is EVEN with zero degree, so the same machine gives `0`, which does not force `>0`. This is the exact
> topological statement of why LRC is hard (THM-582's two-index split: `x=2` odd/forced vs `x=−1`
> even/un-forced, seen here as a concrete walk on `Q_{n-1}`).

## The back-and-forth payoff: B's failure is A's justification
Working A↔B: **because topological forcing fails (B), the existence of a lonely point cannot be certified
by an alternating count — it must be certified by SIZE, i.e. a second-moment/measure bound (A).** The
metagraph is forced *combinatorially* (odd ⇒ ≥1, no measure needed); the LRC, lacking the odd certificate,
falls back on `Var`/`S₂`/the `ζ(2)` floor — exactly the moment machinery whose finite rehearsal is
Thread A. The two threads are the two horns of the two-index dichotomy:
| | certificate | mechanism | tool |
|---|---|---|---|
| metagraph (Rédei, `x=2`) | `H` ODD = `P_n(−1)` alternating count ≠ 0 | **topological forcing** (Ky-Fan) | parity |
| LRC floor (`x=−1`) | lonely count EVEN, degree `0` — no forcing | **size/concentration** | `Var(H)` rehearsal → `S₂`/`ζ(2)` |

So the same `H = Σ_π 1[π HP]` gives BOTH: its *parity* (odd, Rédei — the forcing side) and its *variance*
(Thread A's closed form — the size side). The LRC inherits only the size side; that is the precise
content of "one dimension past Littlewood" at the level of certificates.

## Status
- **PROVED (opus):** `Var(H)=\frac{n!}{4^{n-1}}(G(n)-n!)`; `E[H^2]/E[H]^2=G(n)/n!`; `G(n)=Σ_{odd-part
  comps} k!2^{#≥3}` = succession count, OGF `Σ k!W(x)^k`, `W=x(1+x^2)/(1-x^2)`; verified n≤7.
- **Computed (B):** LRC lonely count even, signed count 0, half-parity mixed (5/12) ⟹ **no Ky-Fan
  forcing**; the precise topological reason the LRC needs the moment floor.
- **Synthesis:** B (no forcing) justifies A (the moment machinery); metagraph `H` carries both the
  forcing certificate (parity, Rédei) and the size object (variance) — the LRC inherits only the latter.
- **Open:** asymptotics of `G(n)/n!−1` (the concentration rate); whether a *non-parity* equivariant
  invariant (e.g. a `Z_7`/apex-7 secondary obstruction, not `Z_2`) could force the LRC where `Z_2` cannot.

Related: the metagraph-covariance/FKG-fails reflection (the `Cov` formula this closes); HYP-2823
(variance extremality), THM-582 (two-index Rédei/lonely), klein THM-587 (`P_n(±1)`), mac-mini HYP-3544
(Ky-Fan on `Q_d`), kyfan-tucker-ham-sandwich-the-even-parity-of-LRC, the Siegel–Rogers moment hierarchy,
Riordan successions ([A002464](https://oeis.org/A002464)), OPEN-Q-108.
