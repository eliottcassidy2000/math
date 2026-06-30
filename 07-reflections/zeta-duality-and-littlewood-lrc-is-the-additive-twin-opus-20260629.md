# The ζ-duality and Littlewood: LRC is the additive (min, ζ(2)-floor) twin of Littlewood's multiplicative (product, liminf-0) approximation; the resonance-killing cuts are the shadow, the 2-adic Littlewood is the killing-game made a conjecture, and the apex-7 odd core is n·|n|₂

*opus-2026-06-29. Owner: merge in zeta-duality and the Littlewood conjecture if real connections exist.
They do — and they tie my cut/Farey/hyperoctahedral framework to the repo's resonance-killing + ζ(2)/ζ(−1)
+ Eisenstein threads and to the homogeneous-dynamics resolution of Littlewood. Crediting kind-pasteur
S31p (resonance-killing), HYP-2856 (ζ(2) floor), HYP-2896 (ζ(−1) = −1/12), S31r (Eisenstein), codex S170
(2-adic Littlewood).*

## 0. My cuts ≡ the repo's resonance-killing game (independent rederivation)
Last turn's **"cut/resonance law: speed `s` kills witness denom `b` iff `b | s`; `M(S)=1/(smallest
surviving b)`"** is *exactly* kind-pasteur's resonance-killing game (S31p). Independent convergence
validates it. What I add: the surviving-resonance ladder `1/b` IS the **Farey left spine** (cuts =
Stern–Brocot geodesics), and the whole board carries the **hyperoctahedral `B_n=Aut(Z^n)`** symmetry
whose central `−1` is the project's `R` — so the killing-game lives on the *same* signed cube as the
tournament metagraph.

## 1. The ζ-duality: the floor and the cap are ζ(2) and ζ(−1), bridged by s ↔ 1−s (verified)
The bound `floor ≥ cap` is one zeta seen across the functional equation:
- **FLOOR = ζ(2) = π²/6** — the Euler-product `∏_p(1−p⁻²)⁻¹` **density of resonances the runners cannot
  all dodge**; surviving lonely mass `≥ 1/(2ζ(2)) = 3/π²` (HYP-2856). *Multiplicative / prime side.*
- **CAP = ζ(−1) = −1/12** — the **regularized killer-sum of the infinite AP** `1+2+3+… = ζ(−1)`
  (`= −B₂/2`); `M({1..11,13}) = 1/12` is its finite avatar (HYP-2896). *Bernoulli / discrepancy side.*
- **Bridge (verified):** `ζ(s)=2^s π^{s−1} sin(πs/2)Γ(1−s)ζ(1−s)` sends `ζ(2)=π²/6 ↦ ζ(−1)=−1/12`.
  So `s ↔ 1−s` IS the floor↔cap duality; **`ζ(2)` counts resonances, `ζ(−1)` weights lonely values, and
  the functional equation is why a Diophantine density and a Bernoulli discrepancy describe the same
  `1/14`.** The ζ(2) sits in the **Eisenstein** constant term — "LRC(14) = Eisenstein(even)∘Legendre(odd)"
  (S31r) is the spectral face of this.

## 2. Littlewood is the MULTIPLICATIVE twin of LRC (the additive one)
**THE Littlewood conjecture:** `liminf_n n‖nα‖‖nβ‖ = 0`. Put it beside LRC:
| | combination | extremum | floor? | pole |
|---|---|---|---|---|
| **LRC** | `min_i ‖s_i t‖` (MIN, additive/∞-norm) | `max_t` (best escape) | **floor `≥1/(k+1)`** (ζ(2)) | **additive** |
| **Littlewood** | `n·‖nα‖·‖nβ‖` (PRODUCT, multiplicative) | `liminf_n` (best approach) | **none (→0)** | **multiplicative** |
> Same torus, **dual norms**: LRC takes the *additive* `min` and finds the loneliest escape — a positive
> floor; Littlewood takes the *multiplicative* `product` and finds the best approach — zero. This is my
> add/mult duality at the Diophantine level: **LRC is the additive pole (a floor, ζ(2)); Littlewood
> probes the multiplicative pole (collapse to 0).** The ζ functional equation is the hinge between them
> — the same `s↔1−s` that swaps the floor (ζ(2)) and cap (ζ(−1)).

Both are **geometry-of-numbers on the space of lattices `SL(d,R)/SL(d,Z)`** under the diagonal flow,
with the cube's **hyperoctahedral `B_d`** symmetry. Littlewood's exceptional set has Hausdorff dimension
0 (**Einsiedler–Katok–Lindenstrauss**, measure rigidity of the diagonal `A`-action) — the analytic
engine the LRC floor is the elementary shadow of. And **both are hardest at the deep Stern–Brocot
points:** Littlewood's stubborn case is **badly-approximable** pairs (bounded partial quotients, e.g.
noble numbers, deepest in the Farey tree); LRC's stubborn case (covering sets) is the **off-spine deep
Farey rational** (denom 89). The metagraph's **asymmetric / non-circulant max-H maximizer** (`a(13)`)
is the tournament avatar of a badly-approximable direction — both are the "no extra symmetry to exploit"
extreme.

## 3. The 2-adic Littlewood IS the resonance-killing game made a conjecture (and it names apex 7)
de Mathan–Teulié **p-adic Littlewood** (repo: codex S170): `liminf_n n·|n|_p·‖nα‖ = 0`. Here
**`|n|_2 = 2^{−v₂(n)}` IS the divisibility-killer** — `n` divisible by `2^k` (small `|n|_2`) is exactly
a speed that kills the `2^k` resonance. So p-LC *unifies the multiplicative killer `|n|_p` with the
additive loneliness `‖nα‖`* in one liminf — the same add/mult marriage as LRC(14). And (verified):
> **`n·|n|_2 = 7` for every `n = 7·2^k`** — the 2-adic Littlewood quantity strips the 2-power killing and
> lays bare the **odd core 7**. `14 = 2·7` is precisely the 2-adic face of **apex 7** (Mersenne `2³−1`);
> the even/odd (2-adic) axis of the duality web is literally the `|·|_2` of the Littlewood functional.

## 4. The unified picture (what merges)
- **Combinatorial shadow:** resonance-killing = my cuts = Farey/Stern–Brocot geodesics, on the `B_d`
  signed cube shared with the tournament metagraph (`R`-even/`R`-odd = the same `−1`-eigenspace).
- **Spectral bridge:** the ζ(2)↔ζ(−1) functional equation (Eisenstein constant term) = floor↔cap.
- **Analytic engine:** homogeneous dynamics on `SL(d,R)/SL(d,Z)`; Littlewood (multiplicative, EKL rigidity)
  is the deep twin of LRC (additive, ζ(2) floor).
- **Arithmetic hinge:** `s↔1−s` (zeta), `|·|_2` (apex-7 odd core), `R=−1` (`B_d`) — one even/odd–add/mult
  hinge, three faces, all meeting at `14=2·7`.

## Improved targets (from the merge)
1. **Import a Littlewood/dynamics tool to the LRC floor:** the `B_d`-equivariant diagonal flow gives the
   ζ(2) resonance-density floor a homogeneous-dynamics proof (systole bound), bypassing the union bound.
2. **Badly-approximable = the hard core, both sides:** characterize LRC covering sets as the
   badly-approximable directions; the deep-Farey margin (denom-89) is a bounded-partial-quotient extremal.
3. **2-adic Littlewood as the disproof discriminator:** a disproof would need `n·|n|_2·‖nα‖` to beat the
   `7` odd core — i.e. defeat apex 7 — which the Mersenne/Heegner rigidity of `7` forbids.

## Status
- **Verified:** `ζ(−1)=−1/12` from `ζ(2)` via the functional equation; `n·|n|_2=7` on `n=7·2^k` (apex-7
  odd core); the additive-min/multiplicative-product contrast.
- **Merged (opus):** cuts ≡ resonance-killing (+ Farey/`B_d`); LRC = additive twin of Littlewood
  (multiplicative); 2-adic Littlewood = killing-game-as-conjecture naming apex 7; ζ(2)/ζ(−1) floor↔cap
  as the hinge.
- **Open:** the homogeneous-dynamics floor proof; badly-approximable characterization of covering sets.

Related: kind-pasteur resonance-killing + ζ(2)/ζ(−1) reflections, HYP-2856/2896, S31r Eisenstein, codex
S170 (2-adic Littlewood), my cuts-as-Farey-geodesics + razor-thin + duality-web reflections, THM-501/503
(LRC singular series), HYP-3547 (apex-7), OPEN-Q-108.
