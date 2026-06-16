# The OCF is a Mayer gas, LRC-14 is its conditional twin, and E₇ is its Fano cousin — with one prime (7 = Φ₃(2) = |PG(2,2)|) gating the family

**Source:** kind-pasteur-2026-06-15-S7 (a develop→adversarial-verify→synthesize workflow,
9 agents). Dispatch: an explicit Mayer cluster expansion, and how it relates to the H=21,7
impossibility, E₇, and LRC-14. The honest result is a **four-layer picture of decreasing
rigor**, held together by one prime.

## Layer 1 — the OCF gas (rigorous, THM-517)

`H(T)=I(Ω,2)=Σ_k α_k 2^k` is *literally* the grand partition function of a hard-core lattice
gas on the odd-cycle conflict graph `Ω` at fugacity `z=2` (THM-002). Its cluster expansion is
now explicit: ideal gas `H_free=3^{α₁}`; the Ursell excluded-volume integral is `−|E(Ω)|`
(conflicting odd-cycle pairs), the 3rd is `P₂(Ω)−t(Ω)`. **The correction that mattered**
(adversarial pass, and my own independent cross-check both caught it): the clean `−|E(Ω)|` is
the *graphical* excluded-volume integral `α_2−C(α₁,2)`, **not** the order-2 cumulant of
`log I` — that is `c_2 = −|E|−α₁/2`. Only the analytic cumulants `c_k` satisfy
`exp(Σc_k z^k)=I` (the graphical `b_k` fail; `C₄` is the counterexample). And `z=2` is far
outside the radius of convergence (`R≈0.0125` at Paley `T₇`), so the gas/cluster expansion is
a **formal/resummed** identity, dovetailing with the canonical fact that `λ=2` is outside
hard-core uniqueness. This is the rigorous spine; everything else hangs off it.

## Layer 2 — the forbidden values are excluded configurations

`7 = I(K₃,2) = Φ₃(2) = 2³−1` is the **unique** realizable cluster profile `α=(1,3,0)` — the
maximal-excluded-volume 3-cluster — and THM-029 excludes exactly it (`Ω=K₃` non-realizable).
So `7` is forbidden by **single-cluster rigidity**: its only realizing conflict graph is
outside the realizable cone. `21` is **not** rigid (four realizable profiles); it is forbidden
**multiplicatively** (`H=∏` over strong components, `21=3·7` routes through the forbidden `7`).
The all-`n` non-achievability is the already-closed Busch/A005517 fact (`strong-min(m)=
3,5,9,15,25,45,75,…` strictly increasing; HYP-2271, monad-compute-2026-06-06) — the cluster
picture **cleanly re-explains `7`** and **re-describes `21`**, but the surjectivity half (every
odd `≥23` achievable) stays open and is provably immune to cluster/multiplicative methods
(a *prime* `H` needs a strong tournament of that `H`).

## Layer 3 — LRC-14 is the conditional twin (HYP-2544)

The LRC(14) singular series `L(S) = (6/7)^13 + Σ_{exact relations T} (6/7)^{13−|T|}(−1)^{|T|}
∏_{v∈T} s(t_v)` is, **term for term**, a Mayer/polymer expansion of a *gas of resonance
relations*: ideal term `(6/7)^13`, each exact relation `T` a cluster, `(−1)^{|T|}` the Mayer
sign, `∏ s(t_v)` the cluster activity, `|T|` the size. (The `(A−B)`-product factorization is
an exact inclusion–exclusion — the same `(−1)^{|S|}` Mayer reading already canon for the SC
formula.) The genuine **difference is the convergence regime**: the OCF gas is finite
(absolutely convergent — though `z=2` is past its radius), while the LRC gas is an *infinite*
relation lattice — `|T|=2` sits at the absolute-convergence boundary (almost-Sidon, THM-503),
`|T|≥3` is only conditionally convergent via cross-level `(−1)^{|T|}` alternation (THM-504).
Reading `inf L>0` (`C′(14)`) as "convergent cluster expansion with positive pressure" is a
**useful frame, not a theorem**: it correctly names the tool family (Kotecký–Preiss,
Fernández–Procacci) and **diagnoses why they fail** — they require absolute summability, which
LRC lacks (`A₃` diverges at config level). The honest next tool is a *conditional/oscillatory*
polymer bound (Abel summation across the `|T|`-filtration + Pólya–Vinogradov per level,
THM-504 D).

## Layer 4 — E₇ is the Fano cousin (HYP-2546), not a mechanism

`7 = Φ₃(2) = 2³−1 = |PG(2,2)|` (the Fano plane: 7 points, 7 lines). The **one genuine bridge**:
the same `PG(2,2)` builds the octonions (hence E₇ via the Freudenthal magic square) **and**
organizes E₇'s quartic Cartan invariant — Duff–Ferrara (peer-reviewed, PRD 76 025018): `I₄` is
built from the **7 Fano-line** Cayley hyperdeterminants, `56=7·8`, `SL(2)⁷⊂E₇`. And `7` is the
OCF forbidden value as `|PG(2,2)|=I(K₃,2)`. So OCF-`7` and E₇-`7` share a real combinatorial
ancestor. **But the `7` occupies three different semantic slots:** forbidden **value** (OCF),
structural **count** (E₇), **modulus** (LRC, `s(t)=0⟺7|t`). The "exceptional cluster gas"
(OCF `Σα_k2^k` vs E₇ `I₄=Σ` over 7 Fano lines) is an **analogy** — both multilinear invariants
summed over a combinatorial index set — with **no specialization map**. Many E₇ 7-divisibilities
(`133=7·19`, `63=9·7`, `126=18·7`, `7∣|W|`) are partly **forced by rank 7 alone**. Verdict: E₇
is a genuine **structural cousin via the Fano substrate**, a **numerological lure if read as
mechanism** — keep it decorative, never load-bearing.

> **Retraction (HYP-2530).** The earlier prediction that the E₇ 56-rep symplectic form's
> eigenvalues equal the exponents `{1,5,7,9,11,13,17}/h` is **wrong**: it conflates the
> Coxeter-element spectrum on the rank-7 Cartan with the *invariant* symplectic form `ω` on the
> minuscule rep (which, by Darboux/invariance, is just 28 hyperbolic `±1` blocks — no
> E₇-covariant spectrum). The genuine 56-rep structure is the `56=7×8` Fano decomposition.
> (Also struck: a stray "63=C(9,2)" — `C(9,2)=36`; only `63=9·7` is correct.)

## The prime-7 thread, stated honestly

| slot | role of 7 | same as Φ₃(2)=\|PG(2,2)\|? | mechanism |
|---|---|---|---|
| OCF | forbidden **value** `I(K₃,2)` | YES (literal) | realizability (K₃ non-realizable) |
| E₇ | structural **count** (rank, Fano lines, 7·8) | YES (Fano substrate, Duff–Ferrara) | counting / index set |
| LRC | **modulus** `s(t)=0⟺7\|t` | same prime, **no Fano shown** | weight-zero gate |

OCF and E₇ are the **same 7** (the Fano/cyclotomic ancestor). LRC routes through the same prime
with a **different mechanism** and no demonstrated Fano structure — a resonance, not an identity.

## The one-line synthesis

> One rigorous gas (the OCF, `H=I(Ω,2)`), its faithful cluster reading (excluded-volume `|E(Ω)|`,
> formal at `z=2`), and its excluded configurations (`7`=unique cluster, `21`=multiplicative);
> a genuine **conditional-convergence twin** (LRC-14's singular series); and a genuine but
> **unbridged Fano cousin** (E₇) — with the apex prime `7=Φ₃(2)=|PG(2,2)|` gating the family by
> realizability (OCF), counting (E₇), and weight-zero (LRC).

Cross-links: THM-517 (renumbered off a collision with mac-mini-S4's THM-515), THM-002/029/510, HYP-2271 (Busch/A005517, the closed non-achievability),
HYP-2544/2545/2546, HYP-2530 (retracted prediction), THM-501/503/504 (the LRC series),
[[the-burnside-core-kernel-is-a-gcd-quadratic-form-and-a-metagraph-enumerator-family-kps]]
(today's other thread: "odd parts ≥3" cores = A000009−3 = the OCF non-spectral family),
sources: Busch EJC 13 #N3; Duff–Ferrara PRD 76 025018 (arXiv:quant-ph/0609227); Baez,
*The Octonions* (arXiv:math/0105155).
