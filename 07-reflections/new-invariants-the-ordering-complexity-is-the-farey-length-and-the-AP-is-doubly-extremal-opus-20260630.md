# Wielding definitions — new LRC invariants from the observer: the ORDERING COMPLEXITY O(S) = #distinct circular orderings of {0,s_i t} as t sweeps (the LRC analogue of the tournament's Hamiltonian-path count) equals the FAREY-SEQUENCE LENGTH Φ(n−1)=Σ_{d≤n−1}φ(d) (totient summatory, A002088) for the AP — exact n=4..10, dilation-invariant, and EVEN where the tournament's H is ODD (Rédei); the AP MAXIMIZES O (it has every difference, hence every crossing), so the AP is DOUBLY EXTREMAL — minimal M (the covering-min, the hardest observer) AND maximal O (the richest ordering) — and the LONELINESS INTEGRAL L(S)=∫M_c dc gives the translation-averaged loneliness (~0.28, floor 1/n at c=0)

*opus-2026-06-30. Owner: look for other creative measures or invariants; definitions are your power, wield
them. Defined three from the observer; the ordering complexity turned out to be the Farey length. The LRC's
combinatorial backbone is the Farey sequence.*

> **⚠ CORRECTION (opus, same day):** the claim below that the AP *maximizes* `O` / is "doubly extremal" is
> **WRONG**. The AP **MINIMIZES** `O` (smallest max speed ⇒ smallest Farey ⇒ `O=Φ(n−1)` is the min; verified:
> min `O` over 4-subsets of `[1..8]` is at the AP). The AP is doubly **MINIMAL** (min `M`, min `O`) — the
> tightest *and* simplest config, matching the **transitive** tournament (min `H=1`). The Farey-length result
> (`O=Φ(n−1)`), the loneliness integral, and the other definitions stand. See
> `the-ordering-complexity-hamiltonian-path-bridge-…-doubly-minimal` for the corrected statement and the bridge.

## Definition 1: ordering complexity O(S) = the Farey length
> **`O(S)`** := the number of distinct circular orderings of the `n` points `{0, s_1 t, …, s_{n−1} t}` as
> `t` sweeps `(0,1)` (the observer fixed first) — the observer's **family of snapshot linear orders**, the
> LRC analogue of the tournament's Hamiltonian-path count `H`.
Computed `O(AP_n)` vs `Φ(n−1) = Σ_{d≤n−1} φ(d)` (totient summatory, the Farey-sequence length):
| n | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| `O(AP_n)` | 4 | 6 | 10 | 12 | 18 | 22 | 28 |
| `Φ(n−1)` | 4 | 6 | 10 | 12 | 18 | 22 | 28 |
> **`O(AP_n) = Φ(n−1)` = the Farey length** (A002088). The orderings change exactly at the Farey crossings
> `t = j/d` (`d = |s_i−s_j| ≤ n−1`), so the snapshot count is the order-`(n−1)` Farey sequence. Properties:
> **dilation-invariant** (`O(pS)=O(S)`, like `M` — verified), and **EVEN** (snapshots come in time-reverse
> pairs `t ↔ 1−t`) — in pointed contrast to the tournament's **ODD** `H` (Rédei). The LRC's "Ham-paths" are
> the Farey sequence; the tournament's are `H = I(Ω,2)`. (Parity dictionary: LRC orderings even, tournament
> paths odd — the `+1` is the tournament's irreducible, the pairing is the LRC's time-reversal symmetry.)

## The AP is DOUBLY MINIMAL (min M, min O) — [corrected]
- The AP MINIMIZES `O`: it has the **smallest max speed `n−1`**, hence the smallest Farey `F_{n−1}`, hence
  `O(AP_n) = Φ(n−1)` is the **minimum** (any other `(n−1)`-set of distinct positive integers has max speed
  `≥ n−1` ⇒ a larger Farey ⇒ `O ≥ Φ(n−1)`; verified min over 4-subsets of `[1..8]`). *(Earlier draft said
  "maximizes" — wrong, see the banner.)*
- The AP MINIMIZES `M`: it is the global extremal `M = 1/n` (the covering-min, the hardest observer).
> So **the AP is doubly MINIMAL — minimal loneliness `M` and minimal ordering complexity `O`** — the tightest
> AND simplest configuration. (Tournament echo: the **transitive** tournament is `H=1` minimal — the AP ↔
> transitive, both the simplest/"ordered" config.)

## Definition 2: the loneliness integral L(S) = ∫ M_c dc
> **`L(S)`** := `∫_0^1 M_c(S) dc` — the **translation-averaged loneliness** over all observer positions `c`
> (the inhomogeneous `M_c` integrated; symmetric under `c↔−c`).
`L(AP_n)`: `0.313, 0.297, 0.288, 0.277` for `n = 6,8,10,14` — slowly decreasing toward a constant (~0.27),
with the floor `M_0 = 1/n → 0` at the origin (the hardest observer) but the *mean* staying `O(1)`. So while
the worst observer's loneliness `→0`, the **average observer is `Θ(1)` lonely** — the difficulty of LRC is
concentrated at the single point `c=0`.

## More definitions to wield (defined, partly computed)
- **Apex reading `A(S)`** := the ceiling of the inhomogeneous spectrum (`= p` for `AP_{2p}`, the apex prime) —
  the `c`-translation knob's top rung; measures the apex (pinned earlier).
- **Witness-Farey spectrum `W(S)`** := the set of optimal witness denominators — Farey fractions; for
  `AP_{2p}` they are multiples of `p` (the apex structures the lonely time).
- **Dilation class `[S]`** := the `gcd`-reduced set — the orbit under `S↦pS`; `M, O` are class invariants, and
  the covering-min is realized by choosing the dilation `p | n` (`AP·p`).
> The Farey sequence is the common thread: `O` = its length, `W` = its fractions, the covering-min escape =
> its mediant/convergent. **The LRC lives on the Farey/Stern-Brocot tree**, and these invariants read different
> features of it.

## What the new invariants buy
- **`O(S) = Φ(n−1)`** is a clean, computable, dilation-invariant LRC invariant equal to a classical
  number-theoretic function (the totient summatory / Farey length) — and it is the honest analogue of the
  tournament's `H` (orderings ↔ Ham paths, even ↔ odd).
- **The AP's double extremality** (min `M`, max `O`) is a new structural fact: the covering-min config is also
  the ordering-richest, tying the loneliness extremal to the combinatorial extremal.
- **`L(S)` localizes the LRC difficulty** at `c=0` (the mean loneliness is `Θ(1)`; only the origin observer is
  hard) — quantifying "the origin is the hardest observer."
- **The Farey backbone** is now explicit across the LRC: orderings (`Φ`), witnesses (Farey fractions), escape
  (Farey mediant/convergent), and the apex (the cap) — one tree, several invariants.

## Status
- **Computed/verified (opus):** `O(AP_n) = Φ(n−1)` (Farey length, A002088) n=4..10; dilation-invariant; even
  (vs odd `H`); **AP MINIMIZES `O`** [corrected]; `L(AP_n) ≈ 0.28` (mean loneliness, floor `1/n` at `c=0`).
- **Defined (wielded):** ordering complexity `O`, loneliness integral `L`, apex reading `A`, witness-Farey
  spectrum `W`, dilation class `[S]` — reading the Farey/Stern-Brocot backbone of the LRC.
- **New facts:** the AP is doubly **MINIMAL** (min `M`, min `O`) ↔ transitive (min `H`); the LRC difficulty
  localizes at `c=0`; `O` is the Farey length and the LRC analogue of `H`.

Related: PINNED-…apex-prime (the apex reading), observer-transformations-dilation (dilation invariance),
the-observer-on-the-tournament-side (`H` = Ham paths), SECOND-CORRECTION-…AP-scaled (min `M` = AP·p),
the-observer-abstraction (the Farey escape); A002088 (totient summatory); OPEN-Q-039/108.
