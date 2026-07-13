# The worst |core|=1 body {1..11,13,84} is the crack between opus's Arguments A and B — flag + refinement

*kind-pasteur-2026-07-11-S127 cont.70. Owner: "calibrate the tight large-sieve against {1..11,13,84}." opus-S268
correctly localized the irreducible core to `|core|=1` (the runner-1 lemma, S265), reporting the `|core|=1`
energy `≤ 0.328`. But my cont.69 worst case is `0.60`, at `{1..11,13,84}` — a body opus's sample **missed**.
Calibrating there shows it is the exact crack between the runner-1 lemma's two arguments: Argument A (measure)
fails and Argument B (equidistribution) is stressed, yet LRC holds by a thin margin. This flags a concrete
verification gap and points at the refinement.*

---

## opus's runner-1 lemma (S265), and the missed body

For a `|core|=1` (speed-1) covering family, LRC(14) ⟺ `S_rest ⊄ D_1` (the rest-safe set is not inside runner-1's
arc `D_1 = {‖t‖ < 1/14}`), proved by two complementary arguments:

- **Argument A (measure):** the rest has a small even speed `s`, `|S_rest ∩ D_1| ≤ |S_s ∩ D_1| = (s−1)/(7s)`; so
  `|S_rest| > (s−1)/(7s)` ⟹ `S_rest ⊄ D_1`. Best at `s=2`: needs **`|S_rest| > 1/14`**. (Covers near-AP rests.)
- **Argument B (equidistribution):** `S_rest ⊄ D_1 ⟺ ε_1 < 6/7`, and `ε_1` is governed by the count of
  **consecutive-difference pairs** `1 = w_i − w_j`; a **spread** rest (few pairs) ⟹ `ε_1` small. (Covers spread
  rests.)

opus-S268's `|core|=1` energy maxed at `0.328` — but the true worst is `core·Σε² = 0.60` at `{1..11,13,84}`
(cont.69), so **that family was not in opus's sample, and A∪B was never checked there.**

## Calibration: {1..11,13,84} fails A and stresses B

Rest `= {2,…,11, 13, 84}` (the speeds `≠ 1`). Measured (`grid = 2·10⁶`):

| quantity | value | verdict |
|---|---|---|
| `|S_rest| = |G'|` | **0.0666** | `< 1/14 = 0.0714` ⟹ **Argument A FAILS** |
| consecutive pairs `1 = w_i−w_j` | **9** (the run `2..11`) | long run ⟹ **Argument B STRESSED** (not spread) |
| `coreCover` | 0.920 | |
| `ε_1 = coreCover − 1/7` | **0.777** | `< 6/7 = 0.857` ⟹ **LRC HOLDS** (margin 0.08) |
| surplus `|S_rest| − |S_rest ∩ D_1|` | **0.0054** | `> 0` ⟹ `S_rest ⊄ D_1` ✓ |

So `{1..11,13,84}` sits exactly in the crack: it is **too concentrated for A** (`|S_rest| = 0.0666 < 1/14`, the
`s=2` threshold) and **too clustered for B** (a length-10 run `2..11`, the opposite of spread). It is neither
near-AP-enough for A nor spread-enough for B — yet LRC holds, thinly. This is why it is the worst `|core|=1`
body and why opus's two-costume split does not obviously cover it as stated.

## Why LRC still holds — and the refinement direction

The truth is milder than A's bound: `|S_rest ∩ D_1| = coreCover·|S_rest| = 0.920·0.0666 = 0.0613 < |S_rest| =
0.0666`. So the *tight* measure statement `|S_rest ∩ D_1| < |S_rest|` **does** hold — Argument A is merely
**too loose**, bounding `|S_rest ∩ D_1| ≤ 1/14` when it is actually `0.0613`. The fix is a **multi-runner
measure bound**: instead of `|S_rest ∩ D_1| ≤ |S_2 ∩ D_1| = 1/14`, use two (or more) small runners,

> `|S_rest ∩ D_1| ≤ |S_2 ∩ S_4 ∩ D_1| < 1/14`,

which is strictly smaller (each extra runner carves the arc further) and — for a rest with a long low run like
`2,…,11` — drops below `|S_rest|`. Concretely: the run `2..11` that stresses Argument B is exactly the surplus
of *small even/near-even runners* that a **refined Argument A** can exploit. So the two arguments are not just
complementary at the endpoints; on this crack body they **swap roles** — the long run (bad for B) is the raw
material for a tightened A. That is the natural closure: Argument A with the `k` smallest runners of the rest,
threshold `|S_{2}∩…∩S_{k}∩D_1|`, which shrinks with `k`.

## Net — the concrete ask for opus

`{1..11,13,84}` is the true worst `|core|=1` body (`core·Σε² = 0.60`, not `0.328`), and it is the crack between
Arguments A and B: A's `|S_rest| > 1/14` fails and B's spread hypothesis is violated by the run `2..11`. LRC
holds (`ε_1 = 0.777 < 6/7`, margin `0.08`), so the lemma is *true* here — but opus's A∪B, as stated and as
sampled, does not verifiably cover it. The refinement that closes it is a **multi-runner Argument A** (carve
`D_1` with the several smallest runners, not just `s=2`), which turns the very long-run structure that defeats B
into the tool that rescues A. Calibrating the runner-1 lemma against this body — not the deep well — is the
concrete step.

*Files: lrc14_runner1_crack_kps_S127.py-style checks in lrc14_runner1_crack_kps_S127.out. Supports opus-S265
(runner-1 lemma A/B), S268 (|core|=1 localization); builds on kps cont.69 (worst body), cont.62/63 (coreCover).
HYP-6248.*
