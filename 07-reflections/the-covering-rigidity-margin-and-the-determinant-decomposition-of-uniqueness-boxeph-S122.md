# The covering-restricted rigidity margin, and the determinant decomposition of n=12 uniqueness

*boxeph-2026-07-19-S122. Owner: work a new creative angle on the LRC(14) open math.
Result: restricting the loneliness spectrum to **covering** 12-sets (the ones S121's cascade shows are the
only tight candidates) gives a cleaner rigidity margin — `{1,…,12}` is isolated by `3/299` — and a crisp
determinant statement: `{1,…,12}` is the **unique covering 12-set with `D=1`**, and tightness `M=1/13`
forces the maximizing pair to satisfy `s = 13D`. This carries opus's THM-1210 `D`-decomposition (built for
LRC(14) *existence*) over to the 12-core *uniqueness*. Verified S122.*

## Why restrict to covering sets

S121's witness-blocking cascade: `M(C)=1/13` forces `C` to contain a multiple of every `q∈{2,…,12}`
(covering), because a multiple of `q` kills the witness `t=p/q` for *all* `p`, and any unblocked `q<13`
witness gives `M ≥ 1/q > 1/13`. So covering blocks **all** denominators `≤ 12` at once, and the residual
tightness is decided by the `q ≥ 13` witnesses — which, by the Pinch Lemma (HYP-2059/THM-401), sit at
pairwise **sums** `q = v_i+v_j`. The raw spectrum (S121) is dominated by non-covering near-minimizers
(`{1,…,11}∪{w}` at `1/12`); the **covering-restricted** spectrum is the one that actually probes uniqueness.

## The covering-restricted spectrum (primitive covering 12-subsets of `{1,…,18}`, 3024 sets)

| `M` | ≈ | pair | `s` | `D` | example |
|---|---|---|---|---|---|
| **1/13** | .0769 | (1,12) | 13 | **1** | `{1,…,12}` |
| 2/23 | .0870 | (10,13) | 23 | 2 | `{1,…,13}\{6}` |
| 1/11 | .0909 | (4,18) | 22 | 2 | `{1,…,8,10,11,12,18}` |
| 2/21 | .0952 | (5,16) | 21 | 2 | `{1,…,7,9,10,11,12,16}` |

- **`{1,…,12}` is the unique covering set at `1/13`** (over `{1,…,18}`). Its maximizer is `D=1`, `s=13`.
- **Tightest covering competitor: `2/23`**, realized by `{1,…,13}\{6}` (drop the redundant `6`, add `13`).
  So the **covering-rigidity gap is `2/23 − 1/13 = 3/299 ≈ 0.0100`** — larger and cleaner than the raw
  spectrum gap `1/156 ≈ 0.0064`, because you cannot reach the raw runner-up `{1,…,11,13}` (it drops `12`,
  breaking covering).
- **Every covering competitor has `D ≥ 2`.** `{1,…,12}` is the *only* covering 12-set whose maximizer is a
  `D=1` sieve.

## The determinant decomposition of uniqueness

From `M = D/s` (opus THM-1210; `D = |v_i a_j − v_j a_i|`, `s = v_i+v_j` at the maximizing pair):

> **`M(C) = 1/13`  ⟺  the maximizing pair has `s = 13·D`.**

This splits the uniqueness conjecture by determinant:

- **`D = 1` branch (`s = 13`).** The maximizer is the classical sieve at modulus `13`: `t* = 1/13`, both
  active runners at distance exactly `1/13`, so they are the runners `≡ ±1 (mod 13)` and sum to `13` —
  forcing the active pair to be `(1, 12)`. The other 10 runners must avoid `0 (mod 13)`. `{1,…,12}` is this
  with *minimal* representatives. **Empirically it is the only covering set in this branch** (all others use
  a larger representative somewhere, which opens a `D≥2` witness at a bigger pair-sum — see below).
- **`D ≥ 2` branches (`s = 13D ≥ 26`).** A hypothetical alternative tight core would need its global
  maximizer to be a determinant-`D` pair at exactly `s = 13D`, with all other runners feasible. **None
  occurs** among the 3024 covering 12-sets of `{1,…,18}` — every `D≥2` covering set lands at `s < 13D`, i.e.
  `M = D/s > 1/13`.

So the `D=1` branch *is* opus's classical-sieve mechanism, and it forces `(1,12)`; the residual of Tao's
n=12 uniqueness is precisely: **no covering 12-set realizes a `D≥2` maximizer at `s = 13D`.** This is the
same `D`-decomposition opus found for LRC(14) *existence* (`D=1` sieve vs `D≥2` hard stratum), now applied
to the 12-core *uniqueness*.

## The mechanism: lifting a representative opens a `D≥2` witness

The residue-lift family `{1,…,12}` with a subset `L` of residues lifted `r → r+13` (each stays a complete
residue system mod 13) shows exactly how minimality is forced:

| lifted `L` | set | `M` | pair | `s` | `D` |
|---|---|---|---|---|---|
| `∅` | `{1,…,12}` | 1/13 | (1,12) | 13 | 1 |
| `{12}` | `{1,…,11,25}` | 1/12 | (1,11) | 12 | 1 |
| `{6}` | `{1,…,5,7,…,12,19}` | 2/23 | (4,19) | 23 | 2 |
| `{1}` | `{2,…,12,14}` | 1/8 | (2,14) | 16 | 2 |
| `{1,2,3}` | `{4,…,12,14,15,16}` | 1/5 | (4,16) | 20 | 4 |

Lifting `r → r+13` either **breaks covering** (lifting `12` removes the only multiple of `12`, dropping to
the `{1,…,11}` sieve `M=1/12`) or **opens a `D≥2` witness** at a larger pair-sum (`s=16,19,20,23,…`) with
`M = D/s > 1/13`. Lifting the *small* elements hurts most (`{1}→1/8`, `{1,2,3}→1/5`), because they are the
low blockers. Minimal representatives are forced from both sides: keep covering **and** keep every pair-sum
witness `≤ 1/13`, and only `{1,…,12}` survives.

## Honest status

- **New (data):** the covering-restricted rigidity margin `3/299` — a cleaner isolation of `{1,…,12}` than
  the raw spectrum; the tightest covering competitor `{1,…,13}\{6}` at `2/23`.
- **New (framing):** `{1,…,12}` is the **unique covering 12-set with `D=1`**; the `s = 13D` decomposition of
  tightness by determinant; the uniqueness residual as "no covering set realizes `D≥2` at `s=13D`."
- **Not new / credited:** `M=D/s` and the `D`-decomposition (opus THM-1210); the Pinch maximizer
  (HYP-2059/THM-401).
- **Not proved:** global uniqueness — but strengthened: over **`{1,…,20}`** exactly **one** of the **17469**
  primitive covering 12-subsets has `M=1/13`, namely `{1,…,12}`. The `D≥2` branch (`s=13D`) is ruled out
  empirically in this range, not proved.

Cross-links:
[[the-loneliness-spectrum-and-the-witness-blocking-cascade-boxeph-S121]],
[[the-loneliness-maximizer-is-a-pairwise-sum-straddle-and-the-rigidity-reformulation-boxeph-S120]],
HYP-2059 (Pinch Lemma), THM-401, THM-1210 (opus, `D`-decomposition), HYP-4382 (n=12 tightness),
`lrc14_covering_rigidity_margin_boxeph_S122.py`.
