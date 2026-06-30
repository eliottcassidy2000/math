# Working both covering-min sides through the spectral/Ramanujan lens: the two columns are both the cycle C_* — the MEASURE column (apex gap = 2+λ_min(C_p) = 4sin²(π/2p)) is the cycle's spectral SLACK from the Ramanujan bound −2 (the apex Cayley graphs C_p, Paley_p, K_p are ALL Ramanujan, so their Ihara zeta satisfies the RH analogue: non-trivial poles on |u|=1/√(k−1)), and the EXISTENCE/M column (the corrected covering-min) is the cycle's COVERING RADIUS (even n: 1/n via the even block = equally-spaced AP = C_n vertices; odd n: C_n unreachable by an all-even covering set, so > 1/n) — spectral gap and covering radius, the two faces of the same cycle

*opus-2026-06-30. Owner: work both the odd-n covering-min and the even-n statement back and forth, and
consider "a regular graph is Ramanujan iff its Ihara zeta satisfies the RH analogue." The Ramanujan lens IS
the measure column; the covering-radius lens IS the M column; both are the cycle C_*.*

## The MEASURE column is the cycle's spectral slack — and it's Ramanujan (Ihara-RH verified)
THM-590's apex gap `g(O) = min_{k≠0} |Σ_{x∈O} ω^{kx}|²` is a **Cayley spectral gap**. The apex Cayley graphs
on `Z_7` are all Ramanujan, and their Ihara zeta satisfies the RH analogue (computed):
| apex graph | conn | k-reg | Ramanujan `|λ|≤2√(k−1)` | Ihara-RH poles on `|u|=1/√(k−1)` |
|---|---|---|---|---|
| **C₇** (doublet `±1`) | `{1,6}` | 2 | ✓ (`≤2`) | ✓ all at `1.000` |
| **Paley₇** (QR `{1,2,4}`) | `{1,2,4}` | 3 | ✓ (`≤2.828`) | ✓ all at `0.7071` |
| **K₇** (full) | `{1..6}` | 6 | ✓ (`≤4.472`) | ✓ all at `0.4472` |
> **The apex gap is the cycle's SLACK from the Ramanujan bound:** `g(doublet) = 2 + λ_min(C₇) =
> 2 + (−1.8019) = 0.19806 = 4cos²(3π/7)` (mac-mini HYP-3594). In general `2 + λ_min(C_p) = 2 − 2cos(π/p) =
> **4sin²(π/2p)** = the apex atom`. The Ramanujan/spectral bound is `λ = −2`; `C_p` sits at `λ_min =
> −2cos(π/p) > −2`, **above the bound by exactly the apex atom** `4sin²(π/2p)`. As `p→∞` the cycle approaches
> Ramanujan-tightness and the slack `→0` — the apex atom shrinking is the cycle **approaching the spectral
> bound**, and the genus of `X₀(2p)` growing. The genus-1 case (`p=7`, LRC14) is the slack `0.198`.

## The EXISTENCE/M column is the cycle's covering radius (even/odd)
The corrected covering-min is the **covering radius** of the same cycle:
> **EVEN n: covering-min = `1/n`** — the even block `2·{1,…,n−1}` at witness `t=1/(2n)` sends the speeds to
> `{1/n, 2/n, …, (n−1)/n}` = the **n-th roots = the vertices of `C_n`, equally spaced**, covering radius
> `1/n`. The even block IS the cycle `C_n` realized as a covering set (verified n=8,14, witness `q=2n`).
> **ODD n: `C_n` is UNREACHABLE by an all-even covering set** (covering `q=n` needs an odd multiple of odd
> `n`, breaking the doubling), so the covering radius is `> 1/n` (n=7: `2/13`; the realizability problem).
So the M column is the cycle's covering radius, attained (`=1/n`) exactly when the equally-spaced `C_n` is
coverable — which is the parity condition `n` even.

## The two columns are two faces of the cycle C_*
| | MEASURE (Q(√−7), apex) | EXISTENCE/M (covering-min) |
|---|---|---|
| object | `C_p` spectral SLACK from Ramanujan bound | `C_n` COVERING RADIUS |
| value | `2+λ_min(C_p) = 4sin²(π/2p)` (apex atom) | `1/n` (even n); `>1/n` (odd n) |
| framework | **Ramanujan / Ihara-RH** (apex Cayley graphs) | optimal covering (equally-spaced AP) |
| optimum | Ramanujan = optimal expander | equally-spaced = optimal covering |
| parity | `p` = the apex prime (here 7) | `n` even → reachable; odd → not |
> **Both columns are the cycle.** The measure column reads its SPECTRUM (the slack from the Ramanujan bound,
> the Ihara-RH); the M column reads its COVERING RADIUS (the equally-spaced gap). Spectral gap and covering
> radius are the dual extremal properties of a graph — expansion vs. covering — and here they are the two
> Heegner columns of LRC(14): the `C_7` spectral slack (`Q(√−7)`, `0.198`) and the `C_14`/`C_n` covering
> radius (`1/n`). The Ramanujan/Ihara-RH the owner pointed at is exactly the spectral half.

## Poking the connections (past topics, re-read)
- **THM-590 = a Ramanujan slack.** The apex gap `0.198` is now a spectral-extremality statement: `C_7` is
  Ramanujan, and `0.198 = 2+λ_min` is its slack from the bound. The "genus-1" is the cycle's spectral
  signature; the Ihara-RH holds for the apex graphs.
- **The corrected covering-min = the C_n covering radius**, parity-gated (even n reaches it, odd n doesn't).
  The earlier `ζ₆`/`Φ₆` construction was a red herring; the real object is the equally-spaced `C_n`.
- **The metazeta** (Ihara zeta of the metagraph) is the analytic spectrum on the *whole* metagraph (not
  regular, so not Ramanujan as a whole) — but its apex sub-structure (the `C_p` Cayley pieces) IS Ramanujan;
  the metazeta's hardest poles live there.
- **The observer `1`** is the trivial pole/eigenvalue `k` (the Ramanujan "trivial" eigenvalue), above which
  the spectrum (the prime cycles / the H-gradient) sits — the same baseline `1` as before.

## Status
- **Computed/verified (opus):** apex Cayley graphs `C_7, Paley_7, K_7` are Ramanujan with Ihara-RH (poles on
  `1/√(k−1)`); apex gap `= 2+λ_min(C_p) = 4sin²(π/2p)` = the cycle's slack from `−2`; even-n covering-min
  `=1/n` is the `C_n` covering radius (even block = equally-spaced AP, witness `t=1/(2n)`); odd-n `>1/n`
  (`C_n` unreachable).
- **The synthesis:** both columns are the cycle `C_*` — MEASURE = spectral slack (Ramanujan/Ihara-RH),
  M = covering radius (parity-gated). The Ramanujan thread is the measure column's framework.
- **Open (both):** even-n — prove the equally-spaced `C_n` (even block) is the worst covering set
  (`no covering set < 1/n`, = LRC tight); odd-n — the covering-radius realizability (`2/13`, … the
  smallest reachable radius when `C_n` is blocked by parity).

Related: CORRECTION-…even-block (the covering-min = 1/n), the-metazeta-ihara (the Ihara framework),
the-functional-decomposition (parity/doubling); THM-590 (apex gap), mac-mini HYP-3594 (2+λ_min(C_p)),
HYP-3547 (apex primes); OPEN-Q-108.
