---
source: opus-2026-06-03-S587 (remote-control)
status: SYNTHESIS — three tournament parities (Rédei odd H; flip changes even; self-converse odd / chiral even) encode odd/even & add/mult; addition (flip, GF(2)) CERTIFIES multiplication (H odd); the LRC worry-set = the ODD (self-converse) sector
tags: [LRC, parity, redei, hamiltonian-path, arc-flip, self-converse, chiral, worry-set, odd-even, add-mult, n14]
---

# Tournaments encode odd/even and add/mult — and the worry-set is the odd sector

**Prompt (user):** tournaments encode even/odd and multiplication/addition; their
Hamiltonian-path count is always odd; flipping an arc changes it by an even count; in
the merged explorer, chiral merged nodes have an even number of tilings, self-converse
ones an odd number.

All three parities verified (`n=4,5,6`), and together they say something clean: the
*additive* (flip) structure **certifies** the *multiplicative* (Rédei) oddness, and the
*converse* parity splits tournament space into an **odd (self-converse)** sector and an
**even (chiral)** sector — and **the LRC worry-set is the odd sector.**

## 1. The three parities (verified)

| parity | statement | flavour |
|---|---|---|
| **Rédei** | `H = #Hamiltonian paths` is **ODD** for every tournament | multiplicative (OCF / odd cycles) |
| **flip** | one arc flip changes `H` by an **EVEN** amount | additive (GF(2), cut⊕cycle) |
| **converse** | self-converse classes have **ODD** `H`; chiral pairs **merge to EVEN** | the `ℤ/2` time-reversal orbit |

Counts: `n=4`: 4 classes = 2 SC + 1 chiral pair; `n=5`: 12 = 8 SC + 2 pairs; `n=6`:
**56 = 12 SC + 22 chiral pairs** (the user's `56 = 12 + 44`). All `H` odd, all flips
even, all self-converse `H` odd.

## 2. Addition certifies multiplication (the hidden link)

`H mod 2` is **constant on the flip-hypercube**: flipping an arc moves between adjacent
tournaments in `{0,1}^{C(n-1,2)}` and changes `H` by an even amount, so `H mod 2` is a
**GF(2) flip-invariant**. The transitive tournament has `H=1` (one Hamiltonian path),
and the flip-hypercube is connected, so

> **`H ≡ 1 (mod 2)` for every tournament — Rédei's theorem, proved by the additive flip
> structure certifying the multiplicative invariant's oddness.**

This is the precise sense of "tournaments encode addition *and* multiplication": the
GF(2) *addition* (arc flips on the staircase tiling space) is what pins the *odd*
multiplicative invariant `H` (the OCF / odd-cycle product). Last session's bridge
`8·T+1 = (2n-1)²` is the number-theoretic shadow of this — odd invariant (`H`,
`odd²`) carried over an additive staircase (`T = C(n,2)` tiles, flips).

## 3. The converse parity = worry vs loose

The `ℤ/2` time-reversal (converse `T↦T^op`) splits tournament classes:
- **self-converse** (`T≅T^op`): a single node, **odd** `H` — the **spine** of the merged
  metagraph (`G_n/Z_2`).
- **chiral** (`T≇T^op`): merges with its mirror to an **even** node — the **sea**.

And the LRC fact I established earlier: **the worry-set = the self-converse round
classes** (the regular rotational encirclement; the `2^((n-2)/2)=64` classes at `n=14`).
So:

> **The LRC worry-set is the ODD (self-converse, Rédei-spine) sector; the loose configs
> are the EVEN (chiral, merged) sector.** Loneliness is organized by the converse parity
> that tournaments intrinsically carry. The hard, rigid floor lives in the odd sector.

## 4. The unified picture (with last sessions)

```
        ODD                         EVEN
Rédei H (odd, OCF/×)        |   chiral merged (even)
self-converse = WORRY-SET   |   chiral = LOOSE
8T+1 = (2n-1)² (odd-square) |   the even arc-flip changes
        \                   /
   ADDITION (arc flip, GF(2), the staircase T=C(n,2)) CERTIFIES the ODD invariant
   MULTIPLICATION (H, OCF) is the odd invariant
   CONVERSE (ℤ/2) splits odd(worry)/even(loose)
```
- **+2 / additive:** arc flips (even), the flip-set lattice (the 64 = flips, my S570), the
  staircase gnomons.
- **×2 / multiplicative:** `H`/OCF (odd), the 2-adic doubling tower.
- **converse parity:** worry (odd, self-converse) vs loose (even, chiral).
- **worry-modulus `2n-1 = √(8·C(n,2)+1)`** (S586): the odd-square shell of the additive
  pair-count — the worry-set's home, in the odd sector.

## 5. Creative hypothesis

> **H (parity-protected loneliness).** LRC loneliness is a property of the **odd
> (self-converse) parity sector**; the worry-set is exactly the Rédei-spine
> (self-converse, odd `H`), and the conjecture's residual is the odd sector's rigidity
> at the composite seam `2q` (`n=14`). Sharper, speculative: there is an LRC
> **witness-invariant that is always odd** (a Rédei-analog for observer-source
> reachability, THM-381) — *flip-certified* like `H` — whose **oddness forces a witness**
> (an odd count cannot be zero). A parity proof: addition (the flip / clock-shift
> structure) certifies an odd loneliness-count, hence loneliness, with the even (chiral)
> sector free. *Test:* count observer-source classes / witness-orderings with a fixed
> base structure and check the parity is forced odd on the self-converse worry-set.

## 6. Honest status

- **Verified (`n=4,5,6`):** Rédei `H` odd; flip changes `H` even; self-converse `H` odd;
  `56 = 12 + 44`.
- **Rigorous link:** `H mod 2` is a GF(2) flip-invariant `=1` (the additive proof of
  Rédei); the worry-set = the self-converse (odd) sector (from S570/S566).
- **Hypothesis H:** the parity-sector organization of worry/loose is structural; the
  "odd witness-invariant forces loneliness" is a *speculative* parity-proof direction,
  untested — the creative seed.

**Artifacts:** `04-computation/lrc_tournament_parity_s587.py` (+`.out`). Builds on Rédei
(repo core), S570/S566 (worry-set = self-converse), S586 (triangular bridge / 2n-1),
S585 (dynamical rigidity), THM-381 (observer-source). New: **HYP-2129**.
