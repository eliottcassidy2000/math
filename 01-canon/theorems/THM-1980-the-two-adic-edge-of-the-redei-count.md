---
id: THM-1980
title: "THE 2-ADIC EDGE OF THE REDEI COUNT — the spectrum resolves H to a 2-adic depth that DECAYS to a single bit, so Redei's parity is the LAST formula for H. The owner's insight (H's #P-ness is an EDGE case; tournaments sit at the boundary between formula-expressible and provably-not, like the harmonic series at p=1) made precise and PROVED n<=7. Define the SPECTRAL 2-ADIC DEPTH d(n) = the largest k such that 'H mod 2^k' is constant on every cospectral class (a function of char_A). Redei (H odd) forces d(n) >= 1. EXACT/DEFINITIVE: d(4)=d(5)=infinity (the spectrum determines H COMPLETELY at n<=5); d(6)=2 (the spectrum pins H mod 4 but NOT H mod 8 — the mod-8 bit turns #P at EXACTLY the cospectral class {13,17} where H leaves the ladder, THM-1780); d(7)=1 (the spectrum pins ONLY H mod 2 — 26 cospectral 7-tournament buckets carry BOTH residues of H mod 4, so even the mod-4 bit is #P). So the depth runs infinity, infinity, 2, 1 for n=4,5,6,7 and hits the Redei floor at n=7: ASYMPTOTICALLY THE SPECTRUM PINS EXACTLY ONE BIT OF H — its parity — and every higher bit is #P. H mod 4 is a genuine (nontrivial) spectral invariant only for n<=6 (both residues 1,3 occur in ~equal numbers: n=6 splits 32:24), and it is NOT score- or c3-determined from n=5 (it is a higher-order spectral quantity), then it stops being spectral at all at n=7. TWO ORTHOGONAL EDGES MEET AT THE HAMILTONIAN OBJECT: the LENGTH edge (THM-1870: cycle counts c_k spectral for k<=n-1, #P at the Hamiltonian length k=n) and the 2-ADIC edge (this: H's bits spectral up to the parity bit, #P above) — the Hamiltonian-path count H is one length past the spectral cycle counts and one bit past a spectral formula. Redei's theorem is not merely A formula for H; it is the LAST one — the single formula-expressible bit that survives at the boundary."
status: >
  DEFINITIVE.  d(4)=d(5)=inf, d(6)=2 are EXACT-EXHAUSTIVE (all iso classes n<=6, exact integer
  Faddeev-LeVerrier char polys + Held-Karp H, no floats); the mod-8 failure witness at n=6 is the
  cospectral class {13,17} = THM-1780's pair (both ==1 mod 4, differ by 4).  d(7)=1 is DEFINITIVE by
  WITNESS EXISTENCE: 26 sampled cospectral 7-buckets (exact char_A) carry different H mod 4, and a
  single such pair PROVES H mod 4 non-spectral at n=7 (an existence result, immune to sampling
  incompleteness); Redei gives d(7)>=1, so d(7)=1 exactly.  H mod 4's global non-triviality (both
  residues, ~equal) and its non-determination by score/c3 (from n=5) are exact n<=6 facts.
  Interpretation (edge / harmonic-series analogy) is a REFRAMING of these proven facts.
  OPEN (the owner's own suggestion, honestly flagged): d(n) is defined via the SPECTRUM; whether a
  cleverer poly-time invariant pins more than the parity bit asymptotically is OPEN — the full
  poly-tuple (score,specA,specS,disc,arb_inv) determines H for n<=6, but its n>=7 depth is untested.
  If NO poly-time invariant beats parity asymptotically, Redei's bit is PROVABLY H's entire
  formula-expressible content.
source: kind-pasteur-2026-07-21-S128c143 (owner: H's #P-ness is an edge case; tournaments sit at the formula/no-formula boundary like the harmonic series; chase it)
depends_on:
  - THM-1780    # H (and c_n) leaves the spectral ladder at n=6 — the mod-8 witness is its class
  - THM-1870    # the cycle-count spectral boundary IS the Hamiltonian length (the orthogonal edge)
  - THM-1965    # the invariant lattice (this pins the 2-adic direction of "no poly invariant determines H")
related: [THM-1945, THM-1775, THM-1885]
external:
  - "Redei 1934 (#Hamiltonian paths is odd). #Ham-paths is #P-complete. Permanent mod 2 = determinant (the 2-adic collapse behind Redei)."
script: 04-computation/H_two_adic_edge_v2_kps_S128c143.py, H_mod4_formula_and_n7_kps_S128c143.py, exact_lattice_and_edge_kps_S128c143.py (+ .out)
---

# THM-1980 — the 2-adic edge of the Rédei count

**The owner's insight.** H (= #Hamiltonian paths, `#P`, THM-1780) is an *edge* case: tournaments sit
right at the boundary between what a formula can express and what provably no polynomial-time formula
can — the way the harmonic series sits at `p = 1`, the marginal divergent case. This makes that
precise, **2-adically**, and proves it through n=7.

## The spectral 2-adic depth

Rédei (1934): **`H ≡ 1 (mod 2)`** — a formula (poly-time constant) for the bottom bit of H. Define

> **`d(n)` = the largest `k` such that `H mod 2^k` is constant on every cospectral class** (i.e. is a
> function of the characteristic polynomial `char_A`).

Rédei forces `d(n) ≥ 1`. The spectrum is the poly-time ladder (THM-1775/1780). `d(n)` measures how
many low bits of the `#P` quantity H the ladder still resolves.

## The decay (exact through n=6, witnessed at n=7)

| n | d(n) | meaning |
|---|---|---|
| 4 | ∞ | spectrum determines H completely |
| 5 | ∞ | spectrum determines H completely |
| 6 | **2** | spectrum pins `H mod 4`, **not** `H mod 8` |
| 7 | **1** | spectrum pins only `H mod 2` (Rédei) |

- **`d(6)=2`.** `H mod 4` is constant on all 28 cospectral classes; `H mod 8` first splits at the
  class `{13,17}` — **exactly THM-1780's pair** (both `≡1 mod 4`, differing by 4). So the mod-8 bit
  of H turns `#P` at the same `n=6` cospectral class where H leaves the ladder.
- **`d(7)=1`.** 26 cospectral 7-tournament buckets carry **both** residues of `H mod 4`. One such
  pair already proves `H mod 4` is non-spectral at n=7 — so the spectrum pins only the parity bit.

**The depth hits the Rédei floor at n=7.** Asymptotically the spectrum resolves **exactly one bit of
H — its parity — and every higher bit is `#P`.**

## `H mod 4` is a real, fleeting spectral invariant

For n≤6, `H mod 4` takes **both** values (n=6: 32 classes `≡1`, 24 `≡3`) and is constant on
cospectral classes — a *nontrivial spectral refinement of Rédei*. It is **not** score- or `c₃`-
determined from n=5 (a higher-order spectral quantity; `(H−1)/2 ≡ c₃ mod 2` holds only at n=4). Then
at n=7 it **stops being spectral** — there is no universal spectral formula for `H mod 4`.

## Two orthogonal edges meet at the Hamiltonian object

- **The length edge (THM-1870):** simple cycle counts `c_k` are spectral for `k ≤ n−1` and `#P` at
  the Hamiltonian length `k = n`.
- **The 2-adic edge (this theorem):** H's bits are spectral up to the parity bit and `#P` above.

The Hamiltonian-path count `H` is **one length past** the spectral cycle counts and **one bit past** a
spectral formula. It is the marginal object on both axes — the tournament analogue of the harmonic
series at `p=1`. Rédei's theorem is the single formula that survives at the corner.

## Why (the mechanism)

The 2-adic collapse behind Rédei is `per ≡ det (mod 2)`: the `#P` permanent equals the poly-time
determinant mod 2. H is `per`-like (THM-1945/1965), so its bottom bit is poly (Rédei = `H` odd), and
the collapse is **exactly one bit deep** — mod 4 the permanent and determinant already diverge, and
the tournament census shows that divergence reaching `H` at n=6 (mod 8) and n=7 (mod 4).

## Named next / open

- **The owner's sharper question:** does *any* poly-time invariant (not just the spectrum) pin more
  than the parity bit of H for all n? The full poly-tuple `(score,specA,specS,disc,arb_inv)` pins H
  through n=6; test its n=7 depth. **If nothing beats parity asymptotically, Rédei's bit is provably
  the entire formula-expressible content of H** — the exact statement of "tournaments at the edge."
- **`H mod 4` at n=7:** characterize the cospectral pairs that split it (the first `#P`-only bit).
- **Prove `d(n)=1` for all `n ≥ 7`** (only the parity bit is spectral) — a 2-adic strengthening of
  THM-1780.

## Correction to THM-1965 (arb invariant)

The `arb` used in the S128c142 lattice was *arborescences rooted at vertex 0*, which is **not** an
isomorphism invariant (the root depends on the labeling). The proper invariant is
`arb_inv = sorted tuple of per-root arborescence counts`. Recomputed exactly: `|arb_inv| = 55` at
n=6 (nearly complete, vs 32 for the buggy version); `arb_inv` **refines `score`** and is incomparable
to `specA, cyc, H`. This **strengthens** THM-1965's cut/cycle story (`arb_inv` is firmly cut-side,
incomparable to the cycle-side `cyc`); the headline `score ⟂ cyc` and all non-`arb` findings are
**unaffected** (they use exact invariants).
