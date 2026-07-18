# BSG/PFR attack the wrong half: `M<1/13` supplies no additive energy — the crux is the Diophantine→energy bridge, not energy→structure

*boxeph-2026-07-18-S104. Owner: attack the inverse theorem head-on with BSG/PFR. Outcome: BSG/PFR
(and all Freiman-type tools) require **additive energy / small doubling as input** and produce structure as
output. But `M<1/13` supplies **no such input** — verified: band-avoiding 12-sets can have Sidon-like low
energy. So the crux **factors**: Half 1 (open, the real crux) `M<1/13 ⟹` high core energy is a
**Diophantine→additive bridge**, not additive combinatorics; Half 2 `energy ⟹ AP` is where BSG/PFR live.
BSG/PFR attack Half 2 and cannot touch Half 1. This corrects the project's recurring "BSG/PFR is the tool"
hope (S89/S92). Crux not proved. Verified S104 computation.*

## The factorization the crux actually has

INV: `M<1/13 ⟹` the 12-core `C` is a dilated AP. Since AP ⟺ minimal doubling `|C−C|=2·12−1=23` ⟺ maximal
additive energy `E(C)` (S89, verified below), INV factors:

> **Half 1 (Diophantine → additive):** `M<1/13 ⟹ E(C) ≥ |C|³/K` (the core has bounded additive doubling).
> **Half 2 (additive → structure):** small doubling `⟹` dilated AP (Freiman `3k−4` / BSG / PFR + the band).

**BSG/PFR are Half 2 tools.** BSG: *high energy ⟹* a large subset of small doubling. Freiman/PFR: *small
doubling ⟹* AP/coset-progression. Every one of them takes an energy/doubling hypothesis as **input**. They
are the energy→structure direction. Half 1 — producing the energy — is not what they do.

## Why Half 1 is not supplied by `M<1/13` (verified)

`M<1/13` at the maximizer `t*=a/q` means the residues avoid the band `B_val=(−val,val)` — i.e. they sit in
`[val,12val+1]`, an interval of length `11val+1 ≈ (11/13)q`. **Band-avoidance is a weak, local constraint**
(avoid a `2/13` fraction of the circle) and gives **no** energy lower bound:

| 12-set inside the band `[14,169]` (val=14) | `|C−C|` | `E(C)` | `E/n³` |
|---|---|---|---|
| AP residues `14·{1..12}` | **23** (min) | **1156** | 0.669 |
| Sidon-like in `[14,169]` | 117 | 316 | 0.183 |
| irregular in `[14,169]` | 121 | 300 | 0.174 |

All three **avoid the danger band** (the local content of `M<1/13`), yet the Sidon-like/irregular ones have
`E ≈ 300` versus the AP's `1156`. So band-avoidance is compatible with near-minimal energy. **`M<1/13`'s
local content supplies no additive energy**, hence no input for BSG/PFR.

And the *full* residue set has large doubling anyway: `R = 14·{1..12}∪{169}` has `|R−R|=47 > 3·13−4=35` —
**dimension 2, Freiman fails on `R`** (S92). The only place small doubling appears is the core `C` *once it
is (near-)AP* — which is the conclusion, not a hypothesis.

## Where the energy must come from — and why it's the open crux

The energy cannot come from the local band-avoidance (just shown). It must come from the **global**
maximality — that `t*` beats every rational, especially the **medium moduli `13<q′<q`** (S102). That is
exactly the untouched regime, and the map from "`t*` is globally optimal" to "the core has additive energy"
is the Diophantine→additive bridge — **Half 1, open**. It is:

- **not** additive combinatorics (BSG/PFR assume the additive input Half 1 must produce);
- **not** the sieve (S102: sieve-complete families exist with low-energy cores);
- **not** maximality-local or CF-descent (S101/S103: reach only far-element divisibility).

So the systematic picture over S101–S104: the crux is a single **Diophantine→additive-energy** implication,
and *both* the elementary toolkit (maximality/sieve/CF) *and* the additive toolkit (BSG/PFR) sit on the
wrong side of it — the former never reaches the additive core, the latter presupposes it.

## What a real attack on Half 1 would need

Half 1 asks: *why does global Diophantine optimality force additive energy on the residues?* The honest
candidates, none elementary:

1. **A Fourier/`L⁴` energy bound from global optimality.** `E(C)=Σ_ℓ|Ŝ(ℓ)|⁴`. Band-avoidance controls
   `Ŝ` weakly; one would need the *medium-modulus* optimality to force a *large* Fourier coefficient
   (concentration on the dual AP `≈ (q/val)ℤ`). This is the S95 "concentration" object — and I proved Weyl
   cannot force concentration (S95). A concentration-producing tool is what's missing.
2. **A transference / dual formulation** turning "no better rational at any `q′`" into an energy inequality
   — a genuinely new lemma, not in the standard additive-combinatorics kit.

Either is a research-level result. This is Tao's `n=12` optimistic conjecture (S89/S94), open.

## Net (honest)

- **Attempted BSG/PFR head-on; found they are Half-2 tools with no Half-1 input.** Verified: `M<1/13`'s
  band-avoidance is compatible with Sidon-like low energy (`E≈300` vs AP `1156`), and the residue set is
  dimension 2 — so no doubling/energy hypothesis is available to feed them.
- **The crux is the Diophantine→energy bridge (Half 1)** — producing additive energy from global rational
  optimality — which is not additive combinatorics, not sieve, not maximality, not CF descent.
- **Correction to the project:** the recurring "BSG→Freiman/PFR is the one right tool" (S89/S92) is only the
  *second* half; the open content is the *first* half, and it is Diophantine.

LRC(14) not proved. I did not fabricate a proof. This session pins the crux to a single missing implication
(global optimality ⟹ core additive energy) and shows the standard additive toolkit is on the wrong side of
it — so closing it needs a new concentration/transference lemma, i.e. real progress on Tao `n=12`.

Cross-links:
[[cf-descent-explains-the-far-element-lcm-13-14-and-stops-three-elementary-tools-converge-the-AP-core-is-open-boxeph-S103]],
[[weyl-is-the-wrong-tool-the-one-line-form-is-concentration-not-equidistribution-boxeph-S95]],
[[the-lrc14-crux-is-sharp-freiman-additive-energy-and-a-discrete-markoff-spectrum-boxeph-S89]],
[[the-abandoned-attempts-decoded-the-crux-is-additive-dimension-two-not-any-scalar-boxeph-S92]].
