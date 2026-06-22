---
id: HYP-2870
title: LRC14 complement-even low-frequency floor -- the spectrum-sum crux is a joint signed covariance, not a Paley imbalance or fixed-denominator witness
status: OPEN / proof order sharpened by exact bounded-core atlas
source: codex-2026-06-22-S95; renumbered after incoming KPS HYP-2868 abstract-object entry
depends_on:
  - HYP-2860
  - HYP-2867
  - THM-567
  - THM-565
  - HYP-2840
  - HYP-2862
related:
  - HYP-2868
  - HYP-2635
  - HYP-2676
  - HYP-2690
  - HYP-2865
  - HYP-+2866
  - HYP-+2869
  - OPEN-Q-108
results:
  - 05-knowledge/results/lrc14_spectrum_consec_channel_atlas_codex_s95.out
---

# HYP-2870 -- Complement-Even Low-Frequency Floor

## Claim

The remaining Node-3 spectrum-sum constant should be proved by splitting the
finite low-frequency covariance into a complement-even structured packet and a
complement-odd cancellative packet.

In the notation of HYP-2860,

```text
R' = 1 + SPEC / baseline,
SPEC = sum_{n != 0} chat(n) conj(ghat(n)).
```

The high-frequency part is already controlled by L2/Parseval.  The low part is
finite.  HYP-2870 proposes that the low part is not a scalar Paley or
fixed-denominator obstruction; it is a joint signed covariance whose dangerous
branch is complement-even and high-additive-energy.

The intended proof split is:

```text
low packet = complement-even AP/Freiman packet + complement-odd fluctuation,
complement-odd fluctuation -> L2/Parseval,
complement-even packet -> finite consecutive/AP/GAP atlas + THM-531 dilation.
```

## Exact Atlas Evidence

Script:

```text
04-computation/lrc14_spectrum_consec_channel_atlas_codex_s95.py
```

Output:

```text
05-knowledge/results/lrc14_spectrum_consec_channel_atlas_codex_s95.out
```

The script scanned all consecutive bounded cores

```text
E = {0,...,k-1},  k=8,...,13,
P subset {1,...,13},  |P|+k=13,
```

for `2380` exact real-space rows.  For each row it computed `R'` exactly as a
Fraction and then used the existing Fourier engine for the first `H=21`
sign-paired low channels.

Key exact minima:

```text
global exact R' minimum:
  k=9, P={1,3,4,5}, E={0,...,8}
  R' = 416640/779291 = 0.534640...
  SPEC/base = -0.465360...
  low/base = -0.438190...
  residue0/base = -0.145523...
  nonzero/base = -0.292667...
  Paley/base = +0.069271...

global low-channel minimum:
  k=10, P={1,3,4}, E={0,...,9}
  R' = 258804/480755 = 0.538328...
  low/base = -0.462829...
  residue0/base = -0.079899...
  nonzero/base = -0.382930...
  Paley/base = +0.085467...

global residue-0 trunk minimum:
  k=8, P={1,3,4,5,11}, E={0,...,7}
  R' = 269955/478676 = 0.563962...
  residue0/base = -0.222025...

global negative Paley-cut minimum:
  k=10, P={1,2,5}, E={0,...,9}
  R' = 276920/410827 = 0.674055...
  Paley/base = -0.233873...
```

The important sign is the first two rows: the worst exact floor and the worst
low-channel floor both have a **positive** Paley cut.  Therefore the floor is
not threatened by QR/NQR imbalance alone.  The deficit is a joint

```text
residue-0 trunk + nonzero-shell mean
```

effect, exactly as HYP-2867's first scout suggested.

Selected probes reinforce the same split:

```text
wide d>=2 routed row:
  R' = 1037554/1240359 = 0.836495...
  low/base = -0.153344...

scale-separated AP row:
  R' = 255750/253489 = 1.008920...
  low/base = +0.015148...
```

So moving away from the bounded AP core pushes the covariance back toward
independence.  The hard branch is not random; it is the structured,
reflection-symmetric, high-energy AP branch.

## Maturity Map From Earlier Routes

The repository history has shaved away several false scalarizations:

1. **Fixed denominator witnesses fail.**  HYP-2865/THM-566 and the incoming
   HYP-2866 refutation show that divisor-loaded covering rows force witness
   denominators up an unbounded AP-core-good ladder.  The proof cannot bypass
   analysis with a bounded rational atlas.

2. **Pure Paley-7 fails as the floor carrier.**  HYP-+2866 first over-specialized
   the floor to apex-7, then corrected it: only a small share of the floor
   covariance sits at multiples of 7.  THM-567 remains useful as a QR/NQR
   diagnostic, not as the whole mechanism.

3. **Absolute Fourier tails fail, signed L2 survives.**  HYP-2860 records that
   crude absolute `1/n` tails are too lossy, while L2 Cauchy-Schwarz is sharp
   enough.  This points to a signed covariance object, not absolute mass.

4. **Per-cell and per-band monotonicity fail.**  HYP-2635 and the Angle-A
   history show that LRC14 extremality is aggregate.  Breaking the object too
   early destroys the cross-band compensation that keeps the cap/floor true.

5. **Freiman/GAP structure keeps reappearing.**  HYP-2635 and HYP-2676 both
   converge on the same inverse principle: coherent signed error means high
   additive energy and should route to AP/GAP normal forms; incoherent error
   should cancel before taking absolute values.

The surviving object is therefore an involutive, signed, low-frequency
covariance quotient: enough structure is retained to see AP/Freiman coherence,
but exact runner identities are discarded until a finite atlas is reached.

## Complement-Even Interpretation

Incoming S35 observes that LRC binding pairs sum to `14` and compares
`v -> 14-v` with tournament reversal/complement.  HYP-2870 uses this as a
coordinate, with one warning:

```text
do not assume v -> 14-v is a literal pointwise symmetry of G_P on [0,1)
before formalizing the finite-ruler/apex grid action.
```

The actionable version is weaker and safer:

```text
decompose the low covariance after the exact finite-ruler quotient into
the part invariant under complement-pair reversal and the residual odd part.
```

The S95 atlas supports this direction rather than proving it.  Consecutive
cores have reflection defect `0` and carry the worst bounded-core covariance.
Rows with asymmetric `P` choose which low channels are active, but they do not
create a vanishing `R'` family in the bounded consecutive atlas.  This is
consistent with:

```text
complement-even core = structured AP/Freiman finite atlas,
complement-odd part = mean-zero fluctuation controlled by L2.
```

## Proof Order

1. **Define the exact complement action.**  Work on the finite-ruler/apex grid
   where the `v -> 14-v` binding-pair symmetry is honest.  Identify its action
   on low channel packets and on the residue-0 trunk.

2. **Project the low packet.**  For fixed low cutoff `H`, write

   ```text
   SPEC_low(H) = SPEC_even(H) + SPEC_odd(H)
   ```

   under the complement action.  Prove `SPEC_odd` has zero aggregate or a
   Parseval-sized bound after summing over complement pairs.

3. **Bound the residue-0 trunk together with the nonzero mean.**  The finite
   target suggested by the atlas is not `PaleyCut >= -eta`; it is

   ```text
   residue0_trunk + nonzero_shell_mean >= -eta * baseline.
   ```

4. **Route coherent even packets to inverse additive structure.**  If many
   low channels are sign-coherent, use additive energy/Ruzsa modeling to send
   the row to an AP/GAP normal form.  Then use THM-531 dilation and the
   bounded-core atlas.

5. **Route incoherent packets to L2.**  If the low-channel sign tournament has
   many sign changes, cycles, or weak coherent projection, recycle the HYP-2860
   L2 tail argument on the finite low packet.

6. **Cash the floor through THM-565.**  Once `R' >= c`, the Node-1
   three-gap/equally-spaced-ruler lemma converts the positive measure floor
   into an actual finite witness.

## Tournament Analysis

Alternate vertex sets considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
residue classes, cover arcs, Fourier modes, complement-paired low channels,
finite-ruler proof obligations, Freiman model classes, and half-tiling mirror
orbits.
```

Chosen vertices:

```text
sign-paired low Fourier channels {n,-n}, optionally quotiented by complement
pairs and by residue mod 7.
```

Pairwise observable:

```text
sign(s_i+s_j), where s_n = 2 Re(chat(n)conj(ghat(n))).
```

Switch/gauge:

```text
i -> j when the two-channel signed contribution is nonnegative in the fixed
Hamiltonian order; reverse otherwise.
```

Tie Hamiltonian path:

```text
n = 1,2,...,H, refined by residue order 1,2,4,3,5,6,0 and then by the
finite-ruler complement quotient.
```

The quotient preserves the predicate

```text
SPEC_low + SPEC_high > -baseline.
```

It destroys exact runner identity and exact pointwise geometry.  That loss is
intentional: exact geometry belongs either to the finite AP/GAP atlas or to
the high-frequency L2 estimate.

Challenged assumption:

```text
the useful tournament vertices must be runners or arcs.
```

For the spectrum floor, the useful vertices are information channels and proof
obligations.  The geometry is visible only after quotienting by the signed
covariance structure.

## Status

Incoming HYP-+2869 gives the stronger assembly claim: universal Farey floor,
uniform `R'`, small-part LRC(<=13), and THM-565 finite-ruler conversion should
close the structure once the uniform bounds and finite check are made rigorous.
HYP-2870 is compatible with that claim but more local: it supplies a
bounded-core low-channel atlas and a proof-order lens for the `R' >= c` piece.

No LRC14 proof is claimed here.  HYP-2870 refines the next sharp target: prove the
uniform floor by a complement-even low-packet theorem plus an L2 complement-odd
bound, rather than by another scalar denominator, Paley, or absolute-tail
shortcut.
