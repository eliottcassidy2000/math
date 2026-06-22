---
id: HYP-2867
title: LRC14 resonance-channel tournament floor -- use low-frequency channel control, QR/NQR diagnostics, and signed tournament structure to close the R' decorrelation constant
status: OPEN / proof route sharpened; exact QR/NQR pairing proved as THM-567
source: codex-2026-06-22-S94
depends_on:
  - THM-567
  - THM-565
  - HYP-+2866
  - HYP-2861
  - HYP-2857
  - HYP-2840
related:
  - HYP-2657
  - HYP-2839
  - HYP-2854
  - HYP-2864
  - HYP-2865
  - OPEN-Q-108
results:
  - 05-knowledge/results/lrc14_resonance_channel_tournament_codex_s94.out
---

# HYP-2867 -- Resonance-Channel Tournament Floor

## Target

Close the remaining LRC14 decorrelation floor by proving a uniform positive
lower bound

```text
R' = meas(coverSet(E)^c cap G_P) / (meas(G_P) * (1-p0(E))) >= c > 0.
```

HYP-2861 gives the exact spectral identity

```text
R' = 1 + SPEC / baseline,
SPEC = sum_{n != 0} chat(n) * conj(ghat(n)),
baseline = meas(G_P) * (1-p0(E)).
```

The high-frequency part is already routed by L2 Cauchy-Schwarz.  The low
part is finite and sits on resonance channels.  This hypothesis proposes that
the low part should be controlled as a signed tournament quotient on those
channels, not as a runner tournament.

Incoming HYP-+2866 is the complementary signal: it asserts that `R'=1` exactly
off resonance, then corrects its own first Paley-7 overclaim.  The floor is a
finite low-frequency resonance problem, not a diffuse decorrelation problem and
not purely apex-7 dominated.  This note takes that finite low-channel carrier
as the right object and asks how to bound its signed contribution.

## Exact Sublemma Now Available

THM-567 proves the shell-balance atom:

```text
q == 3 mod 4 and F(r)=F(-r)
  => sum_{r in QR(q)} F(r) = sum_{r in NQR(q)} F(r).
```

For LRC14, `q=7`, so every even absolute/energy residue mass on `F_7^*`
splits exactly between

```text
QR(7)  = {1,2,4},
NQR(7) = {3,5,6}.
```

This proves the arithmetic half of the QR/NQR residue diagnostic.  It should
not be read as proving that the floor covariance is Paley-7 dominated.  What
remains is the signed low-frequency phase problem: even mass balance does not
prevent all large low channels from contributing negative real parts at the
same time.

HYP-+2866's correction is important: HYP-2657's Paley-7/D7 kernel reality is
about a sector-correction kernel, not the dominant floor covariance.  THM-567
is still the abstract residue-pairing mechanism, but HYP-2867 uses it as a
diagnostic inside a broader low-frequency channel split.

## The Tournament Quotient

Alternate vertex sets considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues mod 7, low Fourier modes, sign-paired Fourier channels
{n,-n}, cover-state atoms in the 64-state Boolean lattice, sheet quotients
b/g, and proof obligations.
```

Chosen vertices:

```text
low resonance channels C_H = {{n,-n}: 0<n<=H and chat(n)conj(ghat(n)) is active},
optionally compressed by residue n mod 7 and by the gcd(P)-support lattice.
```

For real indicators,

```text
chat(-n)conj(ghat(-n)) = conj(chat(n)conj(ghat(n))),
```

so each channel has a real signed contribution

```text
s_n = 2 Re(chat(n)conj(ghat(n))).
```

Pairwise observable:

```text
opposition(c,d) = sign(s_c + s_d), with magnitude |s_c|+|s_d|
used only as a tie-breaker.
```

Switch/gauge:

```text
c -> d if the ordered pair (c,d) is nonnegative after QR/NQR pairing,
and d -> c if the pair is genuinely negative.
```

Tie Hamiltonian path:

```text
increasing |n|, then residue order 1,2,4,3,5,6, then gcd(P)-sheet quotient.
```

The quotient preserves the predicate needed for the proof:

```text
SPEC_low + SPEC_high > -baseline.
```

It destroys exact phase geometry and individual runner identity.  That loss is
intentional: exact phases belong either to the finite low-channel atlas or to
the high-tail L2 bound.  The tournament only records whether signed resonance
channels can align into a globally negative obstruction.

Challenged assumption:

```text
tournament vertices must be runners or arcs.
```

For this proof obligation, runners are too coarse and arcs are too geometric.
The decisive vertices are resonance channels and their residue-shell pairings.

## Proposed Proof Split

1. **Pair shell mass by THM-567.**  On every symmetric low packet, QR and NQR
   carry equal absolute/energy mass.  Thus a negative low sum must come from
   signed phase alignment, not from a missing residue shell.

2. **Separate the signed low sum into three ledgers.**  The first scout shows
   the signed low sum should be decomposed as

   ```text
   SPEC_low = residue0_trunk + nonzero_shell_mean + PaleyCut_component.
   ```

   The Paley cut uses

   ```text
   chi_7 = + on {1,2,4},  - on {3,5,6}.
   ```

   This is the same six-sign object that earlier appeared as K4-edge signs,
   six sector vertices, and the 64-state Boolean lattice aggregate cut.  But
   the S94 scout shows it is not the main negative term in the current AP-like
   rows: the Paley cut is positive there.  The actual damage splits between the
   nonzero-shell mean and the residue-0 trunk (`n=7,14,21,...`).  The new
   target is:

   ```text
   residue0_trunk + nonzero_shell_mean >= -eta * baseline,
   ```

   with the Paley cut used as a balance witness and diagnostic for signed
   phase coherence.

3. **If the signed channel tournament is incoherent, use L2.**  Many sign
   changes or directed cycles in the channel tournament mean the signed vector
   has small coherent projection.  Then Cauchy-Schwarz/Parseval gives the
   existing high-tail type bound with no new arithmetic.

4. **If the signed channel tournament is coherent, route to Freiman/GAP
   structure.**  A nearly transitive all-negative channel tournament means many
   low Fourier channels agree in sign.  That is high additive energy in the
   underlying endpoint and wall-crossing sets.  Use Ruzsa/Freiman modeling to
   send this branch to the already isolated bounded AP/dilate/residue atlases:
   THM-531 scale invariance, HYP-2864 sheet-gcd quotients, and the HYP-2840
   Vitali/rate-V low-resonance patch.

5. **Cash the floor through THM-565.**  Once `R' >= c`, THM-565 converts the
   slow-time floor into a finite apex-ruler witness with

   ```text
   #good >= V * meas(G) - arcCount.
   ```

## Why This May Be the Missing Constant

The previous bounded-denominator route failed globally because divisor-loaded
rows can kill any fixed rational atlas (THM-566/HYP-2865).  The spectrum route
does not require a fixed denominator.  It asks whether the signed resonance
channels can all point against the proof at once.

THM-567 says the raw nonzero residue mass at `q=7` is balanced whenever the
packet is even under `n -> -n`.  Therefore any global negative obstruction in
that quotient must be a signed phase-coherence obstruction, while the residue-0
trunk has to be bounded separately.  Phase coherence is precisely what
additive-combinatorial structure detects: low rank, high energy, bounded
quotient, or AP/dilate shape.  Those are the same branches the current LRC14
proof stack already knows how to route.

So the remaining constant should not be searched for as another real-number
epsilon.  It should be searched for as a finite signed-channel inequality: the
residue-0 trunk plus the nonzero-shell mean must be bounded below, with the
Paley cut certifying whether any QR/NQR imbalance is genuinely hostile.  The
high-coherence failure branch should be classified by Freiman structure.

## First Scout: S94 Low-Channel Atlas

Script:

```text
04-computation/lrc14_resonance_channel_tournament_codex_s94.py
```

Output:

```text
05-knowledge/results/lrc14_resonance_channel_tournament_codex_s94.out
```

At cutoff `H=21`, the AP-like negative rows have positive Paley cut:

```text
k=8 consec:      low/base=-0.2948, Paley/base=+0.0239, r0/base=-0.1223
k=9 consec:      low/base=-0.3772, Paley/base=+0.0561, r0/base=-0.0875
R' floor row:    low/base=-0.4382, Paley/base=+0.0693, r0/base=-0.1455
k=10 consec:     low/base=-0.4092, Paley/base=+0.0217, r0/base=-0.0750
```

The independence-favourable coprime row reverses the trunk:

```text
P={5,7,11}:      low/base=+0.2369, Paley/base=-0.0360, r0/base=+0.1385.
```

So the Paley cut is not the whole obstruction.  It is a diagnostic.  The
binding negative mass is carried by the residue-0 trunk plus the common
nonzero-shell mean.  The scout also records tournament fingerprints: the
21-channel tournaments are strongly connected in most rows, while the exact
`R'` floor row already decomposes as one SCC of size `14` plus seven singleton
components.  This supports the split: incoherent strongly connected sign words
should be L2-controlled, while coherent decompositions should route to bounded
Freiman/GAP atlases.

## Next Sharp Task

Extend the finite `H=21` low-channel atlas with vertices `{n,-n}` and residue
labels in `F_7^*`.  For each admissible bounded core, record:

```text
score histogram of the opposition tournament,
directed 3-cycles,
SCCs,
Hamiltonian-path count under the tie order,
PaleyCut_low / baseline,
nonzero_shell_mean / baseline,
residue0_trunk / baseline,
and the Freiman/additive-energy profile of the endpoint set.
```

The target falsifiable statement is:

```text
Either residue0_trunk + nonzero_shell_mean >= -eta * baseline with eta safely
below the L2 margin, or the endpoint set has bounded-rank GAP structure and
belongs to an existing finite AP/dilate/sheet-gcd atlas.
```

This is narrower than "prove decorrelation" and broader than another local
finite check.  It is the signed low-resonance constant that HYP-+2866 and
HYP-2861 leave open.
