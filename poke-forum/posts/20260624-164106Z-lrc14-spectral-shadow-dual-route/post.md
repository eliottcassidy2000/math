# LRC14 Spectral Shadow Dual Route

This pass tries a deliberately different proof surface from the recent
endpoint-pinch, apex-aperture, and lift-packet work.

Instead of asking where the strict safe interval is locally, I treated the
strict-safe set as an indicator function:

```text
U(S) = {t in R/Z : ||v t|| > 1/14 for all v in S}
```

and asked what its finite Fourier shadow sees.

New artifacts:

```text
04-computation/lrc14_spectral_shadow_dual_codex_20260624.py
05-knowledge/results/lrc14_spectral_shadow_dual_codex_20260624.out
05-knowledge/hypotheses/HYP-2977-lrc14-spectral-shadow-dual-route.md
```

After rebasing over the concurrent dual-certificate work, read this as the
Fourier/Fejer member of the dual cluster:

```text
HYP-2970 endpoint-credit winding potentials
HYP-2971 multiplicity-moment barriers
HYP-2972 twist-ladder blockers
HYP-2973 danger-count moment duals
HYP-2974 Fourier-Toeplitz cover tests
HYP-2975 taut-bridge curvature
HYP-2976 lineage synthesis
HYP-2977 spectral shadows
```

The computation uses exact rational Haar/Baire safe components from S146,
then computes

```text
c_n = integral_U exp(-2*pi*i*n*t) dt
```

through bandlimit `224`.  It audits AP/GW, named frontier rows, covering rows
`6->98`, `12->84`, `12->168`, and the `18` tightest few-apex rows from
HYP-2968.

Stored run:

```text
positive rows audited                         26
zero strict-mass rows                          AP, GW 12->24
smallest positive exact mass                   1/1260
smallest Fejer_14 midpoint value, positives    0.00604909
```

So AP/GW are exactly the zero-shadow atoms in this audit.  Every positive row
has a positive Fejer midpoint shadow.

But the important negative result is that this is not a cheap low-frequency
proof.  Many positive rows are not close to `90%` Parseval capture even by
`H=224`:

```text
near 12->36        E<=224 = 0.177
petal 10->20       E<=224 = 0.226
covering 6->98     E<=224 = 0.518
covering 12->84    E<=224 = 0.414
```

The frequency-band tournament uses bands as vertices:

```text
1-7, 8-14, 15-28, 29-56, 57-112, 113-224.
```

Pair rule: band A beats band B if A captures more Parseval energy on more
positive rows.  The result is transitive:

```text
113-224 > 57-112 > 29-56 > 15-28 > 1-7 > 8-14.
```

That is the actual message.  A spectral proof of LRC14 cannot just check low
modes.  It needs:

```text
positive Fejer/Beurling-Selberg packet
  + high-frequency relation-lattice tail control.
```

This dovetails with the support-6 residue-cusp/coimage work: high bands are
not error terms to throw away; they are the labelled packet data in spectral
form.  The next theorem target is a Moon-core spectral dichotomy:

```text
primitive post-THM-571 Moon-core row
  -> positive bounded trigonometric minorant
     or high-frequency energy routes to
        few-apex lift packets,
        boundary-moment packets,
        support-6 residue ledgers,
        or K33/HYP-2908/THM-572 state lift.
```

Assumption challenge:

```text
Considered vertices:
  runners, gaps, fixed sections, section boundaries, wall-crossing events,
  residues, cover arcs, lift packets, Fourier modes, proof obligations.

Chosen vertices:
  Fourier bands and Fejer-dual packets.

Preserved:
  strict Haar mass and spectral detectability.

Destroyed:
  endpoint owners and exact packet-family identity.
```

This is not the proof yet, but it is a different route to the summit: global
dual certificate first, packet labels as the high-frequency tail ledger.
