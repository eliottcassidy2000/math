---
id: HYP-2704
title: LRC14 survival middle mass is a tail-weighted quotient of the seven-packet recursion
status: OPEN; exact algebraic bridge stored
source: codex-2026-06-20-S65
tangent: T939
depends_on:
  - HYP-2701
  - HYP-2681
  - HYP-2689
  - THM-549
  - THM-551
  - THM-556
  - THM-558
related:
  - HYP-2703
  - HYP-2700
  - HYP-2702
  - HYP-2699
  - HYP-2680
  - HYP-2679
  - HYP-2684
---

# HYP-2704 - Survival Seven-Packet Quotient

## Claim

The HYP-2701 survival middle-mass inequality is genuinely reminiscent of the
`A+B-C+D-E-F+G` / `A+B+C-D-E-F+G` seven-packet recursions, but the correct
identification is not equality of sign patterns.  It is a quotient:

```text
seven packet / tiling side:  labeled singleton-pair-triple packets
survival side:              missed-count depth cuts G1,G5,G6
```

The seven-packet recursions use an inclusion-exclusion/Euler character on
three generators.  The survival currency uses a tail-weighted cut character:

```text
C = p1+p2+p3+p4-4p6
  = G1 - G5 - 4G6,

U4 = 1-C = 1 - G1 + G5 + 4G6.
```

Thus survival middle mass is the **Stokes dual** of the Bonferroni transfer
tax: crossing the live missed-count boundaries gives

```text
dC(t)=C(t-1)-C(t), t=1..6 = (-1,0,0,0,+1,+4).
```

The live transitions are exactly depths `1`, `5`, and `6`.  This is the same
three-live-layer phenomenon as THM-558, with signs reversed because
`C=1-U4`.

## Exact Scout

Script:

```text
04-computation/lrc14_survival_seven_packet_bridge_codex_s65.py
```

Output:

```text
05-knowledge/results/lrc14_survival_seven_packet_bridge_codex_s65.out
```

The script keeps exact `Fraction` arithmetic and verifies:

```text
C coefficients on p0..p6:   0,1,1,1,1,0,-4
U4 coefficients on p0..p6:  1,0,0,0,0,1,5
C = G1 - G5 - 4G6
U4 = 1 - G1 + G5 + 4G6
```

It also compares the relevant packet characters:

```text
actual far correction H(1):       +++++++
pair-tax A+B+C-D-E-F+G:           +++---+
odd half-tiling addressed signs:  ++-+--+
survival cuts on (G1,G5,G6):      + - -4
transfer dC across (1,5,6):       - + +4
```

The conclusion is a useful guardrail: the survival inequality should not be
proved by forcing it into the pair-tax sign vector.  The pair-tax vector remains
a phase/rank coordinate for HYP-2681.  The cap-floor proof needs the survival
cut character.

## Death-Chain Recurrence

For `r` decorrelated far runners, the missed-count depth has the exact
recurrence

```text
K_{r+1}(t) = (1-t/7)K_r(t) + (t/7)K_r(t-1),
K_0(t) = C(t),
```

equivalently

```text
Pr(t -> s after r hits)
  = binom(t,s) * sum_{j=0}^{t-s} (-1)^j binom(t-s,j) ((7-s-j)/7)^r.
```

The exact coefficient table begins:

```text
       t:         0         1         2         3         4         5         6
    r=0:         0         1         1         1         1         0        -4
    r=1:         0       6/7         1         1         1       5/7      -4/7
    r=2:         0     36/49     47/49         1         1     45/49     26/49
    r=3:         0   216/343   307/343   337/343         1   335/343   296/343
```

This is the main new proof signal.  A fully missed state starts as a large
debt:

```text
C(6) = -4.
```

After two independent far hits, its expected boundary coefficient is already
positive:

```text
K_2(6) = 26/49,
K_2(6)-C(6) = 222/49.
```

So the HYP-2701 two-far decorrelated boundary margin is not accidental
smoothing.  It is the exact high-tail transfer tax converting the `p6` debt
into middle survival mass.  The true difficulty is the actual resonant
deviation from this death-chain boundary.

## Relation To The Seven-Term Tiling Recurrences

The full and half tiling recurrences keep labelled packet ownership:

```text
full:      A+B+C-D-E-F+G
odd half: A+B-C+D-E-F+G
```

The survival quotient forgets the labels and keeps only the depth cuts whose
coefficients survive THM-556.  In that quotient the packet labels are replaced
by live boundary depths:

```text
one live low boundary:      G1
one live high-tail edge:    G5
one weighted apex tail:     G6 with multiplicity 4
```

This explains why the survival middle-mass inequality looks like a seven-term
recursion without literally having seven visible terms: the labelled `3+3+1`
packet geometry has been scalarized through missed-count depth.  The scalar
keeps the cap predicate exactly, but destroys sector labels, cyclic order,
relation-lattice phase, and far-runner ownership.  Near equality must therefore
be lifted back to HYP-2681/HYP-2679/HYP-2684 coordinates before bounding the
signed deviation.

## Proof Consequences

1. The two-far proof target should be stated as

   ```text
   actual C(B union {u,v}) = death-chain boundary C_boundary(B,2)
                              + signed resonant deviation.
   ```

   The boundary already flips the `p6` debt into positive currency; only the
   low relation-distance deviation must be controlled.

2. The `far_count>=3` branch should be easier in this quotient.  The table
   gives `K_3(6)=296/343` and all high-tail depths are already near `1`, so a
   third decorrelated carrier almost completely removes the survival apex debt.
   The remaining risk should be a signed multi-far phase deviation, not a raw
   middle-mass shortage.

3. HYP-2703's seven slope bands give the same warning from another direction:
   per-band dominance fails, but the signed aggregate wins.  HYP-2704 says the
   true-wide cap route should also keep the signed operator until the final
   margin comparison.

## Incoming S6 Signal

After the first HYP-2704 draft, pulling `origin/main` brought the mac-mini S6
coverage batch:

```text
07-reflections/the-cover-is-irreducibly-aggregate-seven-bands-not-one-mac-mini-S6.md
05-knowledge/results/lrc14_angleC_pairfloor_opus_0620.out
05-knowledge/results/lrc14_angleE_bonferroni_prefix_macmini_0620s8.out
05-knowledge/results/lrc14_angleG_coverage_majorization_opus_0620.out
05-knowledge/results/lrc14_fourier_z7_resonance_macmini_0620.out
```

This is direct reinforcement.  S6 shows that the seven-sector cover is
irreducibly aggregate: per-band extremality fails, Bonferroni prefix/block
certificates fail in the middle, and coverage-count stochastic dominance fails
at every intermediate threshold even though the top event `|colors|=7` remains
the verified winner in the bounded boxes.  The Fourier scout adds that the
mod-7 relation lattice size is always `7^(k-1)`; the discriminator is the
signed weight profile / which resonant packets appear, not raw relation count.

For HYP-2704 this means: the death-chain boundary is the right scalar quotient,
but the deviation estimate cannot be a raw relation-rank or per-layer bound.
It should keep phase/support labels until it has compared the signed aggregate
against the exact boundary margin.

## Tournament Analysis

The scout uses missed-count depths `t=0..6` as vertices, not runners, arcs, or
sector colors.  The pairwise observable is

```text
K2(t)=E[C after two decorrelated far hits | N=t].
```

The switch/gauge orients toward larger `K2`, tie-broken by smaller depth.  The
fingerprint is transitive:

```text
K2 order:
t3(1) > t4(1) > t2(47/49) > t5(45/49) > t1(36/49) > t6(26/49) > t0(0)

score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
hamiltonian_path_count=1
```

This quotient preserves the exact survival currency predicate and the
decorrelated two-far boundary.  It destroys labelled sector and relation-phase
data, so it is a gate coordinate rather than a complete proof state.

## Status

No LRC14 proof is claimed.  HYP-2704 sharpens the interpretation of HYP-2701:
the survival middle-mass gate is the tail-weighted depth quotient of the
seven-packet architecture.  The next sharp obligation remains the two-far
boundary-margin-minus-deviation lemma, followed by the `far_count>=3` signed
multi-far margin lemma.
