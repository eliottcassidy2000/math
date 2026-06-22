---
id: HYP-2876
status: SUPPORTED guardrail / scaled-atlas proof target
source: codex-2026-06-22-S98
tags: [LRC14, residue-basis, character-sums, bounded-denominator, covering-sets, OPEN-Q-108]
related:
  - THM-523
  - THM-566
  - HYP-2865
  - HYP-2866
  - HYP-2870
  - HYP-2871
  - HYP-2875
  - OPEN-Q-108
results:
  - 05-knowledge/results/lrc14_residue_basis_character_sum_codex_s98.out
---

# HYP-2876: finite residue bases are atlases, not closures

The prompt basis `{83,89,21}` is useful as a sampled residue atlas, but it is
not a global finite certificate basis for covering LRC(14) rows.

The exact obstruction is stronger than the interval-cap obstruction in
THM-566.  For any finite denominator list `B` with all entries at least `2`,
set

```text
S_B = {1,2,...,11,13,84*lcm(B)}.
```

Then `S_B` is primitive and covering.  It is primitive because it contains
`1`.  It is covering because `{1,...,11,13}` covers `2..11,13`, while the tail
`84*lcm(B)` covers `12` and `14`.  For every `D in B` and every numerator `a`,
the tail speed is divisible by `D`, hence

```text
84*lcm(B) * a / D is an integer.
```

So that runner is exactly at the observer and `a/D` cannot be lonely.  Thus
`N(S_B,D)=0` for every `D in B`.  In particular, the prompt basis is killed by

```text
{1,2,...,11,13,13030668}
```

because `13030668 = 84*lcm(83,89,21)`, and the exact counts are

```text
N(S,83)=N(S,89)=N(S,21)=0.
```

So a fixed finite basis cannot close LRC(14).  The surviving target is an
adaptive or scaled residue atlas: choose a bounded collection of moduli after
factoring out the denominators killed by the speed set, or prove that the
needed first unblocked modulus has controlled analytic cost.

## Post-pull collision with HYP-+2876

During this session, incoming `HYP-+2876` asserted the stronger finite
rational-witness claim that every `13`-set has a witness with `D<=41` and that
the basis `{83,89,21}` certifies the covering family.  The construction above
refutes that global reading.  Applying it to `B={2,3,...,41}` kills every
`D<=41`; applying it to `B={83,89,21}` kills the prompt basis itself.

The shared surviving content is the `N(S,D)` character-count formulation and
the apex obstruction.  The corrected status is: finite residue bases are
sample/scaled atlases, not universal closures.

## Sample audit of the prompt basis

The S98 scout tested a deterministic broad sample of `602` primitive covering
rows with speeds up to `10000`.  The prompt basis certified `591/602` rows:

```text
first-positive denominator histogram: {83: 523, 89: 67, 21: 1}
sample failures: 11
first failure least D<=160: (19, 4)
```

This does not contradict the prompt's narrower `602/602` observation.  It says
the observation is distribution-sensitive: `{83,89,21}` is an excellent atlas
for that sampled region, but not a theorem-level closure.

## Strengthened apex floor

The clean apex-7 fragment also strengthens.  If a row is covering, then for
every reduced denominator `D in {2,...,14}` some speed is divisible by `D`.
That speed is at the observer for every numerator, so no denominator
`D<=14` can certify the row.  The obstruction is not only `D=14`; it is the
whole covering floor `2..14`.

This explains why the non-covering complement is easy and the covering core is
hard.  THM-523 gives `tau=1/q` whenever some `q<=14` divides no speed.  Once
every such `q` is present, all those small rational witnesses are dead at the
same time.

## Character-sum reading

For a fixed denominator `D`, define the exact unit count

```text
N(S,D) = #{a mod D : gcd(a,D)=1 and ||s*a/D|| >= 1/14 for all s in S}.
```

Expanding the safe-window indicators on the finite group `(Z/DZ)^*` gives

```text
N(S,D) = main term + resonance error,
```

where the resonances are the additive relations

```text
sum_s k_s s == 0 mod D.
```

The scout confirms the structural picture.  For the champion tower
`{1,...,11,13,84}`, `D=21` is zero because a speed is divisible by `21`, while
`D=41,83,89` have small positive counts.  For the prompt-basis killer, `D=83`
and `D=89` also become zero by divisibility, while `D=41` remains positive.

So a single denominator can dip to zero from resonance or divisibility even
when the main-term heuristic is positive.  A finite basis helps locally because
different moduli avoid different dips; a fixed finite basis still fails
globally because one divisor-loaded tail can align with all of them.

## Tournament and minor-order guardrail

The useful tournament vertices here are not runners.  The S98 Tournament
Analysis used proof carriers and residue atlases:

```text
THM566_divisor_loaded_no_go
scaled_residue_basis
character_sum_resonance_count
finite_sample_basis
single_denominator_certificate
speed_subset_minor_order
```

The selected Hamiltonian path is transitive, with `THM566_divisor_loaded_no_go`
leading and `speed_subset_minor_order` last.  This records the challenged
assumption: loneliness is not minor-closed under runner deletion, so any
Kuratowski-style or odd-hole analogy must live on residue/even-graph addresses
or proof obligations, not on speed-subset deletion.

The `{7,21}` / `E_7` odd-hole analogy remains a useful address diagnostic, but
HYP-2872's warning still applies: no even-graph quotient becomes a proof
carrier unless it preserves the LRC witness predicate or the exact `H` value
being studied.

## Next sharp target

Replace "finite certificate basis" by a two-level statement:

1. Divisor-loaded coherent rows are routed to the AP-core good-denominator
   ladder and treated as a scaled finite atlas.
2. Incoherent rows are handled by HYP-2875's bandlimited spectrum certificate:
   finite low resonances exactly, high tail by L2/Dirichlet cancellation.

That is the residue version of the current spectrum-sum discipline: coherent
low packets become finite atlases; incoherent packets become energy bounds.
