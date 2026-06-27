# HYP-3104: LRC/Tournament Maximizer Signal Atlas

**Status:** SYNTHESIS + exact scout evidence; not a proof.

**Prompt thread:** understand what maximizes values in the LRC and how this
relates to tournaments, while searching prior work for hidden/niche
perspectives and creating new signals.

## Claim

There is no single "maximizer" object at the LRC14 frontier.  The recurring
structure is a multi-currency extremal packet:

```text
maximizer_currency_vector =
  (pair_Pascal_value_shadow,
   higher_order_interaction_debt,
   exchange_local_trap_profile,
   boundary_or_middle_mass_charge,
   far_block_coherence,
   scissors_forgetting_cost,
   tournament_rigidity_or_transfer_guard)
```

The same row can maximize one coordinate and be misleading for another.  The
pair-normalized Pascal cap is a value shadow.  Inclusion-exclusion orders record
the live interaction anatomy.  One-swap exchange tournaments record local proof
geometry.  Hamiltonian-path tournament spectra record rigidity/transfer risk,
not LRC cap value.

The most important missed point is that `k=10`/`j=3` already has live order-3
anatomy even though the final pair-Pascal value is exact:

```text
P=(1,12,13)
O0=1, O1=-3/7, O2=4/91, O3=-1/91
true lonely = 55/91 = C(11,2)/91
pair_excess = 1/91, high_order_net = -1/91, final_dip = 0
```

So "no final dip" is not the same as "pure order-2 mechanism."  The pair shadow
can be value-correct while the proof anatomy is already order-3.

## Evidence

Added exact scout:

- `04-computation/lrc_tournament_maximizer_signal_atlas_codex.py`
- `05-knowledge/results/lrc_tournament_maximizer_signal_atlas_codex.out`

Key readouts:

1. Cap anatomy for the THM-576 minimizers:
   - `j=1`: pure order 1, final dip `0`.
   - `j=2`: pure pair value, final dip `0`.
   - `j=3`: final dip `0`, but `pair_excess=1/91` and `O3=-1/91`.
   - `j=4`: final dip `1/4004`, first live order `3`.
   - `j=5`: final dip `1081/76440`, first live order `3`, and the minimizer is no longer the top endpoint cluster.

2. One-swap improvement graphs have nonglobal local sinks already at `j=2` and
   `j=3`.  For `N=13`, `j=3` has `4` local minima and `3` nonglobal traps:
   `(3,10,11)`, `(2,9,11)`, and `(5,8,9)` in addition to the global
   `(1,12,13)`.  Thus exchange optimality is not the same invariant as scalar
   cap optimality.

3. Small tournament Hamiltonian-path spectra confirm that `H` maximization is a
   balance/rigidity phenomenon.  For `n=6`, `max_H=45`, achieved only by score
   sequence `(2,2,2,3,3,3)` in the labeled scan; missing odd `H` values below
   max include `7` and `21`, matching the forbidden-H warning in HYP-3099.  This
   is a transfer guard, not a direct LRC value formula.

4. Tournament Analysis on signal families, with vertices as signals/proof
   obligations rather than runners, ranks:

```text
high_order_anatomy
> local_exchange_traps
> far_block_coherence
> H_spectrum_rigidity / scissors_forgetting_cost
> deep_tail_corner_signature
> boundary_mass_charge / decorrelation_inward_flux
> pair_pascal_shadow
> raw_scalar_value
```

It has `3` directed 3-cycles, a nontrivial SCC among the transfer/boundary
signals, and `9` Hamiltonian paths.  This says the signal layer is itself not a
single linear hierarchy; the live cycle is exactly the competition between
boundary mass, inward smoothing, scissors cost, and tournament rigidity.

## Prior-Work Synthesis

HYP-3090/HYP-3092/HYP-3097 isolate the pair-normalized Pascal shadow:
`C(14,2)=91`, `1001=11*91`, `2002=22*91`, `3003=33*91`, and `4004=44*91`.
This captures the order-2 volume face, including the one-unit `1/4004` defect at
`k=9`.  HYP-3093 adds the correct caution: equidistribution, equinumerosity,
and equidecomposability are separate shadows and cannot be freely substituted.

HYP-2738 says consecutive rows maximize deep-miss corners `C_a` in the tested
range, but signed `L_y` data is still required for the cap closure.  HYP-2603
and HYP-2675 show consecutive/AP rows as cap-ridge or sector-ridge optimizers,
yet HYP-2797 shows a genuine-wide p0 maximizer of the form "consecutive base
plus tight far doublet."  The recurring structure is not "consecutive always."
It is "boundary concentration until a coherent far block becomes a better
currency."

The gK8 reflection note adds another currency: boundary-weighted miss mass
`10 q0 + q3 + 10 q6`.  Decorrelation moves mass inward and lowers extremes.
This suggests an invariant pair:

```text
boundary_mass_charge  versus  decorrelation_inward_flux
```

The first explains why endpoint clusters look extremal; the second explains why
proof functionals can become easier after smoothing.

Incoming mac-mini S66 makes this proposed signal concrete through the
miss-count PGF

```text
G_N(z)=sum_t q_t z^t.
```

For `k=8`, `consec_8` and even-AP have `#real roots = 0/6`, high extreme mass
`q0+q6=0.3476`, and maximal displayed `L_yK8=3.5823`, while random/spread
rows move into `#real=2` or `#real=4` strata with lower mean extreme mass.  The
reported correlation `corr(#real-roots, extreme-mass)=-0.366` says root
non-realness is a measurable proxy for the same boundary/extreme mass
concentration that gK8 sees.  Add `miss_count_PGF_root_stratum` to the atlas:
it is a coefficient-to-root transform of `boundary_mass_charge`, not a separate
numerology.

HYP-3099 shows that tournament methods are powerful as proof engines but
dangerous as value analogies.  The exchange tournament is bounded but
non-transitive, and the forbidden-H bridge can be a false coincidence.  This
fits the HYP-3104 atlas: tournament signals should measure transfer risk,
local traps, and rigidity of quotient carriers.

## Bold Hypotheses

1. **First-live-order principle.**  The cap proof should not be organized by
   final Pascal dip.  It should be organized by the first interaction order
   whose inclusion-exclusion anatomy is nonzero.  Under this principle, `j=3`
   already belongs to the order-3 regime despite exact pair-Pascal value.

2. **Boundary-to-middle conservation law.**  LRC extremality repeatedly trades
   endpoint/boundary charge against middle survival mass.  Consecutive/AP rows
   concentrate boundary charge; decorrelation and far-block coherence move the
   same obstruction into middle mass or curvature tails.

3. **Exchange cycles are higher-order shadows.**  Nonglobal one-swap sinks are
   the graph-theoretic face of the same higher-order cancellations seen in
   inclusion-exclusion.  A row can be locally unbeatable because every one-swap
   destroys a third-order compensation, even when the row is not globally
   extremal.

4. **Hamiltonian-path spectra are quotient alarms.**  Missing `H` values such
   as `7` and `21` are best treated as alarms that a quotient has forgotten
   orientation or balance data.  They should trigger a scissors/equivalence
   audit, not be imported as scalar LRC inequalities.

5. **The genuine-wide doublet is a curvature maximizer, not an exception.**
   HYP-2797's "consec base plus tight far doublet" is the first visible member
   of a far-block-coherence family.  The right signal is not far speed size but
   the signed curvature tail that survives after Dedekind-style cancellation.

6. **4004 is a portal denominator.**  The number `4004=44*91` is not merely the
   `k=9` dip denominator.  It marks the first place where a pair-Pascal shadow
   must carry an explicit scissors correction.  Wild guess: `4004` is the
   smallest denominator where the LRC cap chart, Pascal mass chart, and
   exchange-trap chart all need a named residual sidecar.

7. **PGF roots are the spectral face of extreme miss mass.**  The S66 signal
   suggests that maximizers live in low-real-root strata of the miss-count PGF.
   Wild guess: the gK8 Delsarte certificate is a disguised root-confinement
   theorem for `G_N(z)`, with consecutive/even-AP rows occupying the all-complex
   stratum because their sector-empty events are maximally correlated.

## New Signals To Measure

- `maximizer_currency_vector`: tuple of the active extremal currencies for a
  row or proof obligation.
- `pascal_anatomy_residual`: pair-Pascal value minus the full
  inclusion-exclusion anatomy, with the order profile retained.
- `first_live_interaction_order`: first `r>=3` with nonzero order contribution.
- `pair_shadow_false_positive`: true when final pair-Pascal value is exact but
  higher-order anatomy is nonzero.
- `exchange_trap_index`: number and depth of nonglobal one-swap local sinks.
- `deep_tail_signature`: the `C_a` corner profile from HYP-2738 together with
  signed `L_y` debt.
- `boundary_mass_charge`: endpoint/edge concentration, boundary-weighted miss
  mass, and extreme moment weight such as `10 q0 + q3 + 10 q6`.
- `miss_count_PGF_root_stratum`: number of real roots, nearest root modulus,
  and log-concavity/root-confinement data for `G_N(z)=sum q_t z^t`.
- `decorrelation_inward_flux`: how much mass moves from extremes to middle bins
  under decorrelation or smoothing.
- `far_block_coherence`: signed curvature / doublet coherence after subtracting
  one-far effects.
- `curvature_plateau_tail`: the HYP-2797/HYP-2679 style residual after the
  plateau cap is removed.
- `scissors_split_count`: minimum retained sidecars needed to separate rows
  equal under a proposed quotient.
- `coarse_to_fine_transfer_risk`: warning score when a coarse invariant is
  value-correct but anatomy-wrong.
- `H_spectrum_transfer_risk`: tournament Hamiltonian-path gap data used only as
  a quotient-transfer alarm.

## Tournament Analysis

The useful tournament vertices in this session are not runners.  We explicitly
considered runners, speed subsets, gaps, complements, fixed circle sections,
section boundaries, miss-distribution moments, inclusion-exclusion orders,
local exchange sinks, far-block curvature packets, observer charts, tournament
Hamiltonian spectra, and proof obligations.

Chosen quotient: vertices are signal families.  The pairwise observable is
whether one signal retains more proof-relevant coordinates than another under
the axes:

```text
predicate preservation, value precision, high-order detection,
local-trap detection, quotient guard, tournament bridge,
computability, hypothesis yield
```

The gauge is majority vote across those axes, with a fixed tie Hamiltonian path.
This preserves the LRC predicate only at the level of "what information would a
proof need next."  It destroys individual runner identities, exact interval
geometry, and actual cap inequalities unless those are included as separate
signals.

Fingerprint from the scout:

```text
score_hist = {0:1, 1:1, 3:2, 4:1, 5:2, 7:1, 8:1, 9:1}
directed_3cycles = 3
Hamiltonian_path_count = 9
top score = high_order_anatomy
```

The challenged assumption is that a tournament method should preserve the
original vertex set.  Here, the useful tournament is on invariants and proof
obligations.  It records which perspectives deserve to dominate a proof search,
not which speed wins a race.

## Next Tests

1. Extend the exact cap anatomy from named minimizers to full bounded banks and
   record `first_live_interaction_order` and `pair_shadow_false_positive` for
   every local minimum.
2. Attach gK8 boundary/middle miss bins to the same rows and test whether
   `boundary_mass_charge - decorrelation_inward_flux` predicts exchange traps.
3. Join the S66 miss-count PGF root stratum to this atlas: measure whether
   `#real roots`, nearest root modulus, and log-concavity defect predict
   `exchange_trap_index`, `first_live_interaction_order`, or gK8 extremality.
4. For HYP-2797/HYP-2679 wide rows, compute `far_block_coherence` against
   `curvature_plateau_tail` and ask whether the doublet maximizer is the first
   coherent far-block family.
5. Build a transfer-risk table: when pair-Pascal, `H` spectrum, exchange traps,
   and scissors sidecars disagree, which proof obligation fails first?
6. Try a small Lean-facing interface where a cap certificate must declare its
   value shadow, first live order, and named residual sidecar before it can be
   used in the LRC14 proof route.

## Dependencies

Builds on HYP-3100, HYP-3101, HYP-3102, the incoming HYP-3103 PGF-zero and
perspective-groupoid notes, HYP-3099, HYP-3098, HYP-3097, HYP-3093,
HYP-3092, HYP-3090, HYP-2797, HYP-2738, HYP-2675, HYP-2603, HYP-+2852,
THM-576, mac-mini S66 miss-count PGF scouts, and OPEN-Q-108.
