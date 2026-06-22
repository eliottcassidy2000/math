---
id: HYP-2905
status: SYNTHESIS / proof-order theorem target
source: codex-2026-06-22
tags: [lrc14, tournaments, induction, boundary-state, strong-ear, half-tiling, finite-comb, open-q-108]
depends_on:
  - HYP-2879
  - HYP-2685
  - HYP-2901
  - HYP-2902
  - HYP-2903
  - HYP-2904
  - HYP-2900
  - HYP-2899
  - THM-082
  - THM-079
  - THM-454
  - THM-523
  - THM-560
  - THM-565
  - THM-566
  - KPS-S31v
  - KPS-S31w
  - KPS-S31x
  - KPS-S31y
  - mac-mini-S47
related:
  - HYP-2906
results:
  - 04-computation/lrc_tournament_induction_switchboard_codex.py
  - 05-knowledge/results/lrc_tournament_induction_switchboard_codex.out
---

# HYP-2905: LRC14 should be finished by boundary-state induction

## Claim

The productive tournament-induction attempts in this repo do not transfer to
LRC as raw vertex deletion, runner deletion, graph minors, or scalar recurrence
relations.  They transfer as **boundary-state induction**:

```text
recursive step = scalar main term + retained boundary/resonance ledger.
```

For tournaments, the clean model is the strong-ear formula from HYP-2879.  When
a new vertex `x` is added to a tournament `T`, with `sig[v]=1` iff `x -> v`,

```text
H(T+x) =
    sum_{b: sig[b]=1} start_T[b]
  + sum_{a: sig[a]=0} end_T[a]
  + sum_{a: sig[a]=0, b: sig[b]=1} Q_T[a,b].
```

The induction step is exact only after retaining the boundary vector
`(start,end,Q)`.

For LRC14, the corresponding large-speed step is HYP-2904.  If
`A=Safe_n(B)` has measure `mu` and `c` interval components, then

```text
meas(A cap Safe_n({v})) >= (1-2/n)mu - 2c/v.
```

The induction step is exact enough only after retaining the boundary state
`(mu,c)` or its multi-large refinement `(core floor, arc count, resonant-pair
graph)`.

Thus the LRC14 proof should be organized as a boundary-state induction tree,
not as a raw size induction.

## Exact audit

The switchboard script verifies both sides.

Tournament side, labelled strong parents through `n=5`:

```text
n  strong_parents  nonconstant_ears  formula_failures  strong_failures
3               2                12                 0                0
4              24               336                 0                0
5             544             16320                 0                0
deletion_min_by_n={4: 2, 5: 2}
```

So the strong-ear formula is exact in the audit and every nonconstant strong
ear remains strong.  Strong deletion reducibility has minimum `2` in the
audited range, matching the larger HYP-2879 audit through `n=8`.

LRC side, AP-core seed at `n=14`:

```text
seed={1,...,11,13}
measure=426/35035
components=4
least comb-certified v=768

v=30030   exact_after=313/30030   comb_lower=7472/735735
v=60060   exact_after=337/32340   comb_lower=1514/147147
v=510510  exact_after=191/18326   comb_lower=26032/2501499
```

So the large-speed induction atom has a positive exact and certified measure
for the radical/lcm committed speeds.

## The flushed-out proof tree

The current LRC14 proof should be stated as:

1. **Omit-prime branch.**  If some prime `p<=13` divides no speed, then
   `t=1/p` is a direct `>1/14` witness.  This is the direct-witness branch of
   THM-523 and the radical handle.
2. **Remove-large branch.**  If a speed is scale-separated from the rest, peel
   it by HYP-2904 and descend to a smaller LRC seed.  Iterating this is
   KPS-S31w's scale-hierarchy descent.
3. **Multi-large branch, `r<=6`.**  KPS-S31v uses the same comb-teeth estimate
   plus union bound: the leading coefficient `1-r/7` stays positive.
4. **Multi-large branch, `r>=7`.**  The union bound becomes vacuous.  The live
   analytic target is a second-moment lower bound with a bounded
   resonant-pair/divisibility defect among the large speeds.
5. **Bounded covering core.**  After omit-prime, dilation, and remove-large
   reductions, the non-descending base is a primitive bounded covering core.
   It must be closed by finite structural extremality: AP/Goddyn-Wong tight
   locus, three-gap/AP-hull rigidity, labelled Fejer/additive-energy packets,
   and the HYP-2903 missing-depth parity guard.

This is the sharp limit of induction.  It reduces every unbounded or
non-covering row to smaller proven cases or direct witnesses, but it does not
remove the bounded covering core.  That core is not a failure of induction; it
is the finite boundary atom, analogous to the exceptional single strong atoms
that survive tournament strong-ear reduction.

## S31x multiplicative atom integration

Incoming KPS-S31x adds the missing product-language around the same proof
tree.  Its scale-separated measurements support

```text
meas(Safe(S1 union S2)) ~ meas(Safe(S1)) meas(Safe(S2))
```

when the two speed clusters live at very different scales.  This is the LRC
analogue of THM-454 tournament multiplicativity over strong components:
separated clusters factor, and the irreducible objects are single-scale atoms.

This does not replace the finite-comb boundary state.  It explains what the
boundary state is measuring: the error term in turning an approximate cluster
product into an effective induction step.  HYP-2904 gives the one-cluster
effective inequality with `(mu,c)`; S31v/S31w give the multi-large/scale-tree
version; S31x identifies the scalar product that should appear once the
boundary error is controlled.

The other S31x suggestion, deletion-contraction on the resonance lattice, also
fits the switchboard thesis.  It should not be run on runners.  It should be
run on a relation graph whose vertices/edges retain the labelled resonance,
depth, and AP/GW tight-locus data.  If the bounded covering core can be
realized as an independence-polynomial or Bonferroni polynomial on that
relation graph, THM-082 becomes the right finite-core reducer.

## S47 Mode-A peel integration

Incoming mac-mini S47 sharpens the same induction map from the tournament
Mode-A side.  In that language:

```text
R1 remove-large  = Mode-A peel / vertex deletion after boundary data
R2 omit-prime    = direct resonance witness
R3 dilation      = primitive normalization
```

The point is not that raw deletion proves LRC14.  The point is that after R1,
R2, and R3 have fired wherever they are allowed, the remaining object is the
irreducible bounded tight-locus atom.  S47 names that atom as the
`{consec, GW}` core: the consecutive AP/dilate boundary and the Goddyn-Wong
sporadic boundary.

This is the LRC analogue of THM-079's `H=21` proof template.  Multiplicativity
reduces to one strong atom; the final proof is then an atom bound.  For LRC14
the atom bound is the "Moon step" analogue:

```text
bounded tight locus = {consec, GW}
and every bounded covering non-tight row has positive cap slack.
```

So HYP-2905's bounded-core obligation should now be read in two layers:
first prove the R1/R2/R3 reduction exhausts all non-atomic rows effectively;
then prove the single-atom extremality theorem, with three-gap/consec
maximality and the GW atlas as the finite boundary cases.

## S117b local interval-peeler sublemma

HYP-2906 sharpens the one-speed remove-large branch inside this switchboard.
It does not need the full seed safe-set component budget when the goal is only
existence.  If a seed with max speed `m` has one witness margin `alpha>1/n`,
then it is safe on a connected interval of length `2(alpha-1/n)/m`; one added
speed has danger teeth of length `2/(nv)`.  Therefore

```text
v > m/(n(alpha-1/n))
```

already produces a witness.  Using `alpha=1/(n-1)` from `LRC(n-1)` gives
`v>(n-1)m`, so in LRC14 a top speed larger than `13` times the second-largest
speed peels immediately by LRC13.

Thus the preferred proof order is: use HYP-2906 first for one locally huge
speed; use HYP-2904/S31v for positive-density and multi-large budgets after
that local gate fails; then send the top-balanced residue to the bounded atom.

## S31y forbidden-clique bounded atom

Incoming KPS-S31y gives a sharper possible description of the bounded atom.
The numerical identity is:

```text
14 = 2 * 7 = (binary arc states) * I(K_3,2),
```

where `I(K_r,2)=1+2r`, so `I(K_3,2)=7`.  Since `H=7` / the `K_3` conflict
graph is forbidden by the tournament obstruction, a fully rigorous LRC14
counterexample would have to pass through the exact correspondence

```text
apex-7 over-cover  <=>  forbidden K_3 conflict graph.
```

This is not yet a proof of LRC14; S31y explicitly leaves that correspondence
as the crux.  It is, however, a strong candidate for the finite atom theorem
after HYP-2906/HYP-2904 remove the locally large and scale-separated branches:
the bounded core should either be `{consec,GW}` with exact-wall witnesses, or
else its over-cover shadow would realize the forbidden `H=7` object.

## Tournament Analysis

The script chooses proof carriers as tournament vertices:

```text
finite_comb_peel
strong_ear_boundary_Q
half_tiling_address
depth_parity_newton
multi_large_union_bound
bounded_covering_core
prime_omission_witness
dilation_normalization
raw_runner_deletion
raw_tournament_minor
```

The pairwise observable is:

```text
(retains boundary state, number of proof predicates retained, size descent).
```

The resulting carrier tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
Hamiltonian path =
finite_comb_peel
> strong_ear_boundary_Q
> half_tiling_address
> depth_parity_newton
> multi_large_union_bound
> bounded_covering_core
> prime_omission_witness
> dilation_normalization
> raw_runner_deletion
> raw_tournament_minor
```

The ranking is not a theorem of importance.  It is a proof-order guardrail:
raw deletion/minor operations sit last because they destroy the boundary state
that actually controls the induction error.

## Consequence

To finish LRC14, the next sharp obligations are now explicit:

1. Make the scale-cluster product effective: prove a uniform boundary-error
   theorem interpolating the HYP-2904 comb bound and the S31x multiplicative
   atom measurement.
2. Prove the `r>=7` second-moment bound with resonant-pair/divisibility defect.
3. Prove the bounded covering core theorem, sharpened by S47 to the atom
   statement:

   ```text
   primitive bounded prime-covering 13-set
   -> {consec, GW} tight locus or positive cap slack.
   ```

4. In that bounded core theorem, use the missing-depth parity Newton ledger
   from HYP-2903 and the AP/three-gap/Legendre-Venn address sheaf from
   HYP-2901/HYP-2902, and test whether resonance-lattice
   deletion-contraction can reduce the finite atom.  Do not try to prove it by
   runner deletion.

This is the tournament-induction lesson applied to LRC: keep the boundary
carrier until the last comparison.
