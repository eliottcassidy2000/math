# HYP-6780 -- The capped-envelope residual is scale-quotient, not a raw finite band

**Status:** exact scaling lemma proved; raw finite-band and `f>=4 => M>=0.097`
claims refuted; quotient classification open.

For every finite positive core `P` with positive-measure good set and integer
`c >= 1`, the `1/14`-good set of
`cP` is the inverse image of the good set of `P` under the degree-`c` circle
cover. Consequently its measure is unchanged, its number of interval
components is multiplied by `c`, and the THM-755 cutoff satisfies

`v*(cP) = c v*(P)`.

Thus a per-core finite interval `(max(P), v*(P)]` does not imply a uniform raw
speed bound. The primitive covering ray

`V_c = {c, 2c, ..., 12c, 13c+1}`

with `13|c` and (`gcd(c,14)>1` or `c=1 mod 14`) is unbounded, satisfies the
top-peel band condition for every such `c`, and is already closed by THM-757 with
`M(V_c)=1/13`. This refutes only the *raw finite-enumeration interpretation*,
not LRC(14) or the capped-envelope theorem.

The open structural target is a finite or recursively controlled
classification after quotienting by core scale, normalized shape, and the
extra runner's residue/offset.

## Proof of the scaling lemma

Write `T_c(t)=ct mod 1`.  The good set is

`G'_{cP} = {t : ||cp t|| >= 1/14 for every p in P} = T_c^{-1}(G'_P)`.

The degree-`c` circle cover preserves normalized Haar measure, so
`|G'_{cP}|=|G'_P|`.  Each positive-length interval component of `G'_P` has
exactly `c` disjoint lifts, so the component count used by THM-755 satisfies
`r_{cP}=c r_P`.  Substitution in THM-755 gives

`v*(cP)=r_{cP}/(pi |G'_{cP}|)=c r_P/(pi |G'_P|)=c v*(P)`.

This also identifies the scale-invariant endpoint complexity

`chi(P)=r_P/(max(P)|G'_P|)`, with `v*(P)/max(P)=chi(P)/pi`.

## Explicit obstruction to the raw cutoff

For `P={1,...,12}`, exact interval arithmetic gives

`|G'_P|=6617/194040` and `r_P=12`, hence `v*(P)=112.011...`.

For `V_c=cP union {13c+1}`, the top speed is only `13c+1`, so it lies below
`v*(cP)=112.011...c` for every `c>=1`.  This comparison is rigorous without a
floating-point value of `pi`: `pi<22/7` and

`(22/7)(13c+1)(6617/194040) < 12c`.

The exact covering condition is

`13|c` and (`gcd(c,14)>1` or `c=1 mod 14`).

Indeed `qc` carries each `q<=12`, the block carries `13` exactly when `13|c`,
and modulus `14` is carried either by some `ic` (`i<=12`) or by `13c+1`.
The first covering scale is `c=26`, not `32760`.  At
`t=(c+1)/(13c)`, every runner has clearance at least `1/13`, while the dilated
12-block gives the matching upper bound.  Thus `M(V_c)=1/13`.

For every covering `c>14`, all 13 speeds are above 14, so `f=13`; nevertheless
`M(V_c)=1/13<0.097`.  This directly refutes THM-758's sampled statement that
all `f>=4` families have `M>=0.097`.  The target `M>=1/14` remains true on this
ray, already by THM-757.

## What the correction does and does not say

- THM-755 remains valid: its residual interval is finite for each fixed core.
- THM-738 remains the strongest complete exact census, closing `f<=3`.
- The `497` terminal cutoff is a statistic from 120 sampled bodies, not a
  uniform theorem.
- THM-746's `339/513` thresholds apply to two named tail shapes, not arbitrary
  large-diameter `f>=4` families.
- The 8260 S105 families are capped, generator-restricted interval-core cases,
  not an exhaustive band.
- THM-757 closes the displayed infinite ray; it does not make the entire
  scale-quotient residual finite.

The corrected frontier is therefore a scale-normal structure theorem:

1. dispatch coherent dilation packs by THM-668/737;
2. dispatch translated/hierarchical clusters by THM-739/740 and their exact
   certificates;
3. apply THM-755/731/732 to genuinely incoherent normalized shapes;
4. prove the remaining `(shape, scale, residue/offset)` atlas finite or
   recursively decreasing.

The concurrent HYP-6785 endogenous pair-sum blocker complex is complementary:
it is finite for each family and retains exact adaptive denominators, while
HYP-6780 identifies the group action that must be quotiented before a
fleet-wide finiteness claim is possible.

## Assumption challenge and Tournament Analysis

Runner vertices are the wrong quotient here: dilation moves every runner at
once, while the LRC value is unchanged and the covering fiber can change.
Alternatives considered were runners, gaps, fixed sections, section
boundaries, wall events, residues, cover arcs, Fourier modes, matroid circuits,
and proof obligations.  The audit uses **quotient carriers** as tournament
vertices: raw set, gcd normalization, normalized core shape, shape plus residue,
and certificate-only.

The pairwise observable is weighted advantage in retaining the labels
`{LRC predicate, covering, primitivity, capped ratio, witness, scale orbit,
finite state}`.  The switch/gauge changes from predicate-first to
compression-first; the fixed tie Hamiltonian path is

`shape+residue > raw set > gcd-normalized > core shape > certificate-only`.

Both gauges are transitive: score histogram `{0:1,1:1,2:1,3:1,4:1}`, zero
directed 3-cycles, singleton SCCs, one Hamiltonian path, and two edge flips
between gauges.  Shape plus residue is first in both.  It preserves the LRC,
covering, primitivity, band-ratio, and witness labels while collapsing the
infinite scale ray.  Core shape alone preserves the band ratio but destroys the
covering and primitivity fibers; a raw runner set preserves truth but destroys
enumerative compactness.

**Artifacts:**
`04-computation/lrc14_band_scale_quotient_codex_S1.py`,
`05-knowledge/results/lrc14_band_scale_quotient_codex_S1.out`.

**Reserved by:** codex-2026-07-14-S1.
