---
id: THM-2941
title: "Critical seven-slot scalar wall, projected aligned-sector closure, and the A6 boundary"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.  On all 3,003
  literal six-body carriers the seven-slot pair-Hunter scalar first crosses
  exactly at h/7, while a literal zero-excess cover is impossible and its
  hypothetical Gram is (h/49)(7I-J).  Aligned safe surplus bounds the first
  drift absolutely for k>=2.  The lossless projected residual closes the
  five-aligned/two-drift face independently of THM-2928, and exact suffix
  filters make k=2,3,4 uniformly finite-reducible.  Exact residue-ray and
  common-status, exact-packet, forced-high-ray, antipodal-phase,
  finite-low-pair torsion, torsion-density pigeonhole, scalar-splice, and
  status-descent addenda improve the k=2/k=3 first-drift caps from 2142/380
  to 1732/275.
  THM-2928's later
  divisor-status/local-needle chain empties k=4, so only k=2,3 remain
  finite-but-uncensused; k=5,6,7 are also empty.  A lossless projection-mass
  addendum closes the common-level reflected-stalk k=1 diagonal
  `Z_(E,q)={qL-e:e in E}` for every `q>=1` and a finite heterogeneous
  two-coordinate `q<=8` box, but not arbitrary k=1.  The zero/one-aligned
  sector outside that family, the remaining finite censuses, the full
  six-body/seven-tail rung, and LRC(14) remain open.  Verification is
  internal exact computation and proof audit; there is no Lean or external
  peer-review claim.
source: root-lrc14-j7-critical-wall-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-1094-exact-two-comb-component-theorem
  - THM-1234-sharp-five-comb-compatibility-floor
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2923-complete-seven-body-six-slot-recursive-pair-hunter-closure
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - LRC(<=13)
related:
  - THM-856-hunter-tree-seven-comb-crossing
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-2184-two-scale-tail-continuation-profile
verification:
  - 04-computation/lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.py
  - 05-knowledge/results/lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.out
  - 04-computation/lrc14_j7_critical_scalar_wall_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_critical_scalar_wall_independent_thm2941.out
  - 04-computation/lrc14_j7_top7_overlap_graph_scout_thm2941.py
  - 05-knowledge/results/lrc14_j7_top7_overlap_graph_scout_thm2941.out
  - 04-computation/lrc14_j7_aligned_projected_arc_suffix_thm2941.py
  - 05-knowledge/results/lrc14_j7_aligned_projected_arc_suffix_thm2941.out
  - 04-computation/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_frontier_fibre_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_frontier_fibre_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_next_frontier_scalar_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_next_frontier_scalar_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_next_frontier_scalar_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_next_frontier_scalar_independent_thm2941.out
  - 04-computation/lrc14_j7_k3_uniform_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_uniform_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z378_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_frontier_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_frontier_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_component_residue_ray_cone_no_go_thm2941.py
  - 05-knowledge/results/lrc14_j7_component_residue_ray_cone_no_go_thm2941.out
  - 04-computation/lrc14_j7_reflected_stalk_k1_mass_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_reflected_stalk_k1_mass_closure_thm2941.out
  - 04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.out
  - 04-computation/lrc14_j7_reflected_two_coordinate_q8_mass_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_reflected_two_coordinate_q8_mass_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_projected_scalar_atlas_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_projected_scalar_atlas_thm2941.out
  - 04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out
  - 04-computation/lrc14_j7_k3_z364_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z364_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z350_scalar_slice_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z350_scalar_slice_thm2941.out
  - 04-computation/lrc14_j7_k3_z350_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z350_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z336_scalar_slice_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z336_scalar_slice_thm2941.out
  - 04-computation/lrc14_j7_k3_z336_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z336_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z330_scalar_slice_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z330_scalar_slice_thm2941.out
  - 04-computation/lrc14_j7_k3_z330_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z330_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z328_scalar_slice_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z328_scalar_slice_thm2941.out
  - 04-computation/lrc14_j7_k3_z328_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z328_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_physical_denominator_reconciliation_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_physical_denominator_reconciliation_thm2941.out
  - 04-computation/lrc14_j7_k3_physical_denominator_reconciliation_audit_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_physical_denominator_reconciliation_audit_thm2941.out
  - 04-computation/lrc14_j7_k3_z324_ray_status_frontier_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z324_ray_status_frontier_thm2941.out
  - 04-computation/lrc14_j7_k3_z324_antipodal_h_drift_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z324_antipodal_h_drift_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z312_ray_status_frontier_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z312_ray_status_frontier_thm2941.out
  - 04-computation/lrc14_j7_k3_z312_finite_low_pair_torsion_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z312_finite_low_pair_torsion_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z312_finite_low_pair_torsion_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z312_finite_low_pair_torsion_independent_thm2941.out
  - 04-computation/lrc14_j7_k3_z306_z302_z298_ray_status_descent_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z306_z302_z298_ray_status_descent_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_z294_to_z276_ray_status_torsion_descent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z294_to_z276_ray_status_torsion_descent_thm2941.out
  - 04-computation/lrc14_j7_k3_z324_located_two_cell_direct_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_z324_located_two_cell_direct_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_2004_2142_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_2004_2142_thm2941.out
  - 04-computation/lrc14_j7_k2_z2060_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z2060_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z2060_ray_status_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z2060_ray_status_independent_thm2941.out
  - 04-computation/lrc14_j7_k2_z2004_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z2004_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z2004_ray_status_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z2004_ray_status_independent_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1992_2003_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1992_2003_thm2941.out
  - 04-computation/lrc14_j7_k2_z1992_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1992_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z1992_ray_status_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1992_ray_status_independent_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1940_1991_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1940_1991_thm2941.out
  - 04-computation/lrc14_j7_k2_z1940_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1940_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z1940_ray_status_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1940_ray_status_independent_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1932_1939_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1932_1939_thm2941.out
  - 04-computation/lrc14_j7_k2_z1932_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1932_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z1932_ray_status_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1932_ray_status_independent_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1837_1931_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1837_1931_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1836_1836_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1836_1836_thm2941.out
  - 04-computation/lrc14_j7_k2_z1836_exact_body_projected_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1836_exact_body_projected_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z1836_high_wall_exact_ray_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1836_high_wall_exact_ray_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1810_1835_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1810_1835_thm2941.out
  - 04-computation/lrc14_j7_k2_z1824_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1824_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z1812_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1812_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1800_1809_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1800_1809_thm2941.out
  - 04-computation/lrc14_j7_k2_high_wall_descent_1800_1810_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_high_wall_descent_1800_1810_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_exact_descent_1800_1824_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_exact_descent_1800_1824_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z1750_literal_packet_projected_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1750_literal_packet_projected_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_band_1743_1749_literal_packet_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_band_1743_1749_literal_packet_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1790_1799_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1790_1799_thm2941.out
  - 04-computation/lrc14_j7_k2_z1790_exact_descent_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1790_exact_descent_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1780_1789_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1780_1789_thm2941.out
  - 04-computation/lrc14_j7_k2_1780_1789_exact_descent_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_1780_1789_exact_descent_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1770_1779_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1770_1779_thm2941.out
  - 04-computation/lrc14_j7_k2_z1776_exact_descent_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1776_exact_descent_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_scalar_band_1760_1769_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_scalar_band_1760_1769_thm2941.out
  - 04-computation/lrc14_j7_k2_z1768_high_wall_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1768_high_wall_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_z1758_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1758_ray_status_closure_thm2941.out
  - 04-computation/lrc14_j7_k2_cap1742_splice_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_cap1742_splice_thm2941.out
  - 04-computation/lrc14_j7_k2_cap1736_scalar_splice_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_cap1736_scalar_splice_thm2941.out
  - 04-computation/lrc14_j7_k2_z1736_hybrid_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k2_z1736_hybrid_closure_thm2941.out
---

# THM-2941 -- critical scalar wall, projected aligned closure, and A6 boundary

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.**

This theorem identifies what the successful THM-2923 recursion loses one rung
earlier.  The complete scalar state `(q_1,...,q_7,B_2)` stops at its exact
critical density on every six-body root.  The critical boundary is not a
counterexample: pointwise topology makes it empty.  What can still escape is a
sequence of positive-excess covers approaching that empty boundary.  Its
natural state is a `7 x 7` restricted-overlap Gram deformation of the
rank-six boundary, joined to the endpoint-grid owner/transition word.

Nothing here proves the six-body/seven-tail rung, the sector with fewer body
speeds, or unrestricted LRC(14).

## 1. Literal six-body carrier

On `T=R/Z`, put

```text
D_w={t in T: ||wt||<1/14}.
```

Fix

```text
E in C({1,...,14},6),
C_E=T minus union_(e in E) D_e,
h=mu(C_E)>0.                                             (1)
```

The allowed external labels are the distinct integers `w>=15`.  Write

```text
c(w)=mu(C_E intersect D_w)
```

and let

```text
q_1>=...>=q_7
```

be the exact global top seven allowed singleton coverages.  Let

```text
B_2=max_(15<=u<v) mu(C_E intersect (D_u union D_v)).       (2)
```

The computations below prove that these are maxima over the infinite label
set, not maxima inside an arbitrary search box.

## 2. General critical scalar lemma

For any mass `h`, nonincreasing ranks `q_1,...,q_7`, and pair cap `B_2`,
define

```text
G_7(a)=a+sum_(j=2)^7 min(a,q_j,B_2-a),
0<=a<=min(q_1,B_2).                                      (3)
```

Suppose

```text
q_7>=h/7,                     B_2>=2h/7.                  (4)
```

Then

```text
G_7(h/7)=h.                                               (5)
```

Indeed, all seven terms in `(3)` equal `h/7`.  If `a<h/7`, every term is at
most `a`, so

```text
G_7(a)<=7a<h.                                             (6)
```

Consequently the first hostile level is exactly

```text
lambda_7=min{a:G_7(a)>=h}=h/7.                            (7)
```

This is a no-go for the scalar certificate, not evidence for a cover.
THM-735 only bounds

```text
c(w)<h/7+gamma/w.
```

At a threshold strictly above `h/7` this gives a finite core; at `(7)` it
does not.

## 3. Infinite-tail seals

If `C_E` has `r` interval components, THM-735 gives the explicit strict tail

```text
c(w)<h/7+gamma/w,                   gamma=99r/490.         (8)
```

For a finite head through `M`, put

```text
tau=h/7+gamma/(M+1).                                      (9)
```

If its seventh rank satisfies `q_7>=tau`, every omitted singleton lies
strictly below the retained top seven.

For the pair cap, let `H_2` be the largest pair union in the finite head and
`q_1` its largest singleton.  The verifier accepts a pair seal only when

```text
H_2>=q_1+tau.                                             (10)
```

A head-tail pair is then at most `q_1+tau`.  Subadditivity also gives

```text
H_2<=2q_1,
```

so `(10)` implies `tau<=q_1`.  Therefore a pair with two omitted endpoints is
at most

```text
2tau<=q_1+tau<=H_2.                                      (11)
```

Thus `(10)` seals the actual global maximum `(2)`, including the
two-omitted-label case.

## 4. Exact all-root scalar census

The primary rational-interval verifier and an independent integer-ruler
implementation both reconstruct all

```text
binom(14,6)=3,003
```

carriers.  Every root satisfies the strict form of `(4)` and hence `(7)`.
There are no scalar closures and no positive hostile-level gaps.

The exact extrema are

```text
min(q_7-h/7) = 11/2548
  at E=(1,2,3,4,6,13),
max(q_7-h/7) = 26077/1284192
  at E=(4,6,9,12,13,14);

min(B_2-2h/7) = 1756/77175
  at E=(1,3,5,7,11,14),
max(B_2-2h/7) = 203624/2207205
  at E=(2,4,8,10,13,14);

min(max_a G_7(a)-h) = 6847/120120
  at E=(1,2,3,4,6,13),
max(max_a G_7(a)-h) = 200138/945945
  at E=(2,4,6,10,12,14).                                (12)
```

The primary verifier's smallest **finite-horizon** pair seal margin

```text
H_2-(q_1+tau)
```

is `1/21351330`.  The independent verifier separately records the smallest
**asymptotic** margin

```text
H_2-q_1-h/7=1613/3783780.
```

These are different quantities.  The maximum pair horizon is `6,635`.  The
primary implementation starts every top-seven scan at `2,000`; the
independent dynamic implementation proves that `717` is the largest actually
required top-seven horizon.

## 5. Multiplicity excess and the exact Gram deformation

Now assume, only in this section, that seven allowed labels give a literal
pointwise cover of `C_E`.  Put

```text
A_i=C_E intersect D_(w_i),
c_i=mu(A_i),
delta_i=c_i-h/7,
m(t)=sum_i 1_(A_i)(t),
Delta=sum_i delta_i.                                     (13)
```

Because the cover is pointwise, `m>=1` on `C_E`, and

```text
Delta=integral_(C_E)(m-1)>=0.                            (14)
```

If every `c_i<=h/7+epsilon`, then

```text
0<=Delta<=7epsilon,
-6epsilon<=delta_i<=epsilon,
mu{m>=2}<=Delta,
sum_(i<j)mu(A_i intersect A_j)
  =integral binom(m,2)<=(7/2)Delta,
integral_(C_E)(m-1)^2<=6Delta.                           (15)
```

Thus a near-critical cover is quantitatively an approximate measurable
partition.

There is an exact covariance form.  In `L^2(C_E)`, put

```text
f_i=1_(D_(w_i))-1/7,
p_ij=mu(C_E intersect D_(w_i) intersect D_(w_j)).
```

Direct expansion gives

```text
<f_i,f_i> = 6h/49+5delta_i/7,
<f_i,f_j> =-h/49+p_ij-(delta_i+delta_j)/7,     i!=j.     (16)
```

For positive excess this `7 x 7` Gram matrix may have rank seven:
`sum_i f_i=m-1` need not vanish.  Rank six belongs to the balanced partition
boundary below, or to an explicitly projected `1^perp` quotient; it is not a
property of every near state.

If a hypothetical cover were both exact and balanced, then

```text
delta_i=0,             p_ij=0,
Gram(f_1,...,f_7)=(h/49)(7I-J).                          (17)
```

This is the Laplacian of `K_7` scaled by `h/49`: it has rank six, kernel the
all-ones line, and realizes the regular `A_6` simplex.  Equation `(17)` is a
description of a hypothetical boundary.  It does not assert that balance
alone gives a cover, and the next section proves that no literal allowed
cover can attain this boundary.

## 6. The zero-excess boundary is empty

The danger combs are open, so `C_E` is a finite union of closed intervals.
All body danger combs contain a neighborhood of zero; hence every
positive-length component has a nonwrapping closed representative.

Suppose the pointwise cover in Section 5 had `Delta=0`.  Equation `(14)` and
the integer-valued inequality `m>=1` imply

```text
m=1 almost everywhere.                                  (18)
```

On a positive-length carrier component `I`, every pair

```text
I intersect D_(w_i) intersect D_(w_j)
```

is relatively open.  By `(18)` it has measure zero, so it is empty.  The
seven relatively open owner sets are therefore a disjoint pointwise cover of
the connected set `I`.  One owner must contain all of `I`.

The word **pointwise** is essential.  An almost-everywhere partition of
strict-open teeth may have an uncovered exit/entry seam, but an LRC
counterexample may not: such a seam belongs to neither open tooth and is a
literal safe point unless a further tooth covers it.  Thus the zero-excess
case has no endpoint handoff or gcd-capacity switching graph.  Its connected
carrier components are chamber-locked to single owners.

If `I` has length `ell` and lies in `D_w`, connectedness puts it inside one
open tooth of length `1/(7w)`.  Since `I` is closed,

```text
ell<1/(7w),                     w<1/(7ell).              (19)
```

The exact all-root component census gives

```text
min_E longest_component(C_E)=23/1092
  at E=(1,6,7,8,10,13).                                  (20)
```

Thus the owner of a longest component would obey

```text
w<1092/161<7,
```

so `w<=6`, contradicting the allowed-label condition `w>=15`.  The
per-root integral owner-bound histogram is

```text
bound   1     2      3     4    5   6
roots  39  1,601  1,049  283   28   3.                  (21)
```

Therefore every hypothetical literal six-body/seven-tail cover has

```text
Delta>0,                         max_i c_i>h/7.           (22)
```

If `w_*` realizes the largest singleton excess, then

```text
delta_*>=Delta/7.
```

Combining this with the strict tail `(8)` gives the packetwise bound

```text
w_*<7gamma/Delta.                                        (22a)
```

There is a complementary transition-width form of the same obstruction.
On a longest component `I`, no one allowed owner contains all of `I`, by
`(19)`--`(20)`.  The relatively open owner sets cover the connected interval
`I`; if all pair intersections were empty, they would disconnect it.
Therefore some pair intersection contains a nonempty relatively open
subinterval `J`.  Put `kappa=mu(J)>0`.  Since `m>=2` on `J`,

```text
Delta>=kappa,             max_i delta_i>=kappa/7,
w_*<7gamma/kappa.                                      (22b)
```

Thus every counterexample packet carries a positive owner-overlap width and a
correspondingly bounded excess owner.  This is still not a uniform bound:
`kappa` depends on the packet and may tend to zero along an escaping sequence.
That limiting degeneration is an overlap-width collapse, not a lawful
strict-open handoff.  At every actual packet an endpoint seam must remain
protected by another open tooth; an unprotected limiting seam is itself a
safe point.

The strictness in `(22)` is packetwise, not yet uniform.  A sequence of
putative covers could still have `Delta` tend to zero while its labels
escape.  Quantifying or classifying that escape is the open problem.

## 7. What pair compatibility recovers

As a scoped control, the second verifier retains the complete labelled
restricted-overlap graph on each root's seven actual top-ranked singleton
labels.  The Hunter lower bound is

```text
h-sum_i c_i+max_T sum_(ij in T)p_ij.                     (23)
```

Among the `3,003` controls:

```text
positive three-edge matching margin          47,
positive best-star margin                   256,
positive maximum-spanning-tree margin     2,200,
nonpositive maximum-tree controls            803.        (24)
```

Every one of the `3,003` literal top-seven unions nevertheless has positive
complement.  Hence compatible pair incidence is real information—the tree
recovers `1,944` controls missed by every star—but maximum-tree Hunter remains
incomplete on `803` of these actual-top-seven controls.  This is only a test
on the actual top-seven labels; it is not a reduction of arbitrary
seven-label packets to those labels, nor a proof that every possible
pair-graph functional fails on literal LRC packets.

There is also an abstract hostile control for any attempt to make first and
second moments equivalent to coverage.  On the even- and odd-parity laws on
`{0,1}^3`, the three coordinate events have the same singleton masses `1/2`
and pair masses `1/4`.  Their union masses are respectively `3/4` and `1`.
Padding by four empty labelled events gives the same seven-event Gram data
with different union predicates.  This is an information-loss example, not
a claim that either law is realized by an LRC carrier.

## 8. The multi-slope Gram/address frontier and a uniform finite sector

Let `L` be the endpoint-grid denominator of the body carrier and write a tail
speed as

```text
w=La+b,                       0<=b<L.                    (25)
```

THM-2928 proves that a literal six-body carrier cannot be covered when all
seven residues `b` vanish, or when exactly one is nonzero.  Combining it with
`(22)` yields:

> Every remaining hypothetical literal six-body/seven-tail cover has
> positive multiplicity excess and at least two nonzero endpoint-grid drift
> residues.

For exactly two drifts, THM-2928 proves that every chart with fixed body
carrier and aligned multiplier set reduces to a finite exact pair-clock
quotient.  Its pointwise argument bounds the smaller slope by `585/154` and
the larger by `13 max A`; its load-bearing extra coordinate is the
carrier-local pair endpoint current, not the global full-circle overlap.
By itself this does not prove those charts empty; the lossless projected
residual below supplies the missing terminal.

The multiplicity deformation in fact gives an **absolute** first apex
throughout the `k>=2` aligned sector.  Suppose there are `k` aligned tails,
`d=7-k` drifts, and let `m_A` and `m` count the active aligned and total tail
combs.  Pointwise on a cover,

```text
m-1 >= (m_A-1)_+.                                      (25a)
```

If `u_A` is the normalized safe mass of the `k` aligned multipliers, Boolean
grid tensorization gives

```text
Delta=integral_(C_E)(m-1)
 >=h integral_T(m_A-1)_+
 =h[k/7-mu(union_(a in A)D_a)]
 =h[u_A-d/7].                                          (25b)
```

The safe floors in THM-2928 therefore give

```text
k                 2       3         4          5          6
d                 5       4         3          2          1
u_k             66/91   55/91    558/1183   478/1365   61/273
eta_k:=u_k-d/7   1/91    3/91     51/1183    88/1365   22/273. (25c)
```

Here `u_A>=u_k`, so `(25b)` gives `Delta>=h eta_k`.  The positivity begins
exactly at `k=2`: for `k=1`, `u_A=6/7=d/7`, so this excess mechanism
vanishes.  This explains, rather than merely observes, why the zero/one-
aligned sector is the remaining infinite frontier.

Each aligned singleton has restricted mass `h/7`.  Ordering the drifts
`z_1<...<z_d` and writing `delta(z)=c(z)-h/7`,

```text
Delta=sum_(q=1)^d delta(z_q)
 <=(6r_E/49)sum_(q=1)^d 1/z_q
 <=6d r_E/(49z_1).                                    (25d)
```

Here `r_E` counts positive-length carrier components, and THM-1094 supplies
the componentwise discrepancy inequality.  The two independent carrier
reconstructions give

```text
max_E r_E/h_E=3993990/32029
```

at `E=(1,10,11,12,13,14)`, where
`h_E=32029/105105` and `r_E=38`.  Combining `(25c)`--`(25d)` gives the
body-uniform integral caps

```text
k                         2      3      4      5     6
smallest drift z_1 <=   6,947  1,852  1,062    473   189.       (25e)
```

There is a complementary lower bound on the largest drift scale.  Put

```text
M=max(E union {z_1,...,z_d}).
```

The `6+d=13-k` nonaligned/body speeds have, by settled `LRC(<=13)`, a point
with clearance at least `1/(14-k)`.  Their distance functions are
`M`-Lipschitz, so the closed arc `I` of radius

```text
R=k/[14(14-k)M]
```

around that point is safe at level `1/14`.  Let

```text
phi_L(t)=Lt mod 1,                 P=phi_L(I).
```

The body and drift teeth miss `I`, while each aligned tooth pulls back from
its normalized danger set.  A cover therefore forces

```text
P subset U_A=union_(a in A)D_a.
```

The interval-image formula gives `mu(P)=min(1,2LR)`.  Moreover `P` is
compact and `U_A` is open and proper, so the inclusion has strict measure:

```text
min(1,2LR)<mu(U_A)=1-u_A<=1-u_k.
```

Consequently every such cover must satisfy the **projected-safe-arc wall**

```text
M > alpha_k L,              alpha_k=k/[7(14-k)(1-u_k)],

k                         2        3          4           5          6
alpha_k               13/150   13/132   2366/21875   2275/18627   819/5936.
                                                               (25f)
```

The earlier whole-cell test is the weaker consequence
`M>kL/[14(14-k)]`.  Thus the aligned sector is squeezed between a bounded
first drift and a quantitatively forced second scale; it is not merely
placed in an undifferentiated finite box.

The arc is only a lower bound for a more faithful quotient.  For any fixed
drift packet `Z`, put

```text
S_(E,Z)=C_E minus union_(z in Z)D_z,
P_(E,Z)=phi_L(S_(E,Z)).
```

For a fixed aligned multiplier set `A`, the original tail-cover predicate is
equivalent to

```text
{D_z:z in Z} union {D_(La):a in A} covers C_E
 iff P_(E,Z) subset U_A.                                      (25g)
```

Indeed `D_(La)=phi_L^(-1)(D_a)`, and the points outside `S_(E,Z)` are
already covered by a drift.  Since `S_(E,Z)` and its image are compact while
`U_A` is open and proper, a completion in particular forces

```text
mu(P_(E,Z))<1-u_A<=1-u_k.                                 (25h)
```

Thus `P_(E,Z)` is a multiplier-free, lossless reduction of the aligned
completion problem, not merely a necessary statistic.  Equation `(25f)`
uses only one guaranteed arc inside `S_(E,Z)`, whereas a finite drift census
can retain the whole projected set, its components, and its endpoint
address.

There is also a finite clause representation.  Let `J_E` be the body-safe
`1/L` cells and, in the normalized coordinate `u` on cell `j`, put

```text
E_z(j)={u in [0,1]: ||z(j+u)/L||<1/14}.
```

De Morgan's law gives, modulo finitely many endpoints and therefore exactly
in Lebesgue mass, the interval identity

```text
T minus P_(E,Z)
 = intersection_(j in J_E) union_(z in Z) E_z(j).        (25i)
```

Closures may be used when evaluating this displayed Lebesgue mass, but the
pointwise predicate always retains the open `D_z`: every drift-tooth seam is
uncovered by that drift and therefore remains in `S_(E,Z)` and in its
projection.  Any isolated body-safe grid point omitted by the positive-cell
ledger projects to zero, which belongs to every aligned `D_a`; hence the
pointwise equivalence `(25g)` is unchanged.

Thus projection changes the drift problem into a finite
cell-by-drift intersection-of-unions with rational interval literals.  The
choice of a drift owner in each cell is precisely the address sidecar that
singleton excess and the Gram forget.

Unlike the safe-surplus bound, the projected statement still has content at
`k=1`.  Here `u_1=6/7`, so the same proof gives

```text
max(E union Z)>L/13,
P_(E,Z) subset D_a
```

for the one aligned multiplier `a`.  This locates the scale and forces a
single-comb shape on the projected residual, but it does not bound the first
of the six drifts.  It is therefore a genuine constraint on the remaining
infinite sector, not a finite reduction of it.

There is nevertheless one canonical six-drift family on which this
constraint is decisive at every body.  Put

```text
Z_E={L-e:e in E}.
```

On the body-safe cell `t=(j+u)/L`, let

```text
U_j=union_(e in E){u:||((L-e)u-ej)/L||<1/14}.
```

Equation `(25i)` gives `T minus P_(E,Z_E)=intersection_j U_j` in mass, so
one cell with `mu(U_j)<6/7` forces `mu(P_(E,Z_E))>1/7`.  This contradicts
`P_(E,Z_E) subset D_a` for every possible final aligned multiplier `a`.
An exact all-`3,003`-body census finds such a cell by a deterministic
two-step selector: the leftmost cell of the first positive carrier
component works for `2,997` bodies, and the corresponding cell of the
second component works for the remaining six.  The worst selected union is

```text
10028643748/12527514945 < 6/7,
```

leaving the uniform residual margin
`4964583434/87692604615`.  Direct multiplier pullback, the reflected
identity, and an independent endpoint-slab union sweep agree on all
`18,054` clauses.  Hence the canonical reflected-stalk `k=1` family is
empty.  This is a scoped family theorem, not a finite reduction or closure
of the arbitrary six-drift `k=1` sector.

The same selector closes the entire **common-level reflected diagonal**.
For every integer `q>=1`, put

```text
Z_(E,q)={qL-e:e in E}.
```

On the selected body-safe cell, write `u=(r+x)/q`.  The toothpick identity

```text
(qL-e)(j+u)/L = x-e(qj+r+x)/(qL) mod 1
```

expresses its union mass as the average of `q` level-one fine-cell masses at
ruler `qL`.  Exact interval clauses close `q=1,2,3,4`.  For `q>=5`, freezing
the base intervals and moving their diagonal endpoints gives

```text
mu(U_j^q) <= 265/336 + sum_(e in E) 2e/(qL-e)
            <= 265/336 + 14/[3(14q-1)]
            <= 6/7-19/23184.
```

Thus `mu(P_(E,Z_(E,q)))>1/7` for all `3,003` bodies and every common level
`q>=1`, so the final aligned comb cannot contain the projection.  A redundant
exact census through `q=30`, the analytic tail, and independent normal and
optimized replays agree.  The hardened LF source/output hashes are
`2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31` and
`22c078474377f0b14297497271d16426b46d2017bc8d838e88b7f7e8ba83275b`;
the semantic hash is
`c5898923940813859f7ce7401e227cdfb5c1d223b322974ceee76142b10f25fd`.
Levels varying independently with `e`, other residue packets, and arbitrary
`k=1` remain open; this diagonal extension removes no additional row from a
finite ledger.

The common-level selector is not monotone under heterogeneous levels, but its
first finite neighbourhood still closes after restoring cell location.  Put

```text
Z_(E,q_vec)={q_e L-e:e in E},   1<=q_e<=8,
```

and allow at most two coordinates of `q_vec` to differ from one.  There are
exactly `3,003*778=2,336,334` labelled packets.  The level-one selected cell
closes all but fourteen; each of those fourteen has union mass strictly above
`6/7`, so fixed-selector monotonicity is genuinely false.  An exact scan of
every body-safe cell repairs all fourteen.  Thirteen have only two bad cells;
the remaining packet has `66` bad cells among `2,260`.  Direct and reflected
arc formulae agree in `144,144` singleton checks, and merge and endpoint-slab
union routes agree in all `84,742` exceptional-cell checks.  Serial, parallel,
and optimized replays have semantic hash
`8c104172b30a3bceaec3fb7f24a48f92a785cf573ab69931f8f1345258409d05`.
This is a finite heterogeneous box, not a common-scale lifting theorem or a
closure of arbitrary `k=1`.

There is an exact all-scale functional form behind the discrepancy tail.
Write the carrier components as

```text
C_E=union_s [a_s/L,b_s/L]
```

with integral endpoints, and put

```text
P(x)=integral_0^x 1_(D_1)(u)du,       delta_E(z)=c(z)-h/7.
```

The integrand in `P` is one-periodic with mean `1/7`, so
`P(x+m)=P(x)+m/7` for every integer `m`.  Direct substitution into

```text
z c(z)=sum_s[P(zb_s/L)-P(za_s/L)]
```

therefore gives the exact recurrence

```text
(z+L)delta_E(z+L)=z delta_E(z).                         (25l)
```

Equivalently, there is a residue amplitude

```text
A_E(b)=z delta_E(z),              z=b mod L,
delta_E(La+b)=A_E(b)/(La+b).                            (25m)
```

Here `A_E(0)=0`.  Extending the even coverage function `c(z)=c(-z)` to
negative labels also gives

```text
A_E(-b)=-A_E(b).                                        (25n)
```

The amplitude itself is a boundary object.  Center the primitive by

```text
g(x)=P(x)-x/7.
```

Then `g` is one-periodic, odd, and piecewise linear, and cancellation of
the linear part gives the finite endpoint formula

```text
A_E(b)=sum_s[g(b*b_s/L)-g(b*a_s/L)].                     (25o)
```

Thus the apparently analytic drift tail is the pairing of a fixed
toothpick wave with the signed endpoint boundary of `C_E` on `Z/LZ`.
Passing from `b` to its denominator remembers the cyclic quotient on which
this pairing lives but forgets the unit direction that evaluates it.  On
one period `g` ranges exactly from `-3/49` to `3/49`; hence each carrier
component contributes at most `6/49` to `A_E(b)`.  This recovers the
component constant in `(25d)` directly from the toothpick amplitude and
identifies the sharp oscillation scale behind that bound.

The ray law is in fact componentwise, not merely true after summing the
carrier.  For a positive-length component

```text
I_s=[a_s/L,b_s/L],
c_s(z)=mu(I_s intersect D_z),
delta_s(z)=c_s(z)-|I_s|/7,
```

the same primitive calculation before summing over `s` gives

```text
(z+L)delta_s(z+L)=z delta_s(z),
A_s(b)=z delta_s(z)=g(b*b_s/L)-g(b*a_s/L),
A_E(b)=sum_s A_s(b).                                    (25o1)
```

Since `L=14 ell`, the restriction to the fourteen-grid is therefore an
exact smaller tooth ladder.  For `z=14m`,

```text
(m+ell)delta_s(14(m+ell))=m delta_s(14m),
14m delta_s(14m)=g(m*b_s/ell)-g(m*a_s/ell).              (25o2)
```

This places the literal spacing by fourteen in the upper `k=3` scalar spikes
below on a visible subladder of the full period-`L` ray system.  It does not
by itself prove that the interstitial heights are empty.

The component vector is a lawful necessary quotient of a tail cover.  If
`Z` is the set of nonaligned tails, then every aligned label `La` has
`delta_s(La)=0`, while pointwise coverage of `I_s` gives

```text
sum_(z in Z) delta_s(z)
 = integral_(I_s)(m(t)-1)dt >=0.                         (25o3)
```

If `|I_s|>=1/105`, the inequality is strict.  Equality would make the
relative-open owner sets a disjoint pointwise cover of the connected closed
interval `I_s`, so one tail tooth would contain all of `I_s`.  Every allowed
tail has `z>=15` and open-tooth length at most `1/105`, which is impossible
even at equality of lengths.

This exact positive-cone refinement has a sharp universal hostile and is not
itself a closure.  For every body speed `e in E`, the carrier is disjoint
from `D_e` up to endpoints, so componentwise

```text
A_s(e)=-e|I_s|/7,
A_s(L-e)=+e|I_s|/7.                                     (25o4)
```

Thus the reflected residue `L-e` supplies a support-one amplitude vector
strictly positive on every component.  Pure component masses consequently
prune no aligned sector, even though `(25o3)` is exact.  They forget the
actual reciprocal height, within-component tooth position, overlap, and
owner order; the denominator/unit, local-needle, and common-status sidecars
used below restore part of precisely that lost information.

Thus every drift residue is an exact hyperbolic ray, not merely an
asymptotic `O(1/z)` tail.  Opposite residue directions have opposite excess;
the self-opposite denominator-two direction has zero excess.  If

```text
d=L/gcd(L,b),             b=(L/d)u,
```

then `u` is a unit modulo `d`, and moving along the ray preserves both `d`
and `u`.  For any fixed denominator multiplicity, its all-label excess
maximum is consequently obtained by merging finitely many decreasing
nonnegative rays.  In fact this maximum is attained: every negative unit
direction has a positive opposite, every zero opposite is already an
attained zero ray, and the denominator-two ray is identically zero.  The
reversal-paired unit directions form the oriented sidecar that a bare
denominator multiset forgets.  This exact ray law will replace finite-horizon
padding in the quotient-first addenda below; it also explains the equal
positive/negative ray counts in their audits.

The exact suffix verifier composes `(25c)`--`(25f)` without treating a
search horizon as exhaustive.  For every root and proposed `z_1`, it
integrates all allowed suffix labels through `H=7,000`, retains the largest
`d-1` distinct exact excesses, and pads omitted labels by
`6r_E/[49(H+1)]`.  When `(25f)` forces a later high label, the corresponding
tail starts at the larger of `H+1` and the strict projected-wall floor.  The
result is the necessary-filter census

```text
k                              2          3        4      5     6
suffix-only max z_1         2,340        432      260    130    44
projected-wall max z_1      2,142        380      182     66  EMPTY
projected surviving rows 2,239,853    376,020   87,975  4,702     0. (25j)
```

A row here is a pair `(E,z_1)` surviving rigorous upper envelopes, not a
realized cover.  The empty `k=6` column independently reproduces the known
one-drift closure; the other columns are finite banks, not emptiness
censuses.

The unique projected `k=3` row at the maximum `z_1=380` can now be removed
exactly.  It has body

```text
E=(1,4,8,10,12,14)
```

and exact suffix dynamic programming through `H=7,000`, with the inherited
omitted-label tail, leaves only

```text
(z_1,z_2,z_3)=(380,410,492),
z_4 in {1164,1220,1358,1500,1836}.
```

All five packets pass scalar capacity, but each has a body-support fibre
larger than the sum of the four sharp divisor-fibre caps: four have margin
`3` at `q=420`, and the `z_4=1358` packet has margin `1` at `q=140`.
Thus, before testing the next integer, the projected `k=3` first-drift cap
was `379`.  This is one closed frontier row, not a census of the remaining
`376,019` necessary rows.  Ordinary and optimized replays are byte-identical; script/output
SHA-256 are

```text
64f98439f677668c82045e7f9107cbfdff467afd8f16975c7e37d8ae5c5c9f26
a1c77b24488240f1ee0295e427ee4583b7d8215caf6615f424bf325350fb56b6.
```

The next integer row is empty already at the scalar suffix level.  At
`z_1=379`, the rigorous all-tail bound leaves `2,579` of the `3,003` body
rows for exact evaluation.  Exact singleton integration through `H=7,000`,
with three distinct later drifts, the projected high-label slot, and
THM-1094 padding for every omitted label, leaves zero rows.  The closest
strict rejection is again

```text
E=(1,4,8,10,12,14),
upper-lower=-4741191283/1316479619000.
```

The live positive control reconstructs the preceding `z_1=380` body and its
positive scalar margin `437649/1736780500`.  An independent rational-interval
carrier and guarded vector primitive reproduce the complete `2,579`-row
semantic digest

```text
48ab29334a93fd0087d9645513be14f884a30bd014c2f05329c1f7d0c295d4ee.
```

Therefore, at this intermediate checkpoint, the projected `k=3` first-drift
cap is `378`.  The source
and output hashes for the integer-ruler referee and independent audit are

```text
88c563a247d59b2d9feb552935d91a2bbc5018beeed56df74c84a37a1174894b
12c1d60a6f1caf7f3a36a9bc890c388b4e44833a474e233038b5f79599715ae3
bd22ce0f86d9f5e359c2a940e0f8133849616e9aa9fc67eb823632ec9371f16d
ea7f6f2c9b189ffa4940fc25c58c74b13af905aed0fc7a6dc02266869775de77.
```

This removes no additional member of the old `376,020` necessary-row ledger:
there was no `z_1=379` row in it.  After the separate `z_1=380` packet
closure, `376,019` necessary rows remain.

The ray law also makes an all-label quotient-first closure possible far
inside that bank.  Fix

```text
E=(1,4,8,10,12,14),        z_1=250,        L=11760.
```

The first denominator is `L/gcd(L,250)=1176`.  For each multiset containing
this denominator and three arbitrary nontrivial divisors of `L`, merge the
first three eligible points on every nonnegative unit ray.  This gives the
exact attained scalar maximum over all later labels, while relaxing their
order and the projected high wall.  Of the `35,990` denominator multisets,
`1,965` survive this exact ray maximum.

The remaining quotient test strengthens THM-2928's common status table
`(37tg+)` by retaining forced pair overlap inside each status.  At an outer
divisor `q`, put `M=D/q`.  If two inner needles have denominators `d,f` and
lengths `e,z`, write

```text
g=gcd(d,f),       e=Ag+r,       z=Bg+s,       0<=r,s<g.
```

Distribution among the `g` common CRT classes and the Fréchet intersection
floor give

```text
I_(d,e;f,z)>=M/lcm(d,f) *
 (gAB+As+Br+max(0,r+s-g)).                              (25p)
```

For any number of needles, Hunter's tree inequality therefore gives the
status-wise union upper bound

```text
|union_i N_i|
 <=sum_i |N_i|-max_(spanning trees T)sum_(ij in T) I_ij.  (25q)
```

One nonnegative `2^p`-cell status table must have the exact activity
marginals and must dominate every tail of the target-load histogram.
Allowing a different table at each load would be unsound; infeasibility of
the common table is instead certified by an exact rational Farkas vector.
The crude all-divisor capacities remove `699` of the `1,965` scalar states,
and the common four-needle Hunter-status test removes all remaining
`1,266`.  Thus this entire `(E,z_1=250)` row is empty uniformly, with no
finite label horizon.  The source/output and semantic SHA-256 values are

```text
34ab29162ed33d90093e6d2bf781def36c420a1cd6596158b5d6579a3a8f3f46
bd40d2bb59946f599026709ed413a4bcb4a9fae5f768bf9a6facc77cbb18df11
c3be93a08ba7315136834ccb8dc15ce68421d3b4f102884ac6383bb0768b3948.
```

Ordinary and optimized transcripts are byte-identical.  Consequently this
separate uniform closure leaves `376,018` rows in the old necessary ledger;
it does not by itself change the current first-drift cap.

Following `MISTAKE-331`, every solver-returned Farkas vector is still checked
over the exact rationals, but the transcript and semantic digest bind only the
deterministic infeasible instances and first witness.  They do not bind the
noncanonical dual-basis representative selected by the LP search.

At the old `k=3` cap, the frontier extraction supplied nine bodies whose
recorded maximal surviving first drift was `z_1=378`.  This was not yet a
complete retest of every earlier frontier body at the lower slice.  A global
cap argument must also retest

```text
E=(1,4,8,10,12,14),
```

because this was the body at the preceding `z_1=380` frontier: closing its
`380` row and proving `379` empty does not logically decide a distinct `378`
row.  The repaired ten-body referee includes it explicitly.  Its exact
all-label ray quotient has zero scalar denominator states at `z_1=378`, so it
dies before any fibre or status transport.  Its sharp suffix maximum uses
denominators `(196,280,588,980)` and label set `{378,380,492,1500}` and misses
the required scalar floor by exactly

```text
11179177/46893073500.
```

This body was not one of the nine rows in the old necessary-row ledger.

On the nine extracted rows, replacing the generic remote tail by the exact
ray merge eliminates four.  Across the other five bodies, the exact all-label
denominator quotient has only `38` states.  The crude all-divisor screen
removes `18`, and the common four-needle Hunter-status table removes the other
`20`, again with exact Farkas certificates.  Hence all nine ledger rows are
empty, the carried-forward body is independently scalar-empty, and

```text
projected k=3 first-drift cap <=377.
```

The ray recurrence was checked at all `670,310` nonzero residues across all
ten candidate bodies (`199,915` across the five bodies with nonempty scalar
state sets); all scalar maxima are attained rather than remote suprema.
Ordinary and optimized transcripts are byte-identical.  The
source/output and semantic hashes are

```text
2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2
51afa1571a99d47b82ec0adbd25bdefb37b2c2d1f7d8ffcb50268208163f4b5c
c0108ee9009ace7ae9270bae2d32aaf9b0de62b6ad05ba859f4c86e25597e0c8.
```

Together with the distinct `z_1=250` closure, this leaves `376,009` rows
from the old `376,020` necessary ledger.

A one-pass scalar atlas over every body and every `200<=z_1<=378` exposes
the descending upper bank exactly.  It contains `6,060` scalar rows, with
successive occupied heights

```text
z_1=378,364,350,336,330,328,324,312
with body counts 9,25,53,8,1,9,45,80,
```

and no other row between these displayed heights.  The all-label ray
quotient closes all `25` `z_1=364` bodies:
`1,109=242+867` quotient states fail by crude fibre capacity and
common-status/Farkas certificates.  At `z_1=350`, an independent global
scalar slice first pins exactly `53` bodies, and their quotient has
`3,200=1,295+1,905` corresponding kills and no survivor.  At `z_1=336`,
the global slice pins eight bodies and the quotient closes all
`109=71+38` states.  The next isolated row, at `z_1=330` and
`E=(1,2,6,8,12,14)`, has one attained state
`(56,84,168,392)`.  At `(D,q,M)=(1176,168,7)`, its load histogram requires
`144` fibres of capacity at least three, while the exact status marginals
supply at most `24+24+24=72`; the independently rebuilt Farkas certificate
has normalized contradiction `-1`.  Ordinary/optimized replays and exact
audits agree.  This first moves the cap to `328`.

The lossless global slice at `z_1=328` has exactly nine bodies.  Their
all-label quotients have `85` states: crude fibre capacity removes `36`, and
independently rebuilt common-status/Farkas certificates remove the remaining
`49`.  Thus every `z_1=328` row is empty.  The next occupied atlas height is
`324`.

The denominator bridge was then rebuilt at the physical-row level before
using the four-drift screens.  It reconciles exactly

```text
raw/support rows       375,913
expected-spike rows    247,566
joint-screen rows      247,565
```

and, with multiplicity,

```text
raw              75,422,968,365
support          62,057,017,675
expected         18,783,903,428
joint            18,778,410,440.
```

This is a reconciliation of necessary denominator states, not a count of
physical covers or numerator assignments.  An independent parser and
arithmetic reconstruction gives the same tables.  At `z_1=324`, the
expected/joint screen leaves `22` of the `45` scalar bodies; the exact
all-label quotient below is therefore still required.

At `z_1=324`, the lossless atlas has `45` bodies.  Exact all-label ray
quotients produce `2,346` states: crude capacity removes `702`, exact
common-status/Farkas certificates remove `1,643`, and precisely one state
survives.  It is

```text
E=(2,8,10,11,12,14),  L=129360,
(d_1,d_2,d_3,d_4)=(3920,4620,10780,10780).
```

Exact scalar thresholding forces the later denominator-`4620` and
denominator-`10780` labels to be `364` and `492`, rules out every low
denominator-`3920` label, and leaves an arbitrary high exact-denominator-
`3920` label `z`.  This tail cannot be cut off by a larger scalar horizon:
the fixed packet already exceeds the scalar lower wall by
`5119/465404940`.

The two cells `j=5880` and `k=19600` are body-safe, and the fixed labels
`324,364,492` miss both cells.  Write every remaining high label as

```text
z=33u+mL,  gcd(u,3920)=1.
```

Then `u` is odd, while `k-j=13720` is congruent to `3920/2` modulo `3920`.
The two `z`-phases therefore differ by `u/2=1/2` modulo one, independently
of the unbounded ray parameter `m`.  Their open danger arcs, each of radius
`1/14`, are disjoint.  Hence the common danger after the two cells is empty,
the lossless projected residual has mass `1`, and it exceeds the three-
aligned open-union cap `36/91`.  Full-cell De Morgan reconstruction and an
independent direct carrier subtraction/projection agree exactly, including
on the representative packet `(324,364,492,12771)`.

Thus all `45` rows at `324` are empty.  These two slices remove `54` further
rows from the old necessary ledger, leaving `375,868`.  Since the atlas's
next occupied height is exactly `312`, this first gave the checkpoint cap

```text
z_1<=312.                                               (25q1)
```

Closing each displayed spike certifies every skipped integer above the next
spike; the descent does not extrapolate from a sparse search.

At the then-occupied endpoint `z_1=312`, the lossless atlas has exactly `80`
bodies.  The exact all-label ray quotient produces `18,249` attained
denominator states.  Crude all-divisor capacity removes `7,481`, and
independently rebuilt common-status/Farkas certificates remove `10,698`.
Exactly `70` states remain, distributed over four bodies as `25+41+2+2`.
Every attained maximizing packet has four distinct labels and begins with
`312`.

The physical-denominator bridge does not shrink this endpoint.  All `70`
states pass its expected-spike predicate, and all `70` pass its one-threshold
support-status predicate.  Among them denominator `3` and denominator `4`
never occur; denominator `2` occurs once.  Thus every small-status allowance
is exactly zero, while the minimum large-capacity support slack is `14,348`.
This bridge is only necessary and forgets unit directions and literal phases.
The terminal closure retains the high-wall coordinate that both quotients
forget.  Since `312` lies below the high floor on all four residual bodies,
the inherited projected scalar theorem requires at least one later high
label.  This hypothesis is essential: the unrestricted zero-high scalar
relaxation passes `69/70` states.  In the inherited high-wall universe, a
duplicate-permitting upper relaxation for two or three high labels falls
strictly below the wall on every body, with positive exact gaps

```text
271403663/168333225060,
22539649297/15003760917600,
295936150144/96567353117235,
41681149/16799037255.
```

Thus every possible packet has exactly one high label and two finite low
labels.  Exact ray maxima reduce the `70` states to `80`
`(state, high denominator, low pair)` cases, or `20` body-distinct low pairs.
For every case, two complete body-safe cells are missed by `312` and both low
labels, and their index difference is `d/p mod d` for the high denominator
`d`: `p=2` in `78` cases and `p=3` in two.  If

```text
z=(L/d)u+mL,  gcd(u,d)=1,
```

then the two high-label phases differ by `u/p mod 1`, independently of the
unbounded height `m`.  This is a half turn for `p=2` and a one-third or
two-thirds turn for `p=3`; in either case the strict radius-`1/14` danger arcs
are disjoint.  The projected drift-safe mass is therefore `1>36/91`.
Primary and independently reconstructed implementations agree on all `80`
cases, including the nonredundant `69/70` zero-high hostile control.  Hence
all `80` atlas rows at `312` are empty.

The next occupied heights are `306`, `302`, and `298`, with `2`, `9`, and
`12` bodies; the intervening atlas levels are empty.  Their exact all-label
ray quotients contain respectively

```text
9=0+9,   32=2+30,   96=15+81
```

states, split into crude fibre-capacity and common-status kills.  The latter
`120` exclusions are all replayed by a second exact rational Farkas checker
against the common 16-cell table; positive feasible and incompatible hostile
controls are both retained.  There are no survivors.  Closing these `23`
further atlas rows leaves `375,765` necessary rows.  The next occupied height
is `297`, with seven bodies.

At that endpoint the exact all-label ray quotient has `1,172` denominator
states.  Crude all-divisor capacity removes `271`; a locally replayed common
16-status screen removes `830`, with all `830` rational Farkas witnesses
independently reconstructed.  Exactly `71` states remain on three bodies.
The inherited `HIGH-TAIL` coordinate is again indispensable: a hostile
zero-high relaxation passes `70/71` states.  Conversely, a
duplicate-permitting two-or-more-high upper lies strictly below the scalar
wall on all three bodies.  Every possible packet therefore has exactly one
high label.  Keeping every finite low label, including negative-amplitude
ones, reduces the residual universe to `73` `(state,high denominator,low
pair)` cases and only seven body-distinct low pairs.

The terminal step uses the following reusable **torsion-density pigeonhole
lemma**.  Let `S` be a set of distinct classes in `Z/dZ`, let `r|d`, and
suppose

```text
2<=r<=7,                 |S|>d/r.
```

The subgroup `(d/r)Z/dZ` has order `r` and only `d/r` cosets, so two distinct
`c,c' in S` lie in one coset.  Their nonzero difference has effective order

```text
s=r/gcd(k,r)<=r,       c'-c=k(d/r) mod d,       0<k<r.
```

For a primitive high ray

```text
z=(L/d)u+mL,             gcd(u,d)=1,
```

the height `m` cancels and the two phases differ by `uk/r`, a nonzero
`s`-torsion fraction.  Multiplication by `u` preserves effective order, so
its circular norm is at least `1/s>=1/r>=1/7`.  Two strict-open danger arcs
of radius `1/14` cannot both contain these phases.  Equality at `s=7` is
allowed: only their excluded endpoints meet.  The lemma requires a set of
distinct residues; counting safe cells with multiplicity would not suffice.

For each of the `73` endpoint cases, `S` is precisely the distinct residue
classes modulo the high denominator represented by complete body cells
missed by `297` and both fixed low labels.  The least qualifying `r`
histogram is

```text
{2:49, 3:14, 4:6, 6:4},
```

and the resulting effective-order histogram is `{2:52, 3:17, 4:4}`.  Thus
every primitive high direction leaves one of the located complete cells
safe; the projected drift-safe mass is `1>36/91`.  All seven atlas rows at
`297` are empty.  Since `295,296` are empty and the next occupied height is
`294` with one body, the necessary ledger falls to `375,758` and the proved
checkpoint cap is

```text
z_1<=294.                                               (25q2)
```

Continuing down the exact atlas, the occupied levels

```text
294, 291, 288, 287, 286, 285, 278, 276
```

contain `45` body rows and `1,549` attained denominator states.  Exact crude
fibre capacity removes `659`; the common 16-status table removes `882`, and
all `882` exclusions are replayed as exact rational Farkas contradictions.
Positive feasible and incompatible hostile controls remain active.  Only
eight states survive, all at `z_1=286` on
`E=(1,8,10,11,12,14)`.

That row is a pinned `HIGH-TAIL` branch.  Its unrestricted zero-high scalar
relaxation passes all `8/8` states, while the duplicate-permitting
two-or-more-high upper misses the wall by exactly

```text
31003267/10140586456 > 0.
```

Thus every possible packet has exactly one high label.  Retaining every
literal low label, including negative-amplitude labels, leaves eight
one-high/two-low cases and only the low pairs `(312,350)` and `(312,364)`.
For their distinct fixed-safe residues, the least qualifying `r` and
effective-order histograms both equal `{2:6,4:2}`.  The torsion-density lemma
therefore gives phase separation at least `1/4>1/7` on every primitive high
ray, uniformly in its height.  All eight states are empty.

The intervening atlas heights `293,292,290,289,284..279,277` are empty.
The next occupied height is `275`, with ten bodies.  Removing the `45` closed
body rows lowers the necessary ledger from `375,758` to `375,713`, and the
proved cap is

```text
z_1<=275.                                               (25q3)
```

The `z_1=312` frontier
source/output/semantic hashes are

```text
24bfd9702d00454782ced222e35d3a003eaea0219c58b34dd9bffacd5e264bd4
e8fa74d4757d4a1947ce93fdc29cf8de00b75d468af0bc9ed33a8d798cfcac85
09fc1f3f6a84fb906a3afdfd5a91268795f9db600e97090650028dca945fb5fd.
```

The finite-low-pair closure source/output/semantic hashes are

```text
6b644fbb4abdbcb9b929b1789a7e73177bcce3116ad1bdc3b4ee4216adb7042c
c94568f5970948e985f920da6bc9873173e85c79e0712b1c84000a939f4612a0
95f9da4ccf85145f974df62d15c79acfa7000e0917e877bb5b32926e4ff6c3e8.
```

The independent closure source/output/semantic hashes are

```text
80e731d139ce97608b3474bdd3d3918d057c205e9ffcff9b78edc54bacacce5e
c0daca33747d968fde46c2cd767e898892dc4fb6762a4db4faa473b7f574c6c7
4667f8f0ef18add7466f849c6f902048d436692aa04c1f9037c2a3623709c91e.
```

The `306/302/298` descent source/output/semantic hashes are

```text
97fbefb8ffb17bf6742948027ec2d18dd9212bf76d35f5511e0b1c8dff64b4fc
316d069c3f560ebc504283ea9f9ab19ce0987f73d8763bb8a27eeda48fcac24a
b589ac0d64fd94468e0971de3b81a3bcffb8e446e506aa210ff85b4aa29dce5b.
```

The `z_1=297` torsion-density source/output/semantic hashes are

```text
f4464e01f0ada1515510a7d59b00582db3d677dd2a91b25407d25702d204e4e5
a0de530aefe273ff74a5494867ca31d29d00d66811173cbfbeaadbfcab99e421
1a4d141dc2dfc49fc0b15b1ec7ba7f4f637827985b3b946f5f367960aa0cbee6.
```

The `z_1=294` through `276` descent source/output/semantic hashes are

```text
f4d7292d5f3ce1b6f10f65706c15457490058a3dc59823e8c176066024307309
f709b8594af9e506f1a6028c9263c743b6340086deb21d41fd28994aa3e2aad7
cb202da83bb1b9b419576c1fb1ed7589d2a25d802f55f1845e0c2023f4f5e644.
```

The physical reconciliation source/output pairs are
`d209433bae4411e1a9597cb60d792b0c18897dffc333e60baaad5abcdaf6cf29 /
6f24829835f1a9eb2f3afcc668ccfa13d21320bb53b27fb89e30911bd20569ab`
and
`ac7d5e5dfda7451f903c1f79bff1422964225b8a2f1400b2f3e53b41b476b776 /
554f73032b2e72f324ee83c2ff751b28d58ccb2a3d367b96ee42166dd2086005`.
The `z_1=324` frontier source/output pair is
`7eaaf551d2bd4ae386e2db4452edac7d30c25f7fc67b71967e48454d688bf78e /
db3c5c68c4aa2f61584ef91dd2171901888270edbf17e860d40f16a64d3a9242`;
the antipodal and independent direct closure pairs are
`e5fbb5166975ddfaf42b30802269c0f76da605df8165f7f2bb05f881d99c5237 /
170523d62180e4779c6eb88b297df75f952b8efa28a58289c5e07dd8d2947d61`
and
`edd44a31aec80daaefb82180e7e26c388d351da914804bcaac73b06387a56d5c /
9f19d0b528606bfd4278b772f6f155eda44ae2c0fa304f5f1718bf356bba944a`.

The same quotient also closes the unique row at the much larger `k=2`
frontier.  For

```text
E=(1,4,8,10,12,14),       z_1=2142,       d_1=280,
```

all `557,845` four-suffix denominator multisets have exact attained
all-label maxima.  Exactly one survives:

```text
denominators=(196,280,840,980,2940),
labels=(2142,2172,2396,2534,3180).
```

At `D=5880`, outer modulus `q=840`, and inner modulus `M=7`, its five
status marginals are

```text
(0,120,120,0,0).
```

The K5 Hunter cap is at most `3` on the zero status.  Consequently every
one of the `300` target fibres of load at least four must activate one of
the only two nonconstant status bits, but their total incidence is only
`120+120=240`.  This is already a transparent one-threshold contradiction;
the full 32-cell simultaneous-threshold LP at loads `4,5,6` independently
returns the exact Farkas certificate.  No monotonicity of the Hunter cap
and no integer-realization claim is used.  Therefore

```text
projected k=2 first-drift cap <=2141.
```

Independent ray, carrier, K5-tree, status, and Farkas reconstructions agree.
The source/output and semantic hashes are

```text
bc4338935603b8971a99905033753458c880845213fddb5c4c19d8d53d6bc95b
b90042cbce51580280d94c84f2b798f00928d92a0dc70d2140dbe55ec3228ff7
c7f8ca7db80c0857a312d1988437d5a938914dd3ca66d054e5e04fe9a1822583.
```

Global scalar bands then make the `k=2` descent exact rather than
row-by-row.  On all `3,003` bodies, `2004<=z_1<=2142` leaves only
`z_1=2142,2060,2004`, all on `E=(1,4,8,10,12,14)`; the latter two
all-label quotients have respectively `16` and `20` states and leave none.
The next bands have unique scalar survivors at `z_1=1992`, `1940`, and
`1932` on the same body.  Their quotients contain `10=3+7`, `1=0+1`, and
`1=0+1` crude/status kills, again with no survivor.

The first four fixed states repeatedly expose the same transparent
quotient: at `(D,q,M)=(5880,840,7)`, `300` target fibres of load at least
four would require an exceptional status, but the status marginals supply
only `240` exceptional events.  At `z_1=1932` the stronger threshold-three
form has `720` heavy fibres but only `360` events.  Exact 32-cell Farkas
certificates and independent no-LP Pruefer-tree audits agree.

The exact all-body band `1837<=z_1<=1931` is scalar-empty.  The next global
slice, at `z_1=1836`, leaves exactly five rows.  For the two
large-period rows, the residue-ray identity
`delta(r+mL)=A(r)/(r+mL)` makes the all-label top-four envelope exact: the
first four positive points on each ray dominate every later point.  Replacing
the fourth unconstrained label by the best label above the forced wall
`floor(13L/150)+1` puts both scalar sums strictly below `h/91`.

The other three rows have `L=11760`.  Their all-denominator quotient contains
`887` scalar states: `180` die by crude fibre capacity, `649` by common K5
status, and `58` remain.  Scalar slack makes the literal ray choices finite,
giving exactly `84` packets.  For every packet the lossless De Morgan
projection satisfies

```text
mu(P_(E,Z))>=25/91,
```

with a prefix of at most two body cells; direct interval subtraction agrees
with the full projected-cell computation.  Since two aligned open combs can
complete only when `mu(P_(E,Z))<25/91`, all `84` packets are impossible.
Ordinary and optimized replays match stored bytes, and an independent audit
checks the ray dominance, finite-slack exhaustiveness, projection direction,
and exact arithmetic.  These five rows lower the old necessary ledger to
`2,239,842`.  At this intermediate checkpoint the proved projected `k=2`
cap was

```text
z_1<=1835.                                              (25q2)
```

Two further all-body bands complete the next descent without interpolating
between sampled heights.  Over all `3,003` bodies and every first label in
`1810..1835` and `1800..1809`, the guarded scalar atlas leaves exactly
sixteen rows, with occupied-height multiplicities

```text
1824:1, 1812:1, 1810:8, 1807:3, 1805:1, 1800:2.
```

Every other height in those bands is scalar-empty.  The atlas rows partition
exactly into ten exact-suffix rows and six forced-high rows; an explicit union
gate verifies that no row is omitted or counted in both branches.

On the ten exact-suffix rows, the all-label residue-ray denominator quotient
has

```text
558 = 137 crude-capacity kills + 328 common-K5-status kills
      + 93 scalar survivors.
```

Exact scalar slack expands the last `93` quotient states into `147` literal
packets.  Direct carrier subtraction and the lossless projected-cell engine
agree on every packet, use prefixes of at most two body cells, and give

```text
mu(P_(E,Z)) >= 25/91
```

with strict minimum margin `1026/16471`.  Thus all `147` packets fail the
necessary two-aligned completion inequality.  As redundant intermediate
checks, the unique `z_1=1824` and `z_1=1812` rows independently close by
`38=15+23` and `11=4+7` crude/status kills, respectively.  Those checks do
not replace or weaken the combined all-row descent.

For each of the six forced-high rows, the exact ray law
`delta(r+mL)=A(r)/(r+mL)` reduces the constrained suffix maximum to the three
largest unrestricted ray values together with the best value above
`floor(13L/150)+1`.  In every case that exact maximum is strictly below
`h/91`, so the high branch is scalar-empty.  Hence all sixteen atlas rows are
empty, the old necessary ledger falls to `2,239,826`, and the current proved
projected `k=2` cap is

```text
z_1<=1799.                                              (25q3)
```

The repaired/replayed source-output SHA-256 evidence pins for the two scalar
bands, the two redundant single-row checks, the high-wall closure, and the
combined exact descent are, respectively,

```text
60a916b6d4cfafc995b9ebc791057dd9afef1a33f312c32ef4da7fdc0151cec4 /
d197eb6179a3f7c7da08d4389fde988c0bd1fbc5db8cfaf8e30435ace3c7d87f
108c55c274c90fbea26131c29110d30d29f29c8133db14869e546ed2c13810b8 /
a652db146760a151572ca2ff8f093cf297cf3a6322df441e530d9da3fb24ba0a
d18157c7e4d074b7d2d2d6081641d801441bb8cd64f38fbc4f8224c597d12e60 /
6f3d7bb75d9b475ba28a005fd796d71f2f910511346c257975df1d1b604107ad
a4120f84a0bab99ccb55596f4a559383b4e5af82b9cf5dcd9d190ed67df0dc21 /
7dcb6602efb3b266f669a3acc0b4c458c6db8f31e8fac66511fbb4cf53184566
853941bc3621ef44e053a2f3382621799c30d89cc1d7ef30c63bf114554270ed /
cf051d65e11743c9357ae361328cebfa6d738a21cec7161d87ad25d4446c393a
e83dfeba64f14c53abcf2a2c67ff000dc4f9e79ac786ba4fe20ab9ad4a76d744 /
c2536ac8100d3dce937c4ca51ca50c0b8a8ef72ce456e339064ae2869ddaab8c.
```

The next two all-body bands are exact.  On `1810<=z_1<=1835`, the scalar
envelope leaves ten rows: eight at `1810`, one at `1812`, and one at `1824`;
the interstitial heights `1811`, `1813..1823`, and `1825..1835` are empty.
On `1800<=z_1<=1809`, it leaves six rows: three at `1807`, one at `1805`,
and two at `1800`; the other six heights are empty.  The sixteen rows split
without overlap into ten exact-suffix rows and six `HIGH-TAIL` rows.

For the exact-suffix class, attained residue-ray maxima, all-divisor crude
capacity, common-`K5` status, scalar-slack literal packets, and the lossless
projected residual give the global pipeline

```text
scalar states          558
crude kills            137
status kills           328
status survivors        93
literal packets        147
projected kills        147
survivors                0.
```

The projected step uses
`P_(E,Z)=phi_L(C_E minus union_z D_z)` and the sharp two-aligned cap
`25/91`.  It is uniform over all distinct later nonaligned labels and has no
finite label horizon.  The separate `z_1=1824` and `1812` packages reproduce
their `38=15+23` and `11=4+7` state closures independently.

For the six `HIGH-TAIL` rows, a later label is forced above
`floor(13L/150)+1`.  On each positive residue ray,
`delta(r+mL)=A(r)/(r+mL)`; hence the first four points dominate the omitted
tail.  Replacing the fourth unrestricted maximum by the best wall-eligible
point puts every constrained top-four sum strictly below `h/91`.  This
closes two rows at `1810`, three at `1807`, and one at `1800` by exact scalar
inequalities, again uniformly over distinct later nonaligned labels.

Thus all sixteen scalar exceptions on `1800..1835` are empty.  They lower
the necessary ledger from `2,239,842` to `2,239,826`, and the current proved
projected `k=2` cap is

```text
z_1<=1799.                                              (25q3)
```

These repeated certificates reveal a quotient family, but do not yet prove
that it closes every lower scalar state: below the isolated upper spikes the
scalar bank becomes rapidly denser.

Four further exact all-body bands descend through `1760`.  On `1790..1799`
the scalar atlas leaves eight rows, all at `1790`.  The canonical combined
engine closes six by the all-label status/projected route and two by strictly
negative forced-high ray maxima.  In particular, the one positive scalar-gap
high row is checked in the larger unrestricted all-label universe: its
`64=1+16+47` states generate `68` packets, all `68` fail the two-aligned cap,
and the minimum margin is `891/5369`.  Dropping the high constraint enlarges
the universe, so this implication has the safe direction.

On `1780..1789` the scalar atlas leaves twelve rows: five at `1780`, six at
`1784`, and one at `1788`.  Ten close by the all-label status/projected engine;
the two forced-high rows have strict gaps

```text
-36718081/292931597140,
-285200240183/3560048135813715.
```

The band `1770..1779` leaves only
`(z_1,E)=(1776,(1,4,8,10,12,14))`; its `77` scalar states split as `29`
crude-capacity and `48` exact common-status/Farkas failures.  The band
`1760..1769` likewise leaves only `(1768,(1,8,10,12,13,14))`, and its
forced-high decreasing-ray invoice has exact gap
`-906233/6582732520`.

Thus all twenty-two newly exposed rows are empty.  Ordinary and optimized
runs match stored bytes, each new script pins its complete profile and census
and fails closed, and an independent replay reconstructed the four bands.
The all-body `1750..1799` atlas certifies every height `1751..1759`
scalar-empty except `1758`; that row has two scalar states and both fail exact
common-`K5` status.  Hence the upper descent reaches `z_1<=1750`.

Two independent lower packages close the contiguous integer block
`1743<=z_1<=1750`.  The all-body slice at `1750` has exactly twelve scalar
rows.  Exact residue-ray high-wall envelopes make four rows scalar-empty.  On
the other eight rows the all-label denominator quotient has `682` scalar
states; common-`K5` status removes `582`, and the last `100` states expand by
positive scalar slack to exactly `229` literal packets.  Every packet has
lossless projected residual at least `25/91`, with minimum strict margin
`4085/54691`, and direct global subtraction agrees with the cellwise De
Morgan projection.

The independent all-body atlas on `1743..1749` checks all `21,021` candidate
rows and leaves only

```text
z_1=1746, E=(1,6,9,10,12,14);
z_1=1743, E=(1,4,8,10,12,14).
```

Thus `1744,1745,1747,1748,1749` are scalar-empty rather than interpolated.
The `1746` row has one exact denominator state and dies by crude fibre
capacity.  At `1743`, `11` scalar states split as four common-status kills
and seven survivors; positive slack gives exactly seven literal packets, all
projected-empty with minimum margin `1026/16471`.  Normal, serial, and
optimized transcripts agree byte for byte, and direct subtraction again
checks the minimum projected packet.

The hash-pinned splice verifier checks the upper cap, the complete `1750`
slice, the complete `1743..1749` band, and every integer in their union; no
height is supplied by interpolation.  The ledger update is

```text
2,239,804 - 1 - 12 - 2 = 2,239,789,
```

and this first splice gives

```text
z_1<=1742.                                              (25q4)
```

The hash-pinned all-body scalar atlas then exhausts every first label in
`1680..1742`.  Its only `58` eligible rows occur at

```text
1683:1, 1694:10, 1702:3, 1708:14,
1722:11, 1724:2, 1732:2, 1736:15.
```

Thus every integer height `1737..1742` is empty without interpolation.  A
self-contained replay against the current canonical uniform-ray and
projected-cell engines partitions the fifteen rows at `1736` into nine
ordinary rows and six explicit `HIGH-TAIL` rows.  On the ordinary rows the
exact ledger is

```text
929 scalar states = 0 crude kills + 774 common-status kills + 155 residuals;
155 residuals -> 286 literal packets -> 286 projected kills -> 0.
```

Five HIGH rows have strictly negative exact forced-wall scalar maxima.  On
the sixth, `E=(1,8,10,12,13,14)`, the scalar maximum is positive by
`91785/20406470812`, so scalar closure alone is genuinely insufficient.
The positive cutoff `4876247/30609706218` leaves `501` literal labels, only
`13260` reaches the forced `13L/150` wall, and exactly one packet remains:

```text
(1736,1836,2004,2340,13260).
```

Its first thirteen projected cells give the hostile lower control zero;
the fourteenth gives residual `14/17`, margin `849/1547` above `25/91`,
and direct full-carrier projection independently gives mass one.  These
three disjoint routes exhaust all fifteen atlas rows over all distinct later
nonaligned labels, with no finite label horizon.  Since `1733..1735` are
also scalar-empty, no interpolation is used.  The updated ledger and cap are

```text
2,239,789 - 15 = 2,239,774,
z_1<=1732,                 ledger=2,239,774.              (25q5)
```

The cap-`1736` splice pins the complete cap-`1742` verifier and scalar atlas.
The hybrid closure pins that splice and the exact atlas partition; normal and
optimized replays match its stored transcript byte for byte.

For `k=5`, there is a second, Gram-facing derivation.  Pointwise

```text
(m_A-1)_+ >=(2/5)binom(m_A,2),
```

while tensorization and THM-1234 give

```text
integral_(C_E)binom(m_A,2)
 =h sum_(i<j)rho(a_i,a_j)>=44h/273.
```

This reproduces `eta_5=88/1365`; the safe-surplus and pair-Gram views are
the same pressure in two coordinate systems.

The `k=5` finite bank is in fact empty.  Its `4,702` rows split exactly
according to the first-drift excess.

In the `4,084` high-excess rows

```text
delta(z_1)>=88h/1365,
```

apply THM-2893's six-tail first-apex gate to
`R=C_E minus D_(z_1)`.  If `h_R` and `r_R` are its mass and component
count, one of the five aligned labels or `z_2` is at most

```text
floor(36r_R/(7h_R)).
```

The aligned labels are at least `L`, while `(25f)` forces
`z_2>2275L/18627`.  These typed lower bounds immediately close `3,827`
rows.  The other `257` rows leave exactly `42,912` integral `z_2`
candidates.  On every candidate, the cell formula `(25i)` gives an exact
rational prefix lower bound

```text
mu(P_(E,{z_1,z_2}))>=887/1365=1-u_5,                    (25k)
```

contradicting the strict containment inequality `(25h)`.

On the subcritical side, `2,290` first rows have a nonempty finite analytic
interval for `z_2`; exactly `618` also survive the projected-suffix predicate.
For such a row put

```text
g=88h/1365-delta(z_1)>0.
```

The component discrepancy bound forces

```text
z_2<=floor((6r_E/49)/g).
```

Across the `2,290` analytic rows there are `7,218,110` row-labelled `z_2`
candidates.  Exact singleton integration together with the suffix predicate
leaves `194,073` admissible drift pairs, supported on `590` `(E,z_1)` rows
inside the `618`-row suffix bank.  Every one again satisfies `(25k)`.  The
smallest certified prefix margin over both banks is `1/378105`; equality at
the five-comb union cap would already be impossible because `P` is compact
and the aligned union is open.

A separate typed recursion independently checks the high-excess bank:
`39,913` of its `42,912` second-drift rows close at the first aligned gate;
the remaining `2,999` close before a one-label terminal, with no multiplier
above four.  Therefore

> No literal six-body/seven-tail cover has five aligned tails and two
> drifts.

THM-2928 now supplies a genuinely independent closure of the same face by
body/divisor projection and relaxed arithmetic-progression address masks.
The two proofs reverse the quantifier order.  THM-2928 retains the aligned
safe set and asks whether two drift masks cover every selected body address;
`(25g)` first removes the drifts and asks whether the aligned union covers
their existential projected residual.  They are dual cell-address
projections, and both succeed precisely because they retain the address
coordinate erased by singleton and Gram statistics.

More explicitly, let `I(j,u)=1` when the two drifts cover phase `u` in body
cell `j`.  THM-2928 fixes `u` in the aligned safe set and obstructs a full
column `I(j,u)=1` for every `j`.  Here

```text
P={u:there exists j with I(j,u)=0}
```

is the set of non-full columns, and its measure is too large for the aligned
union.  The two proofs are arithmetic column obstruction versus measure of
the De Morgan-dual column set on the same body-by-phase incidence object.

The caps in `(25e)` give a uniform **finite reduction** without bounding
the aligned multiplier set.  Delete the bounded `D_(z_1)`.  Six tail labels
remain, with no alignment assumption needed.  At every proper node of their
literal residual tree, the body speeds and already chosen tails total at
most twelve speeds.  Settled `LRC(<=13)` gives a point at clearance at least
`1/13`; the strict margin `1/13-1/14` gives a positive interval in that
residual.  It therefore satisfies the positive-mass interval hypothesis of
the cap-free `p<=6` first-apex recursion THM-2893(7a)--(7b).  Recursing to
depth at most six produces a finite exact decision tree.  This is a
finiteness theorem, not an emptiness census.  At each node the component
estimate inherited from THM-1094 is non-strict with coefficient `6r/49`;
replacing that coefficient by the explicit larger rational `6r/49+1`
supplies the strict form required by THM-2893 without changing finiteness.
The `1/13-1/14` margin is what keeps every proper-node residual positive.

Together with THM-2928's fully aligned and one/two-drift closures, this
proves:

> Every six-body/seven-tail branch with at least five aligned tails is empty.
> Each branch with `k=2,3,4` aligned tails is uniformly reducible to a finite
> exact decision tree.
> Consequently any sector not yet known to admit such a finite reduction
> has at most one aligned tail, hence at least six drifts.

THM-2928's later divisor-status transport and local `98/99` needle terminal
run the `k=4` census and prove it empty.  Thus the current composition
strengthens the first two conclusions to: every `k>=4` branch is empty, and
only `k=2,3` remain finite-but-uncensused.

The next faithful object must therefore join

```text
node data:
  delta_i, residue b_i, normalized slope a_i+b_i/L;
edge data:
  restricted overlap p_ij, gcd/reduced ratio, tooth indices;
component data:
  ordered owner word, transition positions and widths;
quotient data:
  projected residual P_(E,Z), its component word and measure;
hyperedge data:
  the endpoint relation producing the width kappa in (22b).              (26)
```

If tooth endpoints are

```text
(14k+sigma)/(14u),       (14l+tau)/(14v),
sigma,tau in {-1,+1},
```

their signed separation is exactly

```text
((14k+sigma)v-(14l+tau)u)/(14uv).                        (27)
```

Thus vanishing overlap widths are not anonymous analytic errors: they carry
an integer endpoint relation.  Because every `D_w` is open, an exit/entry
coincidence is an uncovered seam, not a two-owner handoff.  Equation `(27)`
records collision address only; it is not a persistence or capacity law.
Any event mesh must use the reduced winding `w/gcd(w,L)`, not raw `w`.
Equations `(16)`, `(25)`, and `(27)` define the multi-slope Gram/address
transition ledger to be classified.

A tournament on the seven labels is not an equivalent quotient.  It can
orient a chosen owner precedence, but it loses symmetric overlap magnitudes,
component order, repeated transitions, ties, endpoint signs, and the
multi-residue clock.  The lawful tournament-related object is a transition
path/graph decorated by the Gram and address sidecars.

Likewise, permuting equal-mass disjoint owners among carrier components
preserves all singleton and pair masses while changing the component owner
word.  This second abstract control isolates the address loss from the
higher-moment loss in the parity example.  Neither control asserts literal
LRC realizability.

## 9. Scope and audit state

The exact scalar census and component bound have two independent
implementations.  The aligned suffix census and five-aligned closure have
independent exact reconstructions and hostile proof audits.  The top-seven
overlap graph is a scoped structural scout.

The consolidated five-aligned closure verifier has LF-normalized
source/output SHA-256
`76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662`
and
`9aecfd75893a537278dcc4e50af7bd45fa2b7925d017748781a18c7163bb716d`;
ordinary and optimized replays are byte-identical.

This theorem does not give a uniform lower bound for `Delta` or `kappa`,
turn the `803` nonpositive actual-top-seven tree margins into certificates,
handle arbitrary packets, finish its `k=2,3` finite decision trees, classify
the zero/one-aligned multi-drift address hypergraph, close the
six-body/seven-tail rung, or prove LRC(14).  The independent THM-2928
divisor-status route now closes `k=4`; `k=2,3` remain open.
