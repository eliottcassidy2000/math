# Sixty clocks, four-dimensional Kakeya, and height-sensitive fibres

**Session synthesis, 2026-08-24.**  The proved promotion is
[THM-4035](../01-canon/theorems/THM-4035-sixty-clock-separation-and-finite-kakeya-spine.md).
This reflection keeps the broader research portfolio and the stopping
boundaries.  Euclidean Kakeya in dimension four and LRC(14) remain **OPEN**;
Sun's 2-4-6-8 conjecture is **REFUTED** by THM-4026.

## Inheritance pass

- **Closest proved mechanism:** THM-4029's twelve rational owners.  Each
  denominator is at most six, so its selector is periodic and the joint tail
  phase is `lcm(1,...,6)=60`.
- **Canonical hostile:** phases `0` and `15` have the same scalar pair
  `(F_r mod 10,T_r mod 30)=(0,0)` but different AP laws, already separated by
  `A_1=37/14` versus `31/14`.
- **Corrected near miss:** MISTAKE-489 records that a finite-table endpoint
  was mistaken for a span crossing; the repaired row is
  `(7,8,10,13,26,infinity)`.
- **Typing boundary:** THM-4029/4035 make the separate distinction that the
  phase templates are 60-periodic, evaluated values retain unbounded height,
  and the owner-denominator lcm—not Sturmian aperiodicity—causes the `60`.
- **Least-used relevant sidecar:** HYP-2235's pinned/concurrency carrier.  A
  direction list is not a Kakeya set; line placement and multiplicity remain
  independent coordinates.
- **External boundary:** corrected Katz--Zahl gives the general `R^4`
  Hausdorff lower bound `3.059849573...`, not dimension four.  The endpoint
  four-linear theorem controls only transverse families.  See the
  [primary-source audit](../05-knowledge/reference/CORE-PAPERS-KAKEYA-4D-2026-08-24.md).

## Live concept board

| object | representation | invariant | operation | quotient loss |
|---|---|---|---|---|
| AP-cover tail | 60 rational functions `L_r(n)` | two-moment fingerprint `(A_1,A_2)` | phase shift and sparse deletion | scalar deficit loses owner, track and gap |
| Fibonacci | state `(F_r,F_(r+1)) mod 10` | Pisano orbit length 60 | companion-matrix advance | `F_r mod 10` has fibres of size 4 or 8 |
| triangular numbers | affine state `(r,T_r) mod 30` | minimal period 60 | triangular cocycle | `T_r mod 30` has only 12 images |
| finite `4D` directions | points of `P^3(F_61)` | ranks of 3- and 4-minors | cyclic diagonal action, translation of lines | projectivization and direction-only data lose scale/basepoint |
| Sun 2-4-6-8 fibres | role-labelled binomial inputs | exact representation multiplicity | residue projection and bounded lifting | local support and average order lose height alignment |

Every new move was compared with all five rows.  The common lesson is not
“everything has period 60.”  It is that a reversible address can coexist with
a badly lossy consumer observable.

## Exact common clock, different evaluators

Three pointed template/state systems are exactly conjugate to successor on `C_60`:

```text
AP phase r             -> full rational law L_r;
r                      -> (F_r,F_(r+1)) mod 10;
r                      -> (r mod 30,T_r mod 30).
```

The AP entry is the cycle of phase templates.  Its evaluated value
`L_(n mod 60)(n)` retains unbounded height `n` and is not periodic.

The causes differ:

```text
AP owner clock:          lcm of denominators 1,...,6;
Fibonacci clock:         matrix order lcm(3,20) modulo 2*5;
triangular clock:        affine cocycle period 2*30.
```

The full AP laws are all distinct.  The first two tail moments recover the
phase, but the scalar Fibonacci and triangular values do not.  Adding the two
scalars still leaves twelve doubletons; parity is exactly the missing bit.
Thus Fibonacci can relabel the AP tail, but it does not derive or explain its
owner coefficients.  Any pointed 60-cycle could perform that relabelling.

## Why `p=61` is the useful Kakeya test bench

In `F_61`, the golden root `phi=44` has order 60 and

```text
phi^r=F_(r+1)+43F_r mod 61.
```

Three independent clocks give the exact nonzero projective chart

```text
(a,b,c) -> [1:phi^a:phi^b:phi^c]
C_60^3  -> P^3(F_61),
```

covering `216000` of `230764` directions.  The missing `14764` are boundary
charts with a zero coordinate.  This dimension count is decisive: one shared
clock among AP, Fibonacci and triangular observables is still one parameter,
not three independent direction coordinates.

The cheapest broad/narrow probe keeps the same clock and changes only its
representation:

```text
B(r)=[1:t:t^2:t^3],       N(r)=[1:t:1:t],       t=phi^r.
```

Every one of the `C(60,4)=487635` four-minors of `B` is nonzero by the
Vandermonde formula; every four-minor of `N` is zero because `N` lies in a
fixed vector two-plane.  Period, phase injectivity and cyclic equivariance
survive both maps.  Four-way transversality changes completely.  Therefore a
Kakeya consumer needs the weight/representation sidecar, not another name for
the phase.

The projective Fibonacci orbit supplies a second hostile: the matrix has
order 60 over `F_61`, yet `A^15=11I`, so its projective orbit has only 15
directions.  The triangular multiplicative shadow
`phi^(2T_r)=phi^(r(r+1))` has minimal period 60 but only 12 values.  Exact
period and direction diversity are different predicates.

## Merge with the Sun counterexample

THM-4026 gives an empty exact integer fibre at

```text
N=896315812331399,
```

while THM-4027 proves every residue modulo every modulus is represented and
THM-4028 proves growing mean representation count `~V X^(1/24)`.  This is the
additive analogue of the direction/basepoint warning:

```text
complete residue support  != exact bounded-height preimage;
complete direction support != a controlled Kakeya union;
periodic scalar address    != owner/transversality data.
```

At `p=61`, the four binomial input roles can be parameterized by the same
nonzero `C_60` clock.  Their image sizes are `(31,24,26,24)`, and the first two
already satisfy `S_2+S_4=F_61`; the target is `21 mod 61`.  Thus the golden
clock is an exact input atlas and an exact hostile to a new modular
explanation.  The remaining information is bounded-lift cost, role-labelled
fibre multiplicity and archimedean height.

This suggests a more faithful shared object: a **phase-coloured fibre**, not a
phase sequence.  Its source is the full input tuple; its target is an output
residue or direction; its fibre retains owners/basepoints and lift heights.
Residue projection preserves local existence but destroys height.  Direction
projection preserves orientation but destroys placement.  AP aggregation
preserves total covered measure but destroys sparse occupancy.  These are
three instances of the same quotient-loss pattern.

## Anchor / Niche / Wildcard portfolio

### Anchor — LRC(14) sparse-owner response

Start from THM-4029's exact twelve-owner law and replace the full AP container
by a sparse deletion mask.  For each phase retain

```text
(owner p/q, side, last track E_s(n), selected time in E?, gap to next point).
```

The immediate target is an exact deletion-response formula or a minimal
witness showing that no bounded phase moment suffices.  Hostiles are
`E={0,...,6,N}` and two masks with the same 60 phase histogram but different
gap placements.  This is the closest route back to LRC(14), because the known
AP theorem loses precisely sparse occupancy and owner placement.

### Niche — cyclic-MDS atlas in `P^3(F_61)`

Classify weight quadruples

```text
W=(w_0,w_1,w_2,w_3) in (Z/60)^4,
v_W(r)=(phi^(w_0 r),...,phi^(w_3 r))
```

up to phase units, translation of all weights, coordinate permutation and
projective scaling.  Record the complete zero spectra of three- and
four-minors.  Consecutive weights give the full-spark control; repeated
weights give the planar hostile.  The open classification asks which
intermediate weight sets are broad on most quartets and which concentrate in
low-degree subspaces.  This is a finite cyclic-MDS problem, not yet a
Euclidean estimate.

### Wildcard — height-conditioned Sun fibres

On the `C_60^4` nonzero input chart modulo 61, stratify every output residue
by phase-orbit type and then attach the least nonnegative lift height in each
role.  Compare the target class `21` with shuffled-phase and equal-local-
density controls.  The goal is a statistic that predicts an empty exact fibre
without contradicting universal modular solubility.  Any candidate must be
tested against nearby represented integers and against THM-4028's positive
fixed-residue mean.

## New bounded attacks

1. **Fourier owner spectrum on `C_60`.**  Decompose the twelve denominator
   selectors and the exact AP moments into characters.  Identify which modes
   are forced by each `q<=6`, then test which survive sparse deletion.  This
   may isolate the smallest phase sidecar beyond `(A_1,A_2)` that retains
   owner response.
2. **General `K`-sector tail law.**  Test and then prove or refute that
   eventual owners have `q<K`, candidate phase `lcm(1,...,K-1)`, and minimal
   period equal to that lcm.  Cancellation can shorten the period, so owner
   periodicity alone is not enough.
3. **Complete the projective atlas.**  Add the `14764` boundary directions to
   the three-clock torus and select one affine line in every direction.
   Compare concurrent, owner-derived and shuffled basepoints through union
   size and multiplicity spectra.
4. **Broad/narrow stability under owner labelling.**  Attach THM-4029's twelve
   owners to both `B` and `N`; measure whether any owner-derived translation
   statistic distinguishes them before determinants are revealed.  A
   positive signal would identify a usable concurrency sidecar; a negative
   one is a clean stopping result.
5. **Arithmetic-Kakeya calibration.**  The current sum-difference exponent
   `1.6751308705...` yields only `1+3/alpha=2.790904850...` in `R^4`, below
   the general `3.059849573...` bound.  To improve that bound through the
   Katz--Tao formula alone requires `alpha<1.456417031...`.  Keep the repo's
   certificate game as its own frontier; a score near `1.675` is not a new
   four-dimensional theorem.
6. **Planebrush interface audit.**  Try to enrich a finite direction spine
   with scale, shading density, two-ends, moving three-plane labels and Wolff
   nonconcentration.  The decisive test is whether the same enrichment can
   still distinguish the broad and narrow controls across nested scales.
   Without those fields, no Euclidean transfer is even typed.
7. **Phase-fibre variance for Sun 2-4-6-8.**  Refine THM-4028's fixed-residue
   average by bounded input-phase classes and measure the second moment of
   representation counts.  The desired theorem is a lower-tail or exceptional
   set estimate; another first-moment asymptotic cannot rule out zeros.

## Stopping boundaries

The synthesis yields an exact algebraic laboratory and several finite
experiments.  It does not yield a Euclidean tube estimate, an induction on
scale, a line-placement theorem, a pointwise additive-basis lower bound, or
an LRC safe time.  Those are the first missing implications, and every future
claim must attach the corresponding sidecar before crossing them.
