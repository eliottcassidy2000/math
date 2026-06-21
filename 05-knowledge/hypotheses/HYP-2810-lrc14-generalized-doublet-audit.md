---
id: HYP-2810
title: LRC14 generalized-doublet audit as a gK8 concentration guardrail
status: OPEN guardrail; exact bounded-span evidence supports gK8 concentration and generalized-doublet fallback; proof not claimed
source: codex-2026-06-22-S78
depends_on:
  - HYP-2816
  - HYP-2815
  - HYP-2814
  - HYP-2812
  - HYP-2811
  - HYP-2807
  - HYP-2809
  - HYP-2806
  - THM-564
  - THM-563
  - THM-557
  - HYP-2799
  - HYP-2684
related:
  - HYP-2808
  - HYP-2805
  - HYP-2775
  - HYP-2788
  - HYP-2708
  - HYP-2694
  - OPEN-Q-108
---

# HYP-2810: Generalized-Doublet Audit As A gK8 Concentration Guardrail

## Claim

After the incoming HYP-2811/HYP-2812 work, the current best LRC14 wide route is
not primarily a `p0` generalized-doublet theorem.  It is the gK8 concentration
route:

```text
L_yK8(E)=10q0(E)+q3(E)+10q6(E) is maximized by bounded/concentrated rows.
max_bounded L_yK8 < 10*cap_k.
Therefore every wide row has p0=q0 <= L_yK8/10 < cap_k.
```

The generalized-doublet route remains the exact structural fallback and sanity
check.  If concentration extremality needs a finite resonant atlas, S78's
evidence says the dangerous genuine-wide `p0` rows are already addressed by
bounded base, far pair `{M,M+g}`, parity/bridge data, full-set primitivity, and
remove-one witnesses, not by a raw slack-floor or sampled maximizer claim:

```text
If E is a cap-dangerous genuine-wide row, then after bounded normalization
its far part has exactly two elements {M, M+g}.

For every bounded base B and gap g, the row B union {M, M+g}
has p0 < cap_k, by a frozen-room plus P/R-tail plus finite-window proof.
```

This is a refinement of incoming HYP-2807.  HYP-2807 is the right conceptual
reframe of the mac-mini k=12 obstruction, but the repo should treat it as an
open theorem until the search domain is exact: full-set primitivity, remove-one
scale witnesses, arbitrary bounded bases, far gaps `g`, and the low-`M` window
must all remain addressed.

## Post-Fetch gK8 Concentration Update

HYP-2812 is the sharpest current target.  Its exact concentration audit reports
that the maximum wide `L_yK8` is below the bounded maximum at k=10,11,12:

```text
k=10: MB=37/7, MW≈4.81347, MW-MB≈-0.47225.
k=11: MB=26603/4410, MW≈5.63223, MW-MB≈-0.40020.
k=12: MB=29287/4410, MW=30494/4851 at E*, MW-MB≈-0.35492.
```

Here `E*=(0,2,4,6,8,9,10,11,12,14,16,18)` is exactly the k=12 odd-bridge row
that broke the old `p0<Q(11)` slack floor.  Under gK8 it is no longer a threat:
`L_yK8(E*)=30494/4851` leaves `10cap-L=11086/4851`.

The proof obligation has therefore moved to a concentration/decorrelation
monotonicity theorem for the missed-sector distribution:

```text
far spreading / decorrelation moves mass away from q0 and q6
and into middle missed-count atoms, so L_yK8 cannot increase.
```

The Mordell-Tornheim side is also sharper after mac-mini S22: the absolute
two-frequency constant in HYP-2808 is

```text
sum_{h,h' != 0, h+h' != 0} 1/(|h| |h'| |h+h'|) = 12*zeta(3).
```

Thus the R-tail route is now a rigorous fallback with a closed-form constant,
while gK8 concentration is the clean route if the monotonicity lemma is proved.

## Post-Rebase q6 Suppression Update

After rebasing over the latest HYP-2814/HYP-2815/HYP-2816 work, the sharp
concentration statement is more specific than "decorrelation moves mass to the
middle."  The current route is:

```text
far insertion on the all-missed atom q6 acts like convolution by Unif(Z/7),
so q6 is multiplied by approximately 1/7 per independent far coordinate;
small far speeds f>=15 have a worse but still strict ratio C(f)<1;
the gK8 weight 10q6 makes this lost all-missed mass pay for any q0 gain.
```

This explains the exact S78 overlay: the k=12 row `E*` can increase `q0=p0`
above the old `Q(11)` floor, but it has already spent enough `q6` that
`L_yK8(E*)` remains below the bounded maximum.  Thus the proof obligation is
now narrower:

```text
Prove a generated-row Krawtchouk/majorization lemma:
for every wide row address (B, far data), the admissible q-profile after
far smoothing is dominated in L_yK8 by a bounded/concentrated profile.
```

The honest small-`f` refinement matters.  The clean `1/7` factor is an
asymptotic law; at the binding boundary `f>=15` the ratio can be larger.
So a complete proof should either:

1. prove the worst-case `q6` contraction ratio for `f>=15` inside the same
   endpoint-period/Dedekind machinery as THM-563; or
2. use the closed-form `12*zeta(3)` R-tail and finite-window atlas to isolate
   the low-`f` exceptions, then apply the asymptotic convolution lemma.

## Exact S78 Evidence

Script:

```text
04-computation/lrc14_genuinewide_generalized_doublet_exact_codex_s78.py
```

Stored output:

```text
05-knowledge/results/lrc14_genuinewide_generalized_doublet_exact_codex_s78.out
```

The first exact run over all normalized rows with `0 in E`, `max(E)<=18`,
`primitive(E)`, and genuine-wide remove-one status gives:

```text
k=10: 27484 genuine-wide rows, over_Q=0
  global best r=2:
  E=(0,1,3,5,7,9,11,13,15,17),
  p0=52685/119119, cap-p0=19310/119119.

k=11: 29724 genuine-wide rows, over_Q=0
  global best r=2:
  E=(0,2,4,6,7,8,10,12,14,15,18),
  p0=2257/4410, cap-p0=12239/57330.

k=12: 24816 genuine-wide rows, over_Q=4
  global best r=2:
  E=(0,2,4,6,8,9,10,11,12,14,16,18),
  p0=238949/388080,
  p0-Q(11)=36613/2716560,
  cap-p0=93691/388080.
```

Thus in this exact finite window every global maximum is a generalized doublet
(`r=2`).  The four k=12 over-`Q(11)` rows are also `r=2` generalized doublets
and are very far below `cap_12`.  This supports the HYP-2807 reframe while
correcting its certainty language: the new exact audit is a finite-window
certificate, not yet the unbounded theorem.

The completed exact run over `max(E)<=20` strengthens the same picture:

```text
k=10: 134502 genuine-wide rows, over_Q=0
  global best unchanged r=2:
  E=(0,1,3,5,7,9,11,13,15,17),
  p0=52685/119119, cap-p0=19310/119119.
  best r=3:
  E=(0,3,5,7,9,11,13,15,17,19),
  p0=2839216/6789783.

k=11: 162633 genuine-wide rows, over_Q=0
  global best r=2:
  E=(0,2,4,6,7,8,10,12,14,17,20),
  p0=25909/49980, cap-p0=134423/649740.

k=12: 156939 genuine-wide rows, over_Q=7
  global best unchanged r=2:
  E=(0,2,4,6,8,9,10,11,12,14,16,18),
  p0=238949/388080,
  p0-Q(11)=36613/2716560,
  cap-p0=93691/388080.
```

At span 20, every displayed over-`Q(11)` row in the k=12 top list is still a
generalized doublet (`r=2`).  The first non-doublet profile is the `r=3` row
`(0,2,4,6,8,9,10,12,14,15,16,18)`, already below `Q(11)`.  The HYP-2807 k=11
sampled witness with max speed `22` is outside this span-20 audit; it should
be rechecked in a later span-22 run before any k=11 global-window language is
used.

## Exact gK8 Overlay

Script:

```text
04-computation/lrc14_genuinewide_gK8_overlay_codex_s78.py
```

Stored output:

```text
05-knowledge/results/lrc14_genuinewide_gK8_overlay_codex_s78.out
```

This overlay recomputes the same `span<=18` normalized primitive genuine-wide
domain by the full missed-sector distribution `q_t`, then ranks rows by
`L_yK8=10q0+q3+10q6`.  It found zero `L_yK8>10cap` violations:

```text
k=10: 27484 genuine-wide rows, L_yK8_violations=0
  max_L row r=2:
  E=(0,1,3,5,7,9,11,13,15,17),
  L_yK8=393157/85085,
  10cap-L=17299/12155.

k=11: 29724 genuine-wide rows, L_yK8_violations=0
  max_L row r=2:
  E=(0,2,4,6,7,8,10,12,14,15,18),
  L_yK8=92777/17640,
  10cap-L=457099/229320.

k=12: 24816 genuine-wide rows, L_yK8_violations=0
  max_L row r=2:
  E=(0,2,4,6,8,9,10,11,12,14,16,18),
  L_yK8=30494/4851,
  10cap-L=11086/4851.
```

The k=12 `p0` breaker is also the span-18 `L_yK8` leader, but its gK8 margin is
large.  This overlay is a finite-window check only; its value is that it ties
the S78 generalized-doublet address audit to the newer HYP-2812 concentration
route using the same exact row domain.

## What The Repo History Says

The LRC14 proof strategy has matured by repeatedly restoring hidden addresses
that scalar summaries forgot:

1. **Consecutive/AP maximizer era.**  Early HYP-2602/HYP-2607 style targets
   tried to make the consecutive row globally extremal.  HYP-2780 killed this
   at k=12, but the cap inequality survived with large room.
2. **Uniform decorrelation-error era.**  HYP-2684 and HYP-2775 showed the right
   Fourier/Weyl architecture, but also showed that independent suprema are false
   currency: large error and small cap room live in different families.
3. **Regime split era.**  HYP-2788/HYP-2793 separated bounded, single-far, and
   genuine-wide rows.  THM-563 then closed the single-far branch exactly by a
   signed endpoint-period/Dedekind ledger.
4. **Frozen-doublet era.**  THM-564, HYP-2799, HYP-2805, and HYP-2806 found the
   correct center for two far elements: the frozen plateau over `y in [0,7)`.
   The key addresses became `(B,g,M)`, full-set primitivity, and the finite
   low-`M` window.
5. **Generalized-doublet era.**  mac-mini S7 found the k=12 over-`Q` row.
   HYP-2807 observed that it is not a new regime but a spread doublet with
   gap `g=2`.  S78 adds exact finite-window evidence and the guardrail that
   generalized-doublet dominance still needs exact-domain proof.

The common lesson is: do not scalarize before the address quotient is stable.
For the current branch the stable address appears to be

```text
(missed distribution q, gK8 moment, bounded base B, far pair {M,M+g},
 full-set gcd, remove-one witnesses, frozen slow coordinate y,
 finite-window row).
```

## External Literature Checkpoint

The current public below-frontier checkpoint is Sungkawichai-Trakulthongchai,
`arXiv:2604.23906`, submitted 2026-04-26.  It proves the LRC cases
`k=10,11,12` in the positive-integer/stationary-runner convention, i.e. up to
thirteen total runners, using Rosenfeld/Trakulthongchai finite checking refined
by new sieves.  The method is not our cap route, but its vocabulary is directly
useful here:

- define proper/improper speed tuples modulo a prime `p`;
- reduce by permutation, sign flip, and unit-multiplier equivalences;
- compute a finite improper set and repeatedly apply lifting/backward
  projection;
- use a polynomial-method proposition for tuples congruent to `(1,2,...,k)`
  modulo `p` when `k+1` and `p>k^2+k` are odd primes.

The paper explicitly identifies efficient computation of the next improper set
as the bottleneck for extending to `k=13` (fourteen total runners).  Its code is
published at `https://github.com/vzsky/13-lonely-runners`.

S78's bridge conjecture is: the cap-dangerous generalized-doublet / even-AP
odd-bridge rows in this repo are the local geometric shadow of the small
surviving improper ansatz families in the public sieve.  A practical finish
could therefore be hybrid:

```text
prime-sieve residual
  -> address quotient by base/gap/bridge/parity profile
  -> far_count r=2 or r>=3 margin split
  -> frozen-room plus P/R tail plus finite cap atlas.
```

This is a bolder, but more checkable, target than a uniform Weyl-error bound:
prove that any k=13 ansatz survivor projects to one of the finite addressed
families already visible in the cap audit.

## Bold Finish Guess

The shortest plausible finish has changed:

1. **Concentration extremality for gK8.**  Prove that far spreading and wide
   decorrelation cannot increase `L_yK8=10q0+q3+10q6`.  The strongest current
   form is a `q6`-suppression theorem: far coordinates convolve the all-missed
   atom toward a `1/7` factor, while generated-word/Krawtchouk constraints keep
   the possible `q0` gain from overtaking the lost `10q6` mass.  This is a
   convex-order or Schur-type statement on generated missed-sector profiles,
   not on arbitrary probability vectors on `q_0,...,q_6`.
2. **Bounded maximum certificate.**  Use the existing bounded gK8 finite check
   and Lean-facing Delsarte certificate: the bounded maximum `MB` is already
   below `10cap_k`.
3. **Finite near-boundary / R-tail bridge.**  If concentration is first proved
   only asymptotically, use HYP-2811 plus the closed-form `12*zeta(3)` R-tail
   constant to reduce the remaining low-`M` window.
4. **Generalized-doublet fallback atlas.**  If the monotonicity proof has
   resonant exceptions, route them through the S78 address quotient
   `(B,g,M,bridge parity, primitive(FULL E), remove-one witnesses)`.
5. **Sieve residual projection.**  Compare any surviving finite atlas with the
   below-frontier lifting/backward-projection sieve, whose next bottleneck is
   computing the improper ansatz set for the fourteen-runner frontier.

If this route works, the generalized-doublet work is not discarded.  It becomes
the finite-address layer beneath a stronger gK8 concentration theorem.

## Assumption Challenge

The useful tournament vertices in this session are not runners or arcs.  The
candidate vertex sets considered were:

- row profiles `(far_count, gap bucket, parity profile)`;
- exact bounded-span row witnesses;
- proof-era strategies from the repo history;
- hidden addresses `(B,g,M,primitive(FULL E), remove-one scaffold)`;
- frozen-law terms and P/R-tail obligations.
- missed-sector distributions `q_t` and the gK8 moment `10q0+q3+10q6`;
- concentration/decorrelation proof states.

The quotient preserves the direct predicate `p0<cap_k` and the identity of the
current proof obligation.  It destroys the stronger false predicate
`p0<Q(k-1)` and the false consecutive-maximizer target.

## Tournament Analysis

The S78 script builds a finite tournament over profile buckets.  The observable
is exact profile maximum `p0`; the switch orients the edge toward the higher
maximum, with lexicographic tie path.  In the span-18 audit the Hamiltonian
path always begins with an `r=2` profile for k=10,11,12, and the k=12 over-`Q`
profiles are still direct-cap safe.
The span-20 audit keeps that order: the k=12 profile path begins with
`('r2','g2','even10_odd2_farodd0')`, and every over-`Q` top row is `r=2`.

The gK8 overlay builds the same tournament with observable max `L_yK8`.  In
the span-18 overlay the Hamiltonian path still begins with an `r=2` profile for
k=10,11,12, with zero `L_yK8` cap violations.

At the proof-strategy level the matured tournament is:

```text
gK8_concentration_extremality
> decorrelation_middle_mass_majorization
> bounded_LyK8_certificate
> addressed_generalized_doublet_fallback
> frozen_room_plus_PR_tail
> survival_live_depth_r>=3_margin
> global_BV_Fourier_error
> uniform_Q_floor
> consecutive_maximizer
> scalar_residue_or_conductance_majorization
```

The edge flips encode the repo's history: `uniform_Q_floor` beat the old raw
decorrelation story until mac-mini S7; `addressed_generalized_doublet` now beats
`uniform_Q_floor` because it preserves the finite k=12 odd-bridge address.

## Status

No LRC14 proof is claimed.  HYP-2810 contributes:

- an exact finite-window audit supporting generalized-doublet dominance at
  `span<=18` and `span<=20`;
- an exact `span<=18` gK8 overlay showing the same genuine-wide window has
  zero `L_yK8` cap violations, with r=2 leaders;
- a post-rebase proof-order refinement: the live concentration lemma should
  prove `q6` contraction under far smoothing, with a small-`f` endpoint-period
  guard and the `12*zeta(3)` R-tail as fallback;
- a proof-history synthesis explaining why the current route should retain
  both the gK8 missed-sector moment and `(B,g,M)` addresses;
- a bridge from the 2026 below-frontier lifting/backward-projection literature
  to the repo's cap/genuine-wide address quotient;
- a sharper next theorem statement: prove gK8 concentration/decorrelation
  extremality, with generalized doublets as the finite resonant fallback.
