---
id: HYP-2810
title: LRC14 generalized-doublet audit, proof-history maturity map, and sieve-frontier bridge
status: OPEN rigor guardrail; exact bounded-span evidence supports the reframe; proof not claimed
source: codex-2026-06-22-S78
depends_on:
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

# HYP-2810: Generalized Doublet As The Mature LRC14 Wide Route

## Claim

The current best LRC14 genuine-wide route should be stated as an addressed
generalized-doublet theorem, not as a raw slack-floor or sampled maximizer
claim:

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
(bounded base B, far pair {M,M+g}, full-set gcd, remove-one witnesses,
 frozen slow coordinate y, finite-window row).
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

The shortest plausible finish is:

1. **Far-count domination lemma.**  Prove cap-dangerous genuine-wide rows have
   `r=2`, or prove the weaker statement that every `r>=3` row has direct cap
   slack by a survival/live-depth margin.  HYP-2708's live-depth kernel is the
   likely human proof currency: adding the third far element makes too many
   high-tail debts nonnegative for the row to remain cap-dangerous.
2. **Generalized frozen-room lemma.**  Extend HYP-2806 from adjacent `g=1` to
   arbitrary fixed `g`:
   `s_{M+g}(y)=s_M(y)+floor(g*y + frac(M*y)) mod 7`, giving an exact
   one-fast-phase frozen law.  Then prove the frozen plateau is below `cap_k`
   for all bounded bases and gaps, with worst gaps small.
3. **Uniform P/R or TV tail.**  Generalize THM-564's
   `M*(p0-Phi)=P+R` to gap `g`.  The periodic part is still built from
   THM-563 endpoint packets; the residual is one-fast-phase Koksma/variation.
4. **Finite bridge atlas.**  Exact-check the remaining small `M,g` window,
   including even-AP plus odd-bridge rows.  This is where HYP-2807's k=12 row
   lives; it is cap-safe and should be treated as a named finite resonant
   family, not as a slack-floor failure.
5. **Sieve residual projection.**  Import the below-frontier sieve lens: for
   the actual `k=13` frontier, characterize the finite improper/projection
   survivors by the same `(B,g,M,bridge parity)` address quotient.  If the
   public lifting/backward-projection machinery leaves only generalized
   doublets plus a finite bridge atlas, the local cap proof and the global
   finite-check proof meet at a single certificate.

If this route works, HYP-2684's broad BV/Fourier lemma is no longer needed in
full generality.  It becomes the conceptual explanation for why the P/R tail
and the survival/live-depth margins are small.

## Assumption Challenge

The useful tournament vertices in this session are not runners or arcs.  The
candidate vertex sets considered were:

- row profiles `(far_count, gap bucket, parity profile)`;
- exact bounded-span row witnesses;
- proof-era strategies from the repo history;
- hidden addresses `(B,g,M,primitive(FULL E), remove-one scaffold)`;
- frozen-law terms and P/R-tail obligations.

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

At the proof-strategy level the matured tournament is:

```text
addressed_generalized_doublet
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
- a proof-history synthesis explaining why the current route should retain
  `(B,g,M)` addresses;
- a bridge from the 2026 below-frontier lifting/backward-projection literature
  to the repo's cap/genuine-wide address quotient;
- a sharper next theorem statement for finishing the genuine-wide branch.
