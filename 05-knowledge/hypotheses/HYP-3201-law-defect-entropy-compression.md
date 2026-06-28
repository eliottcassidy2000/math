---
id: HYP-3201
title: Law-defect entropy extends function compression beyond commutativity
status: SYNTHESIS / exact finite scout; not an LRC14 proof
source: codex-2026-06-28-S279
tangent: T1301
technique: LTI-301
tournament_technique: LTT-201
script: 04-computation/lrc14_law_defect_entropy_compression_codex_20260628.py
result: 05-knowledge/results/lrc14_law_defect_entropy_compression_codex_20260628.out
reflection: 07-reflections/law-defect-entropy-function-compression-codex-s279.md
related:
  - HYP-3200
  - HYP-3162
  - HYP-3199
  - HYP-3161
  - HYP-3160
  - HYP-3154
  - HYP-3153
  - HYP-3152
  - HYP-3151
  - HYP-3150
  - HYP-3147
  - HYP-3146
  - HYP-3142
  - HYP-3140
  - HYP-3132
  - HYP-3122
  - HYP-3109
  - HYP-3092
  - THM-577
  - OPEN-Q-108
---

# HYP-3201: Law-Defect Entropy Compression

## Claim

HYP-3150/HYP-3151 say that a compression `q:X->Y` is proof-legal for a target
function `f:X->Z` only when `f` is constant on the fibers of `q`, or when the
destroyed coordinate is retained as a sidecar.

This continuation turns algebraic laws into the same object:

```text
law L == target function f factors through the quotient q_L
failure(L) == H(f | q_L) > 0
```

The residual conditional entropy is not a metaphor.  It is the number of
sidecar bits left in the fiber after the quotient has been applied.  Thus
commutativity failure is ordered-sidecar debt, associativity failure is
bracketing-sidecar debt, idempotence failure is multiplicity debt,
distributivity failure is context debt, and the K4 flip quotient carries
action/representative debt because it is a transformation/relational monoid,
not a Klein-four group action.

## Incoming Mainline Integration

While this lane was closing, mainline claimed nearby namespace:

```text
HYP-3152: Lee-Yang circle coverage radius / Galois S4 web
HYP-3153: Lee-Yang/Worpitzky/quartic compression packet scout
HYP-3154: Joukowski/De Moivre bridge from LRC circle zeros to
          tournament real-rooted coordinates
HYP-3160: k=8 variance/covariance extremality; entropy ruled out as the
          scalar maximizer, but non-associative/Worpitzky sidecar confirmed
HYP-3161: exact covariance robustness and the correction that 1/7 is not a
          theorem, only a near-coincidence at the crossing
HYP-3162: Joukowski/cyclotomic ideal; cap is the 7th-cyclotomic target and
          the dip is a rational-approximation / real-rootedness defect
HYP-3199: n=4 Einheit/minimality chart; fixed-path score table is a cover,
          exact `x,y` chart is the section, and `c` is filler/canary
HYP-3200: exact k=8 cumulant check; no exact 1/7 law, covariance survives as
          primitive degree-two target, odd associator remains sidecar debt
```

This file was therefore renumbered to HYP-3201.  The connection is positive,
not cosmetic: HYP-3160 says plain row entropy is the wrong k=8 extremal
scalar because consecutive/AP has high entropy, while variance/covariance and
bimodality are sharper target functions.  HYP-3201 is compatible with that
correction.  The information-theory object here is not "maximize Shannon
entropy"; it is `H(target | quotient)`, the residual entropy left by an
attempted compression.  Thus "entropy ruled out" for row extremality and
"conditional entropy detects illegal quotienting" can both be true.

Two guardrails are sharpened by the incoming work.  First, do not turn the
`1/7` associativity-defect smell into a theorem; HYP-3161 gives broad
countervalues and HYP-3200 makes the bounded-bank refutation exact.  Second,
use HYP-3199's exact n=4 cover-versus-section chart
for the monoid discussion: the law-defect packet measures the entropy of
forgetting the `S` representative, while the Einheit chart supplies the legal
section when `x,y` are retained.

Third, HYP-3162 sharpens the root-circle sidecar.  A legal Lee-Yang circle
compression is not just "one radius"; it preserves the 7th-cyclotomic
Joukowski/de Moivre target.  The off-circle dip is the residual
real-rootedness / rational-approximation defect that must remain named until a
proof consumes it.

## Exact Scout

Artifact:

```text
04-computation/lrc14_law_defect_entropy_compression_codex_20260628.py
05-knowledge/results/lrc14_law_defect_entropy_compression_codex_20260628.out
```

Finite readout:

```text
Commutativity quotient q(a,b)={a,b}, domain 1..7:
  a+b          residual 0
  a*b          residual 0
  a^b          residual 0.816327 bits
  {a^b,b^a}    residual 0

Associativity quotient q(tree)=ordered leaves (a,b,c):
  (a+b)+c      residual 0
  (a*b)*c      residual 0
  (a-b)-c      residual 0.800000 bits
  (a^b)^c      residual 0.515625 bits

Idempotence/multiplicity quotient q(word)=support(word):
  max(word)    residual 0
  sum(word)    residual 0.666667 bits

Distributivity rewrite quotient over Z5:
  a*(b+c)      residual 0
  a+(b*c)      residual 0.640000 bits
```

The K4 fixed-path class quotient has a concrete action defect.  Collapsing the
eight exact states to `T,+,-,S`, then flipping one live bit, gives:

```text
H(next_class | class,generator) = 0.701205 bits

flip a: T->{+}, +->{T}, (-)->{S}, S->{-,S}
flip b: T->{-}, (-)->{T}, +->{S}, S->{+,S}
flip c: T->{S}, +->{S}, (-)->{S}, S->{T,+,-,S}
```

So the prompt's correction is exactly right: the flip action is not literally
`V4` after quotienting.  The quotient action is deterministic only with
representative/canary sidecars; otherwise the visible class dynamics are a
transformation/relational monoid packet.

## Lee-Yang / Pascal / Phi4 Reading

The same quotient law explains the Lee-Yang and cap language:

```text
all roots on |z|=R  =>  q0 = q6 * R^6
```

This is a lawful compression of the root multiset to one radius.  The scout's
toy radius packet records:

```text
circle:          product_radius=1.000000, log_radius_variance=0
phi4_off_circle: product_radius=1.045476, log_radius_variance=0.05192603
```

Thus root circularity is a zero-defect quotient, while off-circle spread is a
sidecar.  In the LRC language, the binomial/Pascal cap is the exchangeable
pair-mass compression; de Moivre-Laplace is the Gaussian shadow of that lawful
exchangeable quotient; the k=8,9 dip is the higher-order residual that the
pairwise quotient cannot erase.  HYP-3122's phi4 stabilizer and HYP-3109's
root-curve packet are the same rule in analytic clothing: do not compress the
whole PGF to `q0`, `cap`, or a radius unless the root/cumulant sidecar has
zero defect or is explicitly carried.

## Proof-Route Synthesis

This law-defect view merges the user's listed lanes:

- `L_y = p0 + p6 + (1/10)p3` at k=8 is a target function.  Any quotient used
  to prove consecutive/AP extremality must show that this function, not just a
  neighboring scalar, is fiber-constant or sidecar-restored.
- The k=8 dual remains in the solvable window because the even side descends
  to `u^4-5u^2+4`, hence to a quadratic in `u^2`; the odd Worpitzky side is
  the ordered/bracket-like residual, not something the even fold erases.
- The Galois statement should be used precisely: degree <=4 gives a solvable
  algebraic page.  It does not by itself prove the quotient legal; the missing
  odd/canary/root sidecar still has to be discharged.
- Ear decompositions fit the same typing.  Strong connectivity compressed to
  "has an ear decomposition" preserves reachability; factor-critical compressed
  to "has an odd ear decomposition" preserves an odd-parity sidecar.  The
  LRC/Omega bridge should ask which odd-cycle parity data survives the ear
  quotient and which collision/root event it destroys.
- Newton/Maclaurin quartic moment inequalities become moment-compression
  laws.  The AP/consecutive extremal claim is the statement that the fourth
  cumulant/root-spread sidecar has the right sign and no hidden residual on
  the live fiber.

## Packet Fields

Extend HYP-3151's packet schema with:

```text
law_id
law_quotient_map
target_function
residual_entropy_bits
sidecar_type
sidecar_zero_defect_status
law_failure_examples
root_radius_variance
action_determinism_status
monoid_or_group_status
moment_cumulant_sidecar
terminal_discharge_or_named_debt
```

The useful rule is:

```text
residual_entropy_bits = 0    -> quotient may be legal for this target
residual_entropy_bits > 0    -> keep the named sidecar or stop the route
```

## Tournament Analysis

Vertices were proof carriers, not runners or arcs.  The session explicitly
considered functions, law quotients, bracket trees, supports, distributive
contexts, K4 class-action fibers, PGF roots, phi4 modes, resolvent variables,
Fourier/moment modes, ear events, and proof obligations.

Pairwise observable: retained proof payload under quotient legality, sidecar
explicitness, LRC transfer, root/function retention, and degree control.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
Hamiltonian_path_count=1
priority_path =
  typed_factor_through_schema
  -> law_defect_entropy_bits
  -> root_curve_phi4_radius_sidecar
  -> k8_odd_worpitzky_orientation_sidecar
  -> k4_transformation_monoid_sidecar
  -> associativity_bracket_sidecar
  -> distributivity_context_sidecar
  -> idempotence_multiplicity_sidecar
  -> raw_scalar_value
```

The quotient preserves the LRC predicate only when the target proof function
is constant on its visible fiber or a named sidecar reconstructs the hidden
coordinate.  It destroys orientation, bracketing, multiplicity, context,
representative action, root spread, and odd parity unless those are named.

## Next

1. Put `residual_entropy_bits` and `sidecar_type` on one live HYP-3140
   fiber-PGF row, one HYP-3141 edge-witness row, and the HYP-3142 k=8 moment
   packet.
2. For k=8, separate the even biquadratic page from the odd Worpitzky page by
   the target function `L_y`, not by scalar convenience.
3. Test whether AP/consecutive rows minimize root-radius variance and
   law-defect entropy simultaneously in the bounded bank.  A positive result
   would make "AP is extremal" a typed zero-defect/least-defect theorem rather
   than a bundle of analogies.
4. Build the root-collision ear graph with odd-ear parity labels and measure
   which parity sidecar survives contraction to the zero-real component.
