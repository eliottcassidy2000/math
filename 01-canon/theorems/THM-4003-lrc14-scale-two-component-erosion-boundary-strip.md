---
id: THM-4003
title: "LRC(14) scale-two component erosion and the two-boundary body gain"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. In the THM-3878/3995 conditional scale-two `(2,1,9)` slice, the
  eleven distinct body speeds and THM-3818's internal pair-height cap force a
  safe component of length at least `(2H-1)/(84U(H-1))`, where
  `H=min(U,91^6)`. Combining it with THM-3995's two separately trimmed
  obstruction components gives a new necessary scale inequality. Its
  distinctness corollary is `3t(2U-1)<=8(U-1)^2`; for odd t this is exactly
  `U>=floor(3t/4)+2+1_(t=1 mod 4)`. Exact directed owner residues close 77
  further `(t,U)` cells through `t<=1001`. Orthogonally, if V is the
  second-largest body speed, a top-balance selector closes all 17 conditional
  certificate types when `U/V>=11`, and closes scale two already at
  `U/V>=1001/189`. This does not prove LRC(14).
source: root + body_decoder_geometry / LRC14 long session, 2026-08-24
audit: >
  INDEPENDENT MATHEMATICAL AUDIT PASS. The body-component case split,
  pair-height upgrade, compact/open containment, equality convention, odd
  threshold and overlap with THM-4002 were independently rederived. A
  separate exact path reconstructed the four oriented owner-wall grids and
  audited the finite `(t,U)` census through 1001. Normal and optimized primary
  streams and the independent stored stream byte-match their frozen outputs.
depends_on:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-3995-scale-two-parity-hole-support-and-integer-variance-tariff
  - LRCUpTo13
related:
  - THM-4002-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family
  - THM-1002-quantitative-top-gap-margin-deathstar-S56
  - THM-759-tight-instance-ratio-bound
  - HYP-2906-lrc-one-large-speed-interval-peeler
script: 04-computation/lrc14_scale_two_component_erosion_boundary_strip_thm4003.py
output: 05-knowledge/results/lrc14_scale_two_component_erosion_boundary_strip_thm4003.out
independent_audit_script: 04-computation/lrc14_scale_two_component_erosion_boundary_strip_independent_audit_thm4003.py
independent_audit_output: 05-knowledge/results/lrc14_scale_two_component_erosion_boundary_strip_independent_audit_thm4003.out
script_sha256: f5ed4db4484be206cc17808592654ca2bcb1944e2e0b8c673dc2eaac91123f35
output_sha256: 7e4d2a4fb2ac9d0b841ea8ed6cf396d7da16b49d6224dfb9e95455c40bc5a190
independent_audit_script_sha256: 46825ac964d4154d0519a1acf700595a7f88c9e1b3b42caa87d65d6d4f377b8e
independent_audit_output_sha256: 8ca1c3e33ee3125cd5648c74db516ad56fbccff775cd7360631dcf8f65f27c1e
hash_basis: raw LF bytes
---

# THM-4003 -- scale-two component erosion and the two-boundary body gain

**PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** This is an arbitrary-body refinement inside one conditional
certificate type. It is not a proof of LRC(14).

## 1. Typed setting

Retain THM-3818/3878's exact rank-eleven, two-component `11+2` branch and
consider only the conditional scale-two survivor

```text
n=2(u_1,...,u_11) direct-sum t(1,9),
U=max_i u_i,                 U<=t<4U/3.                (1)
```

The body coordinates are eleven distinct positive integers. Coprimality of
the two component scales makes `t` odd. Put

```text
G(u)={y in R/Z:min_i ||u_i y||>=1/14},
N_t(w)=sum_(a=0)^(t-1) 1_G((w+a)/t).                  (2)
```

For a failure, THM-3878 puts the positive sheet count inside the exact open
quotient obstruction

```text
C=(2/21,8/63) union (55/63,19/21).                    (3)
```

THM-3995 strengthens this eventwise. Modulo null owner walls, each positive
component lies in one of the closed retained cores

```text
K_1=[2/21+1/(42U), 8/63-1/(126U)],
K_2=[55/63+1/(126U), 19/21-1/(42U)].                  (4)
```

Both cores have length

```text
beta_U=2(U-1)/(63U).                                  (5)
```

The separation into two components, rather than only their total support,
is the first load-bearing coordinate.

## 2. Two body boundaries force a longer safe component

Cited `LRCUpTo13`, applied to the eleven body speeds, supplies a point `y_0`
at which every body clearance is at least `1/12`. Let `J` be the connected
component of `G(u)` containing `y_0`. Each endpoint of `J` is a `1/14` wall
of at least one body owner. If owner `w` supplies an endpoint, the
`w`-Lipschitz bound gives distance at least

```text
(1/12-1/14)/w=1/(84w)                                 (6)
```

from `y_0` to that endpoint.

There are two exhaustive cases.

1. Distinct owners `w_L,w_R` can be chosen at the two endpoints. Put
   `W=max(w_L,w_R)`, `V=min(w_L,w_R)`, `Q=91^6`, and `H=min(U,Q)`.
   THM-3818 proves that every same-component primitive pair has height at
   most `Q`. Since `W<=U`, the endpoint pair satisfies

   ```text
   U/V>=H/(H-1).                                       (7)
   ```

   For `U<=Q`, integer distinctness gives `V<=W-1<=U-1`. For `U>Q`, put
   `g=gcd(W,V)`. The pair-height bound `W/g<=Q`, together with
   `W-V>=g`, gives

   ```text
   V<=W-g<=W(Q-1)/Q<=U(Q-1)/Q.
   ```

   These are exactly the two cases of `(7)`. Moreover `1/W>=1/U`.
   Therefore

   ```text
   |J|>=1/(84W)+1/(84V)
       >=1/(84U)+1/(84V)
       >=(2H-1)/(84U(H-1)).                            (8)
   ```

2. No distinct cross-endpoint choice is possible. Then both endpoint-owner
   sets are the same singleton `{w}`. The connected arc `J` is the whole
   safe tooth of owner `w`, so

   ```text
   |J|=6/(7w)>=6/(7U),                                 (9)
   ```

   which is strictly larger than the right side of `(8)`.

Consequently every body in this exact branch obeys

```text
lambda_+(u)>=L_(U,Q):=(2H-1)/(84U(H-1)),
U lambda_+(u)>=1/42+1/(84(H-1)),   H=min(U,Q).         (10)
```

The weaker distinctness-only form, valid for every `U`, is

```text
lambda_+(u)>=(2U-1)/(84U(U-1)).                       (11)
```

This improves the single symmetric interval of length `1/(42U)`: the second
boundary cannot reuse the maximum speed unless it pays for an entire tooth.
Against THM-3910's exact sixteen-type beta table, `(11)` closes two small
scale-one height layers: type `(9,11)` for `U<=17` and type `(8,21)` for
`U<=29`. No other surviving scale-one type is reached for physical `U>=11`.

## 3. An orthogonal top-balance selector

Let `V<U` be the second-largest body speed and put `R=U/V`. Remove the maximum
owner. Cited LRC for the remaining ten speeds supplies a phase of clearance
at least `1/11`; the ten-speed core therefore remains `1/14`-safe on an
interval of length

```text
2(1/11-1/14)/V=3/(77V).                               (11a)
```

In the scaled coordinate `z=Uy`, this interval has length `3R/77`. If
`R<22`, it meets at most one open danger tooth of the maximum owner: distinct
teeth are separated by a safe gap of length `6/7`. Removing one tooth costs
at most `1/7` and splits the interval into at most two pieces. The larger
piece is a body-safe component, so

```text
U lambda_+(u)>=max(0,(3R-11)/154),       R<22.         (11b)
```

For `R>=22`, apply the same argument to subintervals whose scaled lengths
increase to `6/7`; compactness gives

```text
U lambda_+(u)>=5/14.                                  (11c)
```

This is a component-length refinement of the usual one-large-speed interval
peeler. THM-1002 retains a quantitative global margin but not the two
surviving pieces, while `(11b)` keeps exactly the component needed by
THM-3910.

In the conditional lane `t>=U`, THM-3910 gives
`t lambda_+(u)<beta` for every failure. Every one of its 57 scale-one input
types has `beta_1<=1/7`, hence all sixteen surviving scale-one types close
when `R>=11`. Scale two has `beta_2(1,9)=2/63`, and solving

```text
(3R-11)/154>=2/63
```

gives the sharper threshold

```text
R=U/V>=1001/189.                                      (11d)
```

Thus all seventeen conditional certificate types close for `U/V>=11`, and
the scale-two type alone closes already at `U/V>=1001/189`. Equality closes:
a failure needs the strict compact-to-open inequality. This prunes body
shapes rather than pair types; the top-balanced residual remains infinite.

## 4. The new infinite boundary strip

Choose a closed subarc of `J` of length `L_(U,Q)`. Its open interior is
strictly body-safe. Under a hypothetical full-row failure, its image by
`w=ty` lies in the essential positive support of `N_t`; its closure must fit
inside one of the two retained cores `(4)`. Choosing the stated subarc makes
its image nonwrapping in the old strip. Hence

```text
t L_(U,Q)<=beta_U,
3t(2H-1)<=8(H-1)(U-1),          H=min(U,Q).            (12)
```

Using only `(11)` gives the universal, slightly weaker corollary

```text
3t(2U-1)<=8(U-1)^2.                                   (13)
```

For odd `t`, `(13)` is equivalent in the old strip to

```text
U>=floor(3t/4)+2+1_(t congruent 1 mod 4).              (14)
```

For `t=4k+1`, the quadratic right-minus-left side changes sign between
`U=3k+2` and `U=3k+3`; for `t=4k+3`, it changes sign between `U=3k+3` and
`U=3k+4`. Its discrete increment is positive through the old admissible
strip. Thus `(14)` removes one infinite old boundary layer for
`t=3 mod 4` and two for `t=1 mod 4`; the second layer is

```text
(t,U)=(4k+1,3k+2).                                    (15)
```

When `U>Q`, the first inequality in `(12)` retains the fixed strict gain
`1/(84(Q-1))` that distinctness alone would let decay with `U`.

Equality in `(12)--(14)` remains allowed. A body-safe wall is safe at exact
clearance `1/14`, and a first-entering or last-exiting event may occur at an
inner endpoint of a retained closed core. The theorem therefore uses `<=`,
not `<`.

## 5. Exact directed owner residues

For a positive modulus `M`, write `[x]^+_M` for the least strictly positive
residue of `x` modulo `M`. For one hypothetical body owner `u`, put
`d=gcd(t,u)`. The four exact inward distances from the successive boundaries
of `(3)` to a correctly oriented owner event are

```text
e_1(t,u)=[3t-4u]^+_(42d)/(42u),
e_2(t,u)=[9t+16u]^+_(126d)/(126u),
e_3(t,u)=[9t-110u]^+_(126d)/(126u),
e_4(t,u)=[3t+38u]^+_(42d)/(42u).                       (16)
```

For example, the first distance comes from

```text
3t(14a+1)-4u-42uM
 =3t-4u+42(ta-uM),                                    (17)
```

and `ta-uM` runs through `d Z`. Reversing orientation at the first right
endpoint gives base numerator `9t+16u`; the other two formulas are the same
calculation at `55/63,19/21`. Since `t` is odd, all four base numerators are
odd while their moduli are even, so the residues never vanish.

The two obstruction components are exact reflections at this level:

```text
e_4(t,u)=e_1(t,u),             e_3(t,u)=e_2(t,u).      (18)
```

Indeed their base numerators differ by `42u` and `126u`, respectively,
which are multiples of the corresponding gcd moduli. Thus both components
receive the same universal owner-relaxed erosion.

Minimize over every hypothetical owner up to the body maximum:

```text
epsilon_i(t,U)=min_(1<=u<=U) e_i(t,u),
B(t,U)=max(0,
           2/63-epsilon_1-epsilon_2,
           2/63-epsilon_3-epsilon_4).                 (19)
```

The actual body owns only eleven of these integers, so `(19)` can only
underestimate its collars and overestimate its largest retained component.
In the exact finite range `U<Q`, every failure therefore satisfies

```text
t(2U-1)/(84U(U-1))<=B(t,U).                           (20)
```

Replacing the range in `(19)` by the actual owner set gives a stronger
body-specific certificate. Formula `(20)` retains orientation, endpoint,
owner maximum and multiplier; it still forgets which four event owners can
coexist in one decoder body.

## 6. Finite exact census and THM-4002

The primary companion evaluates all `62,989` old-strip cells

```text
11<=t<=1001 odd,
max(11,floor(3t/4)+1)<=U<=t.                           (21)
```

The symbolic gate `(13)` closes `742` cells. The exact residue gate `(20)`
closes `77` more, leaving `62,170` owner-relaxed cells. The additional list
starts

```text
(11,11), (13,12), (13,13), (15,13), (15,14), (15,15)
```

and ends at `(185,141)`; the frozen output records all 77. These are
universal `(t,U)` parameter cells, not an enumeration of body shapes.

The combination is genuinely stronger than taking a union of the old
mechanisms. Relative to THM-3995's original 494 simple closures, the
two-flank bound adds 248 cells and the old residue-only gate adds 92; those
sets overlap in 34 cells. Their separate union adds 306, while applying the
longer body component and the exact residual core in the same inequality adds
325. The remaining 19 cells are exact **synergy rows**: neither refinement
closes them alone.

THM-4002 uses the same obstruction `C` but retains signed endpoint
cross-phase and exhausts fixed bodies `E subset {1,...,21}`. The present
theorem instead combines componentwise parity collars with a generic body
boundary lemma, so it applies to arbitrary bodies. It adds no scope inside
THM-4002's already closed fixed family, and neither theorem closes the
arbitrary `(2,1,9)` certificate type.

Oddness is load-bearing. If even `t` were admitted, `(t,u)=(12,9)` makes the
first and fourth directed residues vanish. Masked or simultaneous events do
not weaken the proof: a first positive cell still needs an entering event and
a last positive cell still needs an exiting event; cancellation can only
enlarge a collar.

## 7. Boundary and replay

The component-erosion/residue strip says nothing further about scale one;
the top-balance selector prunes the body-shape region of all sixteen
scale-one survivors but closes no entire arbitrary-body type. The quadratic
strip is automatically satisfied in most of `t<U`. The theorem constructs no
counterexample and does not prove LRC(14). The next missing coordinate is the
joint assignment of four endpoint events to one eleven-owner decoder body,
with deletion gcd and crossing-code constraints retained.

Reproduce from the repository root:

```text
python3 -B 04-computation/lrc14_scale_two_component_erosion_boundary_strip_thm4003.py
python3 -B -O 04-computation/lrc14_scale_two_component_erosion_boundary_strip_thm4003.py
python3 -B 04-computation/lrc14_scale_two_component_erosion_boundary_strip_independent_audit_thm4003.py
python3 -B -O 04-computation/lrc14_scale_two_component_erosion_boundary_strip_independent_audit_thm4003.py
```

The normal and optimized primary outputs and the independent audit match the
frozen raw-LF streams. **QED.**
