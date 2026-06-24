---
id: HYP-2938
status: PROOF-CARRIER / unit-preservation guardrail, not a proof
source: codex-2026-06-23
tags: [lrc14, pi, unital, flower, unit-distance, codes, tournament-analysis, unit-guardrail]
related:
  - HYP-2937
  - HYP-2936
  - HYP-2935
  - HYP-2932
  - HYP-2894
  - HYP-2892
  - HYP-2445
  - OPEN-Q-108
results:
  - 04-computation/pi_unital_flower_carrier_codex.py
  - 05-knowledge/results/pi_unital_flower_carrier_codex.out
---

# HYP-2938: pi/unital/flower carriers are useful only when the unit is preserved

The prompt linked four families:

1. `22/7` and `cuberoot(31)` as low-complexity approximants to `pi`;
2. a flower/petal process with turn `1/pi`;
3. geometric unitals `2-(q^3+1,q+1,1)`;
4. algebraic/functional-analytic unital structures and maps preserving
   identity.

The executable atlas finds the common proof lesson:

> The word "unital" is useful here as a guardrail: preserve the unit of
> measurement, identity, pair-incidence, or theorem scale before quotienting.

It is not a scalar theorem by itself.

## Numerical Anchors

The atlas records:

```text
pi              = 3.141592653589793
22/7            = 3.142857142857143   error +1.264489267350e-3
cuberoot(31)    = 3.141380652391393   error -2.120011984004e-4
355/113         = 3.141592920353983   error +2.667641894050e-7
pi^3 - 31       = 0.006276680300
```

Thus `cuberoot(31)` is a real low-complexity improvement over `22/7`,
because `pi^3` is just above `31`.

The flower statement needs a unit correction:

```text
1/pi turns      has convergent 7/22, drift after 22 petals about 1.014 degrees.
1/pi radians    has turn fraction 1/(2*pi^2), first good denominators 20 and 79.
22/7 inserted into the literal-radian normalization gives 49/968 turns.
```

So the claimed `22` petal families belong to `1/pi` as a **turn fraction**,
not to literal `1/pi` radians.  This is exactly the same kind of warning as a
non-unital homomorphism: changing the unit changes the structure.

## Unital Parameter Rows

For a formal unital parameter row,

```text
v = q^3+1 = (q+1)(q^2-q+1),
k = q+1,
r = q^2,
b = q^2(q^2-q+1).
```

Let `h=q^2-q+1`.  The notable rows are:

```text
q=3:  v=28=C(8,2),  k=4,  h=7
q=4:  v=65,         k=5,  h=13,  and 65+7=72
q=5:  v=126,        k=6,  h=21
q=6:  v=217=7*31,   k=7,  h=31
q=7:  v=344,        k=8,  h=43
q=10: v=1001,       k=11, h=91=C(14,2)
```

The `q=3` row is already active in HYP-2892/HYP-2894: it gives a real
`28`-point pair-frame but not a canonical AP8 pair-slot design.

The C27-facing next experiment is to try lifting the HYP-2937 marked shell
transfers into q=3 unital 4-point blocks after category-1 AP/Goddyn-Wong
labels are attached.  A positive lift would make the unital's `lambda=1`
balance a transfer-packet smoother before HYP-2908; a negative lift would
confirm that the unital remains only a weighted labelled frame.

The `q=6` row is the new prompt bridge.  As finite geometry, it should be
treated as a formal parameter row unless an actual design is supplied; but as
a carrier it is clean:

```text
block size = 7,
h = 31,
points = 7*31,
cuberoot(31) ~= pi.
```

That makes it a compact address for the LRC seven-sector world plus the
`pi^3 ~= 31` approximation.

The `q=4` row is a code hint rather than a construction: `65+7=72`, suggesting
a possible way to think about self-dual `[72,36,16]` code coordinates as a
`65`-point unital-like body plus seven sector/unit coordinates.  HYP-2894 is
the guardrail: the labelling must be observer-relative; the symmetric design
alone is not enough.

The `q=10` row has `h=91=C(14,2)`, echoing the LRC14 pair-slot shell.

## Relation to LRC14

This does not change the current LRC14 proof state.  It supplies labels for
the Binary Relational Exploration mandate:

- `unit_preserving_guardrail`: keep turn/radian, identity, measure-one,
  exact denominator, or pair-incidence units explicit;
- `C27_shell_transfer`: current S136 low-gap quotient;
- `q3_unital_pair_frame`: known weighted frame, after HYP-2894's guardrail;
- `q6_formal_unital_31`: new `7*31`/`pi^3` carrier;
- `q4_unital_65_plus_7_code_hint`: possible `[72,36,16]` coordinate prompt;
- `pi_turn_22_flower`: meaningful only after the turn unit is fixed.

The most useful transfer to LRC14 is methodological: any quotient that changes
the unit can create false structure.  The same issue appears in recent work as
`q` versus `p+q` versus `p*q`, apex phase versus optimal phase, and unital
pair-frame versus observer-relative exact tiling.

## Tournament Analysis

Vertices were proof carriers rather than runners:

```text
unit_preserving_guardrail,
C27_shell_transfer,
q3_unital_pair_frame,
q6_formal_unital_31,
pi_turn_22_flower,
literal_radian_correction,
q4_unital_65_plus_7_code_hint,
raw_pi_numerology.
```

The observable was a lexicographic tuple:

```text
(unit fidelity, LRC relevance, incidence strength, novelty, anti-scalar guard).
```

The tournament was transitive:

```text
unit_preserving_guardrail
> literal_radian_correction
> C27_shell_transfer
> q3_unital_pair_frame
> q6_formal_unital_31
> pi_turn_22_flower
> q4_unital_65_plus_7_code_hint
> raw_pi_numerology.
```

Fingerprints:

- score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`;
- directed `3`-cycles: `0`;
- SCC sizes all `1`;
- Hamiltonian path count `1`.

## Next Experiments

1. For `[72,36,16]`: test whether any known length-72 coordinate splitting
   naturally looks like `65+7`, and if so whether the seven coordinates behave
   as sector/unit labels or merely as parity padding.
2. For LRC14: attach the unit-preserving guardrail to every new quotient:
   specify whether the quotient preserves exact `q`, phase unit, sector unit,
   or pair-incidence unit.
3. For unital design carriers: extend HYP-2894 from q=3 to formal q=4/q=6
   parameter rows only as labelled incidence frames, not as symmetric tilers.
4. For flowers/phyllotaxis: compare turn fraction `1/pi`, literal radian
   fraction `1/(2*pi^2)`, and golden-angle fraction by their continued-fraction
   denominator families and see which relation packets survive quotienting.
