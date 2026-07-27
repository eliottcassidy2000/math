---
id: THM-2566
title: "Two-chart saturated cusp atlas and parasitic-plane ledger for the sporadic Keller map"
status: >
  PROVED + INDEPENDENTLY AUDITED over C (indeed all polynomial identities and
  ideal calculations are over Q).  The two raw cusp pullbacks, their divisor
  multiplicities and saturations, the punctured two-chart cover, both hostile
  controls, both origin images, and both cusp-point fibres are exact.  The
  generic fibre interpretation is inherited separately from THM-2473/2546.
  This concerns only the fixed sporadic three-variable Keller map and proves
  nothing about JC(2) or a general Keller-monoid tower.
source: codex-2026-07-27
audit: root-2026-07-27 (independent algebra and normal/-O/stored/hash/docs replay)
depends_on:
  - THM-1335-trisection-modulus-master-identities-trace-polynomiality
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
related:
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
  - MISTAKE-287 (raw cuspidal pullback lost a coefficient-zero plane)
script: 04-computation/keller_two_cusp_atlas_thm2566.py
output: 05-knowledge/results/keller_two_cusp_atlas_thm2566.out
script_sha256: 6758420fb5e599f3f86e3adca73c8b601c01da9973c4628c4f131c424748f6e9
output_sha256: 39b229cc575eea92f2df30efed4051b9a61e83d1ad219a9ed26abf940a6af5c5
---

# THM-2566 -- the Jelonek hypersurface is a saturated cusp pullback

Work in

```text
R = Q[a,b,c],
L = 27a^2c^2 - 18abc + 16a + b^3c - b^2.
```

For the fixed sporadic Keller map `F`, THM-2473 proves that the Jelonek
non-properness set is the irreducible hypersurface `A_F=V(L)`.  This theorem
audits two polynomial maps from target space to the same cuspidal cubic

```text
C = V(S^2-T^3) in A^2.
```

The central distinction is between a raw inverse image and its saturation.

## 1. Two exact cusp identities

Define the escape-coordinate chart

```text
T_c = 4 - 3bc,
S_c = 27ac^2 - 9bc + 8,
Phi_c(a,b,c) = (S_c,T_c),
```

and the companion `u=1+xy` chart inherited from THM-1335,

```text
T_a = b^2 - 12a,
S_a = 54a^2c - 18ab + b^3,
Phi_a(a,b,c) = (S_a,T_a).
```

Direct expansion gives

```text
S_c^2 - T_c^3 =  27 c^2 L,                         (1)
S_a^2 - T_a^3 = 108 a^2 L.                         (2)
```

Since `c` and `a` are prime and neither divides the irreducible polynomial
`L`, the scheme-theoretic inverse images and their effective divisors are

```text
Phi_c^{-1}(C) = V(c^2L),       div = 2V(c)+V(L),    (3)
Phi_a^{-1}(C) = V(a^2L),       div = 2V(a)+V(L).    (4)
```

Consequently their underlying sets are

```text
|Phi_c^{-1}(C)| = V(c) union V(L),
|Phi_a^{-1}(C)| = V(a) union V(L).                  (5)
```

This is the exact repair of HYP-9033's former raw-pullback sentence; see
MISTAKE-287.  The coefficient squares cannot be discarded globally.

The same correction applies to the associated Weierstrass presentation.  If

```text
p = 3(12a-b^2),             q = 2S_a,
```

then

```text
-4p^3 - 27q^2 = -2^4 3^6 a^2 L.                   (6)
```

Thus the raw singular-pencil divisor is `2V(a)+V(L)`, not just the Jelonek
divisor.  It becomes exactly `V(L)` on `D(a)` or after saturation by `a`.

## 2. Saturation recovers exactly the Jelonek hypersurface

In the UFD `R`, unique factorization in (1)-(2) gives the principal-ideal
equalities

```text
((S_c^2-T_c^3) : c^infinity) = (L),                (7)
((S_a^2-T_a^3) : a^infinity) = (L).                (8)
```

For example, if `c^n f` is divisible by `c^2L`, then primality of `c` and
`c not| L` force `L|f`; the converse follows after multiplying `Lf` by a
sufficient power of `c`.  The same proof works for `a`.  The companion
verifies (7)-(8) independently by the standard elimination ideals

```text
(S_c^2-T_c^3, 1-uc) intersect R,
(S_a^2-T_a^3, 1-ua) intersect R.                   (9)
```

Equivalently, each localized chart is a literal cusp pullback:

```text
D(c) intersect Phi_c^{-1}(C) = D(c) intersect V(L),
D(a) intersect Phi_a^{-1}(C) = D(a) intersect V(L). (10)
```

These charts cover the punctured Jelonek hypersurface.  Indeed,

```text
L|_(a=c=0) = -b^2,                                 (11)
```

so a point of `V(L)` with `a=c=0` must be the origin.  Therefore

```text
V(L) = {0}
       union (D(c) intersect Phi_c^{-1}(C))
       union (D(a) intersect Phi_a^{-1}(C)).        (12)
```

At the exceptional origin the two charts remember different cusp strata:

```text
Phi_c(0,0,0) = (8,4),       a smooth point of C,
Phi_a(0,0,0) = (0,0),       the cusp point of C.    (13)
```

Thus (12) is genuinely a two-chart atlas with one explicitly retained point,
not an equality obtained by intersecting the two raw pullbacks.

## 3. The parasitic planes are genuine affine target loci

On the first coefficient plane,

```text
Phi_c(a,b,0) = (8,4)                                 (14)
```

for every `(a,b)`.  The cusp is smooth there because the gradient of
`S^2-T^3` at `(8,4)` is `(16,-48)`.  This plane is not a boundary or a
compactification chart: it is an ordinary target plane and has the global
polynomial section

```text
F(0,b,a-4b^2) = (a,b,0).                           (15)
```

THM-2473 separately proves the generic `1+2` fibre law on this plane: away
from `L|_(c=0)=16a-b^2`, the section is accompanied by two finite residual
points, so the complex fibre has its full three points.

The clean hostile

```text
(a,b,c)=(1,0,0):       Phi_c=(8,4),       L=16      (16)
```

lies on the raw cusp pullback but off the Jelonek set and has a full fibre.
This is the first failed point of the global claim.

The second coefficient plane behaves differently but is equally real:

```text
Phi_a(0,b,c) = (b^3,b^2),                           (17)
```

so `V(a)` maps through the normalization parametrization of the cusp.  The
hostile `(0,1,2)` has `Phi_a=(1,1)` and `L=1`.  More sharply,

```text
(a,b,c)=(0,1,0):
  L=-1,  Phi_c=(8,4),  Phi_a=(1,1),                 (18)
```

so even the intersection of the two **raw** inverse images contains the
extra line `V(a,c)`.  Localization or saturation is load-bearing.

## 4. The cusp-point fibres and the omitted curve

For the `c`-chart,

```text
T_c=0  iff  bc=4/3,
S_c=(27ac^2-4)+3T_c.
```

Hence

```text
Phi_c^{-1}(0,0)
 = {(4/(27t^2), 4/(3t), t) : t in C^*}.            (19)
```

In particular `c` is automatically nonzero, so no parasitic `c=0` point
enters this fibre.  THM-2546 proves that (19) is exactly the empty-fibre curve
`E` of `F`.  Thus the statement “the cusp point detects total escape” survives
globally even though the raw cusp-curve pullback did not.

The `a`-chart retains its own coefficient-zero component at the cusp point.
Modulo `T_a=0`, exact division gives

```text
8S_a = b^3(3bc-4).                                  (20)
```

Therefore, set-theoretically,

```text
Phi_a^{-1}(0,0) = V(a,b) union E,                   (21)
```

the target `c`-axis together with the same omitted curve.  On `D(a)`, only
`E` remains.  Equations (19)-(21) explain the different origin images in
(13) rather than hiding them in one scalar discriminant.

On `D(a)`, these identities and THM-2473/2546 give the exact local Kodaira
dictionary:

| target stratum | Weierstrass cubic | affine `F`-fibre count |
|---|---|---:|
| `L!=0` | smooth | `3` |
| `L=0` off `E` | nodal (`Delta=0`, `(p,q)!=(0,0)`) | `1` |
| `E` | cuspidal (`p=q=0`) | `0` |

The same table is false globally: the parasitic axis `V(a,b)` is cuspidal,
but THM-2546 gives the target `(0,0,c)` its one affine preimage `(c/2,0,0)`.

## 5. Mechanism, boundary, and loss ledger

- **Mechanism:** a cubic discriminant records both degeneration of the actual
  escape divisor and degeneration caused by a selected coefficient becoming
  zero.  Here those pieces separate exactly as `a^2L` or `c^2L`; saturation
  removes the coefficient plane without discarding the Jelonek closure.
- **Equality boundary:** on the stated charts, raw pullback equals Jelonek
  pullback after restricting to `D(c)` in the escape chart or `D(a)` in the
  companion chart.  Globally the extra planes, their multiplicity two, and the
  exceptional origin are unavoidable.
- **Preserved predicate:** the saturated maps preserve membership in the
  fixed map's non-properness set.  The raw maps preserve only membership in
  “Jelonek or coefficient-zero.”
- **Destroyed information:** the cusp coordinate alone forgets which factor
  caused the discriminant to vanish.  The coefficient `a` or `c` is the
  required sidecar.
- **Scope:** no statement here concerns a general Keller map, composition of
  Jelonek divisors, atomhood, `JC(2)`, or the conjectural genus axis.  Any
  tower law must carry and saturate every leading/constant-coefficient factor
  before reading a discriminant zero set as an asymptotic set.

## Reproduction

```bash
python3 04-computation/keller_two_cusp_atlas_thm2566.py
python3 -O 04-computation/keller_two_cusp_atlas_thm2566.py
diff -u 05-knowledge/results/keller_two_cusp_atlas_thm2566.out \
  <(python3 04-computation/keller_two_cusp_atlas_thm2566.py)
sha256sum 04-computation/keller_two_cusp_atlas_thm2566.py \
  05-knowledge/results/keller_two_cusp_atlas_thm2566.out
```

The ordinary and optimized transcripts agree exactly.  The script uses no
`assert`, performs both saturation eliminations over `QQ`, verifies both
hostiles, the affine section, the open-cover boundary, the origin images, and
both cusp-point fibres, and matches the frozen output byte-for-byte.
