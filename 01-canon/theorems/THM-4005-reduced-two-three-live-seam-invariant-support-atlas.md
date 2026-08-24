---
id: THM-4005
title: "The live reduced 2:3 seam excludes the 3x4 and 4x3 support cells"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the oriented
  reduced pole-depth 2:3 cell, after the specific linear target normalization
  of THM-3992 and on the continued THM-3997 seam gamma=-a^3/2, the first two
  source-normal diagonals descend to the residual-mu5 coordinates
  A5=a^5, b=[y](R/gamma), d=[py](R/gamma). They force both outputs to have at
  least four retained nonconstant source weights. Thus the 3x4 and 4x3 cells
  are empty in this live seam. The exclusion transfers across exactly the
  diagonal scalings, translations, and C-by-A linear shear recorded in
  THM-3992. In the fixed gauge the sharper strata are 5x6 for b!=0 or d!=0
  and 4x5 for b=d=0; the 4x5 floor itself does not transfer across the shear.
  The known companion has an exact Weierstrass packet modulo t^3, and the
  first missing residual row is
  t^3*(c40+c21*x+c02*x^2). No atlas-point existence, other reduced cell,
  arbitrary target-automorphism invariance, global factorization, ownership,
  or JC(2) conclusion follows.
source: root + jc2_live_34_invariant_atlas + independent hostile audit, 2026-08-24
audit: >
  PASS (root + independent no-import audit, 2026-08-24). The primary exact
  companion rewrites the canonical THM-3997 rows, reconstructs every retained
  weight, derives the finite Weierstrass packet, and checks endpoint and
  target-shear hostiles. The independent verifier starts from the raw
  canonical rows, groups actual x^j*t^n monomials, retains special coefficient
  cancellations, and independently derives the companion by Euclidean
  division. It proves that only the seven-piece exclusion transfers across
  the recorded gauge and exhibits why the stronger 4x5 floor does not. Both
  normal and optimized primary streams match; both independent streams match
  the frozen transcript across 111 exact gates.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3987-gwozdziewicz-every-line-height-two-three-weight-floor
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
related:
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-3998-reduced-two-three-three-by-at-most-three-source-weight-support-obstruction
  - THM-3999-companion-divisor-boundary-endpoint-and-class-ledger
script: 04-computation/jc2_live_23_support_atlas_thm4005.py
output: 05-knowledge/results/jc2_live_23_support_atlas_thm4005.out
script_sha256: 6e83e7aa534c88cd6c50860cdda0dec862485b325e9c62d3fd03abad2afd63d4
output_sha256: e4378e7c3e2b2d55b561a0f09d9114464030f15bfed8ef46e1eb2b5b94c6e6c5
independent_audit_script: 04-computation/jc2_live_23_support_atlas_thm4005_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_live_23_support_atlas_thm4005_independent_audit.out
independent_audit_script_sha256: 94783fe667e66f8865c6296de133ef4b6f4a5ff597a5ec351a75d65cba4f48c9
independent_audit_output_sha256: 15c2bc6b34a1a2cc255313799a6f97fe429df79ca9483d7c4d5b92dea4681c6a
independent_audit_semantic_sha256: 432395a4b9ae2e17325e93294f04fcb156fd756d1e86bd63c6fe533860be7fb7
hash_basis: raw LF bytes
---

# THM-4005 -- invariant support atlas on the live reduced 2:3 seam

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Let `(A,C) in B_2^2`
be a hypothetical polynomial Darboux pair with

```text
J_(x,t)(A,C)=1.
```

Assume the oriented reduced pole-depth `(2,3)` cell, the specific target gauge
fixed in THM-3992, and the continued live seam of THM-3997:

```text
a!=0,                 gamma=-a^3/2,
A5=a^5!=0,            b=[y](R/gamma),
d=[py](R/gamma).                                      (1)
```

The theorem gives necessary identities for every such hypothetical pair. It
does not assert that any triple `(A5,b,d)` lifts through later rows or comes
from a Darboux pair.

## 1. Complete invariant source-normal jet through `t^2`

Divide `A` by the residual character `a` and `C` by `a^4`. The first two
THM-3997 source-normal diagonals become

```text
A/a = 1+(A5/4)x^2
    + t*(4/(3A5)+(A5 b/4)x+2x^2)
    + t^2*((9A5^3 b^2-512)/(144A5^2)
           +(A5 d+5b)x/4-4x^2/(5A5)+(3A5 b/8)x^3)
    + O(t^3),                                             (2)

C/a^4 = -3x/4-(A5/8)x^3
      + t*(-3b/8-4x/A5-(3A5 b/16)x^2-(3/2)x^3)
      + t^2*(-(3A5 d+7b)/(8A5)
             +(2816-45A5^3 b^2)x/(480A5^2)
             -3(A5 d+12b)x^2/16-12x^3/(5A5)
             -(9A5 b/32)x^4)
      + O(t^3).                                           (3)
```

These formulas are invariant under THM-3997's residual fifth-root action and
depend only on `(A5,b,d)`. They are exact rewritings of the canonical rows,
not a bounded-jet existence ansatz.

For a source-normal monomial `x^j t^n`, put

```text
weight(x^j t^n)=j-2n.                                    (4)
```

Target constants are deleted from retained support, while a nonconstant
weight-zero term such as `x^2t` remains. The pair `(weight,n)` determines
`(j,n)`, so a nonzero coefficient at one `x^2t` order cannot be canceled by
a later diagonal carrying the same weight.

## 2. Support exclusion and exact strata

Equations `(2)--(3)` already force

```text
{2,0,-2} subset supp(A),
{3,1,-1} subset supp(C).                                 (5)
```

THM-3987 independently forces each output to contain a surviving even
negative tail weight at most `-4`. That weight is new in both sets in `(5)`.
Consequently

```text
|supp(A)|>=4,                 |supp(C)|>=4,              (6)
```

and neither output can have retained support three. This excludes both
`3x4` and `4x3` in the fixed live-seam gauge.

The second diagonal gives the sharper complete three-stratum invoice.

1. If `b!=0`, then `A` also has weight `-1`, while `C` also has weights
   `-2` and `0`. Together with the compulsory tail,

   ```text
   |supp(A)|>=5,              |supp(C)|>=6.              (7)
   ```

2. If `b=0,d!=0`, the `t^2` row forces `A)-weights `-4,-3` and
   `C)-weights `-4,-3,-2`. Hence the same `5x6` bound holds. The proof
   retains the special algebraic loci on which other displayed coefficients
   cancel and does not use those unsafe terms.

3. If `b=d=0`, the exact forced sets are

   ```text
   {2,0,-2,-4} subset supp(A),
   {3,1,-1,-3} subset supp(C).                           (8)
   ```

   The required even `C)-tail is distinct from all four odd weights in
   `(8)`, so

   ```text
   |supp(A)|>=4,              |supp(C)|>=5.              (9)
   ```

Thus the first cell not rejected by these exact invoices in the fixed gauge
is `4x5`, and it occurs only on the stratum `b=d=0`. This is a necessary
candidate floor, not a construction.

## 3. What transfers across the THM-3992 target gauge

The target operations actually used by THM-3992 are nonzero diagonal
scalings, translations, and one triangular linear shear of the form

```text
C_old=c*C_new+kappa*A_new+constant,        c!=0.         (10)
```

They preserve the nonconstant `A)-support. If `A` had exactly four retained
weights, the strata above force `b=d=0` and

```text
supp(A)={2,0,-2,-4}.                                     (11)
```

The four forced normalized `C)-weights `{3,1,-1,-3}` are disjoint from
`(11)`; therefore undoing `(10)` cannot cancel them. The `3x4/4x3`
exclusion consequently transfers across the recorded THM-3992 normalization.

The stronger fixed-gauge `4x5` floor does **not** transfer. The invariant
atlas permits additional normalized `C)-weights forming a full copy of the
four `A)-components; a shear can cancel that copy and leave the four forced
odd `C)-weights. Thus this proof does not exclude a pre-gauge `4x4` packet.
No invariance under arbitrary nonlinear target automorphisms is claimed.

## 4. The known companion packet and its first missing row

Put `Rtilde=R/gamma`. On the live seam, THM-3997 gives

```text
Rtilde = -16p^2/(3A5^2)+b*y+d*p*y
       +(b^2/4+2752/(135A5^3))*p^3+O_source(t^4).        (12)
```

With `z=1+x^2t`, `p=tz`, `y=xtp`, and
`Qtilde=Q/gamma=G/(gamma*t)`, every coefficient through `t^2` is fixed:

```text
Qtilde = x^2+6/A5
       +t*(b*x+6x^2/A5-16/(3A5^2))
       +t^2*(b^2/4+b*x^3+d*x-32x^2/(3A5^2)
              +2752/(135A5^3))+O(t^3).                 (13)
```

Modulo `t^3`, exact Weierstrass division gives `Qtilde=U*W`, where

```text
U = 1+6t/A5+t^2*(b*x-32/(3A5^2)),                       (14)

W = x^2+6/A5+t*(b*x-124/(3A5^2))
  + t^2*((d-12b/A5)*x+b^2/4+44872/(135A5^3)).           (15)
```

The discriminant of `W` is

```text
-24/A5+(496/(3A5^2))*t-(179488/(135A5^3))*t^2+O(t^3),  (16)
```

independent of `b,d` through this order and nonzero at `t=0`. If
`e_+^2=e_-^2=-6/A5` with `e_-=-e_+`, the two finite normalization germs
satisfy

```text
x_+(t)+x_-(t)=-b*t+(-d+12b/A5)*t^2+O(t^3).              (17)
```

These are two normalization addresses, not a component or owner census.

THM-3999's normalized endpoint polynomial is

```text
E(y)=1-Rtilde(0,y)=1-b*y+O(y^2).                        (18)
```

Thus `b!=0` forces a boundary endpoint, but `b=0` does not imply boundary
disjointness. The boundary-disjoint hostile `Rtilde in (p^2,py)` retains the
mandatory coefficient `[p^2]Rtilde=-16/(3A5^2)`.

The first residual layer not determined by the available rows is exactly

```text
([p^4]Rtilde,[p^2y]Rtilde,[y^2]Rtilde)=(c40,c21,c02).    (19)
```

After division by `t`, it contributes

```text
t^3*(c40+c21*x+c02*x^2)+O(t^4)                         (20)
```

to `Qtilde`. Here `c02` is the first missing endpoint coefficient after
`b`, while `c21` is the first missing source-odd sidecar on `b=d=0`.
Computing that next Hasse/source-normal row is the cheapest exact continuation.

## 5. Scope and reproduction

This theorem proves only necessary statements in the oriented reduced
`(2,3)` live seam and the recorded linear target gauge. It does not prove:

1. existence or all-row consistency of any `(A5,b,d)` atlas point;
2. the reversed depth orientation or any other reduced cell;
3. a target-automorphism-invariant `4x5` floor;
4. global factorization or irreducibility of the companion;
5. completeness or ownership of the two visible normalization addresses;
6. emptiness of the entire reduced `(2,3)` problem or `JC(2)`.

Reproduce from the repository root:

```text
python -B 04-computation/jc2_live_23_support_atlas_thm4005.py
python -B -O 04-computation/jc2_live_23_support_atlas_thm4005.py
python -B 04-computation/jc2_live_23_support_atlas_thm4005_independent_audit.py
python -B -O 04-computation/jc2_live_23_support_atlas_thm4005_independent_audit.py
```

The two primary streams agree exactly. The two independent streams agree
exactly with the frozen audit transcript, pass `111` exact gates, and have
semantic SHA-256
`432395a4b9ae2e17325e93294f04fcb156fd756d1e86bd63c6fe533860be7fb7`.
