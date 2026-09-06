# THM-4453 inert-sum/five-ray disjointness and entry composition

**Status: PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY LOGIC-AUDITED;
LRC(14) OPEN.** The strict-topology pre-audit is `VERIFIED-EXACT`; the new
inert-atlas separation follows from THM-4153/3818/4450 plus exact integer
arithmetic.

## 1. Strict versus essential topology

For the physical clock-two failure set

```text
F_T=(union_(n in T) D_n) intersect ((union_(n in T) D_n)+1/2),
T=(1,7,13), D_n={x: ||nx||<1/14},
```

an independent rational-wall reconstruction gives

```text
actual strict topology:    longest 1/91, 12 components,
                           widths 4*(1/91)+8*(1/98);
a.e.-merged representative: longest 1/49, 8 essential components,
                           widths 4*(1/49)+4*(1/91).
```

The joining point `x=1/14` is absent from the strict set. Thus `1/49` is not
an actual connected-component width at this shape.

The only unrepaired positive assertion found is scratch history:

```text
.scratch/lrc14_composite_residual_next_sep06/
q2_tail_network_probe_height199.out:23
```

Its source merges touching intervals at
`q2_tail_network_probe.py:33-36` before reporting `longest` and `components`
at lines `92-95`. The downstream scratch scout
`.scratch/lrc14_q2_general_odd_sep06/general_odd_probe.py:14-18,25,44-53`
imports that routine. The parent report's displayed winner `(1,7,11)` at
lines `541-556` is strict-correct, so its finite mass conclusion is harmless;
only nonleader topology rows are contaminated.

Maintained canon, results, scripts, and `05-knowledge/results/INDEX.md` now use
the `(1,7,13)` occurrence only as correction prose. In particular THM-4451
states the distinction explicitly at theorem lines `38-41,170-181`. THM-4449
now calls its BV count `essential circle components` and explicitly disclaims
actual strict topology at lines `200-210`; this terminology repair has no
effect on its mass or BV cutoff.

## 2. New inert-atlas/five-ray separation

Write `G_H={y: ||hy||>=1/14 for every h in H}`. Let `a,b` be distinct positive
integers coprime to six, put

```text
s=gcd(a,b), a=sp, b=sq, 1<=p<q, gcd(p,q)=1.
```

### Proposition (theorem-ready)

Assume every prime divisor `ell` of `p+q` satisfies

```text
ell == 2 (mod 3), with v_ell(p+q)<=2.                    (A)
```

For every finite positive body `H`, if

```text
mu(G_H) >= 4/91,                                       (B)
```

then the row `2H disjoint-union {a,b}` has a common phase of clearance at
least `1/14`.

This proposition does not require a decoder graph, a finite height box, or a
literal unit. Its LRC14 applications below use `|H|=11` so the row has thirteen
speeds.

### Proof

Under `y=2x`, failure says that the compact set `G_H` lies inside the proper
open pulled two-tail cross-comb at primitive ratio `(p,q)` and odd common scale
`s`. Compact-to-proper-open containment makes the measure inequality strict.
Thus (B) and failure imply primitive cross-comb mass `>4/91`, not merely
`>=4/91`.

THM-4153's exact all-height classification (theorem lines `70-90,137-149`)
says that a primitive odd ratio of mass `>4/91` belongs to its eleven-ratio
set. Because both `p,q` are coprime to six, precisely five remain:

```text
(1,11), (1,23), (5,11), (1,37), (1,25).                (C)
```

Their primitive sums and obstructions to (A) are

```text
(1,11): 12=2^2*3       (prime 3),
(1,23): 24=2^3*3       (prime 3 and v_2=3),
(5,11): 16=2^4         (v_2=4),
(1,37): 38=2*19        (19 == 1 mod 3),
(1,25): 26=2*13        (13 == 1 mod 3).
```

Hence (A) and (C) are disjoint, a contradiction. The `4/91` boundary is
inclusive exactly because failure forced strict `>4/91`.

The six-coprime hypothesis is essential to this five-ray version. In the
larger all-odd exception set, `(1,9)` survives; `p+q=10=2*5` satisfies (A) and
the ratio belongs to the THM-3818 atlas. This is the minimal scope hostile.

## 3. Actual THM-3818 entry consequence

THM-3818 requires its selected primitive pair to have `p+q<=356` and satisfy
(A) (theorem lines `129-145`). Re-enumeration gives exactly `5,855` ratios,
with `p<=177,q<=355`, and empty intersection with (C).

Therefore every actual THM-3818 `11+2` decoder equality row in the odd-3-unit
residual branch is safe whenever its eleven-coordinate component can be typed
as the physical `2H` in the proposition and `mu(G_H)>=4/91`. The two labels in
`{a,b}` must be **exactly the two-vertex decoder component**; it is not enough
that some unrelated pair in the row is atlas-admissible. Under the actual
entry hypotheses that two-vertex component is connected, so its unique pair
is a decoder edge and its primitive ratio is indeed in the 5,855-entry atlas.

At this mass level the conclusion is stronger than incoming commit
`9db26e9f3`: no normalized unit, `W=V_dec`, or finite-box assumption is needed
once the physical pair itself satisfies (A). The incoming result remains
valuable below the mass gate.

## 4. Original-body q=4 and q=2 gates

These are genuine sufficient gates **inside the inert-pair class** (and hence
inside the matching actual THM-3818 branch), not merely preprocessing
conditions. For arbitrary odd-3-unit pairs they only trigger localization and
leave the five rays (C).

### q=4, exactly one tail of valuation one

For a ten-body `C` and odd 3-units `r,a,b`, write

```text
4C union {2r,a,b} = 2H union {a,b},  H=2C union {r}.
```

Assume the thirteen speeds are distinct and the reduced pair `(p,q)` obeys
(A). THM-4450's exact hybrid identity gives

```text
mu(G_H) >= max(mu(G_C)/2, mu(G_C)-1/8).
```

Consequently

```text
mu(G_C) >= 8/91  ==>  mu(G_H) >= 4/91  ==>  row safe.   (D)
```

Equality in (D) is allowed. The half-mass term is the cheapest route; the
other term alone would require `mu(G_C)>=4/91+1/8=123/728`. Thus `8/91` is a
conditional sufficient original-body gate, improving the unconditional
odd-3-unit pair-cap gate `8/77` only because inertness removes all five rays.

### q=2, exactly one even tail

Write

```text
2C union {2r,a,b} = 2H union {a,b},  H=C union {r},
```

where `r,a,b` are odd 3-units, `|H|=11`, and the reduced pair obeys (A).
THM-4450's sharp absorbed-label floor is

```text
mu(G_H) >= mu(G_C)-8/63.
```

Therefore

```text
mu(G_C) >= 4/91+8/63 = 20/117  ==>  row safe.          (E)
```

Again equality is allowed, and `20/117` is a conditional sufficient
original-body gate. It improves `124/693`, which used the full `4/77` pair cap
and did not assume inertness.

There is no corresponding direct q=2 zero-even composition: that signature is
the natural `10+3` split `2C union {a,b,c}`. An odd tail cannot be absorbed as
an integral member of an even core `2H`. THM-4451's strict three-tail component
gate remains the correctly typed mechanism there.

## 5. What incoming `9db26e9f3` does below the gates

For a primitive row in THM-3818's actual finite box, suppose the actual decoder
components are exactly the natural `2H` and `{a,b}`, have sizes `11+2`, and
the full support-at-most-three, height-at-most-`Q` span obeys `W=V_dec`. Put

```text
d=gcd(H), V=H/d, t=2d, g=gcd(a,b), (p,q)=(a/g,b/g).
```

Primitivity gives `gcd(t,g)=1`. If `(p,q)` is the actual atlas pair and
`1 in V`, the incoming theorem closes the row at every body mass.

For q=4, `H=2C union {r}` and `r` is odd. Then `d` is odd, so no even entry
`2c` can equal `d`; hence

```text
1 in H/d  iff  d=r  iff  r divides every c in C.        (F)
```

For q=2 one-even, `H=C union {r}` and the literal unit condition is merely
that some member of `H` equals `gcd(H)`; it need not be `r`.

This is a fixed-row conditional, not a universal body-component lower bound.
Its LRC(11) step gives a normalized component of `G_V` length at least
`1/(42 max V)`; after the `d` pullback this is `1/(42 max H)`, exactly the
height-dependent bound already used by THM-4450. The breakthrough is instead
the bounded-crossing relation forcing a large scale separation and the cyclic
grid gluing. Mass alone, gcd(V)=1 alone, or component cardinalities alone do
not preserve those coordinates.

THM-4451 also cannot be applied directly to q=4's tail triple `(2r,a,b)`.
Pointwise, because `D_(2r)+1/2=D_(2r)`,

```text
F_(2r,a,b) = D_(2r) union
             ((D_a union D_b) intersect ((D_a union D_b)+1/2)).
```

At `r=1`, a `D_2` tooth already has physical width `1/14`, exceeding both
THM-4451 strict odd-tail caps `17/693` and `19/1001`. Oddness is not cosmetic;
the correct q=4 target is THM-4450's two-tail quotient with its sheet/address
sidecar.

## 6. Typed map and cheapest decisive tests

Source for q=4:

```text
(C;r,a,b) |-> H=2C union {r}, d=gcd(H), V=H/d,
              t=2d, g=gcd(a,b), (p,q)=(a/g,b/g).
```

Target: an actual THM-3818 datum `(I,J;V,(p,q);t,g;W,V_dec)` with `I=2H`,
`J={a,b}`.

Preserved: the physical thirteen-speed row, labelled `11+2` partition,
primitive pair ratio, strict failure implication, and scale coprimality.

Lost unless retained as sidecars: the actual decoder-graph identity of the
partition; `W=V_dec`; the physical scales `t,g`; the literal coordinate one;
and, under `y=2x`, the half-lift/sheet address. `mu(G_H)` discards all
component addresses and endpoints.

Cheapest tests, in order:

1. In the high-mass lane, compute `mu(G_H)` exactly and factor `p+q`. If
   `mu(G_H)>=4/91` and (A) holds, the proposition closes immediately; no
   decoder computation is needed.
2. Original-body sufficient prefilters are the exact inequalities (D) and
   (E), followed only by factorization of the physical tail pair.
3. Below those gates, test the literal unit first. In q=4 this is the cheap
   divisibility test (F). Only then is it worthwhile to build the complete
   support-three relation graph and verify the actual `11+2` split and
   `W=V_dec` required by `9db26e9f3`.

## 7. Reproduction

From the isolated worktree root:

```text
python -B 04-computation/lrc14_inert_sum_five_ray_disjointness_thm4453.py
python -O -B 04-computation/lrc14_inert_sum_five_ray_disjointness_thm4453.py
```

Both normal and optimized outputs match their frozen `.out` files and end in
`PASS`. The composition script independently enumerates all 5,855 atlas pairs,
checks the five exclusions, retains `(1,9)` as a hostile scope control, checks
all fraction transfers, verifies 45,045 q=4 unit/divisibility controls, and
tests the even-tail Boolean identity.

Raw-LF SHA-256 values (before this report):

```text
04-computation/lrc14_inert_sum_five_ray_disjointness_thm4453.py
  e456fb062a4e13ed29b7c68c5b926107d4524f6903a0ac1c70b95c6755bfa74d
05-knowledge/results/lrc14_inert_sum_five_ray_disjointness_thm4453.out
  da67b3db4d130994ff7eb643b2ee595ea23950d630a3d3e3b219502800314c8d
```

The independent quantifier audit is
`05-knowledge/results/lrc14_inert_sum_five_ray_disjointness_thm4453_independent_audit.md`.
