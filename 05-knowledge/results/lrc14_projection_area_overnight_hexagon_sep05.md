# Exact continuous network roofs and the arithmetic boundary they omit

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED; FINITE-EXACT named
controls.** The `observer_collision` referee independently checked the
projected owner-density factor, layer-cake strip integral, simplification,
strict ordering, and sharp continuous bounds. Independent rational polygon
integration verifies the formula without using its layer-cake derivation.
This is a new sidecar and a method obstruction, **not** a new universal
discrete network bound or a claim to LRC(14).

The continuous benchmarks of all three LRC14 degree-zero networks have
surprisingly simple formulas. For primitive sorted positive odd ternary-unit
`w=(a,b,c)`, let `E_i` be the exact raw projection sum in
[THM-4414 — six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md).
Replace the owner-live lattice measure by its continuous density and call
the resulting integral `J_i`. Then

```text
J_i = [6c(a+b+c)-6w_jw_k-2ab]/(343c^2),
0<J1<J2<J3<12/343,             J1<10/343.               (1)
```

Both constants are sharp suprema for their specified continuous quantities.
They do not dominate the actual network sums: even `E_i<=2J_i` fails at a
complete three-direction carrier dictionary. The precise missing object is
the owner-coset boundary discrepancy of the projection roof.

## Inheritance, concept board, and what is not being transferred

The closest proved mechanism is the raw roof formula THM-4414. The earlier
[THM-4392 — raw-carrier box-spline Fourier--Poisson duality](../../01-canon/theorems/THM-4392-lrc14-raw-carrier-box-spline-fourier-poisson-duality.md)
provides the lattice/continuous normalization for the **physical** three-
facet minimum. Here the summand is one network projection instead; it can
remain positive on a transverse boundary of the carrier body, so its
physical box-spline formula must not be substituted.

The canonical hostile is `(19,23,29)`; the corrected near miss is uniform
area density implying a discrete bound. The least-used sidecar is the
transverse lattice width and owner offset. The concept board is: live
lattice basis; strict convex carrier body; projected roof; continuous area;
owner discrepancy; selector ordering. The incoming
[midpoint-deficit report](overnight_20260906_lrc_midpoint_deficit.md)
already proves that curvature can miss an affine boundary deficit. The
parallel conclusion here is more limited: continuous density can miss a
persistent boundary discrepancy, and the actual selector can differ from
the continuous one. These are two loss audits, not an asserted equivalence
of curvature and lattice-point discrepancy.

Targeted searches of the recent LRC canon and overnight result notes for
the constants, continuous projection formula, owner density, and named
ratio hostile found no earlier statement of (1). The existing physical
coarea identity and raw network formula are explicitly inherited. No
reserved atlas result is used, and no larger height census is performed.

## 1. Integrate the exact projection roof, not the physical roof

Put

```text
A_i=3(sum(w)-w_i)/14,      q=3/(7c),
B={C in w-perp : |C_i|<A_i for all i},
f_i(C)=min(q,(A_i-|C_i|)/(w_jw_k)) on B.                (2)
```

The ambient kernel lattice projected onto coordinates `C_j,C_k` has
covolume `w_i`, since `w` is primitive. The live set is two cosets of the
index-nine dilation `3 Lambda` inside the kernel lattice, so its continuous
projected density is `2/(9w_i)`. Define

```text
J_i=2/(9w_i) integral_B f_i(C) dC_j dC_k.               (3)
```

This fixes every normalization: there is no additional sign, owner, or
three-sheet factor. Open and closed polygon boundaries have the same
integral; all discrete carriers elsewhere retain the strict inequalities.

By layer cake, at height `0<=t<=q` only the `i`-th strip shrinks:

```text
|C_i|<A_i-t w_jw_k; the other two strips stay fixed.    (4)
```

Use `U=w_j C_j`, `V=w_k C_k`, and abbreviate

```text
P=w_j A_j, Q=w_k A_k,
R(t)=w_i A_i-abc t,      V0=abc.
```

The section is the rectangle `[-P,P] x [-Q,Q]` intersected with
`|U+V|<R(t)`. Its area is

```text
4PQ-(P+Q-R(t))_+^2+(|P-Q|-R(t))_+^2.                  (5)
```

For the actual LRC parameters the first hinge is always positive:

```text
P+Q-R(t)=3w_jw_k/7+V0 t>0.                            (6)
```

The second hinge is identically zero on `[0,q]`. Indeed `R(t)` decreases,
and `R(q)=|P-Q|` for `i=1,2`, while
`R(q)-|P-Q|=3a(c-b)/7>0` for `i=3`.

Let `X=3w_jw_k/7`. Integrating (5), with the Jacobian in (3), gives

```text
J_i=2/(9V0)[4PQ q-X^2 q-X V0 q^2-V0^2 q^3/3].         (7)
```

Now `4PQ-X^2=(9/49)V0(a+b+c)`. Substitution of `q=3/(7c)` proves the
polynomial formula (1). This argument is analytic; parity, residue, and
primitivity are needed for interpreting the chosen density as the LRC
owner density, not for the algebraic strip integration itself.

## 2. Continuous order, sharp bounds, and the exact selector discrepancy

The coordinate differences are

```text
J2-J1=6(b-a)/(343c)>0,
J3-J2=6a(c-b)/(343c^2)>0.                              (8)
```

Write `x=a/c`, `y=b/c`, so `0<x<y<1`. Then

```text
343J1=6+6x-2xy <6+6x-2x^2<10.
```

For `J3`, the numerator is `6+6x+6y-8xy`. If `x<=3/4`, maximize the affine
function of `y` at `y=1`, obtaining `12-2x<12`. If `x>3/4`, maximize it
at the lower endpoint `y=x`, obtaining `6+12x-8x^2<=21/2<12`.
Together with (8), this proves both bounds in (1).

The primitive unit triples `(6n+1,6n+5,6n+7)` show `J1 ->10/343`.
The triples `(1,6n+1,6n+5)` show `J2,J3 ->12/343`. Thus these are genuine
sharp continuous suprema, not equality values of the discrete networks.

Define the exact arithmetic discrepancy `Delta_i=E_i-J_i`. Then

```text
E2-E1=6(b-a)/(343c)+Delta2-Delta1,
E3-E2=6a(c-b)/(343c^2)+Delta3-Delta2.                   (9)
```

Every reversal of the fixed continuous ordering is therefore accounted for
by an explicitly sized difference of arithmetic discrepancies. At the
inherited selector hostile `(1,19,79)`,

```text
(E1,E2,E3)=(184/10507,8/553,12/553),  E2<E1,
```

even though `J1<J2`. Formula (9) specifies the missing sidecar rather than
recommending the wrong fixed projection from a continuous picture.

## 3. Hostile boundaries and the strongest surviving statement

At `(19,23,29)`, all three actual sums exceed their continuous benchmarks.
More sharply, the complete three-direction example `(5,23,47)` has

```text
Omega=+/-{(1,10,-5),(13,-11,4),(14,-1,-1)},
E3=2038/37835,    J3=2890/108241,
E3/J3=335251/166175>2.                                (10)
```

No minimal-height or global-extremum assertion is made. This is a genuine
multi-direction hostile, not a one-ray artifact. The exact three-direction
theorem already closes its network target; (10) refutes the proposed
intermediate domination, not that final target.

There is also a fully analytic abstract hostile to an additive universal
boundary constant. Let `L=Z^2`, `H={(x,y):x=0 mod3}`, take `0<eps<1/2` and
an integer `T>=3`, and set

```text
B=(-1-eps,1+eps) x (-T,T),
f(x,y)=min(1,T-|y|) on B.
```

The complete surviving set is `{(x,y):x=+/-1, y=1-T,...,T-1}`. It spans
the plane and contains a live lattice basis, with an unbounded number of
primitive directions as `T` grows. Its exact weighted sum and owner-density
integral are

```text
E=4T-2,
J=(4/3)(1+eps)(2T-1),
E/J=3/[2(1+eps)]>1.                                  (11)
```

For fixed `eps<1/2`, `E-J` grows linearly in `T`. The roof vanishes at the
two longitudinal ends but remains positive along most of the transverse
boundary, which falls immediately outside the live columns. Thus no bound
`E<=J+K` with an absolute constant `K` follows merely from a symmetric
convex body, an index-three deleted subgroup, a live basis, and this type
of capped projection roof. The example is **abstract**, not asserted to be
realizable by an LRC speed triple.

The first failed implication is replacing the complete owner measure by
its average density without a boundary term. The strongest survivor here
is the exact continuous formula and exact discrepancy identity (9). The
next test for any proposed bound is to control the transverse integral
width and owner offset at every layer (4), including layers whose live
support drops to one direction or becomes empty. An unweighted basis or
the physical three-facet minimum does not supply that missing control.

## 4. Exact verification and remaining scope

The
[independent polygon script](../../04-computation/lrc14_projection_area_overnight_hexagon_sep05.py)
clips the original `(C1,C2)` carrier polygon into a central plateau and two
affine roof pieces, then triangulates and integrates each affine function.
It does not implement the layer-cake formula as its integration engine.
All `220` sorted positive integer triples with `c<=12` test the continuous
algebra, including nonprimitive and nonunit values solely as algebraic
controls. Eight named LRC dictionaries are rebuilt using the inherited
congruence-row engine; no larger census or unproved atlas is imported.

```bash
python3 -B 04-computation/lrc14_projection_area_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_projection_area_overnight_hexagon_sep05.py
```

The [output](lrc14_projection_area_overnight_hexagon_sep05.out) records exact
rational integrals and all named positive/hostile controls. Gates remain
enabled under optimization. Normal and optimized outputs byte-match, with
`1,360` primary checks and `294` inherited row-engine checks. Raw-LF
SHA-256 values are:

```text
source 8d02f14cdd19a744a645024bf053151c1e3819cb5d0d7f22c6cd322b51664aeb
output f13051b8ca103d529fe83b744360fcdde4f52986782f9b27e5023bb35be97a42
```

Arbitrary-support discrepancy control, the
universal network ceiling, entry, synchronization, and LRC(14) remain
**OPEN**; neither the small continuous suprema nor the geometric basis
lemma is presented as a solution to those problems.
