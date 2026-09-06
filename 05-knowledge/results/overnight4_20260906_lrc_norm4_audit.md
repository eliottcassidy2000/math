# Independent audit of the sharp norm-four physical and selected-network ceiling

**Verdict: PASS -- analytical tail, pointwise fixed roof, exact finite base,
and saturated physical coarea.** This audits Sections 3--4 of the
[root native report](overnight4_20260906_lrc_parityfree_native.md).
It does not independently re-audit that report's entire general coefficient
classification or infer an LRC(14) theorem.

The audited all-height result is, for sorted **distinct** positive primitive
ternary-unit triples admitting a relation of absolute pattern `(1,1,2)`,

```text
physical mass = min_i E_i <= 11/140,
equality exactly at (2,11,20).
```

Oddness is not assumed. This is a physical/selected-projection ceiling,
not a lower bound for a convex-body representation. The exact finite base
has 82 distinct-speed cases; permitting weak order adds `(1,1,1)`, which
has mass zero and is retained only as a negative control in this audit.

## 1. The fixed roof is pointwise, not merely an integral identity

Put `q=3/(7c)` and, for `x>=0`,

```text
L_i(x)=[3(w_j+w_k)-14|v_i|x]/(14w_jw_k),
f(x)=max(0,min(q,L_1(x),L_2(x),L_3(x))).
```

The three sorted sectors exhaust the absolute coefficient pattern:
`c=2a+b`, `c=a+2b`, and `2b=a+c`. Relations doubling the smallest or
largest speed alone on one side would force equality of all speeds, and
the remaining sign arrangements violate their strict order.

Let `d` be the doubled coordinate, and let its complementary speeds be
`p,c`. The candidate closed profile is

```text
T(x)=max(0,min(q, 3(p+c)/(14pc)-2x/(pc))),
u0=3(c-p)/28, u1=3(c+p)/28.
```

It equals the full profile for every real `x>=0`. The following table gives
the complete endpoint check; pairs of values refer to the two competing
roofs in their natural coordinate order.

| Sector | Doubled coordinate | `u0,u1` | Competing values at `u0` | Competing values at `u1` |
|---|---:|---|---|---|
| `c=2a+b` | `a` | `3a/14,3(a+b)/14` | `3/(14a),3/(14a)` | `3/(7c),0` |
| `c=a+2b` | `b` | `3b/14,3(a+b)/14` | `3/(14b),3/(14b)` | `3/(7c),0` |
| `2b=a+c` | `b` | `3(b-a)/14,3b/14` | `q,3/(7b)` | `3/(14b),3/(14b)` |

Every competing value at `u0` is at least `q`; in the first two sectors
this uses `c>=2a` or `c>=2b`, respectively. Their monotonicity therefore
puts them above the cap on `[0,u0]`. On `[u0,u1]` their differences from
the doubled roof are affine and nonnegative at both endpoints. Beyond
`u1` the doubled roof is nonpositive. Thus `f=T` on the entire half-line,
and evenness handles negative indices. Selection commutes with summation
here because this **one fixed** roof works everywhere. It does not follow
from equality of integrated areas alone.

At radius `3/14`, the strict norm-four defect bound is `|v.n|<6/7`, so
every integer defect is zero. The full live carriers are exactly `kv`
with `3` not dividing `k` **and** `f(k)>0`. The finite support gate is
essential. The fixed projection indexed by `d` therefore equals physical
mass, while every other projection dominates it.

## 2. Exact integral and the strict residue-counting error

The trapezoid gives

```text
I = integral_R f(|x|) dx = q*(u0+u1)=9/98.
```

For `R>0`, let `N` be the largest integer strictly below `R`. The number
of positive integers below `R` not divisible by three is `N-floor(N/3)`.
For the three possible residues of `N`,

```text
N-floor(N/3) <= 2N/3+2/3 < 2R/3+2/3.              (1)
```

This is the centered specialization of the two-residue open-interval
discrepancy bound already proved in Section 3 of the incoming
[variable-radius coarea report](variable_radius_empty_core_sep06.md).
No new priority is claimed for the arithmetic constant.

For every level `0<y<q`, the superlevel set of this continuous even
profile is an open interval `(-R(y),R(y))`. Doubling (1) and integrating
over `y` yields

```text
mass < (4/3) integral_0^q R(y)dy +(4/3)q
     = (2/3) I +(4/3)q
     = 3/49+4/(7c).                              (2)
```

The strict inequality survives integration: (1) is strict at every
positive superlevel radius, including integer endpoints, and `q>0`.
For `c>=33`, (2) is strictly below `11/140`; at 33 the exact gap is
`1/32340`. Equivalently the real cutoff is `560/17<33`.

## 3. Independent finite base and equality control

The [audit source](../../04-computation/overnight4_20260906_lrc_norm4_audit.py)
generates every `1<=a<b<c<=32` with primitive speeds, no zero residue
modulo three, and one of the three displayed identities. It does not
derive its universe from the comparison TSV. There are exactly 82 rows.

For each row it evaluates the closed trapezoid literally at all live
integer indices. It independently builds the full minimum of all affine
roofs and compares at **every** pairwise intersection, every roof zero,
and both closed-profile breakpoints. Both functions are affine between
those endpoints, so this also supplies a complete pointwise finite
certificate rather than a chosen grid of samples.

All 82 masses agree with the independent native head's physical entries,
all fixed doubled-coordinate projections equal the mass and selected
projection, and every strict live-carrier count agrees. The sole maximum
is `(2,11,20)`, with

```text
f(1)=3/140, f(2)=1/56,
mass=2*(3/140+1/56)=11/140.
```

The tail is strict, so no equality can occur above this finite base.
The older candidate ceiling `6/77` is explicitly violated by this same
triple; its status is not revived by the new tail estimate.

## 4. Independent saturated-plane and coarea audit

Let `w` be primitive positive and `v` primitive integral with `v.w=0`.
Choose integral `z` with `z.w=1` and put `n1=v cross z`. The vector triple
product gives `w cross n1=v`. This also proves saturation without assuming
it from a drawing: for an integer `m` in `v^perp`, write
`m=alpha*w+beta*n1`. Then `w cross m=beta*v` is integral, and primitivity
of `v` forces `beta` integral. Primitivity of `w` then forces `alpha`
integral. Thus `(w,n1)` is a basis of the full plane lattice.

The map `(t,y) -> t*n1-y*w` has constant Euclidean area Jacobian `||v||`.
Fubini therefore gives, for every `r>0`,

```text
integral_R length{y:t*n1-y*w in [-r,r]^3} dt
 = area([-r,r]^3 intersect v^perp)/||v||.           (3)
```

The parameter `y` is the physical line parameter, not Euclidean length
along `w`; no extra `||w||` belongs in (3). Integer saturation fixes the
carrier multiplier spacing. Boundary choices do not change the integral.

The cross-product right inverse is inherited from Section 1 of the
variable-radius report. Its Section 3 integrates a different object:
**carrier width over defect**, giving `4r^2 sum(w)`. Here (3) integrates
**physical interval length along a single carrier line**. The distinct
fiber length is retained; replacing it by a projection width would destroy
the physical mass predicate.

Solving the plane equation for a coordinate of largest absolute coefficient
reduces (3) to planar square area divided by that coefficient. For `(1,1,1)`
the square loses two triangles of total area `r^2`, giving `3r^2`. For
`(1,1,2)` the whole square survives and division by two gives `2r^2`.
For `(1,2,2)`, two triangles of area `r^2/4` each are removed; division
by two gives `7r^2/4`. The checker reconstructs these clipped polygons
with exact rational edge intersections, independently integrates full
physical envelopes, and obtains all three constants at radii `3/14`,
`1/2`, and `2/3`.

## 5. Reproduction and scope

[Normal output](overnight4_20260906_lrc_norm4_audit.out) and
[optimized output](overnight4_20260906_lrc_norm4_audit_optimized.out) accompany
the source. Reproduce from the repository root:

```text
python 04-computation/overnight4_20260906_lrc_norm4_audit.py
python -O 04-computation/overnight4_20260906_lrc_norm4_audit.py
```

All **2,833 explicit gates** pass. The exact head is used only below the
proved cutoff; it is not extrapolated. This audit accepts the root native
report's fixed-roof, norm-four ceiling, coarea, and equality claims. Its
all-parity general-sector classification remains covered by the separate
referee cited there. No shared navigation or canon was edited.

Source SHA-256:

```text
5a2213bdafeb4a7c9829dab683bc325b0e407bed57a05b925a05408c887a3485
```

Compared raw TSV SHA-256:

```text
c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa
```

Normal and optimized outputs are byte-identical after explicit LF
normalization, with SHA-256:

```text
180c431b309b17a5656910d4cb43289d072dabfe51c0d8e6844799a9b62af3b9
```
