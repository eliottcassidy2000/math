# Variable-radius coarea counting: exact owner densities and a sharp bulk constant

**Status: PROVED ELEMENTARY + INDEPENDENTLY AUDITED (root); exact small controls.** This is a standalone
lattice-count statement for every real `r>0`. It does not assume odd speeds.
The count is defined below at every radius; physical network claims below
concern the inherited ternary-unit setting at radius `r=3/14`.

For primitive positive integer `w=(w_1,w_2,w_3)`, put `c=max(w)` and

```text
Lambda_r(w)={C in Z^3: C.w=0, C_i!=0 mod3 for all i,
                         |C_i|<r(sum(w)-w_i) for all i},
N_r(w)=|Lambda_r(w)|.
```

When a speed is divisible by three, this coordinate-unit gate is an algebraic
filter; an identity with the actual distinct-sheet failure comb is not assumed.

Choose **any** primitive nonzero integer relation `v.w=0`; write
`S=||v||_1`, `M=max_i|v_i|`, and choose `i` attaining `M`.
If every speed is a unit modulo three, then

```text
|N_r(w)-(8r^2/9)sum(w)|
 < 4r(sum(w)-w_i)/(3M)+4rS/3+4/3.                    (1)
```

In particular,

```text
N_r(w)/c <8r^2/3+8r/(3M)+4rS/(3c)+4/(3c).           (2)
```

The leading constant `8r^2/3` is sharp, even within distinct positive odd
ternary-unit triples: for `n=5 mod6`,

```text
w_n=(n^2,n(n+2),(n+2)^2),
N_r(w_n)/(n+2)^2 ->8r^2/3                            (3)
```

for every fixed `r>0`. A speed divisible by three changes the answer:
exactly one such speed doubles the bulk density; two such speeds make the
owner-live support empty. The exact residue table and sharper relation-specific
remainder are given below.

## Inheritance and the connection contract

The source is the affine cube-slice mechanism and zonotope quadrature of
[the universal-slope note](lrc14_global_slope_empty_core_certificate_sep06.md),
which feeds
[THM-4434 — universal scale-three network projection bound](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md).
The new target is a two-sided integer-count estimate with variable radius,
arbitrary parity, and all possible speed residues modulo three.

The map is

```text
error cube -> planar zonotope -> scalar defect sections
           -> complete affine integer lines with owner residues.
```

It retains the area integral, central section width, every integer defect,
the primitive scalar lattice spacing, and the number of live owner residues
in each defect class. Replacing individual sections by area plus peak loses
their detailed positions; the peak term and exact defect intercept restore
an explicit arithmetic error bound. Raw multipliers are never discarded.

The concept board is: coarea integral; nearest short relation; residue-fiber
density; two-sided interval discrepancy; balanced large-coefficient families;
strict defect endpoints. The anchor is a reusable counting theorem, the niche
is the exact empty-defect boundary, and the wildcard is the doubled density
when one speed vanishes modulo three.

The inherited hostile is the
[mixed-parity row `(2,5,7)`](lrc14_parity_empty_core_sep06.md). It has a norm-three
owner-unit relation and exceeds the odd-only `6/77` network ceiling. The
present theorem keeps that short-relation correction rather than hiding it
in a bulk asymptotic.

## 1. Every raw carrier and every affine lattice point is retained

Primitivity of `w` makes `n -> w cross n` surjective onto its integer kernel:
if `u.w=1`, then `n=C cross u` satisfies `w cross n=C`. For a lift of `C`,
the strict raw roofs are exactly the pairwise intersection conditions for
the three open intervals

```text
(n_j/w_j-r/w_j, n_j/w_j+r/w_j).
```

The one-dimensional interval intersection property therefore gives
`C=w cross e` for an `e=n-tw` in `(-r,r)^3`, and conversely.

The defect `delta=v.n=v.e` is an integer with `|delta|<rS`, and
`v cross C=delta w`. Because `v` is primitive, its scalar product with
integer vectors takes every integer value. Every defect line is thus an
actual integer affine line `C_delta+Z v`; any two integer points on it
differ by an integer multiple of primitive `v`.

Let `f(delta)` be the scalar length, in units of `v`, of the image of
`{e in [-r,r]^3:v.e=delta}` under `e -> w cross e`. The live points correspond
to the open scalar interval. Set `f=0` at `|delta|=rS` and outside that
support, preserving the strict defect convention.

## 2. Exact owner fibers, including speed residues that were formerly excluded

Modulo three, each affine defect fiber consists of three points. Let
`m_j` be the number having all coordinates nonzero when `delta=j mod3`.
Negation gives `m_1=m_2`. The complete table is:

| Zero speed coordinates modulo three | Relation residue type | `(m_0,m_1,m_2)` | `h=m_0+2m_1` | `kappa=(m_1+|m_0-m_1|)/3` |
|---:|---|---|---:|---:|
| 0 | all relation coordinates nonzero | `(2,0,0)` | 2 | `2/3` |
| 0 | one relation coordinate zero | `(0,1,1)` | 2 | `2/3` |
| 1 | all relation coordinates nonzero | `(2,1,1)` | 4 | `2/3` |
| 1 | some relation coordinate zero | `(0,2,2)` | 4 | `4/3` |
| 2 | every possible relation | `(0,0,0)` | 0 | 0 |

Three zero speed residues contradict primitivity. For zero speed residues
absent, the owner-live kernel contains exactly two opposite points; a
relation direction either contains both or is transverse to their line.
For exactly one zero speed residue, permute and scale unit coordinates to
write the kernel as `(x,y,-y)`. Its owner-live set is `x,y!=0`, four points.
A full-support relation is a diagonal in this two-dimensional residue plane:
its zero fiber has two live points, the other fibers one each. An axis
relation has no live zero-fiber point and two in each nonzero fiber.
With two zero speed residues, the kernel equation forces the remaining
carrier coordinate to be zero, proving emptiness.

This also proves the **exact empty-defect-list criterion**. If

```text
D={delta in Z: |delta|<rS, m_(delta mod3)>0},
```

then

```text
D is empty iff h=0, or (m_0=0 and rS<=1).             (4)
```

The first allowed nonzero defects are `+/-1`; their endpoints are strict.
Emptiness of this complete defect list forces `Lambda_r(w)` empty. A
nonempty list need not realize a carrier, so (4) is not an iff criterion
for carrier-support emptiness.

## 3. Two-sided interval discrepancy and coarea

For an open interval of length `T`, the count `n` in one residue class modulo
three obeys

```text
-1 <=n-T/3 <1.
```

For two selected classes the successive gaps alternate between one and two,
giving

```text
-4/3 <=n-2T/3 <4/3.                                  (5)
```

For the lower bound, if the interval contains `n` selected points, extending
it to the nearest exterior selected points gives length at most `3n/2+2`.
For the upper bound, its interior selected-point span is at least `3n/2-2`.
The same lower bound covers an empty interval since the largest gap is two.
The lower constant is attained, for example, by the open interval `(0,2)`
and residue classes `0,2`. For those same classes, the upper constant is
approached by `(2-epsilon,3+epsilon)` as `epsilon` decreases to zero: the
selected points are `2,3` and the error is `4/3-4epsilon/3`.

Put `beta(0)=0`, `beta(1)=1`, `beta(2)=4/3`, and define

```text
F=sum_(|delta|<rS) (m_(delta mod3)/3) f(delta),
B=sum_(|delta|<rS) beta(m_(delta mod3)).
```

Completeness of the affine lines and (5) imply `|N_r-F|<=B`. When `D` is
empty, `N_r=F=B=0`; no strict inequality `0<0` is used.

The image of the cube under

```text
e -> (v.e,(w cross e)_i/v_i),  v_i!=0,
```

is a centrally symmetric planar zonotope. Its three pairwise generator
determinants have absolute values `w_1,w_2,w_3`, so

```text
I=integral f(t)dt=4r^2 sum(w).                        (6)
```

Its section width `f` is even and nonincreasing on the positive half-line.
Resetting a nonzero endpoint width to zero is a downward jump and preserves
this monotonicity. Rectangle comparison, which needs no endpoint continuity,
gives

```text
|sum_k f(hk)-I/h|<=f(0), h>0.
```

Writing `rho_j=m_j/3`, the exact sampled width is

```text
F=rho_1 sum_k f(k)+(rho_0-rho_1)sum_k f(3k).
```

Consequently the complete two-sided master inequality is

```text
|N_r-(h/9)I| <= kappa f(0)+B.                        (7)
```

This separates geometric volume, relation-dependent sampling, and the exact
finite affine-count intercept. It remains valid in every row of the residue
table and for every `r>0`.

## 4. Explicit relation corrections and the sharp asymptotics

At defect zero, `w cross e=t v`. In a coordinate with magnitude `M`,

```text
f(0)<=2r(sum(w)-w_i)/M.
```

For ternary-unit speeds, `h=2`, `kappa=2/3`. The two possible intercepts obey

```text
B <8rS/9+4/3             if v is owner-unit;
B <4rS/3+4/3             if v has a zero residue.
```

The second upper bound also dominates the first. Substituting in (7) proves
(1)-(2). These bounds hold for odd and mixed speeds alike; the relation norm
is not assumed even. They show exactly why a bounded short relation can
prevent a pure volume approximation.

More generally, along any ternary-unit sequence with relations satisfying
`M->infinity` and `S/c->0`,

```text
N_r/c -(8r^2/9)sum(w)/c ->0.                          (8)
```

Both hypotheses matter: a large coefficient alone does not control the
number of sampled defect lines. For the family in (3), take
`v_n=(n+2,-n,0)`. It is primitive, `M=n+2`, `S=2n+2`, and the three speeds
are primitive, distinct, odd, and ternary units for `n=5 mod6`. Their ratio
sum tends to three, so (8) proves (3). In fact their shortest relation norm
is exactly `2n+2`: reduction modulo `n` and `n+2` writes every relation as

```text
((n+2)s,-ns-(n+2)t,nt), s,t in Z.
```

An axis choice attains `2n+2`; if both parameters are nonzero, the outer
coordinates alone have total magnitude at least `2n+2`.

For exactly one zero speed residue, (7) instead has bulk term
`(16r^2/9)sum(w)`. Here `kappa<=4/3` and the simple universal intercept bound
`B<(8/3)rS+4/3` suffices. The doubled constant is also attained as a limit:

```text
u_n=(n^2+1,n^2+n+1,n^2+2n+2), n=1 mod3,
v_n=(n+1,-2n-1,n),
N_r(u_n)/(n^2+2n+2) ->16r^2/3.                        (9)
```

The first two speeds have gcd one because their difference is `n`; the
middle speed alone is zero modulo three. The displayed relation has
`M=2n+1`, `S=4n+2`, so its correction vanishes after division by `c`.
Thus the speed-unit hypothesis in (2) cannot be silently removed.

## 5. Sharp low-defect and parity controls

For `w=(1,2,5)` and `v=(2,-1,0)`, the complete defect list is empty precisely
when `r<=1/3`. At that boundary the raw support is empty; for every `r>1/3`
the carriers `+/-(-1,-2,1)` are present. At the exact sample `r=3/8`, they
are the complete support. This realizes the strict threshold in (4).

At the inherited radius `r=3/14`, the mixed-parity row `(2,5,7)` instead has
unit relation `(1,1,-1)` and defect zero, so its list never disappears.
Its complete support is `+/- (1,1,-1)`, and the projections and physical mass
are

```text
E=(22/245,6/49,1/10),    physical mass=22/245>6/77.
```

The two-sided count theorem survives; the odd-only network ceiling does not.
The additional tiny certificates `(1,10,11)` and `(2,11,20)` from the parity
note were independently reconstructed by full integer boxes, including their
physical masses `6/55` and `11/140`. No new parity census was run here.

## 6. Quantitative bounded resonance

**Root derivation; independently audited by `three_ray_geometry`.** Let
`w=(a,b,c)` be primitive, positive, distinct, sorted, and ternary-unit, with
no parity assumption. Define the shortest relation norm

```text
lambda_1(w)=min{||v||_1: 0!=v in Z^3, v.w=0}.
```

A minimizing relation is primitive. Put `S=lambda_1(w)` and choose its
maximum coordinate magnitude `M`, so `M>=S/3`. Dividing (1) by `c`, and using
`sum(w)-w_i<2c`, gives a first error term strictly smaller than `8r/S`.
The elementary planar lattice argument in Section 4 of the
[universal-slope note](lrc14_global_slope_empty_core_certificate_sep06.md)
gives `S<4sqrt(c/3)`: its area computation uses positivity, distinctness,
sorting, and primitivity, but not oddness. Hence, for every `r>0`,

```text
|N_r(w)/c -(8r^2/9)sum(w)/c|
 <8r/lambda_1(w)+16r/(3sqrt(3c))+4/(3c).              (10)
```

For **fixed** `r>0`, every sequence with `lambda_1(w)->infinity` therefore
satisfies the bulk convergence (8). Indeed the lattice bound itself forces
`c->infinity`, and all three terms in (10) vanish. No limiting ratios of
the speeds are required for this centered convergence.

Conversely set `e_r(c)=16r/(3sqrt(3c))+4/(3c)`. If the absolute normalized
deviation on the left of (10) is at least `eta>0` and `e_r(c)<eta`, then

```text
lambda_1(w)<8r/(eta-e_r(c)).                         (11)
```

In particular `e_r(c)<=eta/2` forces `lambda_1(w)<16r/eta`. Thus a fixed
positive deviation persisting at large height is confined to finitely many
primitive relation-coefficient patterns. This bounds relation types, not
the heights or the number of speed triples on each relation plane. It
recasts the earlier finite coefficient gate as a quantitative short-relation
principle. The balanced family (3), with exact shortest norm `2n+2`, is the
positive control; bounded-relation families retain the correction term.

This corollary is analytic: the source and frozen exact output below are
unchanged. The referee checked the maximum-coordinate choice, strict roof
bound, parity-free lattice hypotheses, fixed-radius limit, and inversion
with a positive denominator.

## Reproduction

[Source](../../04-computation/variable_radius_empty_core_sep06.py) and
[output](variable_radius_empty_core_sep06.out).

```bash
python3 -B 04-computation/variable_radius_empty_core_sep06.py
python3 -B -O 04-computation/variable_radius_empty_core_sep06.py
```

The verifier exhausts the finite-field fiber table, checks small open-interval
words, and retains complete raw carrier sets for seven named threshold,
parity, and residue controls plus four samples of each limiting family.
Its `997` exact gates pass, and normal and optimized outputs are identical.
The root referee independently audited the residue table, both interval-error
directions, zonotope area and quadrature, exact empty-defect boundary, master
inequality, and both balanced families. The referee repaired the upper-error
sharpness example above to use consecutive selected points; the inequality,
source, and exact output were unchanged.
Raw-byte SHA256:

```text
source a5c6f4d47cbda80fe2d0c6335e093ad76df74d8488f6d7560e8c5a9c512b946e
output 50e9fd2db6dc507d1b20bf8d3d429db52913fb6b1b11d18fb4f79a7f446c3d95
```

The limits in (3) and (9) are proved by (7); the samples are controls, not
evidence extrapolated into those limits. No new theorem ID, speed census,
entry argument, or LRC(14) conclusion is asserted here.
