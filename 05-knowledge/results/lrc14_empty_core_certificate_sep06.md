# Short-relation cube slices certify the LRC14 network through norm twenty

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This note proves a statement about the complete degree-zero network
projections of THM-4414. It does not prove entry, synchronization, or LRC(14).

Let `w=(a,b,c)` be primitive, distinct, positive, odd, nonzero modulo three,
and sorted. Suppose there exists a primitive integer relation `v.w=0` with
every coordinate nonzero and `||v||_1<=20`. Then

```text
min_i E_i(w)<=6/77.
```

Equality occurs only at `(1,5,11)`. If the complete carrier support has at
least two primitive directions, the stronger conclusion holds:

```text
E_i(w)<6/77 for every i.
```

The same stronger conclusion holds for a one-ray support whose primitive
direction has norm other than four. Empty support has three zero projections.
The result includes families with arbitrarily many primitive live directions.

## Inheritance, scope, and the new operation

The closest proved mechanism is the affine defect decomposition from
[THM-4386 — canonical component relation and zero-defect incidence](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md),
with the one-zero residue dichotomy of
[THM-4398 — one-zero relation residue dichotomy and small norm atlas](../../01-canon/theorems/THM-4398-lrc14-one-zero-relation-residue-dichotomy-and-small-norm-atlas.md).
Those physical comb bounds do not directly control the larger network sums.
The consumer here is the exact raw formula in
[THM-4414 — six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md)
and the sufficient count gate from
[THM-4422 — projection deficit and Beatty-row reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md).

The inherited hostile is `(1,5,7)`: two projections exceed `6/77`, so an
every-projection statement needs the norm-four exception. The corrected near
miss is retaining only one live ray in `(17,23,25)`; the omitted ray contributes
real capacity. The least-used sidecar is the *length of the whole projected
error slice*, including its affine defect. The slice, not a chosen live
point, is the legal object to replace by its two endpoints.

The concept board is: raw carrier; affine integer defect; cube slice;
mod-three interval word; convex speed normalization; exact network count.
The anchor is the network ceiling, the niche is affine lattice counting, and
the wildcard is replacing an entire fiber by its extremal endpoint data.
The connection contract is exact: project the cube slice to the carrier line,
preserving all possible raw addresses and the scalar interval width; retain
the defect and owner residue as sidecars; discard individual roof values only
after applying the common upper cap `q=3/(7c)`. This loses sharpness, not mass.

## 1. Exact affine-line decomposition

Put `r=3/14` and let

```text
Lambda(w)={C in Z^3: C.w=0, C_i!=0 mod3,
                         |C_i|<r(sum(w)-w_i), all i}.
```

By THM-4386, every such `C` is `w cross n` for an integer lift `n`, modulo
integer multiples of `w`. Its component has positive length exactly when
some real `y` gives `e=n-yw` in the open cube `(-r,r)^3`. Define

```text
delta=v.n=v.e,
P_delta(v)={e in [-r,r]^3: v.e=delta}.
```

For completeness, the error-cube implication is elementary. The three
allowed `y` intervals are `((n_i-r)/w_i,(n_i+r)/w_i)`. The `j,k` intervals
overlap strictly exactly when
`|w_j n_k-w_k n_j|<r(w_j+w_k)`, the `i`th raw roof condition. Pairwise
intersecting open intervals on a line have a common nonempty intersection:
their maximum left endpoint is strictly below their minimum right endpoint.
Thus all three strict roofs give a real `y` and an interior error vector.
This uses actual positive speeds and does not assume that arbitrary points
of a bounding box lie in the cube image.

The scalar `delta` is an integer independent of the lift representative and
obeys `|delta|<r||v||_1`. The vector identity

```text
v cross C = v cross (w cross n) = delta w
```

shows that fixed-defect carriers lie on an affine line parallel to `v`.
Conversely, since `v` is primitive, `v.n=delta` has an integer solution for
every integer `delta`; any two integer carriers on this line differ by `kv`,
`k in Z`. Thus no lattice index is omitted.

The image `w cross P_delta(v)` is a closed segment on that affine line.
Its length in units of `v` is

```text
T_delta(v,w)=max_(e in P_delta) (w cross e)_i/v_i
             -min_(e in P_delta) (w cross e)_i/v_i,    (1)
```

for any coordinate `i`; every `v_i` is nonzero. It is coordinate independent
because differences of image points are multiples of `v`. Live carriers lie
in the corresponding open interval. Closed slices are used only to compute
the width, so a zero-roof boundary is never counted as live.

## 2. The owner word on each affine interval

There are exactly two possibilities for a primitive full-support relation.
If all `v_i` are units modulo three, its three terms `w_i v_i` are equal in
`F_3`, and pairwise distinct owners imply `delta=0 mod3`. Conversely, on a
fixed such defect line the residues of `C_0+kv` run through the line spanned
by `v` modulo three; precisely two of the three residues have every
coordinate nonzero. These are the live owner residues.

If exactly one `v_i` is zero modulo three, the other two `w_i v_i` are
nonzero and opposite. The owner difference makes `delta!=0 mod3`.
Conversely the affine residue line for a fixed nonzero `delta` meets the
two-element owner set in exactly one point, so exactly one `k mod3` is live.
Two zero coordinates are impossible because `w.v=0 mod3` would force the
third coordinate to vanish. Three zero coordinates contradict primitivity.

On an open interval of length `T`, one residue class modulo three contains
strictly fewer than `T/3+1` integers. Two residue classes contain strictly
fewer than `2T/3+4/3`: if there are `m` selected integers, the span between
the first and last is at least `3m/2-2`, because consecutive gaps alternate
between one and two. Open endpoints make the inequality strict.

Consequently, with `D(v)` the complete allowed defect list,

```text
N=|Lambda(w)| < rho_v sum_(delta in D(v)) T_delta(v,w)+B_v,       (2)

unit relation:     D={delta: |delta|<r||v||_1, delta=0 mod3},
                   rho_v=2/3, B_v=4|D|/3;
one-zero relation: D={delta: |delta|<r||v||_1, delta!=0 mod3},
                   rho_v=1/3, B_v=|D|.
```

Every defect and every affine multiplier is retained. No assumption bounds
the number of primitive directions in the complete support.

## 3. A finite rational polytope certifies the entire speed cone

For fixed `v,delta`, formula (1) is a width of the fixed polygon
`P_delta(v)` under a linear functional in `w`. It is therefore convex,
continuous, and positively homogeneous as a function of `w` on `v.w=0`.
The sum in (2) has the same properties. Normalize `max_i w_i=1` and enlarge
the speed domain to

```text
W_v={w in [0,1]^3: v.w=0}.
```

A convex function on this polygon attains a maximum at some vertex:
express any point as a convex combination of vertices and use convexity.
It is harmless that vertices may have zero or equal speeds; they are used
only for an upper bound on the original positive distinct speed domain.

All vertices of both `P_delta(v)` and `W_v` arise by fixing two cube
coordinates at endpoints and solving for the third. They are rational,
and can be exhaustively computed without numerical optimization. Define

```text
A_v = max_(w a vertex of W_v)
           rho_v sum_(delta in D(v)) T_delta(v,w).
```

Take the maximum also over the three possible isolated coefficient signs;
permutations need no separate speed optimization because `W_v` already
allows every coordinate ordering. The resulting `A_pattern` depends only
on the coefficient magnitude pattern. Equation (2) gives

```text
N < A_pattern c+B_pattern.                             (3)
```

Every projection summand is at most `q`, so `E_i<=qN`. If `A_pattern<2/11`,
then every projection is strictly below `6/77` whenever

```text
c >= B_pattern/(2/11-A_pattern).                       (4)
```

The maximum is an exact finite rational calculation for a fixed relation,
not a census over speed triples and not an asymptotic extrapolation.

## 4. Complete full-support shell through norm twenty

The primitive positive magnitude triples with even sum at most twenty and
at most one coordinate divisible by three give exactly `73` patterns. Odd
sum cannot annihilate three odd speeds. The unique pattern failing
`A_pattern<2/11` is `(1,1,2)`, the already proved norm-four sector of THM-4422.
For all remaining `72` patterns, (4) is a complete finite-head reduction.
The largest required cutoff is

```text
pattern (1,7,8): A=15/98, B=4,
c >=4312/31 <140.                                     (5)
```

Every other cutoff is smaller. The accompanying output records all `72`
rational maxima, maximizing signed relation/speed vertices, defect lists,
intercepts, cutoffs, and finite-head row counts. The verifier generates all
signed coordinate permutations, solves their equations for `c`, and keeps
exactly the eligible primitive sorted rows below the relevant cutoff.
It reconstructs the complete raw support, computes every exact projection,
checks the affine-layer count bound for every relation presentation, and
compares all heads with an independent literal six-sheet contact engine.

The finite heads satisfy the selected target, with equality only at
`(1,5,11)`. Every head row whose support is empty, multi-ray, or on a ray
outside norm four has all three projections strictly below the target.
There are `4,083` relation-pattern/head memberships and `1,944` distinct
head triples: `33` empty, `21` norm-four, `733` other one-ray, and `1,157`
multi-ray supports. The largest individual projection outside the norm-four
class is `12/161` at `(1,19,23)`; among these finite multi-ray heads it is
`18/301` at `(5,37,43)`. These leader values are finite-head statements only.
All infinite tails are strict by (3)-(4). The remaining norm-four sector
has the same unique equality by THM-4422. This proves the stated theorem.

## 5. A transparent new infinite subclass

For signed relation magnitudes `(1,2,3)`, the defects are just `+/-1`.
After changing error signs, the positive-defect slice is

```text
x_i in [-3/14,3/14], x_1+2x_2+3x_3=1.
```

Writing `t_i=3/14-x_i` gives `t_1+2t_2+3t_3=2/7`, a triangle;
the other cube walls are inactive. The three projected vertex differences
have scalar magnitudes `(2/7)w_i/(|v_jv_k|)`, hence

```text
T_1=T_(-1)=(2/7) max_i(|v_i|w_i)/6 <=c/7.
```

There is one owner residue per affine line, so

```text
N <2c/21+2,
E_i <2/49+6/(7c) <6/77 for c>=25.
```

This counts both whole affine lines. For example `(1,401,601)` has `38`
carriers on `19` primitive directions; `(1,1201,1801)` has `116` carriers on
`58` directions. Both satisfy `a-3b+2c=0`, and the complete family is covered
at arbitrary height. Thus the result is not a bounded-direction extension
of the earlier one-ray and two-ray theorems.

## 6. Independent audit of incoming one-ray and two-ray work

The incoming
[one-ray note](lrc14_one_ray_overnight_hexagon_sep05.md) and
[two-ray note](lrc14_two_ray_overnight_hexagon_sep05.md) have sound analytic
reductions. Their strict multiplier cutoffs, raw sign counts, primitive
coefficient bounds, and finite-head thresholds have the correct directions.
The two-ray determinant has the load-bearing factor three: full owner
vectors are parallel modulo three and their cross product is an integer
multiple of primitive `w`. Its reciprocal inequality uses `M_u,M_v>=7`.
The conclusion counts exactly two primitive directions, not the dimension
of their real span. The independent literal engine is a declared shared
verification dependency and is not evidence of producer independence.

The source THM-4386 inverse/address and zero-defect arguments are sound,
including the use of strict error inequalities to force an integer multiple
of three to vanish at norm fourteen. THM-4414's spacing lemma is sufficient
and sharply separates length feasibility from mere forest structure; its
two owner orientations must remain separate before the final raw sum.
No correction was found in these load-bearing arguments.
Both incoming producer replays also reproduce their checked-in output bytes.

## Reproduction and remaining scope

```bash
python3 -B 04-computation/lrc14_empty_core_certificate_sep06.py
python3 -B -O 04-computation/lrc14_empty_core_certificate_sep06.py
```

The complete output is
[lrc14_empty_core_certificate_sep06.out](lrc14_empty_core_certificate_sep06.out).
The source is
[lrc14_empty_core_certificate_sep06.py](../../04-computation/lrc14_empty_core_certificate_sep06.py).
The only imported mathematical implementation is the explicitly independent
literal-sheet audit from the one-ray verifier. Assertions are not used for
checks, so optimization cannot remove the validity gates.

Normal and optimized runs pass `97,780` explicit checks and have byte-identical
output. Frozen raw-byte SHA256 values:

```text
source 9d098b57d8852fb2b152c8120ecd2c4ee38ef78fecb87de7f11554ab59305d59
output 63a8bba2f5d15f647194720328319b84bb856ea74dd912f60cf33af03b433e62
```

The separately written
[independent audit](../../04-computation/lrc14_relation_slice_audit_empty_core_three_ray_sep06.py)
uses rational polygon clipping instead of cube-edge enumeration, a different
cross-product coordinate, every signed permutation, every finite-field owner
fiber, a complete height-first universe of `16,089` eligible triples through
`139`, and independent integer-box carrier enumeration. It reproduces all
coefficient constants, `4,083` memberships, and `1,944` unique proof rows,
including every projection and the physical mass. Its `13,383` explicit
checks survive optimization, and normal/optimized output bytes agree.

```bash
python3 -B 04-computation/lrc14_relation_slice_audit_empty_core_three_ray_sep06.py
python3 -B -O 04-computation/lrc14_relation_slice_audit_empty_core_three_ray_sep06.py
```

```text
independent source cdcc95a63788191135e0d341c5f72b83a61b302cf867c34497653f48649b71a8
independent output 573eebbd23cdab2d567c5e1e6df4e519677c8904cc7eb71f1de9b85b73107903
shared semantic head digest ff9eb3c383b719cd5e2fd767d3ec8858f1b851df4c6f4d97dca0fa390c788168
```

The finite coefficient shell and coefficient-dependent heads are part of the
proof. The two large affine-family examples are controls, not completeness
evidence. Relations with zero integer coordinates are outside this theorem.
No theorem supplies a full-support relation of norm at most twenty for every
eligible triple. Arbitrary long relations, the universal projection target,
entry, synchronization, and LRC(14) remain open.
