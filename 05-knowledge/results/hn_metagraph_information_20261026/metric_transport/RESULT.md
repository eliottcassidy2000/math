# Planar three-pair transport already carries a four-state section

2026-10-26. **PROVED elementary (MT1--MT3, scoped below), independently
[audited](../enumeration/METRIC_AUDIT.md).** Exact computations
are **FINITE-EXACT**, with two independent metric tests. This is a restriction
on a proposed realization mechanism, not a six-colouring theorem for arbitrary
unit-distance graphs.

## Inheritance and the decisive extra coordinate

The [previous palette result](../../hadwiger_six_seven_eight_20261025/fano_eight/RESULT.md)
gives eight endpoint choices `A=F2^3`, parity `p(x)=sum x_i`, and the split
`W=F2^3 semidirect S3 = C2 x S4`. The four states are the antipodal pairs
`{x,x+111}`, equivalently `e=x+p(x)111` in `E=ker p`.
The [previous metric result, M1 and section 3](../../hadwiger_six_seven_eight_20261025/metric_spindle/RESULT.md)
classifies the seven-image unit cube as the regular flower and proves that
only its collision-pair stabilizer `C2 x S3` descends. Those are inherited
results. The extension here is a canonical section for **every** nondegenerate
planar three-pair frame, followed by its consequence for arbitrary networks.

Closest mechanism: the flower's unique antipodal collision fibre. Hostile:
the two abstract parity-zero holonomies with disjoint fixed pairs. Corrected
near miss: an abstract permutation of endpoint choices need not preserve a
planar frame. Least-used sidecar: the **signs of the unique linear dependence**.

## MT1. A canonical balanced antipodal pair

A geometric frame consists of a centre `o` and six distinct points
`o +/- u_1, o +/- u_2, o +/- u_3`, with `|u_i|=1` in the Euclidean plane.
Distinctness excludes parallel direction pairs. The following proof also
works with unequal positive lengths provided the three direction lines remain
distinct. No unit edge is identified with a palette inequality by definition.

The matrix with columns `u_i` has rank two. Its kernel is spanned by

```
c = (det(u_2,u_3), det(u_3,u_1), det(u_1,u_2)).             (MT1)
```

Every `c_i` is nonzero. Let `b_i=1` if `c_i<0`, and `b_i=0` otherwise.
Then `sum |c_i| (-1)^b_i u_i=0`, so the origin lies strictly inside the
triangle obtained by choosing endpoints according to `b`. Conversely, any
endpoint triple whose triangle contains the origin strictly gives a positive
dependence; one-dimensionality of the kernel forces its sign pattern to be
`b` or `b+111`. Thus exactly these two endpoint choices are balanced.
Multiplying `c` by a nonzero scalar leaves the unordered pair unchanged.

An isometry between two whole frames sends their centres to each other,
because each centre is the mean of its six points. Suppose its signed
permutation is `g(x)=pi(x)+a`, using the target-coordinate convention
`R u_i=(-1)^(a_{pi(i)}) v_{pi(i)}`. Applying R to (MT1) sends the kernel
coefficients to `c'_{pi(i)}=(-1)^(a_{pi(i)}) c_i`. Therefore

```
g {b,b+111} = {b',b'+111}.                               (MT2)
```

This is also immediate from preservation of strict containment in a triangle.
It is a geometric **four-state section**, not merely a parity class. In the
regular flower it is the collapsed cube pair. For a generic frame the cube
subset sums can all be distinct, and the balanced pair still exists.

## MT2. Parity is sufficient within the whole-frame isometry class

Let any finite graph carry one geometric frame at each vertex. For each
oriented edge `u->v`, prescribe an isometry of the whole six-point source
frame onto the whole target frame, with endpoint map `g_e(x)=pi_e(x)+a_e`.
Reverse edges use inverse maps. Loops and parallel edges are allowed. We
ask for endpoint sections `x_v=g_e(x_u)`, retaining all eight choices.

Then the full section problem is soluble **if and only if**

```
eta_v = eta_u + p(a_e)                                 (MT3)
```

is soluble. Necessity follows by taking parity. For sufficiency, use MT1:
the balanced pair at each vertex contains exactly one endpoint choice of
each parity, since `p(111)=1`. Choose its unique member of parity `eta_v`.
Equation (MT2), together with (MT3), forces every seam equation. Each
connected component has at least two full sections when parity is soluble;
there may be additional sections. This proof works for arbitrary graph size,
not just the finite triangle census below.

Equivalently, after choosing local four-state coordinates that place the
balanced state at the origin, every metric seam acts through a point
stabilizer `S3 < S4`. The residual `E`-translation cocycle is a coboundary.
The full abstract `S4` obstruction is thus unavailable in this construction
class. The chosen origin is supplied by geometry, not by an illicit global
frame choice on an arbitrary abstract bundle.

## MT3. The inherited two-loop obstruction fails a metric gate

Take `pi=(0 1 2)`, `g_0(x)=pi(x)`, and `g_5(x)=pi(x)+101`.
Both have parity zero. Their fixed sets are `{000,111}` and `{011,100}`,
respectively; they have no common fixed point. They also have no common
invariant antipodal pair. For a parity-zero map, preserving an antipodal pair
forces each of its members to be fixed individually. By MT1, these two maps
cannot both be whole-frame isometries at **any** planar three-pair frame.

Each map can occur separately in a suitably marked regular flower. Thus
testing each loop in its own geometry misses the obstruction: they need a
single shared local frame. A map such as `g_5` can also be an isometry
between differently marked flowers. It is the consistent network, not the
individual signed permutation in isolation, that MT2 constrains.

## Exact controls and failure boundary

`python3 05-knowledge/results/hn_metagraph_information_20261026/metric_transport/check.py`
uses rational arithmetic and tests all 48 signed coordinate permutations
between all 64 ordered sign-marking pairs for each of three frames:

| Frame | Self-isometries | Balanced pair | Cube subset-sum images |
|---|---:|---|---:|
| Regular A2 flower | 12 | 0, 7 | 7 |
| Generic rational unit directions | 2 | 2, 5 | 8 |
| Reflection-symmetric rational directions | 4 | 1, 6 | 8 |

The first metric test compares Gram matrices; an independent test compares
all 15 pairwise squared distances of the six endpoints. The cofactor-kernel
construction is separately checked against triangle containment signs.

All `4^3 * 12^3 = 110592` marked regular-flower triangle networks are tested.
Full eight-state section counts are obtained from the composed maps, without
assuming the parity verdict. The histogram is `0:55296, 2:18432, 4:27648,
8:9216`; parity passes exactly in the nonzero cases. Every passing network
retains its distinguished balanced pair. These finite controls supplement
the all-graph proof, rather than extending a bounded census by inference.

The hypotheses matter. Collinear/duplicated spokes destroy the nonzero
kernel-sign construction. More importantly, arbitrary unit-distance gadgets
need not transport an entire antipodal frame by an isometry: partial boundary
relations, changed pairings, and interior vertices are outside MT2. A
successful six-colour obstruction must therefore supply its full unprescribed
boundary relation and its geometry, rather than identify an abstract seam
permutation with rigid planar transport. No claim of universal failure for
those other construction classes is made.

## Connection contract and next decisive test

Source: three signed planar direction pairs. Target: the four-state palette
bundle. Map: take the sign pair of the rank-two kernel, or the balanced
endpoint triangle. Preserved predicate: a common endpoint section under
whole-frame isometries. Destroyed information: lengths, angles, interior
vertices, and proper-colouring extension counts. Required sidecar: the full
metric frame and the exact signed seam map. Cheapest hostile test: require
both inherited three-cycle holonomies at the same frame; MT3 rejects them.

The next positive search should use genuinely relational seams that escape
the point-stabilizer restriction, while checking their actual six-colouring
tables. This theorem does not replace those tables with geometric labels.
