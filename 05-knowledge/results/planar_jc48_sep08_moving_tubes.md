# Exact moving-centre tubes on one complete frozen braid path

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.** This changes the certificate
representation, not the polynomial curve, its paths, or any existing
frozen engine. It supplies an actual configuration homotopy and a full
word check, rather than endpoint root matching alone.

## 1. The uniform certificate

Let `F(u,z)` be a monic polynomial of degree `n` in `z`. On a straight
base segment and one affine centre path write

    u(theta)=u0+theta du,       z(theta)=z0+theta dz,
    F(u(theta),z(theta)+w)=sum C_jk theta^j w^k,
    0<=theta<=1.

For a rational complex number define `A(z)=|Re z|+|Im z|` and
`B(z)=max(|Re z|,|Im z|)`. A constant radius `r>0` is certified if

    B(C_01) r > sum_(j,k)!=(0,1) A(C_jk) r^k.             (1)

This bounds **every** other coefficient, including higher powers of
`theta` and all varying linear terms. Since `B<=absolute value<=A`,
Rouché's theorem applied against `C_01 w` gives exactly one zero
inside the moving disk for every parameter, counting multiplicity.
In particular that zero is simple. Cancellation of the first residual
alone is not the criterion.

For each pair of affine centres the prototype computes exactly

    min_(0<=theta<=1) max(|Re(d0+theta delta)|,
                          |Im(d0+theta delta)|).         (2)

This is the minimum of a convex piecewise affine function. A minimum
occurs at an endpoint, a coordinate zero, or a solution of
`Re=Im` or `Re=-Im`; a constant interval includes such an endpoint
or breakpoint. Evaluating this finite list therefore gives (2),
not a sampled lower bound. If each pair's radius sum is strictly
less than (2), their Euclidean disks are disjoint for all parameters.
The prototype uses each centre's minimum pair distance divided by
16, a sufficient fixed choice without an optimal-radius claim.

With `n` such disks, all `n` roots are exhausted. Continuity and
the one-root property give labelled root paths. Consecutive pieces
share the identical rational centre at their joint endpoint, so
their endpoint disks are concentric and nested; both contain one
root, which must be the same one. This pays label gluing without
assuming that arbitrary nearby endpoints have the intended labels.

Contract each actual root to its affine centre within its moving
disk. Disjointness makes this a homotopy in the configuration space.
At shared endpoints the contractions agree. If the starting and
ending unordered centre sets coincide, this is a homotopy of loops,
with the base configuration identified by the same endpoint paths.
Consequently the exact centre polygon computes the actual braid.

## 2. The operation and its retained data

The source is the previously audited rational witness for the
complete smooth loop in
[the two-cusp braid bundle](planar_jc48_sep08_two_cusp_braid.md).
The target is a shorter sequence of segments on precisely the same
six straight base edges. The map discards intermediate vertices,
then independently certifies the full moving tubes. It retains the
literal polynomial, all labelled endpoint centres, base loop,
complex projection, and full braid word. The discarded fixed-disk
subdivision is not used as a proof of the new inequalities.

The prototype imports no old continuation or crossing functions.
It verifies the old source hash, safely reads its literal coefficient
table through the AST, and checks equality with its own table and
the complete monic degree-four condition. The input gzip and raw
JSON hashes are also checked. Its expansion and rational crossing
extraction are standalone. Greedy doubling/bisection produces valid
chunks; no monotonicity of acceptability or maximal coarsening is
asserted. Every accepted chunk is checked directly.

## 3. Positive and hostile controls

The closest mechanism is the frozen fixed-centre Rouché tube; its
expansion loses common root motion when bounding the base increment
separately. For `F=product_(a=0)^3(z-u-a)`, taking `dz=du` makes
every positive-`theta` coefficient vanish exactly. The moving
certificate therefore permits arbitrarily large common translations
while the centre separation stays fixed. The literal control uses
`du=10^6(1+i)`. This is a structural cancellation, not an assertion
that every moving tube improves the previous bound.

Two cheap hostiles retain the missing information. Swapping the
endpoints of two constant roots keeps all endpoint separations
positive but makes their affine centres coincide halfway; (2)
rejects it. For

    F=(z-u^2(1-u))(z-2)(z-3)(z-4),

the first root has zero endpoint centres and `C_10=0`, yet at
`u=2/3` it equals `4/27`, outside the proposed radius `1/16`.
The complete higher-term bound (1) rejects the segment. Thus neither
endpoint isolation nor linear residual cancellation alone is a tube.

The live concepts are common root velocity, exact centre separation,
endpoint labels, configuration homotopy, and the full projected word.
The least-used sidecar is the affine-centre separation envelope.

## 4. Named finite result and reproduction

The complete inherited smooth loop has **9,553** fixed-tube segments.
The prototype constructs **92** moving-tube segments in **654** exact
attempts and obtains the identical full reduced word `[3,1,-3]`.
This is a reduction in certificate segments, not a universal runtime
claim. An initial instrumented run took approximately 3.1 seconds
on this host. Timing is omitted from the frozen deterministic output.

Run from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_moving_tubes.py
python3 -B -O 04-computation/planar_jc48_sep08_moving_tubes.py
```

Optional `--progress` reports edge counts and wall time. The output
contains a semantic digest of every chosen endpoint index, radius,
and positive Rouché margin. Its acceptance universe is this one
complete named path plus the declared synthetic controls, not every
previous or incoming certificate.

Normal and optimized runs pass **10,324 always-active gates** and
agree with all 374 bytes of the frozen output.

| Artifact | SHA-256 |
| --- | --- |
| Source, 11,199 bytes | `735e928d0d39bc43d6f1ec3e3d8a4965b8181813cf6094d6e8b50bb8fcf59990` |
| Frozen output and both replays | `14606ab9af4e2b6d25f26b4b6c7c7638c9577c51c4eb8f652bf87cf717d12ea5` |
| Semantic segment certificate | `77c95276801a4096438b0c2c018124b2bf35fb3323b7f572689f99a6752f73e9` |

The [independent method/source audit](planar_jc48_sep08_moving_tubes_audit.md)
accepts the full common-base configuration homotopy and exact replay.
It also independently reconstructs the resultant and checks affine
coefficient expansions. Its frozen SHA-256 is
`8bebfabf1b42227d37214e58b1321bcdc1afb168748aa0a4a6cd9dceccc64629`.
Status promotion did not alter the source or frozen output.

The safe next operation is to measure this independent representation
on more frozen paths after audit. No running infinity-eleven producer
or old frozen certificate was modified.
