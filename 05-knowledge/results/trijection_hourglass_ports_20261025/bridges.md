# The common operation is retaining the obstruction coordinate

2026-10-25. **CITED / INHERITED PROVED** where linked; the connection ledger
below is a methodological synthesis, not an equivalence between open problems.

## The geometric carrier and the root cover are different geometries

The [audited divisor carrier](../reflection_bridge_20261025/geometric_carrier.md)
is U=(P2 x P1 x P1) minus two opposite coordinate sections. Its ten
rank-zero points retain the two divisor reflections. The seven nonempty
coordinate supports belong to its P2 factor, but only three of those seven
are rank-zero points. The divisor order needs a coordinate flag; it is not
the specialization order of ten closed points. Over a finite field of size q,

    #U(F_q)=(1+q+q^2)(1+q)^2-2.

This toric open construction is not the seven-root hourglass surface. The
[root hierarchy](../session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md)
is a complete intersection of n-3 quadrics in P^(n-1), over characteristic
zero with distinct branch lines and no fourfold crossing. Adjunction gives

    K_X=O_X(n-6),       K_X^2=(n-6)^2 2^(n-3).

Its six-root resolutions are K3; seven roots give p_g=7,K^2=16.
The canonical threshold is explained by 2(n-3)-n=n-6. The seven-root
sign group has rank six for another explicit reason: seven independent
coordinate signs modulo the one simultaneous projective sign. Neither
calculation identifies a K3 surface with the preceding toric open variety.

In the recovered [code construction and new four-bit model](codes/README.md), seven branch labels are
assigned the seven nonzero vectors of F2^3. Here Fano incidence is genuine:
three labels are dependent exactly when their vectors sum to zero. This
addition is additional structure absent from a bare seven-support monoid
space. It detects the actual singularity of the quotient above a triple
crossing. This is a precise place where the seven-point object does more
than count seven things.

## The Monster partition has a specified action

The recovered [valuation and modular calculation](../session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md)
and [Sylow fixed-point interpretation](../session32_square_monster_20261023/SYLOW_FIXED_POINT_PARTITION.md)
already explain the user's 4,5,6 partition. It is not necessary to guess
another partition of the fifteen primes.

Take the Monster action on its 2A involutions, whose point stabilizer is
2.B, the double cover of the Baby Monster. The three regimes are:

| Primes | Size | Fixed-point interpretation in this action |
|---|---|---|
| 29,41,59,71 | 4 | No subgroup of order p fixes a point |
| 17,19,23,31,47 | 5 | Every Sylow p-subgroup fixes some point |
| 2,3,5,7,11,13 | 6 | Some subgroup of order p fixes a point, but no full Sylow p-subgroup does |

This follows from the p-adic orders of the group and its stabilizer. The
group orders were rechecked against the primary
[Monster ATLAS](https://brauer.maths.qmul.ac.uk/Atlas/spor/M/) and
[Baby Monster ATLAS](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/) on
2026-10-25. The stabilizer has the extra factor2; omitting it changes its
2-adic depth even though the three sets stay the same. These are inherited
results, not discoveries of the current session.

There is also a rigorous arithmetic bridge: supersingular elliptic curves
with level structure can characterize prime orders in Monster centralizers;
see Aricheta, [Supersingular Elliptic Curves and Moonshine](https://arxiv.org/html/1809.07421),
Theorems1.1 and1.3. The recovered exact computation shows that each of the
four absent Baby primes acquires one Frobenius-conjugate pair of coarse
level-two isogeny points despite rational underlying supersingular j-points.
Forgetting the isogeny identifies away this obstruction. The cited theorem
and inherited finite computation have different scopes.

## Connection ledger and decisive tests

| Source -> target | Map and retained predicate | Information lost / required sidecar | Decisive hostile |
|---|---|---|---|
| Marked tournament -> ordered graph | Backward-edge bits, with the observer retained; every original arc is recoverable | A new path requires its actual joins and reversal, not just a new coordinate graph | A tree-chart switch can change H from5 to1 |
| Root sign group -> F2^3 | Quotient sends branch signs to seven distinct nonzero labels | Actual crossing triples must be retained | A dependent triple gives two 1/4(1,1) singularities instead of one A1 |
| Abstract box -> physical plane support | Labelled support embedding | Physical line partitions, which can merge disjoint abstract lines | At p29,L11 the kernel gains two modes beyond all adjacent cubes |
| Supersingular isogeny data -> j-invariant | Forget level structure; endpoints retained | Frobenius on distinct edges is erased | The four primes have rational vertices but a nonrational parallel-edge pair |
| Divisor lattice -> rank-zero toric points | Exponents to labelled coordinate points; reflections retained | Flags are needed for divisibility, and support points are a different functor | Ten rank-zero points versus61 literal pointed-monoid points |

The strongest research route shared by these examples is to make the
quotient computable together with the missing local constraint. For the
black-core problem this means carrying positional reversal through every
intermediate path. For the hourglass it means carrying the branch matroid.
For the incidence box it means updating line partitions as length grows.
These are concrete proof obligations; no common scalar count settles them.

The newly proved detour and code statements live in the companion packets.
All-order black-core connectivity, an unrestricted small A000568 recurrence,
Collatz and the Riemann hypothesis remain open in this work. No reduction
between these open statements is claimed.
