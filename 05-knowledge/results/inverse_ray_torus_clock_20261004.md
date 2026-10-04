# Torus direction layers and guarded inverse rays share an exact quotient clock

2026-10-04. **PROVED** the group, source-torsor, prime-power geometry and
recursive tournament statements below. **FINITE-EXACT** for the declared
complete coset and cellular controls. No historical priority claim and no
transfer of a finite phase cycle to an integer Collatz cycle or convergence
theorem.

Artifacts: [script](../../04-computation/experiments/inverse_ray_torus_clock_20261004.py)
and [output](inverse_ray_torus_clock_20261004.out).

## 1. Inheritance and the precise new map

Incoming commit `78c99c600`, now integrated, adds
[checked switches and torus layers, section 6](checked_switch_phase19_20261004.md).
It partitions the Paley directions at a prime `p=7 mod12` into cosets of
the order-three subgroup, builds a triangular torus per coset, and computes
the fourth-power layer periods at 19, 5779 and 87211. It warns that the
87211 clock visits only nine of 14535 layers.

The [inverse-ray codec](inverse_ray_ternary_addresses_20261004.md) retains
the actual hub, source row, exponent block and first-hit certificate.
The [mod19 observer](mod19_recursive_observers_20261004.md) proves the
source-aware address selector, including collapsed hubs. The
[resonance-depth theorem](mod19_resonance_depth_20261004.md) identifies
`R(p)=max(0,v3(ord_p(2))-1)` as the shared ternary block depth.

The added bridge is an explicit commuting map, rather than equality of
periods: **cubing a multiplicative direction coset turns the layer action
by 4 into the inverse-source block action by 64**. An anchored affine
change then recovers the source residue itself.

Closest mechanism: a quotient that retains its acting element. Canonical
hostile: nine visited layers do not mean nine total tori. Corrected near
miss: a cube coordinate recovers an H-coset, not a distinguished direction.
Least-used sidecars: the hub normalization, the initial layer, and the
additive fibre at a higher prime power. The board is **direction / layer /
block / hub / valuation / cellular link / exact source**.

## 2. Cubing is a genuine layer isomorphism

Let `p=7 mod12` be prime, let Q be the nonzero squares in `F_p`, and let
`H={1,omega,omega^2}` be the subgroup of order three. Since an element of
order three is a square, H lies in Q. The incoming torus layers are Q/H;
there are `L=(p-1)/6` of them.

**PROVED quotient isomorphism.**

    c: Q/H -> Q^3,          c(qH)=q^3.                 (1)

It is well-defined, multiplicative, onto its stated image, and has trivial
kernel: the kernel of cubing on Q is exactly H. Thus it loses no layer
label. It does forget which of the three directions in that layer was used.
In particular

    c(4qH)=64*c(qH).                                  (2)

Write `d=ord_p(2)`. Every orbit of layer multiplication by 4 has length

    m=ord_p(64)=d/gcd(d,6).                            (3)

Indeed (1) conjugates it to multiplication by 64 in Q^3. A return of any
point is a return of the multiplier itself, so every orbit has this same
length. There are exactly L/m orbits.

| p | Total layers L | Direction clock ord(4) | Layer/block clock m | Layer orbits | Layers missed by the orbit of H |
|---:|---:|---:|---:|---:|---:|
|19|3|9|3|1|0|
|5779|963|2889|963|1|0|
|87211|14535|27|9|1615|14526|

At 87211 a chosen orbit misses 1614 **other orbits**, not merely some
representatives of one cycle. Each missed orbit would require its own
initial coset label if the whole layer space were to be represented.

This cubing occurs in the **multiplicative direction group** `F_p^*`.
It is different from the incoming phase-field operation
`(eta^x)^3=eta^(3x)`, which multiplies an additive vertex label by 3.
Neither operation identifies a cubic map on vertices with a Collatz map.

### The direction coordinate cannot always be recovered multiplicatively

Since Q is cyclic of order 3L, cubing has a multiplicative section
`Q^3 -> Q` iff `3` does not divide L. In that case exponentiation by
`3^(-1) mod L` gives a section. Otherwise the unique order-L subgroup of
Q contains H, so cubing is not injective on a possible section image.

At p=19 the order-three cube 7 has square cube roots 4, 6 and 9, all of
order nine. No group homomorphism can select one while preserving the
order-three relation. A set-theoretic representative can still be chosen;
its multiplication then needs an H-valued correction.

This is a global section obstruction, not a statement that every orbit
loses an additional phase. At p=127, Q has order 63 and the global extension
does not split, but the chosen direction clock and layer clock both have
order seven. Its 21 layers form three seven-cycles.

## 3. An anchored source torsor makes the inverse-ray map literal

Fix a positive odd hub u, with `3` not dividing u, and fix one source row.
The inherited exact inverse channel is

    n_b=(2^(kappa+6b)*u-1)/3,        b>=0,
    U(n_b)=u,                      n_(b+1)=64n_b+21.

Here kappa in `{1,...,6}` is fixed by the hub modulo 9 and the source row.
If u has a first-hit certificate, so do all these sources, with the single
root exception `u=1,kappa=2,b=0` removed.

Assume additionally `p` does not divide u, and put `alpha=2^kappa*u mod p`.
For the whole layer space define

    Phi_alpha(qH)=(alpha*q^3-1)/3 mod p.              (4)

This is a bijection onto the source-address set
`{x:3x+1 in alpha*Q^3}`. That set is a translate and rescale of Q^3; it is
a torsor with the chosen origin `n_0`. Direct substitution gives

    Phi_alpha(4qH)=64*Phi_alpha(qH)+21.               (5)

For the actual inverse channel,

    (3n_b+1)/(3n_0+1)=64^b mod p,
    Phi_alpha(4^b H)=n_b mod p.                      (6)

Thus (4) does more than match clock orders: it maps the marked layer orbit
to the actual source residues and intertwines the two operations. The
fixed inverse channel covers exactly the orbit of H under this marking.
The other layers in (4) are legitimate residue coordinates, but are not
members of that channel merely because they occur in the torsor.

The unit condition concerns the **hub**, not the source. The certified
source `(2^18-1)/3=87381` is zero modulo 19, but `3n+1=1 mod19` is a unit
and normalization is harmless. Conversely, if `p|u`, every source is
`-1/3 mod p` and (6) would divide zero by zero. To recover the phase one
must retain `v=v_p(u)` and source precision at least `p^(v+1)`, then use
the peeled units `(3n+1)/p^v` and `2^kappa*u/p^v`. This is extra information,
not an inverse operation on the collapsed residue.

Advancing b is an **inverse-family parameter change** adding six halvings
to one edge. It is not forward odd time: all n_b map to the same u.
The actual integer, the unbounded b, the hub and its certificate remain
necessary. Residue equality does not establish exact membership in the ray.

The resonance theorem now has a geometric interpretation. By (3),
`v3(m)=R(p)`. A ternary source address at precision `3^a` retains
`b mod3^(a-1)`, so its compatibility with this layer orbit is agreement
modulo `3^min(a-1,R(p))`. At 5779 and 87211 this quotient has order nine
when a>=3, despite their very different full layer spaces. The quotient
does not identify the remaining factor 107 of the 5779 clock.

### Every positive shared depth has an oriented geometric prime

The inherited resonance theorem defines, for r>=2,

    t=2^(3^(r-1)),       C_r=(t^2-t+1)/3,

and proves that every prime divisor p of C_r has
`ord_p(2)=2*3^r`. In particular every such prime is `1 mod6`, and its
shared depth is `R(p)=r-1`. No primality assumption on C_r is needed.

There is always at least one such divisor with `p=7 mod12`, so the
oriented torus and tournament construction applies at every positive
shared depth. Indeed t is divisible by four, hence
`3*C_r=t^2-t+1=1 mod4`, giving `C_r=3 mod4`. In its prime factorization
at least one prime is `3 mod4` with odd multiplicity. Since every prime
factor is `1 mod6`, that prime is `7 mod12`.

This is an existence statement, not a claim that all factors qualify.
For example the inherited exact factorization
`C_4=163*135433*272010961` contains 163 congruent to 7 modulo 12 and the
other two factors congruent to 1 modulo 12. The construction does not
choose a canonical prime or identify the distinct full layer spaces.

## 4. The same construction extends over the residue ring Z/(p^k)

Keep `p=7 mod12` prime and put `m=p^k`; this is the ring `Z/(p^k)`, not
the field with `p^k` elements. Let Q_m be its square units. Each square
unit has exactly two square roots, and -1 is not a square because it is
not one modulo p. Therefore `|Q_m|=phi(m)/2`.

The three roots of `X^3=1` modulo p lift uniquely to all p^k: their
derivative `3X^2` is a unit. Call the lifted subgroup H_m. It lies in Q_m.
Cubing again has exactly this kernel and gives

    Q_m/H_m ~= Q_m^3,          [q] -> q^3,
    [q] -> [4q]               corresponds to z ->64z.   (7)

There are `phi(m)/6` unit-direction layers. Choose a lifted omega and
consider the additive map

    Z[zeta] -> Z/m,       a+b*zeta -> a+b*omega,
    zeta^2+zeta+1=0.

Its kernel has the integer lattice basis `(m,0),(-omega,1)`, hence index
m. This translation lattice quotients the usual triangular tessellation
to a connected orientable torus. Multiplication by a unit q relabels its
vertices; its three positive directions are qH_m and its six directions
are `+/-qH_m`. They are distinct already modulo p. Each layer has

    (vertices,edges,faces)=(m,3m,2m).

These tori partition exactly the pairs whose difference is a unit. They
do **not** yet partition the complete graph when k>1.

### Valuation strata supply all remaining edges

For a pair with difference of valuation `v<k`, retain its common additive
fibre `a mod p^v`. Write its vertices as `a+p^v X` and `a+p^v Y` and use
the unit-direction construction modulo `p^(k-v)`. There are p^v fibres,
each with `phi(p^(k-v))/6` tori. Consequently every valuation stratum has
exactly

    phi(p^k)/6 tori, each with p^(k-v) vertices.        (8)

The strata and direction cosets are disjoint, and every distinct pair has
a unique valuation and fibre. Thus they partition all edges of K_(p^k).
Their total edge count is

    sum_(v=0)^(k-1) [phi(p^k)/6]*3*p^(k-v)
      =p^k*(p^k-1)/2.

At `19^2=361`, this gives 57 tori with 361 vertices and 57 with 19 vertices.
Their edge counts are 61731 and 3249, totaling 64980. Separating all surface
layers gives 114 disjoint tori. At each of the original 361 glued vertices,
the cellular link is a disjoint union of 60 hexagons:

    sum_(v=0)^(k-1) phi(p^(k-v))/6 = (p^k-1)/6.

These are **face-incidence links**, not induced neighbor subgraphs. The
glued one-skeleton is K_361; the glued object is not a single torus.

Reduction modulo `p^(k-1)` gives p-sheet coverings on strata `v<=k-2`.
Each parent direction layer has p possible child direction layers. The top
stratum `v=k-1` instead collapses to vertices. At 361->19, the 57 large
tori cover the three parent tori, while all 57 small tori collapse into
their 19 fibre points. There is no blanket surface-covering statement for
the entire normalized union.

Direction-layer branching is also different from clock-orbit branching.
The order of 64 grows by either 1 or p on a one-level lift. A parent
clock orbit therefore has p child clock orbits in the flat case, and one
in the growing case. At 19->361 there are 19 child direction layers per
parent layer, but the single clock orbit of length three becomes one
clock orbit of length 57. An additive-fibre label must additionally be
retained when discussing nonunit strata or a global vertex dilation.

## 5. The stratified directions define an intrinsic recursive tournament

For distinct x,y in `Z/(p^k)`, let `v=v_p(y-x)<k`. Orient

    x -> y  iff  ((y-x)/p^v mod p) is a nonzero square.  (9)

The normalized residue is independent of the integer lifts of x,y.
Swapping them multiplies it by -1, a nonsquare modulo p, so exactly one
orientation holds. Only the diagonal ties. Every positive direction in
every torus of (8) agrees with (9).

Write vertices in base p, least significant digit first. At the first
place where two digit strings differ, their digit difference is exactly
the normalized residue in (9). Thus this tournament is the iterated
lexicographic product of the prime Paley tournament with itself. This
is a proved digit recursion, not a field identification or an arbitrary
orientation of tied unit residues. It is regular of outdegree
`(p^k-1)/2`. The script checks every pair at p=7,19 and k=1,2, including
all 64980 pairs of the 361-vertex tournament.

The hypothesis `p=7 mod12` is doing two jobs: order-three direction roots
exist, and the square half is asymmetric. At p=13 the group cubing
quotient still exists, but -1 is a square. Its two H-cosets in Q give the
same unoriented six-direction layer. The square rule no longer orients
exactly one edge per pair, so it is not the tournament construction above.

An unoriented survivor exists for every prime `p=1 mod6`, and likewise
for its prime powers. Write U for all units and mu_6 for the six roots of
unity. Then `[x] -> x^6` is an isomorphism `U/mu_6 -> U^6`, taking layer
multiplication by 2 to multiplication by 64. The kernel is exactly mu_6;
these are the correct unoriented six-direction layers. This includes
primes congruent to 1 modulo 12 without asserting a quadratic tournament.
At 13 it gives two distinct unoriented layers, whereas the two Q/H labels
above duplicate just one of them. This algebraic boundary does not alter
the oriented theorem or the declared computational universe.

## 6. Exact audit and transfer boundary

Reproduce:

```text
python -X utf8 -B 04-computation/experiments/inverse_ray_torus_clock_20261004.py
python -O -X utf8 -B 04-computation/experiments/inverse_ray_torus_clock_20261004.py
```

The script enumerates every square-direction coset and every layer-clock
orbit at primes `{7,19,31,127,5779,87211}`. It checks the full torsor
conjugacy, independently computes orders by all prime-divisor exclusions,
and compares one complete marked layer orbit with actual positive root-ray
integers at each prime: 988 exact first-hit edges in total. It does not
enumerate an 87211-vertex triangular mesh.

It independently enumerates every edge, triangle, and cellular vertex link
in all valuation-stratum tori for `(p,k)=(7,1),(19,1),(7,2),(19,2)`, proves
the complete-graph edge partitions by exact equality, and checks the
digit-recursive tournament orientation for every distinct pair. The
19->361 direction fibres and clock cycles are compared separately.
Controls include the nonsplit cube section at 19, the locally unshortened
clock at 127, zero source versus collapsed hub, peeled normalization,
and the p=13 orientation failure. All checks survive optimization.

| Source -> target | Exact map and preserved predicate | Lost information / retained sidecar |
|---|---|---|
| Direction coset -> cube coordinate | `[q] -> q^3`, a group isomorphism | Individual direction; retain an H-coordinate if needed |
| Marked layer orbit -> inverse source residues | (4)–(6), exact block action | Actual magnitude and exponent quotient; retain hub, row, b and certificate |
| Layer clock -> ternary quotient | `b mod3^R`, exact address compatibility | Other clock factors and unvisited initial layer orbits |
| Prime-power vertices -> torus strata | Difference valuation, fibre and unit H-coset | Gluing forgets a surface label and creates disconnected cellular links |
| Prime-power digits -> tournament | First differing digit, equation (9) | Itineraries and certificates are not encoded by the orientation |

The phase is now tied to an actual guarded family by an explicit formula.
It remains a finite observation of that family. Establishing a certificate
for an independently supplied integer requires its exact source identity
and legal route, not only its layer, cube class, or recursive digit address.
