# Prime-field rigidity and the missing global period checks

Date: 2026-10-19. Status: **PROVED** for the statements below, with **FINITE-EXACT** controls. This is a separate addendum to [the unit-distance carry checkpoint](../unit_distance_carry_20261019/README.md), published in commit `068e07831225888397f60fc0640ba64824ee70a6`. It does not change that checkpoint or claim a new global coloring bound.

## 1. The full rigidity rank is exact at every odd prime p>=5

Let ell=F_(p^2), p>=5 prime, and let G be its norm-one Cayley graph. The rigidity matrix R has one row per undirected unit edge and two coordinates per vertex, over F_p.

**Proposition 1.** The kernel of R consists exactly of the two translations and the infinitesimal rotation. Consequently

    rank R = 2p^2-3,
    dim stress = p^2(p+1)/2 - 2p^2 + 3.

Choose two nonparallel unit directions u,v and write positions as xu+yv, with x,y in F_p. For an infinitesimal motion h, set H_u(x,y)=Tr(c(u)h(x,y)) and H_v(x,y)=Tr(c(v)h(x,y)). The u-edge equations make H_u independent of x, so H_u=f(y). The v-edge equations likewise give H_v=g(x). These two projections determine h, since the trace pairing is nondegenerate and u,v are independent.

There is a third unit direction w=a u+b v with a,b nonzero: each one-dimensional F_p-line contains at most two unit vectors, while p+1>=6 unit vectors exist. The w-edge equations give

    a [f(y+b)-f(y)] + b [g(x+a)-g(x)] = 0.

The two summands depend on separate variables, so both are constant. Since any nonzero element generates the additive group of F_p, f and g are affine. Their slopes have opposite signs. There are exactly three free parameters, already supplied by translations and rotation. This proves the proposition.

The field restriction is essential to this proof: for F_(p^f), a nonzero additive step need not generate the whole field. No extension to nonprime q is asserted.

At p=3 there are only two direction pairs. The same two-direction calculation gives a six-dimensional motion space and rank 12, matching the six independent triangle parity checks in the first checkpoint.

For p=29, the graph has 841 vertices and 12,615 edges. Its rigidity rank is **1,679**, and its stress dimension is **10,936**. The 435 disjoint straight-cycle stresses are independent. There are therefore exactly **10,501 additional independent stress conditions** after those cycle checks.

## 2. Equality is a position problem with fixed error values

Fix arbitrary base representatives X_v modulo p^2. For an oriented edge e=(v,w), define

    c_e=(N(X_w-X_v)-1)/p modulo p,
    (Rh)_e=Tr(c(w-v)(h_w-h_v)).

Here c inside the trace denotes conjugation; the vector c_e denotes the carry. A section X_v+p h_v has normalized norm-error vector z=c+Rh. For any section from depth one to depth two, write this error on edge e as z_e in F_p. The sum of these errors around any straight p-cycle is -2. This is independent of the initial vertex representatives and of edge orientation: reversing an edge changes both factors in its rigidity row and leaves its norm unchanged.

Thus equality in the universal edge-loss bound is equivalent to the following requirements:

1. Exactly one edge is selected from every straight p-cycle.
2. Its error is exactly -2; all other errors are zero.
3. The resulting error vector z lies in the affine coset c+im(R).

The first two conditions alone do not imply the third. At p=29, the candidate universe is 435 choices of one position in a 29-edge cycle, rather than arbitrary field-valued sparse error vectors. A basis of the remaining stress space gives 10,501 linear syndrome tests over F_29. This describes a finite problem, not an efficient exhaustive algorithm.

There is also a smaller way to test any selected support. First solve the equations on two nonparallel unit-direction families by integration along their straight cycles. Their cycle sums already vanish, so these equations are always consistent. They determine the two motion projections up to one offset per line: **2p unknown offsets** remain. The other direction families then give an ordinary linear system in those offsets.

The two-family system has rank 2p(p-1); Proposition 1 shows that the remaining system has rank exactly 2p-3. At p=29, each candidate support can therefore be tested using **58 unknowns**, with rank **55**, instead of solving afresh in 1,682 vertex-correction coordinates. The remaining thirteen direction families supply 10,933 equations. This is an exact elimination, not a relaxation.

## 3. Local wheel checks can all pass while a global period fails

For p>=5 with p=2 modulo 3, use triangular coordinates with norm

    Q(x,y)=x^2+xy+y^2.

Its discriminant -3 is nonsquare, so this is the same unramified quadratic norm setting. The three unit directions

    u=(1,0), v=(0,1), w=(-1,1)

form a triangular direction family. Consider this especially simple candidate support, with edge origins reduced modulo p:

- select a u-edge exactly when its origin has x=0;
- select a v-edge exactly when its origin has y=0;
- select a w-edge exactly when its origin has x=0.

It has exactly one selected edge on every straight cycle, hence 3p selected edges in this family. Prescribe error -2 on each selected edge.

At each vertex, the surrounding unit hexagon and its six spokes have a local stress: rim weights +1, spoke weights -1. The candidate support satisfies every such stress equation **over the integers**. Indeed, for each of the three edge directions separately, the two selected-rim indicators equal the two selected-spoke indicators, as translated indicator functions. An exact regular hexagon lifts using a primitive sixth root of unity, so its carry syndrome is zero. Thus all these local wheel checks pass.

Nevertheless, this support is impossible for a section. Here is a short global certificate.

Use vertex representatives (x,y) with 0<=x,y<p. Put

    H_u=2h_x+h_y,   H_v=h_x+2h_y.

Integrating the first two edge families gives

    H_u(x,y)=A(x)+f(y),  A(x)=-2 when x!=0 and A(0)=0,
    H_v(x,y)=B(y)+g(x),  B(y)=-2 when y!=0 and B(0)=0.

The w-edge carry is

    c_w(x,y) = -1_(x=0) - 1_(y=p-1).

Its prescribed error is -2 times 1_(x=0). After subtracting the already determined A and B contributions, the remaining equation is

    g(x-1)-g(x) - f(y+1)+f(y) = T(x,y),
    T(x,y) = -3*1_(x=0) + 2*1_(x=1)
             +2*1_(y=0) - 1_(y=p-1).

Sum over x in F_p. The g terms telescope and the f terms are multiplied by p, so the left side is zero. The right side is **-3+2=-1**, a contradiction. Summing over y similarly gives the opposite defect +1.

This exhibits a missing global period condition even after every disjoint-cycle check and every local hexagonal-wheel check passes. It does **not** prove the minimum section defect exceeds its cycle bound: another support may satisfy the global checks, as the q=5 witness in the first checkpoint does. In particular, the failed seam ansatz must not be promoted to a universal nonattainment claim.

At p=29 the norm-one directions split into five triangular families, since 6 divides 30. One may put this seam support independently in each family. The resulting 435-edge support passes every straight-cycle check and every wheel check within those families, but is rejected by the displayed global-period certificate already inside one family. This is an explicit hostile to a local-only decoder.

## 4. Precisely one period is missing from cycles and wheels

The preceding obstruction admits an exact dimension statement for the
three-direction triangular torus of section 3, with p>=5 prime and
p=2 modulo 3. This subsection concerns
that subgraph, rather than all p+1 unit directions at once.

**Proposition 2.** Its full stress space has dimension p^2+3. The translated
hexagonal-wheel stresses span a space of dimension p^2-3. Adding all straight
cycle stresses increases the dimension to p^2+2. Thus cycles and wheels miss
**exactly one** global compatibility condition.

The full rank is 2p^2-3 by the proof of Proposition 1, which already uses only
these three directions. There are 3p^2 edges, giving the full stress dimension.

For the wheel count, let lambda(x,y) weight the wheel centered at (x,y).
The coefficient on the edge in direction w originating at (x,y) is

    lambda(x-1,y)+lambda(x,y+1)-lambda(x,y)-lambda(x-1,y+1).

If every edge coefficient vanishes, this mixed difference is zero, so
lambda(x,y)=f(x)+g(y). The u-edge coefficient then becomes

    g(y-1)+g(y+1)-2g(y),

and the v-edge coefficient becomes the analogous second difference of f.
They vanish precisely when f and g are affine functions on F_p. The space
of lambda is therefore three-dimensional, proving wheel rank p^2-3.

To determine the missing conditions, let r be any edge-residual vector that
annihilates all cycle and wheel stresses. Integrate its u- and v-edge equations
to obtain a correction h_0. This is possible because those line sums vanish.
Subtract R h_0. The remaining residual is supported only on the w-edges;
call it T(x,y).

The translated wheel equations reduce exactly to

    T(x+1,y)+T(x,y-1)-T(x,y)-T(x+1,y-1)=0.

Therefore T(x,y)=F(x)+G(y). The w-cycle checks say

    sum_x F(x) + sum_y G(y) = 0.

Full solvability asks whether T can be written as

    g(x-1)-g(x) - f(y+1)+f(y).

A cyclic first-difference equation over F_p is solvable exactly when its
right-hand side sums to zero. Thus the final condition is

    P(T)=sum_x T(x,y_0)=sum_x F(x)=0.                 (P)

It is independent of y_0. It is also independent of the choice of h_0:
changing h_0 adds differences of the displayed form, whose period is zero.
Moving a constant from F to G changes neither sum, since p=0 in F_p.
Once P(T)=0, both sums vanish and the required f,g exist. There is exactly
one possible remaining obstruction. The seam of section 3 has P(T)=-1,
so that obstruction is genuinely nonzero for some residuals. Equivalently,
the quotient of the cycle-and-wheel-compatible residual space by im(R) is
one-dimensional. Duality gives the claimed stress span p^2+2.

At p=29, one triangular family has stress dimension 844, wheel rank 838,
and cycle-plus-wheel rank 843. These are characteristic-p counts derived
here, not characteristic-zero torus counts imported by analogy.

The five families at 29 have disjoint edge sets. All stresses supported inside
individual families therefore span dimension 5*844=4,220. The full graph's
stress dimension is 10,936, leaving codimension **6,716** for compatibility
between families. Checking only cycles and wheels has dimension 5*843=4,215:
five missing within-family periods and those 6,716 cross-family conditions
remain. This quantifies why local checks alone cannot certify the full lift.

## 5. Controls and stopping point

Run `python3 verify_followup.py`; it uses only the standard library. It checks the seam support, all cycles, all local wheel equations, the full residual identity, and the two nonzero global period sums at p=5,11,17,23,29. It also directly
checks the three stress dimensions in Proposition 2 at p=5 and p=11. At p=29 it also constructs the full five-family 435-edge hostile and checks all 435 cycles and 4,205 wheel equations exactly.

The first checkpoint's q=3 and q=5 witnesses remain positive controls for the fixed-error formulation. The new result proves a sharp exact rigidity dimension, reduces candidate testing to 58 unknowns at 29, and identifies a concrete global condition that local wheels lose. It leaves the actual minimum defect at q=29 open. A future attempt should search supports with global-period compatibility built in; repeating the naive seam construction or merely adding more local wheel checks cannot solve that gap.

An exploratory mixed-integer feasibility run on the wheel-only subsystem reached its five-second cap at p=29 without a candidate. This is not evidence of infeasibility: the explicit seam support above supplies a candidate for that subsystem. The exploratory solver and its outputs are intentionally excluded from the publication subset. The proved obstruction concerns compatibility with the full rigidity system.
