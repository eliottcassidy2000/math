# One exact order for all future O/Z/J constructions

**PROVED IN NOTE.** An exact finite criterion compares the Hamiltonian-path
counts of two rooted states under every common future word, including an
unbounded number of root changes. The actual operations and their tail
formula are inherited; no literature-priority claim is made. This order
does not assert exact-value coverage or prove Chan–Pak.

## 1. The common frame

For a rooted tournament write its state as `(H,E,U,d)`: total Hamilton paths,
paths ending at the root, ordered spanning two-path covers with the marked
path starting at the root (the other may be empty), and Hamilton paths
after deleting the root. Use the inherited physical maps

    O(H,E,U,d) = (H,0,2U,d),
    Z(H,E,U,d) = (H+U−E,d,2U−d,d),
    J(H,E,U,d) = (H+U−E,H−E,H+U,H).

O adds a sink, Z adds a vertex beating only the retained root, and J makes
that same addition and marks the newborn. See
[Session225 recovery, §1](../../session225_root_descent_recovery_20261108/PROOF.md)
and its cited combinatorial graft proof. Words act from left to right.

Make the invertible integer change of coordinates

    q(H,E,U,d) = (h,u,v,t) = (H,U,U−d,H+U−E).

Its inverse is `(H,E,U,d)=(h,h+u−t,u,u−v)`. Direct substitution gives

    Oq = (h, 2u, u+v, h+2u),
    Zq = (t, u+v, 2v, t+2v),
    Jq = (t, h+u, u, h+2u).                         (1)

All three matrices are nonnegative. These are response coordinates, not
a quotient: no information in the four-coordinate state has been lost.

**Theorem.** For any two rooted states x,y, the following are equivalent:

1. `H(Wx) ≥ H(Wy)` for every finite word W in O,Z,J, including the empty word.
2. The same inequality holds for every finite word using only O,Z.
3. `q(x) ≥ q(y)` coordinatewise.

The equivalence also holds for arbitrary real four-vectors acted on by
these linear maps; only genuine tournament states are called constructions.

## 2. Proof and the extremal words

Let β=C mod2 and use the inherited binary function

    g(0)=g(1)=1,
    g(2a+1)=g(a),
    g(2a)=g(a)+a+(a mod2).

The O/Z word with address C reads its bits from least to most significant,
with 0=O and 1=Z. Its output is

    H_C = H+d(g(C)−1)+(U−d)C+(d−E)β.                (2)

This is the old tail formula, in the present g convention. It also follows
immediately by induction on the lowest bit using the displayed O/Z maps.
Let `f(C)=C+1−β−g(C)`. Induction gives `g(C)≥1` and `f(C)≥0`.
Rearrange (2) as

    H_C = H+U(g(C)−1)+(U−d)f(C)+(U−E)β
        = (1−β)h+βt+u(g(C)−1)+v f(C).             (3)

For differences between x and y, C=0 forces Δh≥0 and C=1 forces Δt≥0.
At C=2^k, k≥1, equation (3) gives Δh+2^kΔu, so nonnegativity for every k
forces Δu≥0. At C=2^k−1 it gives Δt+(2^k−2)Δv, forcing Δv≥0. Thus (2)
in the theorem implies (3). Conversely (1) preserves coordinatewise order,
so (3) implies (1). The remaining implication is immediate. QED.

This proof uses only four explicit families of words. No finite search is
being substituted for the universal quantifier.

At the same starting count H the criterion reduces to

    (U, U−d, U−E)_x ≥ (U, U−d, U−E)_y.            (4)

In particular decreasing U alone is not a sufficient comparison when E or
d also changes. A finite list of sources can be ordered by (4), but deleting
the larger sources need not preserve the set of *exact* future values.

## 3. Application to the actual paired reset

The [native construction](../modular_cycles/PROOF.md) replaces each of the
H13-source twins at H=211,213 by a different tournament with that same H:

    old = (H,13,347,13),     new = (H,25,269,25).

Thus

    q(old)−q(new) = (0,78,90,90).                   (5)

By the theorem, the new source is no more expensive under **every common
O/Z/J word**, not just a single tail. It is strictly cheaper as soon as the
word contains Z or J: a leading run of O keeps the first difference zero
and the fourth positive; the first Z or J makes the first difference
positive, and subsequent O/Z/J preserve strict positivity of both the
first and fourth differences. An all-O word leaves H unchanged.

For a pure O/Z address C the exact gap is

    90C−12(g(C)−1) = 78(g(C)−1)+90f(C)+90β
                   ≥ 78C+12β.

The same-count switch shortens the displayed graph order and raises the
minimum deletion count. Therefore neither graph order nor d alone is the
right future-cost measure.

The separate target-preserving ray rewrite from H13 to H5 changes state by
`old−new=(0,8,8,8)`, hence q-difference `(0,8,0,0)`. Equation (1) proves
that it too weakly lowers every future O/Z/J count. Its exact equality
language is

    Z* O*   or   Z* J O* .                         (6)

Indeed Z fixes `(0,a,0,0)`. A first O or J makes the fourth difference
positive while keeping the first zero. Any later Z or J then makes H
strictly larger forever. This identifies precisely when a local equal-count
replacement can be transported without paying a count difference.

## 4. Why the old dynamics connection mattered

The route here came from
[quadratic-cycle phase, §4](../../quadratic_cycle_phase_20261015b.md),
continuing THM-4139 in
`01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md`.
Its two order-six lifts preserve different hexagonal forms, while their
product is hyperbolic and at a specified shear is Pell-conjugate. The
lesson used here is constructive: find one frame in which **all** source
operations preserve the proposed comparison. Equation (1) supplies exactly
that frame, rather than merely testing each operation in its preferred norm.

The source is a changing-generator matrix problem; the target is the rooted
tournament response module. The map is the explicitly invertible q change
of basis. The preserved predicate is all-future count comparison. Its lost
information is the actual graph (already lost in the inherited response
state), so each new source reset still needs literal graph witnesses. The
cheap decisive test was substitution of J into the same q coordinates,
which produced a nonnegative matrix rather than a frame-change cost.

The Collatz [joint-port note](../../collatz_hexagonal_growth_20261102/JOINT_PORTS.md)
provides a complementary construction lesson: two source branches may share
an actual tail, and the shared endpoint does not determine their two clocks.
Here equal H is the shared endpoint; `(E,U,d)` retains the separate continuation
responses. This is a transfer of method, not a claim that tournament counts
are Collatz iterates or modular parameters.

## 5. Reproduction and the remaining implication

Run `python3 -B 05-knowledge/results/h13_paired_descent_20261111/dynamics/probe.py`.
The controls compare literal O/Z/J recurrences with the changed coordinates,
use independent bit recursion for (2), test the two actual reset differences,
and retain a witness for necessity of each of the four inequalities. The
all-word result is proved above, not inferred from these finite controls.

**OPEN:** a mechanism that turns such a decrease into a target-preserving
descent for an arbitrary missing odd N. Changing a source before a future
word usually changes N. The positive cone precisely measures this cost; it
does not silently restore it. The H13-to-H5 ray provides one exact equality
language, and further equality languages or a count-restoring operation are
needed for a completeness proof.
