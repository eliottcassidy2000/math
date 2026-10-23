# Exact ternary observation of a mod-eight guard

23 October 2026. **PROVED:** the elementary algebraic statements below for
an authenticated complete guard table. **FINITE-EXACT:** the accompanying
controls. This is a compiler for supplied guards, not an all-size algorithm
for deriving native electrical potentials or recognizing tournament wiring.

## 1. The correct coefficient ring

Let `R=(Z/8Z)[w]/(w^2+w+1)`. Retain both coefficients of `a+bw`; there are
64 elements. We have `w^3=1`, `1+w+w^2=0`, and

    (w-1)(w^2-1)=3.

Since 3 is a unit modulo 8, both root differences are units. This last fact is
essential: a nonzero Fourier coefficient may be a zero divisor, but cannot
be killed by either nontrivial root difference.

For a function `f:F3^r -> Z/8Z`, use the natural inclusion into R and put

    fhat(k) = 3^(-r) sum_a f(a) w^(-k dot a).
    f(a)    = sum_k fhat(k) w^(k dot a).

The inverse follows by multiplying the one-coordinate orthogonality sums:
they are 3 for zero frequency and 0 otherwise. This is ordinary finite
Fourier inversion with an explicitly suitable ring, not a claim of a new
general Fourier theorem. It avoids floating-point roots and retains all
three bits of the original guard.

For a tuple of guard residues, take the unionS of its nonzero Fourier
supports. A Boolean validity flag is a legitimate additional residue. If a
phase is defined only on valid paths, extend by the record

    (valid(a), valid(a)*phase(a)).

The separate first component distinguishes invalid entries from legitimate
zero phases. Additional original guards require additional components;
integer class IDs must be retained losslessly, for example in enough base-8
digits, rather than reduced modulo 8 without justification.

## 2. Exact dependence, kernel lifts, and what is lost

An axis coordinate i is essential if and only if some k in S has `k_i !=0`.
Indeed changing only that coordinate multiplies coefficientk by
`w^(k_i)-1`, a unit when nonzero. Inversion proves the converse as well.

More strongly, the full translation stabilizer is exactly

    K={v: f(a+v)=f(a) for everya} = (span_F3 S)^perp.

The forward implication uses the same unit argument; the reverse follows
term by term. If B is any row basis of `span S`, then

    f(a)=g(Ba),     each observationBa has3^(r-rankB) lifts.

This is the smallest-rank **linear** observation through which the entire
record factors. If another linear observation C suffices, then
`ker C <=K`, hence `rank C >= rank B`. It is not necessarily the smallest
arbitrary nonlinear encoding. A record fibre can be a union of several
K-cosets; equal output does not itself give a universal translation move.

The axis footprint and the observer rank differ. For `f(a,b)=h(b-a)` with
nonconstant h, both axes are essential but the rank is 1: the diagonal
threefold kernel survives. A pinned coordinate must be substituted first
and the restricted table analyzed anew; restrictions can create or destroy
translation symmetries relative to the remaining free coordinates.

For a supplied actual frame parametrization, this preserves the complete
guard record and the exact kernel multiplicities. It forgets the individual
owner alignment unless the fibre coordinate is retained as a sidecar.
Native path lists, clocks and endpoint class labels remain external inputs.

## 3. Minimal universal coupled-move budget

Give a parameter translation v cost `wt(v)`, the number of changed F3
coordinates. For b>=0 define

    K_b=span{v inK: wt(v)<=b}.

The graph on a fixed coset `a+K`, with every such translation as an allowed
move, has precisely the cosets of K_b as connected components. One direction
is closure of sums; the other writes each vector of K_b as a linear
combination of generators and executes those additions. Consequently the
sharp uniform budget for connecting that coset is

    min{b: K_b=K}.

Equivalently it is the smallest possible largest weight in a basis of K.
For K=0 the budget is 0. This is a parameter-coordinate cost: it must not be
renamed number of moved vertices, native path length, or electrical energy.

In the prior shared-frame automorphism compiler, a parameter change from a
to a+v is realized by the actual automorphism `g_(a+v) g_a^(-1)`. This remains
true although the automorphism group law is nonabelian. Fixed indexed ports
must already have been compiled into the retained coordinates and pins.
The [original phase analysis](../guards/README.md) proves that the native
electrical phases themselves are constant on this automorphism kernel.

This theorem does not connect distinct K-cosets in the same guard fibre.
The finite hostile is the validity table for four cyclic-cell phases with
endpoints 0 and binary consecutive differences. Its accepted internal states
are `(0,0)` and `(1,2)`, corresponding to 000 and 111. Its universal stabilizer
is zero, but one source-specific two-coordinate move connects the two.
The [block theorem](../moves/COUPLED_PHASE_MOVES.md) supplies exactly such adaptive moves.

## 4. Converse, residual trees, and the modulus boundary

Negating a phase does not change its Fourier support, since -1 is a unit.
For an authenticated chart change `a -> Ma+c` with M invertible over F3,
the frequency support is transformed by `M^T`; nonzero coefficients are
multiplied by units. Rank and kernel dimension therefore survive. An axis
footprint need not survive an arbitrary linear chart change. Applying this
to tournament converse additionally requires transporting the actual node
potential and path chart; bare phase negation is not that whole operation.

There is also a lossless ternary decision tree: at each depth, branch on
one actual coordinate 0,1,2. Two prefixes may be merged exactly when their
complete residual guard tables on every suffix agree. These are the unique
coarsest deterministic states for the fixed coordinate order: any merger
of unequal residuals is separated by one suffix. The checker records these
counts. This may be much smaller than the full tree, but need not be. The
incoming native exponential coset-cover and NP-completeness results forbid
inferring a generally efficient native compiler from the existence of a
ternary representation.

At modulus 9 this particular Fourier inversion fails. In
`(Z/9Z)[w]/(w^2+w+1)`, the nonzero function `(3,3,3)` has unnormalized
transform `(0,0,0)`: the zero-frequency sum is 9 and the other two sums
are `3(1+w+w^2)`. Thus the distinction here is invertibility of the
three-branch averaging operator, not a numerical analogy between an 8/9
spectral gap and moduli 8 or 9. Over any coefficient ring where 3 and the
needed root differences are units, the same proof works.

## 5. Decisive tests and scope

Run `python -B 05-knowledge/results/native_guard_moves_farey_20261023/fourier/check.py`.
It exhausts all 512 one-coordinate mod-8 guards and all 512 two-coordinate
Boolean guards. Every case checks inversion, literal translation symmetries,
literal essential axes, quotient lift counts, and low-budget connectivity
by an independent graph traversal. Positive controls include sum compensation
(two-coordinate moves) and a common three-coordinate phase. Hostiles include
an eight-element value fibre with trivial translation stabilizer, the adaptive
000/111 bridge, and the mod-9 annihilation above.

The native consumer independently loads the original order-eight owner map
and harmonic potential. On the pinned two-cell fixture, the complete record
at relative phases 0,1,2 is respectively `(1,1)`, `(1,5)`, `(0,0)`.
Here the components mean validity and valid times original phase. Its
minimal axis footprint contains both coordinates, but its linear observer
has rank one and kernel `span{(1,1)}`. It retains exactly three actual
paths in each valid phase level, connected with sharp phase budget two.
Including the exact target ID as four base-8 digits does not enlarge the
observer. Pinning an actual first-cell owner removes the diagonal freedom.

The checker authenticates Hamiltonian junctions, masks, source and target
owners, endpoint pins and converse directly. Under tournament converse,
the actual path is reversed and the merged endpoint classes remain the
same. Thus the original merged-class phase is unchanged; this is different
from reversing one oriented electrical port, which negates its phase.

The source-specific native consumer and provenance are recorded in the
parent audit. Full-table input can be exponential. No uniform local bound
for the original electrical guard, and no all-size metagraph wiring
classification, is asserted.
