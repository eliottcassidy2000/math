# Independent boundary and polarity audit

Verdict: PASS, with no corrections required. The frozen producer manifest is
9b52b7d7495d418de0f561b8ae7e95c08736f316b2d0c7926a01a54d08248439.
The producer directory and repository were not edited. The audit imports no
producer Python and uses only the copied, manifest-closed input packet.

## Independently established facts

The Paley-7 arc set and 189 Hamiltonian paths were reconstructed from vertex
permutations. All five- and six-arc deletion sets were tested against these
paths, yielding no five-blocker and exactly the same 63 six-blockers. Thus the
input census was not merely trusted. The phase-10 code, native endpoint parity,
and all 2,646 blocker/action combinations for the 42-element extended group
were reconstructed independently.

The 21 decoder values are distinct and are exactly the native boundaries:
9c maps to V minus c; 9c+7 to the closed incoming neighborhood of c; 9c-7 to
the closed outgoing neighborhood. All 3,969 ordered blocker pairs have equal
boundary if and only if their exponents are congruent modulo 21. A separate
phase-1 reconstruction fails this equivalence, as required.

The field checker uses carryless polynomial convolution and polynomial long
division, not the producer's shift-and-reduce multiplication. With modulus
X^6+X+1, the powers of X enumerate all 63 nonzero elements. The fixed set of
fourth power is the four-element subfield. Every native boundary fiber is
exactly a multiplicative coset of its three-element unit group.

All 63 nonzero derivative directions have the stated 0/2/4 distribution:
33 empty outputs, 30 outputs with two preimages, one with four. The exceptional
output is a^-1 and its preimage is exactly aF4. This is an affine four-element
set, whose three nonzero vectors represent ONE projective point. It is not a
projective line. The audit independently constructed the 21 projective lines
both as two-generator F4 spans and as trace kernels; the constructions agree.
Each line has five points, and none maps to a line under point-only inversion.

The trace relation gives 42 typed base objects and 105 flags, degree five in
both types. All 18,522 base pair/action combinations preserve incidence, with
the 21 type-reversing operations correctly swapping the endpoints. At blocker
level all 39,690 flag/action combinations also preserve incidence. The lifted
structure has 126 typed objects, 945 edges, and degree 15 in both sheets. Equal
neighborhoods in EACH sheet are exactly the three-element native boundary
fibers. This is a vertex blow-up, not a locally bijective graph covering.

The absolute point labels are 3,6,7,12,14. The explicit pairs (63,28679) and
(63,455) have three shared native arcs, but respectively are and are not flags
of the new relation. Thus native intersection count does not determine this
transported incidence.

## Proof and typing review

The centered families all have boundary V minus their center and exponent
shifts by 21. For a source pair with difference d=2^k, its closed-in center
is u+5d, and 10*16^k = 3d+7 modulo 21 gives the claimed decoder. The analogous
sink formula gives the closed-out center u+3d. This proves the commuting square
beyond the exhaustive test.

The projectivization is a genuine THREE-TO-ONE map F64* to F64*/F4*, not a
bijection from 63 blockers to 21 points. F64 has dimension three over F4, so
these cosets are exactly points of PG(2,4). Its lines correspond instead to
two-dimensional vector subspaces.

For relative trace T(z)=z+z^4+z^16, T(1)=1. Hence T(ax) is a nonzero functional
for every nonzero a, by evaluating at x=a^-1. Nondegeneracy identifies scalar
classes of coefficients with all projective lines. A point label e denotes
[g^e]; a line label l denotes kernel T(g^-l*x). The relation T(g^(e-l))=0 is
well-defined modulo 21 by F4-linearity. Simultaneous affine label action is
induced by x -> g^(9b)*x^(16^k), which is F4-linear. Its induced line-coefficient
action is a -> a^(16^k)/g^(9b), exactly the claimed line-label formula.

The typed converse (e,point)->(-e,line) is the trace polarity: the target line
coefficient is g^e. Its square is identity, and exchanging the two types makes
the new incidence exponent unchanged. This is distinct from the invalid
point-only inversion collineation claim. Absolute points satisfy T(x^2)=0,
equivalently T(x)=0, explaining the five-point absolute line.

For the derivative, outside x=0,a the substitution y=x/a gives y^2+y=1/(ab).
The absolute-trace criterion yields the two-root counts. At b=a^-1 the two
nontrivial F4 roots join 0 and 1, giving aF4. No projective-line inference follows.

Scope is strictly this minimum-blocker universe at order seven and the newly
defined typed relation. There is no all-size, all-order, Hamiltonian-count,
native-XOR, native-distance, or native-arc-incidence transfer.

## Reproduction and hostile controls

Run python -B replay.py from any directory. It checks the exact closed audit
manifest, runs independent audit.py in normal and optimized modes with identical
saved outputs, replays the pinned producer packet, and rejects added, missing,
mutated, and changed-manifest files in disposable copies. Mathematical hostile
controls reject phase 1, bare-inversion collinearity, and native-intersection
identification. No external repository, installed package, or network is needed.
