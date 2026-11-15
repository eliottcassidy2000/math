# A literal nonlocal escape: four consecutive odds from independent switches

The [two-site obstruction](TWO_SITE.md) is a scope theorem, not a global law for one apex. There is a literal **additive** escape with an 11-site transitive spine, a one-site control, and a five-site control. In the order from first site to last, use the 11-bit binary presentation of base integer 1158. Flip the bits in either or both disjoint xor masks 32 and 220:

| Controls | Spine integer | Hamiltonian path count |
|---|---:|---:|
| none | 1158 | 721 |
| 32 | 1190 | 723 |
| 220 | 1114 | 725 |
| both | 1146 | 727 |

All four are strong: the high spine bit is 1, the low spine bit is 0, and that marks a directed spanning cycle via the apex. The count can be checked combinatorially by the inherited recurrence (p(0)=0) and (p(2u)=p(u)+u+[u odd], p(2u+1)=p(u)), with count (H=1+p). Our verifier also checks all four by a separate last-vertex subset DP on the actual tournaments. This is an exact repair of width six above a fixed literal baseline, **not** a construction choosing that baseline from an arbitrary odd target.

With nonlinear interaction allowed an earlier example is a 9-site spine, base integer 232 and disjoint xor masks 458 and 16. In that order the four counts are 213,215,219,217. Trim neutral leading 0s and trailing 1s to obtain strong versions of each individual tournament; the shared switch construction itself uses the same nine spine sites.

## How to find the escape without searching larger construction trees

For each spine length, group binary words by exact count. Suppose a rectangle consists of four consecutive odd numbers. Take its **least-count** word as base (reorient each switch if necessary). Take candidate adjacent corners from any two of the other three count fibers. XOR each against the base; the two masks must be disjoint for independent controls. The fourth corner is then forced as the XOR of the base and both masks. Check its count and, for the additive version, check that the adjacent corners have increments 2 and 4. This is a target-side inverse test in three small fibers, not a forward search over every switch-pair. It is complete for this type of rectangle. Replay checks no shorter witness through 8 and 10 spine sites respectively, but minimum length is only a finite certificate, not a premise for the examples.

## An address-neutrality test for reusing this construction

The additive example cannot simply be inserted into an arbitrary surround. For each switch δ, with core positions numbered i=0,…,10, let P=Σ_i δ_i 2^i and Q=Σ_i δ_i 2^(10−i). On fixing the two exterior parities, the coefficient of the two exterior address quotients (U,V) in its half-step is (−P,Q). For mask 32 these numbers are (P=Q=32); for mask 220 they are (P=-56,Q=-44). The resulting two linear constraints have nonzero determinant, so at each exterior parity there is at most one possible pair of exterior quotients preserving both steps. In particular there is no unbounded fixed-core padding extension of this particular unit rectangle. The next check for candidate escapes is to look for collinear but physically positive exterior coefficients, or allow the marked switching pattern itself to change with the prescribed target.

This example exhibits a distinction from mere residue control: a nonlocal correlated switch can cancel all integer mixed counts and leave exact unit steps, but a further independent condition is needed to let those steps travel with a moving baseline.
