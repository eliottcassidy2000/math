# Sources and boundaries

Checked 2026-11-04. No forbidden ledger or startup helper was used.

## Primary mathematical sources checked on the web

1. Pavel Etingof, *Groups Around Us*, Theorem10.2, finite subgroups of SO(3): https://math.mit.edu/~etingof/groups.pdf . Used for the classification only.
2. Donu Arapura, *Algebra*, Chapter14, Theorem14.1, finite subgroups of SO(2) and O(2), including D1 and D2 conventions: https://www.math.purdue.edu/~arapura/algebra/algebra14.pdf . Used for the planar classification only.
3. Etingof et al., *Introduction to Representation Theory*, Section4.1, Theorem4.1.1 and Proposition4.1.2: https://math.mit.edu/~etingof/representationtheorybook.pdf . Used for Maschke and its converse, not for the original computations in this packet.
4. Ethan Akin, *Rock, Paper, Scissors, Etc: Topics in the Theory of Regular Tournaments*, v4 (2020): https://arxiv.org/html/1806.11241v4 . The inherited THM-4644 attributes the odd-multiplier and Fermat rigidity mechanisms to Theorems3.13-3.14 plus the prime-degree affine theorem. No fresh priority claim is made.
5. PrimePages, the five known Fermat primes: https://t5k.org/lists/top_ten/topten.pdf . Cross-check with university source https://sites.lsa.umich.edu/hochster/wp-content/uploads/sites/1337/2024/11/fib20B.pdf . No claim that there can be no further Fermat primes.
6. Czech Academy of Sciences, Gauss-Wantzel theorem and original references: https://www.cs.cas.cz/portal/AlgoMath/Geometry/PlaneGeometry/GeometricConstructions/GaussWantzelTheorem.htm . Used only for the classical constructible-polygon criterion.

## Inherited repository facts

Full filenames and links are listed in `PROOF.md` Section0. Read locally: LEM-003, THM-3447, THM-4644, THM-4553, THM-1370; the Gilbreath/Fermat/Platonic synthesis; and the Session17 metagraph trijections/reset synthesis. The path-complement operation follows that labelled definition. The present proof does not inherit any aggregate-to-original lifting claim.

The proof is self-contained for edge stabilizers, edge-orbit formulas, the affine stabilizer diagram, the sharp real-dimension lower bound, and the binary construction. Finite facts about Paley endpoints and signed path orbits are certified by `verify.py`. The SO(3) classification and the complete prime-degree circulant theorem are explicitly external/inherited inputs.

## Connection contracts

- Platonic solid -> rotation group: keep rotations, forget vertex/face distinction; restore the duality marking. Five solids map to three exceptional groups.
- Directed pair -> unordered pair: quotient by reversal C2; forget orientation; restore a half-set gauge. AGL(1,7) has stabilizers1 andC2 respectively.
- Hamilton path -> metagraph neighbour: complement all off-path arcs; map is injective and Aut(T)-equivariant; no label loss. Forgetting source or labels requires separate orbit bookkeeping.
- Ordinary edge -> signed edge orbit: quotient by global converse; retain the reflection-fixed sidecar to distinguish orbit sizes and residual C2.
- Real Paley permutation module -> binary permutation module: retain the group action, change the coefficient field; do not reduce the skew signs and pretend orientation survived.
- Prime p -> circulant multiplier groups: use the odd part of p-1; a pure2-power forces trivial multipliers. The number of presently known Fermat primes plays no role in the proof.
