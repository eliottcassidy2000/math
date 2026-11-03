# Constructive edge jumps and source-preserving recurrences

**PROVED scoped; FINITE-EXACT controls. Same-agent proof review with
independent counting and move implementations.**

A two-apex cancellation construction realizes every even Hamiltonian
edge-count jump 2d between strong tournaments on bitlength(d)+6 vertices.
A uniform program reverses that edge in two native moves using one temporary
vertex, retaining both physical apex ports.

The derivative has an exact six-coordinate integer recurrence. A unit
6 by 6 Hankel minor proves sharp linear dimension and supplies an integral
decoder for the normalized actual source words. Three infinite unit-jump
families have closed formulas and rational generating functions, including
an order16 pair with counts 305835 and 305837=7*43691.

The same derivative defines an integral alternating lattice. Its
source-augmented form peels into unit symplectic pairs. Removing the source
coordinate gives a complete cokernel classification, including every
(Z/d)^2, measurable by mixed differences of actual strong sources. Oddness
determines whether an affine center exists.

- [THEORY](THEORY.md): actual Boolean constructors and native move clocks.
- [CARRIER](CARRIER.md): six-coordinate recurrence, minimality and source decoder.
- [GEOMETRY](GEOMETRY.md): integral peel, paired torsion, half-grid centers,
  fifteen F2 planes and nonzero hourglass-style curl.
- [RECURRENCES](RECURRENCES.md): infinite unit jumps and the exact 7-adic phase.
- [RECOVERY](RECOVERY.md): inherited mechanisms, neglected ideas and typed bridges.
- [AUDIT](AUDIT.md): universes, independent implementations and inherited controls.
- [HANDOFF](HANDOFF.md): constructive next obligation.

These results control changes and retain sources. They do not yet match
every required starting H to a unit-jump source; the full tournament
Hamiltonian spectrum remains open. The six-coordinate carrier is closed for
common chain append and edge-jump observations, not arbitrary re-marking.

From the repository root, with Python and SymPy:
~~~sh
python -X utf8 -B 05-knowledge/results/apex_wedge_constructor_20261103/reproduce.py
~~~
The runner executes all four certificates normally and with -O, compares
their receipts, checks source pins, and writes reproduction.json.
