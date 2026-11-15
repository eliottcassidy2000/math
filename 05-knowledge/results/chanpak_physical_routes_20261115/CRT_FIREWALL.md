# Even separate prime-power universality does not preserve a joint route

Here is an exact counterpart to Artin-style entanglement, stronger
than the coprime-to-6 example in the analogy note. Start at y=1 and
allow only the positive integral moves

    f(y)=2y,     g(y)=3(y+1).

Every residue modulo **each prime power** is reachable from the
same seed. Nevertheless residue 5 modulo 6 is not reachable.

* For p≠2,3 the commutator displacement is −3, a unit modulo p^k.
  The two maps are permutations, so inverse letters there are finite
  positive powers. The commutator translation fills the entire ring.
* Modulo 3^k, f alone fills the units. The order of 2 is 2·3^(k−1):
  v₃(2²−1)=1 and each cubing raises that valuation by exactly one.
  To reach a target z divisible by 3, induction at 3^(k−1) first
  reaches z/3−1 and one final g reaches z. The induction starts
  at the trivial ring.
* Modulo 2^k, the coordinate w=2y+3 identifies the residue ring
  with the odd residues modulo 2^(k+1), and g becomes multiplication
  by 3 in that coordinate. The two seeds f(1)=2 and f²(1)=4 give
  w=7 and 11, in different cosets of the 3-generated subgroup
  (the subgroup is precisely the classes 1,3 modulo 8 for k≥2:
  the order is 2^(k−1), obtained by squaring from v₂(3²−1)=3).
  Their orbits fill all odd w. For k=1, g already toggles parity.
* Modulo 6, any word containing g lands at a multiple of 3 and
  never leaves that class. The f-only words give 1,2,4. Thus 5
  is missing. This is a terminal-language obstruction, not
  failure of any individual prime-power test.

Even an appropriately chosen positive semigroup with complete local
orbits at each prime can fail to lift simultaneous local witnesses
to a single word. The Chinese remainder theorem assembles coordinate
targets, not admissible words that visit them together. This is not
a tournament construction; it explains why the proposed physical
route must preserve the **joint** marked interface, not merely reach
every prime-power projection. The absolute output set of this same
example is power sparse, by the 2,3 weighted stopping-frontier bound
in the analogy note.
