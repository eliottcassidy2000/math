# The 67/69 pair: from a single arc to a portable skew selector

**Proved physical families, exact obstructions, and finite exact controls.**
Chan–Pak and universal placement remain open. The 67 and 69 in this
packet are literal Hamiltonian-path counts, not removed residue classes
or a claim that 69's factorization controls its neighboring prime.

## The pair is more structural than its original presentation

- **[Painted domino, changed exteriors, and the absolute bound](PAINTED_DOMINO.md).**
  The exterior triangle reversal becomes one arc reversal after a
  relabeling. Contracting the used edge into an ordered domino isolates
  an exact skew of a two-path endpoint matrix. Two full binary contact
  tapes can vary arbitrarily with exact gap two; e.g. five independently
  variable runs start at 1489/1491. The same proof reaches 89/91 and
  2325/2327 with different exterior cores, not merely a rescaling of
  the original triangle. Yet any **finite** bank of fixed exteriors
  gives only O(X^log_3(2) (1+log X)) such counts up to X. More free
  run lengths are not enough to obtain cofiniteness.
  [Literal replay](pair_structural.py).
- **[Domino audit and the original module boundary](TRIANGLE_6769.md).**
  Direct structural cancellations explain why the original uninflated
  third vertex is different from a newly inserted tape. Its independent
  [replay](triangle_6769_controls.py) checks nontransitive inputs too.
  The [sharp third-module formula](THIRD_MODULE.md) classifies all
  exact gap-two cases at the old three sites. Inflating mark 6 to a
  pair changes 67/69 to 135/133, reaching 7·19 with the opposite sign.
  The same note gives the exact ternary parity law for the two
  original marks and the converse asymmetry that explains a port alias.
- **[Exterior growth boundaries](EXTERIOR_GROWTH.md).**
  Keeping an automorphism does not keep the skew. Inflating the fixed
  point of the seven-vertex exterior changes the bracket to
  2H+17C2+12C3, proved for arbitrary nonempty modules.
  Another exact finite test rules out all 252 strong equivariant
  one-orbit extensions of this particular exterior as unit selectors.
  [Fixed-point replay](exterior_fixed_point.py) and
  [extension census](equivariant_extension.cpp).

## Port and orbit combinations

- **[Two independent modules and their roots](TWO_SCALE_AND_PORTS.md).**
  The original 67/69 first two marked modules independently expand
  to any positive transitive sizes. The exact tensor also gives a
  prescribed-target decoder. At the last clone of mark 1, one
  nontrivial seam transports the unit gap; iterating that strategy
  cannot keep increasing the baseline while preserving it. Two
  companion window banks repair additional recorded eleven-block
  omissions by literal constructions. [Replay](replay_two_scale.py).
- **[A separate port-only switch](PORT_SWITCH.md).**
  The same order-five reader with a changed entry and exit gives
  72Q+67h vs 72Q+69h on any nonempty source module; two more
  reader pairs give slopes 73 and 74. The free module can be transitive
  for exact gap two, or cyclic for explicit mixed powers of 2 and 3.
  The same packet certifies 67 outside the unrestricted O/Z/J
  grammar and 69 inside the eleven-block grammar. [Replay](port_switch.py).

The cleanest new obligation is to choose the exterior and its
automorphism-related masks **from the prescribed target**, keeping
their skew bracket 2 while changing the baseline. The fixed
exterior sparsity theorem and the two failed growth steps give
concrete tests of that obligation, rather than another numerical
spectrum extension.
