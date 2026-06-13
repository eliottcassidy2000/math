# Cauldron Complexity And Variant Atlas (S619)

The clean answer to the user's difficulty question is: base fixed-`k`
last-boil is hard, but not hard in the same way as A000568, LRC, or unit
distance.

It is hard because it is a Schur-frontier problem. Even before last-boil
removal, the weak two-term `k=4` first-boil frontier crosses `250000` canonical
safe states at `n=13`. The exact `k=3` last-boil search is already `184948`
states for the weak rule. That is real computational pressure.

But it is not A000568-hard in the base form. A000568 is an orbit/Burnside
problem over all tournaments on `n` vertices. Base cauldrons quotient only by
interchangeable color labels inside a one-dimensional prefix process. No
tournament isomorphism class has to be retained unless we add a hidden-label or
observer quotient on purpose.

It is also not yet LRC-hard or unit-distance-hard. LRC has a continuous
all-orders coverage object, the `p_0` coimage, and finite moments cannot decide
the hard wall. Unit distance has geometric embedding and totally-unfaithful
side channels. Base cauldrons are purely additive color-frontier states.

The good news is that cauldrons are a small laboratory for importing those
missing structures one at a time.

The first imported structure that actually bites is adversarial parity. In the
attack-only game, A places odd numbers into B's cauldrons and B places even
numbers into A's. Under the two-term weak rule, B wins every small case
checked: `1v1` and `1v2` at `n=6`, `2v1` and `2v2` at `n=14`. The reason is
not raw density but residue correlation: odd + odd is even and cannot make an
odd target, while even + even can make an even target. When finite or
three-term sums are allowed, A can use `1+3+5=9`, flipping the `2v1` case.
This is the same kind of lesson as Collatz/two-block and the Vitali wall: the
low-order mass is misleading because the residue stream is correlated with turn
order.

The best creative variants from S619:

- attack-only adversarial parity cauldrons: cleanest finite game and already
  shows correlated residues;
- gift-or-poison cauldrons: players may place `n` into their own or the
  opponent's cauldrons, adding genuine strategic choice;
- modular CRT cauldrons: boil on zero-sum/subset-sum in `Z/mZ`, directly
  importing CRT and Collatz residue channels;
- LRC depth cauldrons: a cauldron is a family of danger arcs and boils when
  `p_0=0`;
- unit-distance geometric cauldrons: place points into cauldrons and boil on
  unit edges, triangles, or unfaithful embedding certificates;
- hidden-quotient A000568 cauldrons: observer sees only orbit type, forcing
  Burnside-style side data;
- OCF odd-cycle cauldrons: ingredients are odd cycles and boil values are
  independence-polynomial obstructions such as the forbidden `H=7,21` face.

My current taste: the adversarial parity variant is the one to push first. It
is small enough to solve exact cases, but it already carries the deep pattern:
when the process itself sorts numbers into correlated residue streams, density
is blind.

Artifacts: `04-computation/cauldron_complexity_variants_s619.py` and
`05-knowledge/results/cauldron_complexity_variants_s619.out`.
