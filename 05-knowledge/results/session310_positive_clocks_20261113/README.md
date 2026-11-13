# Session 310: one observer, two structures

This session follows Session 309's exact-refill work and the nearby odd-cycle, Lonely Runner, Sierpiński and contextual-cover threads. It distinguishes a physical reserve, the polynomial that can certify a neutral move, and the arithmetic guards needed for an exact target. **Chan–Pak remains open.**

## A neutral refill that the entire count polynomial cannot see

The [full-profile pair](NEUTRAL_TAILS.md) gives a strong eight-vertex tournament B and a strong seven-vertex tournament A with exactly Q_B(t)=Q_(A→TT1)(t), where Q counts labeled path slots. Substituting B by A→TT1 preserves the Hamilton count in **every uniform exterior**, not only contexts with two visits. It exposes a real detachable trailing vertex. For r ordinal copies the same equality exposes r sites at once, starting from none, with no change in the dyadic valuation. In the Π clock, all r resulting single-arc transfers are literal. This particular clock family is divisible by 5. Separate exact exteriors exhibit the same neutral rewrite at the prime 2647 and at 231161=7·33023 (plus a larger seven-prime example); we do not claim Π transfers in those exteriors.

[The cycle trade](ODD_CYCLE_TRADE.md) explains the identity. The strong core trades six disjoint triangle pairs for twelve additional pentagons, with triangles and heptagons unchanged. Both kinds cost four powers of t, with respective weights four and two. This coincidence scales: two globally strong tournament families with the same full slot polynomial have triangle-packing ranks 2r+1 and r+1. A full contextual equality can hide an unbounded structural difference.

## A weighted hypergraph, not just a packing number

[The contact polynomial](CONTACT_POLYNOMIAL.md) groups odd-cycle packings by the connected sets of blocks they touch. A mixed pair of triangle blocks costs at least 8 above its independent baseline 9; a pairwise-uniform cyclic triple costs 3132 above its baseline 27. Rather than selecting just a largest matching of interacting blocks, the exact expansion charges all disjoint matchings with their weights. For a pure triangle blowup of a quotient Q, this also gives H≥3^k H(Q), with equality precisely for transitive Q. These are proved structural filters, not an integer coverage argument.

## The same expansion exposes determinant low bits

[The positive clock expansion](POSITIVE_CLOCK.md) explains why the slot polynomial cannot distinguish the pentagon trade, and why its low-degree terms still constrain exact refills. The minimum power counts each singleton or odd strong component once, and each even strong component twice. Independently, the known polarized determinant exits obey Y+e=12P_M(2)(3P_N(2)+2H(N)). The positive coefficients and exact binary valuations force all e=±1 exits into ±1 mod96, including the singleton boundaries. This improves the earlier mod-12 obstruction. A compiler must therefore use additional last-step families rather than merely find more witnesses for this determinant.

The next construction needs a physical return that selects the target, transfers the newly exposed sites in its own exterior, and replaces the spent module with another ready one. Unbounded one-shot capacity, all-context neutral substitution, and exact prime or seven-prime fixtures are now available, but those are not yet such a return.

## Replay

Run `python verify_neutral_tails.py`, `python verify_seven_splice.py`, `python verify_trade.py`, `python verify_contact.py`, and `python verify_positive.py` beside this file. All are standard-library verifiers tested in normal and optimized mode.
