# Collatz addresses, changing tournament sources, and the exact target

**14 November 2026.** Eight scoped all-order results and six independent replay programs. Chan–Pak remains open. No computational spectrum record is claimed. The question here was whether a Collatz itinerary can be made into an **exact physical inverse**, rather than only a residue address.

| Direction | Result | Boundary |
|---|---|---|
| [Induced-source triangle tariff](TRIANGLE_TARIFF.md) | If the old source has r cyclic strongly connected components, the smallest possible positive path-count increase is at least 2·3^(r−1), sharp. Retaining k literal triangle components cannot implement a batched 3^k(H+1)−1 by **any** induced extension, apart from the single step 1→5. | Change the actual core, not merely the later contact mask; the argument cannot prohibit replacement. |
| [The leading trit](TRIANGLE_CLOCK.md) | For one apex over an ordinal triangle train, the first potential surviving base-3 digit is a sum of products of marks on **adjacent** triangles. With no abutting triangles, the stronger digit includes the two actual endpoint marks at each gap. | Generic trit balance does not cover 2 mod 3 after three triangles are retained. |
| [Abstract versus physical trits](ABSTRACT_TRIT.md) | Independent-set graphs implement the shifted ternary rule exactly, but two natural encodings violate distinct necessary cycle-conflict constraints. | Only the exact conflict-graph encodings are ruled out, not all equinumerous tournaments. |
| [Why 63 needs two chronological paths](READOUT_63.md) | No single ordinal time reading, even from arbitrary real numbers rather than Collatz, has H=63. A five-triangle certificate rules it out, and two specified Hamilton paths over one rank order suffice. The minimum path budget composes by maximum over strong components. | Two paths over one rank order do not yet give a two-orbit arithmetic decoder. |
| [Short Collatz clock compiler](WINDOW_COMPILER.md) | Given a generic ordinal permutation, a least-height staircase decides feasibility, and a backwards congruence compiler produces an actual positive seed on either 3x±1 sheet. For strong generic patterns, a seed needs at most (2 log₂3−1)(m−1)+O(1) bits. | The ordinal permutation still has to be chosen from H, and exceptional actual windows need not be generic. |
| [Exact binary suffix decoder](BINARY_SUFFIX_DECODER.md) | At a fixed marked source, the finite-word absolute inverse reduces to an exact address split and O(√Y) candidate tests for target displacement Y. This gives a harmonic word-capacity law. Every prescribed initial mark string still has the *same* residue image as unrestricted suffixes. | Even taking all source cut masks of one fixed cyclic core, this suffix family misses at least 10/21 of the odds in lower relative density. |
| [The sparse itinerary](PARITY_TAPE.md) | One increasing block schedule gives a 2-adic isometry onto all odd residues, but its typical address length is quadratic in the number of required bits; its exact absolute output count is asymptotic to a distinct-partition law in log₂X. | Coherent residue completeness is compatible with a subpolynomial exact spectrum. |
| [Finite exact slope menus](FINITE_SLOPE.md) | In transitive-bit attachment at fixed cycle count, a specified pair of slopes and length fixes the whole physical word. A finite slope menu transported as affine lines gives only O(log X) counts. | Not an obstruction to nonlinear, changing-core, or target-tested isolated transitions. |

## What seems worth trying next

Keep the smallest uncovered test **63** separate from the universal exact decoder. It tells us that a second chronological path is an actual new resource, as opposed to a longer valuation clock on the old resource. Try a two-path seam that reroutes a source while preserving the target equation, with the two paths read against the same rank. Start by looking for moves for which the triangular correction is carried by a local interaction **without** leaving an ordinal bank of independent triangles intact. At each branch require a physical witness, correct absolute output, and an actual decrease in an integer target or auxiliary rank. A modular hit or a local decrease relative to a later high point does not discharge any of those conditions.

The close analogy is with an exact Collatz itinerary: the parity word or valuation clock is good input information, but it does not by itself buy return below a prescribed baseline. The prior repo audit of that failed induction (Session31, `session31-collatz/AUDIT.md`) motivated the separation. For the original questions see [Chan–Pak, §6](https://doi.org/10.1016/J.TCS.2024.114776) and [Bernstein–Lagarias](https://websites.umich.edu/~lagarias/doc/bernstein.pdf). No implication between the two open conjectures is claimed.

## Replay

Run any of the following in this folder (normal or with `python3 -O`):

- `python3 replay_triangle_tariff.py` (literal extensions, subset path DP)
- `python3 replay_triangle_clock.py` (120 marking controls with gaps, subset path DP)
- `python3 replay_readout_63.py` (literal path and rank certificates)
- `python3 replay_windows.py` (both sheets and an independent finite generic alphabet)
- `python3 replay_binary_suffix.py` (binary split, actual words and fixed-core seeds)
- `python3 replay_parity_tape.py` (coherent inverse and literal small tapes)

Replays are finite controls on unbounded symbolic proofs. Historical priorities and inherited results are credited within the notes.
