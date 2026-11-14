# Session 332: local coverage, hidden digits, and exact tournament returns

Odd counts have two natural addresses, \(H=2i+1=2m-1\), with \(m=i+1\). This session treats three different things separately: reindexing, observing a residue, and supplying an absolute count with a ready interface.

## The proof waves

| Prediction tested | Result and what it suggested next |
|---|---|
| The two odd-address origins might just relabel the gaps. | They do, but there is also an exact old arithmetic cycle underneath: the minus-sheet greedy inverse cycle \(\{4,11\}\) becomes \(\{7,21\}\). Forward closure fails because attained \(15,31,63\) all map to 21. Test actual backward moves instead. |
| A simple one-apex growth run should supply the next \(3\cdot2^s+1\) ray. | It cannot. A different marked core of count 13 supplies the entire ray by transitive inflation. A second family now has all four backward branches physically supplied for arbitrarily large lengths, and another continuation exposes the next missing port. |
| The visible 3 bits and hidden 5 bits might be independent groups. | They are a cyclic kernel of order 32 over an eight-element quotient in the actual multiplicative address; the exact sequence does not split. A prime-register odometer can keep the carry separately. |
| Fixed-size products using 3 and 5 should already fill the octet. | They reach one of two cosets. A same-size exchange between multiplier 17 and \(3\cdot5\) supplies the second. |
| Once every phase fits an added size, the next size should also work. | False. Every even-size full condition is exactly the odd-size condition three below. Modulo 512 the first added size is 159, size 160 has a single gap, and 161 onward is full. |
| The large-size cost might stay comparable to the modulus. | A lattice and a forced change of shortest direction prove infinitely many precisions with optimal square-root cost, plus a general \(O(2^n/n)\) conductor. |
| Complete binary and seven-power views should combine freely. | False. Quadratic reciprocity prescribes \((7/H)\) by the added order. One equal-size substitution by \(R_5\) repairs it; a rank calculation adds an independent axis at 13 via \(R_7\). |
| Once every physical register transition is reproduced, PRIMEGAME can output a tournament of the prime count. | Its unary state remembers the exponent in \(Q'(1)\), not in its Hamiltonian count. A uniform decoder either assumes the desired quotient or overshoots exponentially. The smallest counterfeit seed 50 shows why initialization itself matters. |
| A primitive four-state should always admit a primitive move, and an unrestricted protected two-port language should then give many absolute targets. | A strong infinite family of primitive states has no such next move; reached primitive states do have protected old/first choices, but all genuine lengths in that protected grammar still give sublinear absolute coverage. |
| Two different local move orders might truly return to the same future interface. | Yes, uniquely in the elementary old/first comparison at the shortest block. Iteration gives \(2^m\) nonisomorphic strong tournaments at one exact \(H_m\). A two-slot observer opens the final collision, and this sparse ray alone covers every scalar \(7\)-adic residue. Thus modular surjectivity and even huge exact multiplicity still do not imply an absolute cutoff. |

## Results and replays

- [Odd-coordinate predecessors](ODD_COORDINATE_PREDECESSORS.md): transported rings versus ordinal multiplication, the exact cycle and full predecessor menu, three-bit erasure with five-bit retention, the four-way physical slice and the next blocked one-apex ray.
- [Reciprocal cycle rigidity](RECIPROCAL_CYCLE_RIGIDITY.md): an elementary Vieta descent excludes every even exponent in the generalized multiplicative/affine loop; a Pell argument fixes the multiplier 3 uniquely at head 7.
- [Fixed-size Galois phases](FIXED_SIZE_GALOIS.md): the literal construction behind the 3+5 split, exact size holes, all-order size estimates, and the independent quadratic phase repairs.
- [Three-fraction odometer](EIGHT_BIT_ODOMETER.md) and [PRIMEGAME return analysis](PRIMEGAME_RETURN.md): a literal copy of the natural numbers with 3+5 carry bits, exact false-positive and looping initialization, and the divisor-dependent clock. [Tournament registers and their decode wall](TOURNAMENT_REGISTER_WALL.md) gives the physical program lift but also its precise limitation.
- [Commuting diamonds](COMMUTING_DIAMOND.md) and [primitive guards](PRIMITIVE_GUARDS_AND_SCARCITY.md): exact same-state returns with geometric differences, quadratic/valuation conditions for allowed ports, and absolute limitations even after no gcd guard remains.

The next proof obligation is still an absolute inverse: given an odd target, find it in a physical class with the required next state, or prove that a finite set of exact inverse branches reduces every remaining target. These investigations supply new candidate branches, remove several impossible approaches, and identify distinct phase and size obstructions. They do not prove Chan–Pak, ordinary Collatz, or undecidability of either.
