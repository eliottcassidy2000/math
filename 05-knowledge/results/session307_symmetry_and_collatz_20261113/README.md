# Session 307: symmetry budgets under growth

**Result:** proved scoped size-growth obstructions in both problems, not a proof of Chan–Pak, ordinary Collatz, or either of its two open halves.

## Tournaments: global symmetry has a measurable price

[The orbit fold](ORBIT_FOLD.md) gives an exact characteristic-p factorization of the falling-factorial path-cover polynomial for a prime-power automorphism. Each vertex-orbit length has its own weighted directed quotient and its own Frobenius exponent. Only symmetry-invariant covers survive reduction; a directed path in an invariant cover must use vertices from a single length tier, once per orbit, so the quotient weights count genuine phase choices. Reverse-path cancellation makes the quotient skew and forces the parity of all surviving cover lengths. This measures what remains visible as the symmetry tiers change, without pretending that modular coefficients determine the actual path count.

[The integer symmetry budget](SYMMETRY_BUDGET.md) complements this: choose a directed path in the majority tournament on the orbits, then count block-contiguous fillings from the unit increments of each cyclic action. This bounds the number and lengths of occupied symmetry tiers at fixed Hamiltonian count, even if arbitrarily many fixed vertices are present. In particular, every tournament whose Hamiltonian count is a prime at least 5 has trivial automorphism group. A proof of Chan–Pak by symmetrical intermediate carriers must eventually dispose of their global automorphisms on every prime target; arithmetic divisibility alone did not show this.

## Collatz cycles: a finite motif dictionary cannot hide long cycles in few seams

[The cycle motif theorem](FIXED_MOTIFS.md) proves that a primitive positive Collatz cycle of length N built from repetitions of any fixed finite set of parity motifs needs Ω(log N) block seams. A cyclic cut bounds the actual block-boundary value polynomially in N using a two-logarithm bound. Exact 2-adic divisibility then limits each repetition to the logarithmic precision available at its entry. The result handles independently varying repetition counts, not just one-parameter rays. It does not rule out a cycle that uses many seams or a growing dictionary.

## Collatz divergence: necessary constraints on actual tails

[The descent horizon](DESCENT_HORIZON.md) makes the first contracting affine parity threshold polynomial in its length, rather than using an exponential denominator-only estimate. Thus an actual future minimum y of a hypothetical divergent orbit has a growing interval of slope noncontraction. [The factor-complexity note](PARITY_COMPLEXITY.md) bounds how many small integer orbit values a parity factor language can encode; in particular, zero factor entropy is incompatible with a global polynomial upper growth envelope for a nonrepeating orbit. These constraints are on genuine integer orbits, not on the existence of an arbitrary compatible 2-adic sequence.

The common mechanism is precise but not a transfer of proofs: a repeated unit is cheap only when its interface information can be reused. In the tournament the information is the phase-weighted orbit quotient and its integer count cost; in Collatz it is the signed affine slope and 2-adic precision at each seam. The remaining constructive Chan–Pak gap is a symmetry-breaking return that stays ready for the next target; the remaining Collatz gap is excluding orbits outside the bounded-complexity regimes.

## Replay

Run `python verify_orbit.py` or `python -O verify_orbit.py` next to this file. The standard-library exact subset/cover replay checks 95 prime-order controls and six mixed prime-power controls; both modes agree. The cycle and divergence notes give direct proofs rather than finite cycle censuses.
