# Chan–Pak: put the target between two physical counts, then control their phase

[Chan–Pak, §6](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf)
conjectures that every positive odd integer except 7 and 21 is a
Hamiltonian path count of a tournament. **Exact universal placement
remains open.** This packet makes several all-order advances on the
absolute error and tests why a change of contact or alphabet is needed.

1. **[Paired binary tapes](TWO_TAPE_STEERING.md):** from a single
   fixed cyclic carrier, construct strong tournaments whose counts
   bracket every sufficiently large target with error O(N^(3/8)).
   One tape had a first-order checksum problem. Two tapes move along
   a tangent to a product, with their common low bits holding the
   largest digit-dependent term fixed. Freeing the digits but retaining
   their common carry gives infinitely many bands of width of order
   sqrt(N) with O(N^(1/4)) mesh. The band centers are not universal.
2. **[Modular programming](PREFIX_PROGRAMMING.md) and
   [complementary suffixes](PAIR_PROGRAMMING.md):** in any aligned
   binary cylinder one tape can program its address and surcharge
   independently at any fixed modulus. Programming two tapes
   separately destroys the tangent. Complementary suffixes instead
   preserve their exact sum, giving
   **[physical O_m(N^(2/5)) bracketing](CONGRUENT_PAIRED_STEERING.md)**
   at every fixed modulus M, with both counts congruent to the
   prescribed odd target mod 2M. Dyadic precision can grow with N
   without changing the exponent by more than o(1). The exponent is worse than in the
   unconstrained bracket; rounding to an admissible high-sum cylinder
   costs real absolute error. Fixed-modulus control does not permit
   silently increasing the modulus beyond the error.
3. **[Moving the exposed vertex](NESTED_RECONTACT.md):** if the
   same baseline is kept, the [exposed charge](FIXED_BASELINE_CHARGE.md)
   has only finitely many dyadic rays. Moving the contact introduces
   a genuine unimodular operation, not another reset, but the
   three-letter nested grammar remains sparse on the prime and 7p
   lanes. Directed triangle blocks break the reset. Adding only
   either one of the two new cyclic types still gives an all-order
   sublinear bound on prime and 7p outputs. Using the appropriate
   cyclic/acyclic returns together gives a free ternary tree with an
   exact inverse for the **pair** of stored coordinates, not for the
   prescribed scalar. The literal two-contact counterexample shows
   what to retain next: ordered paths between the two contacts, not
   just their individual start data.

The common principle is **preserve the quantity that controls the
next error, not every coordinate of the current construction**.
Freezing both low digit tapes handicaps placement; preserving their
sum or carry allows useful compensation. Conversely, erasing the
joint contact readout makes indistinguishable states produce
different outcomes. A next exact proof obligation is to construct a
physical repair comparable to the remaining absolute debt, or use a
target-directed carry selector whose cost stays bounded when the
precision grows. Neither a rounded real target nor a compatible
residue is an exact witness.

## Replay

Run `replay_steering.py`, `replay_prefix.py`, `replay_pair.py`,
`replay_even_pair.py`, `replay_lift.py` and `replay_cyclic_pressure.py` with Python 3.
Each has a corresponding text receipt; standard and optimized
runs agree. Literal small controls count paths by subset-endpoint
DP. Large-target integer controls check the specified decoder; the
infinite statements follow from the proofs.
