# MRB, Leech glue, and the five-corridor Chan–Pak clock

**Chan–Pak remains open.** The strongest conclusion of this session is
a local/global separation for a single physical family, not a larger
bounded census. The historical order \(4,2,5,3,7,13\) suggests an
active five-step clock and a separate exact reset at 13, not a
six-element arithmetic identification.

- **[The same clock matches every finite congruence](FULL_CONGRUENCE_CLOCK.md).**
  Repeating \(4,2,5,3,7\), using only positive widths, realizes every
  compatible residue of an odd integer modulo *any* modulus. The odd
  part uses rank-one idempotents at 3 and 7 and determinant-one
  groups elsewhere ([proof](ODD_MONODROMY.md)); the binary part uses
  exact bit toggles on a preserved line. These pieces work
  simultaneously, not just on separate witnesses.
- **[Absolute size is the missing variable](CLOCK_HEIGHT.md).**
  The same full clock has only \(o(X^{0.406})\) distinct outputs
  up to \(X\). Almost every residue modulo \(m\) needs a clock total
  above \(m^2\); the largest least height is
  \(\omega(m^{500/203})\). An [independent phase-pressure proof](CLOCK_SPARSITY.md)
  also gives an ordered-word bracket, a physical lower-bound
  injection, and a scoped barrier to combining two such reservoirs.
- **[MRB and Leech sources](MRB_LEECH_SOURCES.md).**
  Revisited the two Burns papers supplied in Session 31 and the
  contour, phase-filter, and Niemeier work. The six labels are all
  genuine pure-\(A\) Coxeter numbers, but no Leech cause of the
  reset follows. The actual rank-two analogue is an integral glue
  condition modulo 5345 at the 13 seam. Its transverse coordinate,
  not the norm alone, distinguishes physical repair candidates.
- **[A reset turned into an exact repair](TRANSVERSE_REPAIR.md).**
  Two explicit terrace prefixes agree after length 13. Their hidden
  difference becomes exactly \(-2\) or \(+2\) after length 10 or 14,
  respectively, for every terminal width. These are unbounded,
  absolute-height repair rails. Choosing the right rail for an
  arbitrary target remains open.
- **[MRB-style residue rounding](RESIDUE_ROUNDING.md).**
  All three Tribonacci pole residues have rank one. Given the two
  reservoir states, the leading term pins down at most one late
  corridor length, with a uniform analytic error bound and a final
  exact cubic-field and integer check.
- **[Unbounded nests](NESTED_WRAPS.md)** and
  **[crossing ports](PARALLEL_GLUE.md).**
  Outer wraps plus ordinal sums are still sparse even at unlimited
  nesting, with wrapped atom counts \(o(X^{0.928})\). The crossing
  construction gives a genuine new three-channel readout and an
  explicit \(7p\) witness. Even arbitrarily wide, order-preserving
  matched-port gluing at unlimited depth, interspersed with ordinal
  sums, has only \(o(X^{0.952})\) distinct outputs. Its matched-port
  scope matters. The [free-matching extension](FREE_MATCHING.md) proves
  the same exponent for arbitrary matching permutations. It is not
  an obstruction for general changing interactions among reservoirs.

**Proof frontier.** By the inherited ordinal-product reduction it is
enough to construct strong witnesses on odd primes other than 7 and
on \(7p\) for primes \(p\ge5,\ p\ne7\). Congruence completeness cannot
supply those witnesses: we now have it in a family provably much too
sparse. Nor is another layer of symmetric wrappers enough. The
transverse pair suggests a more concrete test: find a source pair
with provably reusable \(\pm2\) (or another exact additive correction)
*and* a target-selected common baseline, while keeping sufficient
ports to return to the same readiness state. A rank, a matching
congruence, or a small relative error is not that baseline theorem.

Run the Python replays next to their notes. In particular,
`python clock_replay.py` and `python replay_full_clock.py` verify the
finite pieces of the local theorem, `python replay_height.py` and
`python clock_sparsity_replay.py` verify the rational pressure
bounds, `python replay_transverse.py` checks the exact signed rails,
and `python cross_replay.py` checks the literal crossing counts.
