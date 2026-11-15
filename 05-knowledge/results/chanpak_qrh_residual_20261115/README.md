# Chan–Pak: keep the unit repair, contract the placement error

**Uniform result:** the same physical repair now gives a complete
three-odd window on each side of every sufficiently large target N,
all within **O(N^(1/3))**, improving O(N^(3/8)) for this construction.
For **every fixed odd modulus m**, both windows can also have the
same actual address sum and a designated member exactly congruent
to N modulo m, at error O_m(N^(1/3)). The windows themselves still
have both exact +2 steps. No Riemann hypothesis is assumed.

## The constructive steps

1. [Complementary addresses](COMPLEMENT.md): condition on x+y and
   complement the two low blocks. Their leading checksum terms
   add to a constant, while x−y can move in unit steps. Two sign
   bisections, which require no monotonicity of the checksum, give
   literal bracketing tournaments.
2. [Let N choose the alignment](ADAPTIVE.md): using the actual
   midpoint overshoot instead of its worst-case bound gives
   **O(N^(1/4))** on explicit unbounded target intervals of length
   Θ(sqrt(N)). The [zero rail](RAIL_ZERO.md) adds separate
   aligned intervals where the target height ends in zeros.
3. [Program the suffix without changing the address sum](ALL_ODD.md):
   a finite-state complementary transfer adjusts the checksum and
   residue, preserving the physical center and the actual sum. The
   3- and 5-parts are handled by choosing the center and, for 5,
   one extra constraint on the sum. The shorter guard compiler
   certifies 13, 12, and 13 bits for moduli 7, 11, and 13.

## What old sessions contributed

The [October 17 quasi-Riemann/ROOT session](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/preprint_connections_solfeggio_digits_20261017.md)
read the [October 5 quasi-Riemann manuscript](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex)
as a distinction between a family estimate and a genuine transfer
that makes the remaining scale smaller. Here exact binary
complementation completes a family and then realizes its
cancellation in a **single** tournament, as opposed to averaging
two proposed counts. The older Collatz complementary-routing
work supplied the related requirement that the composed choice
remain legal. These are design transfers, not implications of a
zero-free half-plane.

The [family moments](FAMILY_MOMENTS.md) make the distinction
quantitative: over the whole complementary phase orbit, the
dominant linear variance is real. Pairwise averaging removes it
but does not generally produce a tournament. The [Fourier and
finite-state transfer](QRH_TRANSFER.md) separates coarse modular
coverage, exact modular steering, and absolute placement.

## What is still missing

The [capacity proof](CAPACITY.md) shows that our uniform lane has
Θ(X^(2/3)) outputs at scale X. Radius X^(1/3) is the right power for
this lane, not just a loose estimate. Fixed-modulus steering closes
a local congruence problem inside it; it cannot fill its absolute
gaps. To reduce the uniform radius further or hit every target
exactly, we must add independent physically selectable states
outside that lane. Chan–Pak remains open.

The proof and replay files are colocated. For a concise replay use
`python complement.py`, `python moments.py`,
`python verify_short.py`, and `python full_odd.py`.
