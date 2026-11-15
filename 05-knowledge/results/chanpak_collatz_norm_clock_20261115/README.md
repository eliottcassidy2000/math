# Mod 8 as a gate, not an address: Collatz and the tournament spectrum

**Main constructive advance.** The same physical three-center bank can
now place literal +2,+2 windows near both sides of every sufficiently
large odd target, with one designated count in each window agreeing with
the target modulo **every fixed integer** (not just odd integers) and
within O(N^(1/3)) for that modulus. A new binary excursion construction
achieves agreement modulo 3·2^k with k ≥ c sqrt(log₂ N) while retaining
the O(N^(3/8)) absolute budget. This is expanding residue precision, not
zero residual or a proof of Chan–Pak.

## What the three mod-eight pictures have in common

- **[The table is an exact midpoint chart](TABLES.md).** The central four
  entries of the modulo 4x+ε multiplication table are x and 3x+ε twice
  each. The two independent gates are the sign ε and the parity of x.
  The even-x sheet does not enumerate the “other half” of the odd edges:
  it re-encodes the exact two-halving subfamily. It also identifies a
  symmetric two-sign bridge: 4d+ε and 8d+ε reach 6d+ε in the opposite
  signs, at one and two halvings respectively. This explains both the
  perfect gap and its rigid absolute midpoint. Modulus 9 gives the
  concrete bridge 5 → 7 ← 9.
- **[Diagonal joins see a different character](NORM.md).** If h counts
  Hamiltonian paths and P₂ ordered nonempty two-path covers, the profile
  is (x,y)=(P₂+h,h). Cyclic gluing with a singleton reads
  x_A x_B+2y_A y_B. Diagonal joins give x²+2y², not x²−2y², and both
  physical coordinates are odd, so the actual diagonal counts are
  always 3 mod 8. They exclude not only 7 and 21 but also 63, which is
  attained elsewhere. Off-diagonal gluing already gives the inert prime
  23 using transitive components of orders 2 and 3.
  Every allowed dyadic residue is actually attained by diagonal joins,
  yet their prime outputs have relative density zero in that same class. A norm collar alone
  cannot provide universal +2 repairs: a fixed source switch cannot
  continue producing exact +2 on unbounded host orders.
- **[Clock moves must preserve the carry](DYADIC_CLOCK.md).** The
  complementary address rail freezes one mod-4 bit. A final 1 on each
  address unlocks it without altering the three-center increment.
  Equal-length binary excursions are neutral on the shared address sum
  and can independently set higher bits. The [mixed CRT proof](FIXED_CRT.md)
  keeps the fixed odd residue simultaneously. The [five-bit bypass](TABLES.md)
  independently handles mod 24 by keeping both the address sum and the
  reversed coefficient sum unchanged.

The relevant local-global distinction is the same in all three examples:
a projected label reveals a lane but forgets how to continue inside it.
In the Collatz return equation the missing data are the ordered affine
carry and the integer starting point. In the norm picture they are the
signed pair determinant and the physical set of input profiles. In the
three-center bank they are the absolute baseline and the address cost of
programming more bits. This is a concrete analogue of separating
first-passage transport from population averages, not a logical reduction
between either conjecture or a consequence of the quasi-Riemann hypothesis.
See [Tao's first-passage transport paper](https://arxiv.org/abs/1909.03562)
and the [prior QRH transfer](../chanpak_qrh_residual_20261115/README.md).

## Metagraph interpretation: retain marked path-pair flux

There is a general exact port ledger behind this comparison. Suppose a
labeled metagraph move reverses the single arc u→v. Let C(a,b) count
ordered pairs of nonempty directed paths covering every vertex once,
the first ending at a and the second starting at b. Each pair uses
separate vertices for a and b. Its validity does not depend on the arc
between them. The edge weight of the move, measured by its change in
Hamiltonian count, is exactly C(v,u)−C(u,v): cut each lost or gained
Hamiltonian path at the reversed arc and concatenate to invert the cut.
Its halves are integral as a signed difference divided by two, since
both Hamiltonian counts are odd. Conversely, projection to the sole
scalar H or unmarked P₂ does not supply this signed flux formula's two
endpoint inputs. Merging by converse interchanges the endpoints and
reverses the local choice, but leaves the change between the corresponding
source and target counts unchanged. This is a port-local counterpart to
the [older electrical metagraph ledger](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/electric_metagraph_20261021_hodge_kirchhoff_determination_consecutive_ports.md);
it does not assert that any particular marked projection reconstructs
all metagraph wiring. It says which interface needs to be kept when
searching for a literal sign-controlled unit edge.

**Next proof obligation:** find a physical change of marked source or
reader with enough choices to absorb the *remaining integer quotient*
without leaving the crossing band or destroying the +2 steps. Bracketing
on a growing modulus is real progress, but it does not set that integer
quotient to zero.

Replay `python dyadic_clock.py`, `python verify_crt.py`, and
`python verify_tables.py` from this folder.
