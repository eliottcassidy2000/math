# Five spreads, first-25 quintiles, and target-selected returns

The principle is to **design the allowed continuation before throwing
away information**. A projection that reads the current scalar exactly
may still be useless for choosing the next move. This packet gives
three exact cases: a lossless reset with one retained bit, a periodic
sieve where the quotient can never prove primality, and a physical
tournament recursion where an otherwise relevant coordinate can
actually be discarded for a specific alphabet.

## 1. The 15 Petersen moves are five isotropic triples

Use the [five-label switching atlas](../planar_shadows_petersen_switching_20261116/README.md).
In one parity sheet of the switching observer, subtract a desired
target state from the current state, leaving an even Eulerian graph x
on five labels. Addition is symmetric difference. Write z=K5, the
converse displacement, and q(x)=|x|/2 mod 2. Its polar form b is edge
overlap parity, with radical <z> and q(z)=1. The 15 possible free-frame
moves are exactly the nonzero q=0 states, all 4-cycles.

Each square omits one of five labels. On each remaining K4 there
are three squares. Their pairwise sum is the third and the overlap
parity vanishes (any two share two edges). In the four-dimensional
representation W modulo <z>, the five triples are therefore a
partition of the 15 nonzero points into five two-dimensional
isotropic subspaces, a symplectic spread. Equivalently W is F4²
and these are its five F4 lines. To see the coordinates, choose two
of the subspaces as the two coordinate planes; each other one is
the graph of an invertible binary 2-by-2 matrix. Normalize one to
identity. The last two matrices have no fixed nonzero vector, so
they are the two inverses of order 3 in GL2(F2), exactly the two
nontrivial F4 scalar multipliers.

There is a real Fermat connection here: 3=F_0 is the number of
nonzero scalars on a line, 5=F_1 the number of lines, and 15=F_2−2
the number of nonzero vectors. In general these counts are q−1,
q+1,q²−1 for q=2^(2^j), with q+1=F_j. The A5 action on the five
omitted labels is the projective SL2(F4) action in these coordinates:
a binary linear spread-preserver is semilinear over F4 (fix three
lines first, leaving a simultaneous GL2(F2) map). The A5 subgroup
has no nontrivial homomorphism to the order-two field-automorphism
group or to the determinant group of order three. Its faithful
action therefore lies in SL2(F4), itself of order 60. This is
independent of the archive's Fermat-prime theorem for *circulant*
tournaments (THM-4644). The five known Fermat primes do not
canonically label these five lines; in particular 257 and 65537
are outside the first-25 prime window.

**Deterministic reset.** Every nonzero projected target difference
in W selects its unique square move. If the missing converse bit
is available, the full 32 states have exact distance layers
1,15,15,1. For q(x)=0 and x≠0 use x. For q(x)=1 and x≠z, select
a square u with b(x,u)=1 (squares span W); then x+u is another
square. The state z is not one move, and cannot be a sum of two
because the two lifts of a projected move have opposite q. Any
square u takes it to the preceding distance-2 case. Lexicographic
choice of a distance-decreasing move gives one outward edge per
state, except the base, and reaches the base in at most three.
This is a literal terminating error-corrector in the switching
observer. It does not make cut switching preserve the Hamilton
count of a single tournament.

## 2. Monster, first 25 primes, and first 25 square-sum graphs

[The quintile follow-up](QUINTILES.md) classifies the exact square
reflection rail frame. The first width is 5, i.e. the graph on
[1,25]; subsequent widths are 120,2760,63365,..., by a positive
Pell descent with inverse strictly decreasing to (1,1). Here
25 is a **first possible configuration**. The five consecutive
blocks and the five rails give different decompositions of the
same labels. Along graph order the first path sizes at most 25
are 15,16,17,23,25; the two rail matchings alone are not a
Hamilton-path proof. Pulling the prime marks back to ranks is
legitimate, but replacing ranks by primes breaks nine of the
ten rail edges. That is the concrete loss test.

The shared 15 with Ogg is not an identification. The modular
filtration (levels 1,2,6) gives 6+5+4, not five natural triples.
Fi22's first six primes survive to level six; the Baby layer
adds {17,19,23,31,47}, with unchanged valuations in the Monster;
the Monster adds {29,41,59,71}. These four have a particularly
useful loss mechanism: the earlier level-two calculation found
two parallel coarse isogeny edges exchanged by Frobenius, with
rational vertices. Forgetting the edge makes the obstruction
vanish. This is genuinely analogous to forgetting a continuation
bit, not a map from A5 moves to supersingular primes. Sources:
[Aricheta's theorem](https://arxiv.org/abs/1809.07421),
[Duncan and Ono on Ogg's observation](https://arxiv.org/abs/1411.5354),
and the repo [marked level-two calculation](../session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md).
The earlier [Paley/Monster crossing](../session32_square_monster_20261023/paley/PALEY_MONSTER_BRIDGE.md)
does give a literal group bridge, via the seven-point affine
quotient of a Monster normalizer, not via equality of lists.

The exact exponent sieve for 99·2^k+23 has least period 61,209,720,
accepting 32,353,278 exponents per period. The whole third prime
quintile adds no rejection. The nested Fi22/Baby/Monster wheels
keep 40,806,480 / 34,855,535 / 34,708,960 classes per the same
common period; the first-25 wheel keeps fewer still. These are
different ways to cut the same prime window, by rank or by group
support. At k=16 the number 6,488,087=1301·4987 slips through the
first-25 wheel; in fact every k=16 mod 306,048,600 slips through
and is composite. More generally any nonempty finite sieve of
this form has infinitely many false positives: take any positive
survivor with value >1, pick one of its prime divisors q outside
the sieve, and repeat its phase modulo ord_q(2) as well as the
finite sieve period. The growing value stays divisible by q.
This is the precise obstruction to treating periodic residue
readiness as an absolute certificate.

Archive tour (distinguish exact transfers from similarities):
- [First arrivals, Pell resets, and Monster support](../monster_levels_square_sum_7_21_20261022.md)
  and [Session32 port synthesis](../session32_square_monster_20261023/SYNTHESIS.md).
- [Prime block hostile test](../session32_square_monster_20261023/PRIME_BLOCK_CARDINALITY_HOSTILE.md):
  the square-sum graph on Ogg primes *themselves* is a five-leaf
  tree, not the Hamiltonian path on consecutive ordinal labels.
- [Ogg, triangular numbers and Chowla descent](../ogg_triangular_chowla_20260927.md);
  [Fermat/Sylow symmetry](../gilbreath_fermat_platonic_20260926.md);
  [Fermat and Wieferich order](../fermat_wieferich_order_20261014e.md).
- [Golden route compiler](../golden_primes_carries_route_compiler_20261003.md):
  carry bits, elapsed time and the original descent threshold all
  require marks; a compressed route can hide the condition needed
  at the next step. This is the operational analogy to the sieve.
- [Finite-cut square-sum retention](../square_sum_finite_cut_retention_20261019/FINITE_CUT_RETENTION.md):
  retain endpoint incidences, not merely counts.

## 3. Physical three-way target selection and its safe alphabet

See [TARGETS.md](TARGETS.md). Three different nine-vertex raw seeds
compile every odd target in the bands 756·2^s−33+{0,2,4}, s≥2,
with an identical response for future near-sink blocks, even though
deletion counts are unequal. The inverse is deterministic on these
bands (recover the unique digit and dyadic height). Arbitrary
near-sink blocks enlarge them by an exact formula, with certified
prime and 7p outputs from composite blocks. Arbitrary binary
control does **not** preserve the shared response, so this is a
restricted alphabet, not a proof of Chan–Pak.

Replays: python switch_reset.py, python wheel.py, python quintiles.py,
and python check_tournaments.py. All use the standard library.
