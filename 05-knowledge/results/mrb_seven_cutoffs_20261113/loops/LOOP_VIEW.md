# What the loop viewpoint preserves, and where the cutoff lives

**PROVED elementary statements / inherited source recovery / OPEN transfer.**
The useful version of the loop viewpoint is a family of explicitly named
observations. The new [profile theorem](PROOF.md) realizes it on actual
tournaments: all uniform contexts modulo q see a finite path feeding an
exact cycle. Its transient records decay of nonunit modes; its loop records
unit modes. Changing q changes both quantities in a proved way.

## Incoming priority and complementarity

After the first derivation, commit `a03df47463` supplied
[Session297 CUTOFF](../../session297_parity_cutoff_20261113/CUTOFF.md),
[CYCLIC](../../session297_parity_cutoff_20261113/CYCLIC.md), and
[BRIDGES](../../session297_parity_cutoff_20261113/BRIDGES.md), all read here.
Credit that packet for its residue-plus-size cutoff criterion, its precise
independent-base automaticity criterion, and its actual cyclic-frame
origin action. It also recovers Session199's physical profinite-density
hostile with missing67; this is stronger than a merely artificial sparse
set. No discovery claim for those cutoff principles is made here.

Our new contribution is different: the *complete uniform-context profile*
of TT_n has an exact transient tau(q), exact least period L(q), and a
source-realized finite-box reduction for all transitive inflations of a
fixed quotient. Session297's cyclic-frame census acts on tournament-plus-
Hamilton-cycle markings; our loops act on congruence observations of module
sizes. Neither action identifies the other's physical states.

The pair {4,11} also has an independent source meaning, recovered in
Session297 BRIDGES: it is a literal least-exponent greedy inverse cycle for
the sign-specific3x-1 map. Both requested pairs occur in an older n14 CRT
runner construction. Our shifted H-coordinate below is an exact additional
coordinate relation, not an identification of these different dynamics.

## Three precise meanings of overlapping periods

1. **Finite periods.** A subset A of the integers is periodic with period m
   exactly when its membership factors through Z/mZ. The quotient is a loop
   of m residues. Finitely many such tests combine through the lcm of their
   periods; the CRT further separates coprime coordinates.
2. **All finite congruence observations together.** The map
   n -> (n modulo m)_m is injective: a nonzero integer cannot be divisible
   by every m. Its completion is the inverse system of compatible finite
   cyclic quotients, the profinite integers. Thus the collection of loops
   can retain every integer without equating n with n+m for any one m.
3. **General harmonic phases.** The map n -> n alpha modulo1 into the
   ordinary circle is injective for irrational alpha, since a nonzero
   integer times alpha is not integral. It preserves addition. It does not
   make an arbitrary membership indicator continuous, nor supply a finite
   periodic description. Rational alpha gives a finite cyclic quotient.

The inverse-limit definition and finite factorization of continuous
finite-valued observations are standard; Brian Conrad's
[Math249B notes, §1 and Lemma3.1](https://virtualmath1.stanford.edu/~conrad/249BW09Page/handouts/profcohom.pdf)
provide a primary course-note source. The theorem here is elementary and
does not require cohomology. In particular, a continuous Boolean observation
on the profinite integers factors through a finite quotient: each of its
two clopen fibers is covered by residue cylinders, compactness supplies
finitely many cylinders, and their moduli have a common multiple. Its
restriction to integers is therefore periodic. Finite startup exceptions
must be retained separately.

This recovers [THM-3364, cyclotomic Boolean clocks](../../../../01-canon/theorems/THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase.md):
eventual periodic supports modulo finite difference correspond to clopens,
and the full finite rational Fourier transform retains phase. A scalar
density or Fourier magnitudes alone lose that phase. The least period is
the lcm of the orders of supported characters. No new theorem is claimed
for this inherited equivalence.

## An exact boundary for infinite superpositions

Every set A of integers is a *pointwise* limit of periodic Boolean functions.
For N>=1, choose a function of period2N+1 agreeing with 1_A on [-N,N].
For each fixed integer n this agreement becomes permanent once N>=|n|.
Thus an unrestricted infinite periodic approximation can encode arbitrary
sets. This is a real representation theorem, but its approximation periods
and information content may grow without any useful bound.

Uniform convergence is sharply different. If Boolean periodic functions
converge uniformly to a Boolean function f, then eventually their sup-norm
distance to f is less than1/2. They must equal f everywhere. Thus f itself
is periodic. A continuous Boolean observation on the ordinary connected
circle is even more restrictive: it is constant, since the two fibers
would otherwise disconnect the circle.

An infinite intersection of periodic sets also need not be periodic:

    intersection_(m>=1) mZ = {0}.

The integer squares are infinite but not eventually periodic, because
their gaps 2k+1 are unbounded, whereas any infinite eventual periodic set
has bounded gaps. The primes have arbitrarily long composite gaps
(N+1)!+2,...,(N+1)!+(N+1), and are infinite by Euclid; they likewise cannot
be eventual periodic. These examples do not forbid useful local or
harmonic approximations. They identify the extra regularity or error
estimate such an approximation needs to support a proof.

An arbitrary finite set is eventually empty and hence has eventual period1.
That fact does not explain its particular cutoff: the discarded finite
startup contains all of its membership information. The new tournament
profile theorem improves on this tautology by deriving the *exact* cutoff
from prime-power valuations and the exact surviving loop from unit orders.

## Zero is a movable coordinate only if the operations move with it

For an additive torsor, choosing a different origin is legitimate and
preserves differences. For a group or ring with its operations held fixed,
zero is distinguished as the additive identity, and one as the
multiplicative identity. If new coordinates are x'=x+c, transported
arithmetic is

    x' plus_c y' = x'+y'-c,
    x' times_c y' = (x'-c)(y'-c)+c.

Simply shifting all labels while retaining ordinary addition and
multiplication does not preserve the original arithmetic. In finite loops,
phase or an origin mark is the sidecar that restores a named divisibility
predicate.

This is directly relevant to the two pairs. The recovered source
[microcosm gate note, §3](../../microcosm_gate_doubling_resolvent_20261106.md)
uses q=(H-1)/2, sending the forbidden tournament counts H=7,21 to q=3,10.
Writing instead r=(H+1)/2=q+1 sends the same forbidden counts to r=4,11.
In these defined coordinates the pairs are two affine descriptions of one
exclusion; this does not claim that every historical occurrence had that
meaning.
For ordinal sums, H multiplies. Its operation in these coordinates is

    q_join=q+q'+2qq',             identity q=0,
    r_join=2rr'-r-r'+1,           identity r=1.

This records exactly how the origin shift affects the operation. Modulo14,
translation by7 swaps3 with10 and4 with11; modulo7 each pair collapses.
That cyclic geometry is correct, but forgetting which representative and
which H-coordinate was used erases the distinction between exact counts7
and21. It does not construct either excluded count or a return past them.

## Recovered 77/14 precedents and their boundaries

- [THM-621, fourteen-fold ladder](../../../../01-canon/theorems/THM-621-fourteen-fold-ladder-law.md)
  repairs the old set identification: the surviving lifted set at the last
  row uses168=14*12, not25. Its arithmetic factor14 comes from the stated
  divisibility sieve and bounded lift family. It is not a universal
  fourteenfold law on arbitrary mathematical sets.
- [THM-3363, D14 complement-clock terminal](../../../../01-canon/theorems/THM-3363-d14-complement-clock-small-lrc-terminal.md)
  is a pointwise circle-cover argument on one specified LRC carrier. Its
  useful mechanism is adding the carrier's missing complement clock. A
  periodic residue word alone does not retain the carrier geometry.
- [Session240 source recovery](../../session240_seven_eleven_symmetry_20261110/provenance/README.md)
  distinguishes Paley affine actions, projective/design quotients,
  Collatz labels, binary codes and physical marked tori. Its nearby
  [THM-5102, minimum-order77 symmetry](../../../../01-canon/theorems/THM-5102-minimum-order-seventy-seven-symmetry-and-exact-mixed-path-cost.md)
  has a sharp fixed-order path minimum on19 vertices, explicitly not an
  all-order H floor. The common factorization77=7*11 does not identify
  that action with the 30-state uniform-context loop proved here.
- The [microcosm audit](../../microcosm_audit_H1_20261106.md), MISTAKE-729,
  repairs an asymptotic gate-tower error when the cyclic-triangle
  difference is zero. The theorem says normalized towers can forget a
  core while exact profiles still differ. Equation(9) in the new proof is
  another precise example of an observer forgetting growing magnitude.

The same workflow survives these different settings: specify the observer,
retain its phase and marks, compute its exact kernel, and prove the allowed
operations preserve that kernel. A modular loop is then a valid reusable
quotient. Turning it into an exact target return requires a new theorem
controlling what the observer discarded.

## Next proof obligation toward Chan–Pak

Use the exact profile orbit as a rejection/cache layer for uniform-module
experiments. A proposed replacement at modulus q can be certified in every
uniform exterior by finitely many run coordinates, and transitive-module
sizes reduce to the exact tail-plus-loop classes. Then retain the integer
D_2 difference, root/port data and full Hamiltonian count as explicit lift
obligations. The unproved step is producing positive physical replacements
that kill those differences at the requested exact target; neither CRT nor
an infinite collection of modular loops supplies that positivity.
