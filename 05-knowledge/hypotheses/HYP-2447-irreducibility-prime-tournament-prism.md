# HYP-2447 - Irreducibility and prime values form a tournament witness prism

**Status:** OPEN synthesis; finite carrier atlas verified.
**Source:** codex-2026-06-12.
**Companions:** HYP-2226, HYP-2227, HYP-2446, HYP-2445, HYP-2443,
THM-410.
**Computation:** `04-computation/irreducibility_prime_tournament_prism_codex.py`;
stored output `05-knowledge/results/irreducibility_prime_tournament_prism_codex.out`.
**External anchors:** Singh arXiv:2411.18366; Iravanian arXiv:2410.15880.

## Statement

The two-way relationship between primes and irreducible integer polynomials is
a clean model for the repo's larger scalar/support theme:

```text
prime values can certify irreducibility,
irreducibility plus no fixed local obstruction should force prime values,
prime digits can generate irreducible polynomials,
reducible integer polynomials are recovered by recombining looser-field atoms.
```

This should be treated as a tournament prism.  The vertices are not raw
integers, but witness slots:

```text
value factorization,
fixed divisor / admissibility,
digit place-value address,
root-trace recombination atoms,
Heegner/class-group side channel,
raw prime-run length.
```

Edges point from the side channel that preserves more irreducibility
information to the one it explains after scalar projection.

## Evidence

The script verifies a small typed atlas.

### Bunyakovsky / fixed divisor lane

```text
Euler41:        irreducible=True, fixed_divisor=1, prime_run=41
Rabinowitsch41: irreducible=True, fixed_divisor=1, prime_run=40
BlockedEven:    irreducible=True, fixed_divisor=2, prime_run=1
NonHeegner7:    irreducible=True, fixed_divisor=1, prime_run=1
```

This separates three distinct facts that scalar prime-run length can conflate:

- irreducibility;
- admissibility / no fixed divisor;
- unusually long early prime horizon.

Bunyakovsky's conjecture is the infinite horizon version of the middle two
conditions, not a claim that every admissible irreducible polynomial has a long
initial prime run.

### Singh / Murty value-factor lane

The Singh paper proves upper bounds on the number of irreducible factors from
prime-factor data of sufficiently large values.  The atlas reproduces the
paper-style examples:

```text
x^4 - x^2 - 2x - 1 at m=3:
  f(m)=65=5*13, actual Z-factors=2.

x^6 - x^2 - 2x - 1 at m=3:
  f(m)=713=23*31, actual Z-factors=2.

x^4 + 10x^2 + 1 at m=12:
  f(m)=22177=67*331, actual Z-factors=1.
```

The first two attain the factor-count bound; the quartic shows why a value
factor bound is a support gate, not the whole proof.  The value has two prime
factors, but coefficient symmetry excludes the quadratic split.

The script also records the threshold nuance in the source: the theorem line
prints an `H_f+d+1` condition, while the displayed examples use `H_f+d`.
Future formalization should resolve this from the published text before
canonizing the exact inequality.

### Cohn digit-prime lane

Cohn's criterion sends a prime written in base `b` to an irreducible digit
polynomial.  The atlas samples primes that already occurred in the recent
Frobenius/LRC residue work:

```text
2339 -> 2*x^3 + 3*x^2 + 3*x + 9
1091 -> x^3 + 9*x + 1
6551 -> 6*x^3 + 5*x^2 + 5*x + 1
7643 -> 7*x^3 + 6*x^2 + 4*x + 3
10007 -> x^4 + 7
```

All are irreducible over `Z` in the check.  This is the reverse arrow:

```text
prime number + place-value address -> irreducible polynomial.
```

### Real-factor recombination lane

Iravanian's paper revives factorization by first splitting over `R`, then
recombining linear/quadratic real atoms via subset-sum.  The toy polynomial

```text
(x^4+1)(x^2+1)
```

has real quadratic atoms with traces `+sqrt(2)`, `-sqrt(2)`, and `0`.  Integer
factors are exactly the subsets whose trace obstruction cancels and whose atom
groups are complete:

```text
{+sqrt(2), -sqrt(2)} -> x^4+1
{0}                  -> x^2+1
all atoms            -> whole polynomial
```

This is the recombination face of the same prism.  Reducibility is hidden when
we project to real atoms; it reappears as an integer subset-sum condition.

## Tournament Analysis

The stored proof-route tournament uses the observable

```text
(side_channel_strength, exactness, repo_bridge, computability-risk).
```

It is transitive:

```text
singh_value_factor_bound
> bunyakovsky_admissible_irreducible
> cohn_digit_prime
> heegner_euler_horizon
> real_factor_recombination
> raw_prime_run_length

score histogram: {0:1,1:1,2:1,3:1,4:1,5:1}
directed 3-cycles: 0
Hamiltonian paths: 1
edge flips versus scalar-only ranking: 3
```

The edge flips are the point: scalar exactness alone ranks the lanes
differently from side-channel strength.  This is the same lesson as HYP-2446
for p-curvature and HYP-2443 for LRC14 blockers.

## Transfer To The Repo

The proposed dictionary is:

```text
integer polynomial factor       <-> module / SCC / support packet
prime value                     <-> single irreducible witness
fixed divisor                   <-> unavoidable local obstruction
Bunyakovsky admissibility       <-> no local obstruction at any finite gate
Cohn digit address              <-> scalar prime plus place-value side channel
real-factor recombination       <-> lift to looser atoms, then subset-sum back
tournament Hamiltonian path     <-> recombination order through witness slots
directed cycle                  <-> incompatible side-channel commitments
```

For LRC14, this says the finite witness is analogous to a prime value of the
denominator-cover polynomial.  A counterexample would be an admissible-looking
irreducible object with no prime witness anywhere.  The project has been
trying to prove that no such object exists by retaining the missing side
channels: carry, owner, Pisano class, divisor fiber, and marked support.

## Next Moves

1. Formalize a small "fixed divisor" analogue for LRC denominator ledgers:
   which runner residues block every lift the way `x^2+x+2` is always even?
2. Extend the recombination toy to tournament factors: split a tournament into
   real/spectral/Pfaffian atoms and test which atom subsets recombine to
   integer `H` or `Pf` packets.
3. Add a Cohn-style place-value lane to number-theory tournaments: primes as
   digit polynomials in base `b`, with arcs determined by digit-carry
   comparisons.
4. Compare Bunyakovsky admissibility with the p-curvature ledger HYP-2446:
   both say that local obstructions must vanish everywhere, but only a
   retained side channel can explain why prime/witness values should persist.
