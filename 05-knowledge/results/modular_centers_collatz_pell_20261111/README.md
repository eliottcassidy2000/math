# Modular centers, Pell plateaus, and an endpoint-symmetry obstruction

**11 November 2026. PROVED scoped statements; CITED Dirichlet corollary;
FINITE-EXACT computations. Universal odd seeds, Chan-Pak and Collatz
remain OPEN.** No priority claim or claim of a complete metagraph taxonomy.

The strongest tournament result is a new integration of the existing
qualifying-path compiler and local symmetries: if every set of at most2k
old vertices has a pointwise-fixing subgroup of order d, **every k-vertex
extension has Hamiltonian-path count divisible by d**, even if the
extended tournament itself is asymmetric. This gives rigorous exclusions
before an attachment search, rather than a condition on an average.

Consequences proved in [TOURNAMENTS.md, P11](TOURNAMENTS.md):

- For every nonprime Paley tournament on F_(p^f), q=p^f=3 mod4 and f>1,
  every one-vertex extension has H divisible by f. No selected set gives
  a unit-offset seed; neither branch3H+2 nor3H+4 is possible in one step
  from that actual source. This strengthens the inherited pair-only
  seed exclusion.
- For the recursive family R_(p,e)=P_p[R_(p,e-1)], p prime=3 mod4,
  every k-vertex extension has H divisible by p when2k<p^(e-1).
  At p=3,e=3 this covers every extension by up to four vertices.
- A ten-vertex example has trivial automorphism group and H=15033,
  yet inherits this divisibility from its nine-vertex base. The relevant
  symmetries act on endpoint-labelled coefficients, not on the child.

## The user's modular-table construction is exact

For m=4x+sigma, sigma=+1 or-1, the central entries are two copies each
of x and b=3x+sigma. Each row sums to m. These are the two occurrences
of m if row sums are included; the four individual entries are all <m.

| m mod8 | encoded raw edge |
|---:|---|
| 5 | x odd: x -> b=3x+1 |
| 1 | x even: x -> (b-1)/6=x/2 |
| 3 | x odd: x -> b=3x-1 |
| 7 | x even: x -> (b+1)/6=x/2 |

[PROOFS.md, P1–P3](PROOFS.md) gives the exact conjugate map on moduli,
the inverse gates17 and7 modulo24, and a residue guard plus affine
formula for every finite shortcut parity word. By Dirichlet's theorem,
every such word occurs at infinitely many prime table moduli and at
infinitely many composite ones. Initial-table primality alone therefore
does not rule out a finite growth pattern.

## What {24,25,26} contributes

In the square-sum graph their new-vertex parent sets are

    24: {1,12}; 25: {11,24}; 26: {10,23}.

The common degree2 is a Pell cancellation. Both inherited infinite Pell
families now have explicit formulas for the changing parent sets.
Every odd-square table R=(2h+1)^2 encodes the halving edge
h(h+1)->h(h+1)/2 and has the unique incoming table2R-1.

At R=25 the table path is49->25->13, encoding12->6->3. One Pell family
has square incoming moduli; the other has incoming moduli k^2+1. In the
latter family, writing the center R=r_j^2 gives

    M_j=2r_j^2-1: 1,17,577,19601,665857,...,
    M_(j+2)=34M_(j+1)-M_j.

Only power-of-two indices j can give prime M_j, but even those need not:
M_8=886731088897=257*1409*2448769. The surviving dyadic subsequence is
pairwise coprime, so it has infinitely many distinct prime factors.
Prime values belong to A002496 and are detected exactly by
A135244(2R)=(2R-1)^2. This is a candidate-family reduction, not a proof
of infinitely many prime values. See [P6–P8](PROOFS.md).

## Finite descriptions and tracked aspect sequences

- All multiplication-table row fibres are classified by gcd(a,m).
  Prime moduli maximize reversible nonzero rows. Under an explicit
  uniform-column model a known row loses exactly log2 gcd(a,m) bits.
  CRT composes the prime-power tables, retaining unit coordinates.
- The first singular centered row has radius(p+1)/2, where p is the
  least prime divisor of a composite odd modulus. At25 the radius is3.
  A fixed central unit mask cannot distinguish all primes from composites.
- Every finite Collatz parity word is an infinite arithmetic progression
  with an explicit endpoint formula; the guard must accompany the formula.
- In the recursive tournaments, arcs at valuation v number
  p^e(p-1)p^(e-v-1)/2. For the ternary family, the explicit wreath group
  has exactly2*3^((3^e-1)/2) free subsets, in two orbits with selected
  sizes(3^e-1)/2 and(3^e+1)/2. Rigidity alone still fails the seed equation.

The Jacobi comparison construction has missing edges at nonunit
differences for composites. Its maximum induced-tournament size is the
least prime factor. The recursively completed ring tournament on27
vertices is distinguished from field P27 by common-outneighbor counts
{4,10,12} versus{6}. These are actual named family classifications,
not a replacement for A000568 or its isomorphism quotient.

## Verification and reproduction

Run from this directory:

```text
python3 -B check.py
python3 -B -O check.py
```

Both runs produce identical [verification.json](verification.json) and
[manifest.json](manifest.json). The checks cover16,380 signed word
guards;2,727,100 literal multiplication-table cells;200 Pell terms;
all512 attachments of the recursive nine-vertex tournament; all351
field27 pair fixers; and ten recursive tournament family samples.
Hamiltonian and ordered-cover counts on three and nine vertices are
independently checked by literal permutation enumeration. Universal
claims rest on the proofs, not these finite controls.

[Inheritance and connection contracts](BOARD.md) ·
[Arithmetic proofs P1–P8](PROOFS.md) ·
[Tournament proofs P9–P12](TOURNAMENTS.md)
