# Berggren q-spine cap-seven atom sieve and harmonic bifurcation

**Research synthesis, 2026-08-15.**  The proof source is
[THM-3455](../01-canon/theorems/THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum.md).
This file records how the newly incoming THM-3453 changed the concept board
and which LRC obligations remain.

## Incoming theorem, exact intersection

THM-3453 promoted the all-modulus literal half-twist support antichain

```text
{8,9,10,11,12,13,14,15,23,25,29,38,51,68,148}.
```

The parabolic Berggren branch supplies

```text
q_t=(2t+1)^2+2.
```

Their intersection is not a numerical scan.  It is a local-solubility sieve:
even atoms die because `q_t` is odd, the `5`-atoms die because `-2` is not a
square modulo five, and `13,23,29` die by `(-2|p)=-1`.  Only

```text
9,11,51
```

survive.  Rank priority then gives exact ranks `4,6,7`; rank five is absent.

This was the strongest mathematical signal in the incoming work.  THM-3453
did not itself mention the Berggren branch, while THM-3415 had reached only
the rank-four slice there.  The global atom classification turned the earlier
partial observation into a complete spine spectrum.

## Concept board update

1. **Global atom antichain:** exact divisor support for all literal cap-seven
   moduli.
2. **Quadratic branch:** the norm label `q_t=x^2+2` converts divisibility into
   square-root conditions.
3. **Rank priority:** overlaps select the smallest atom rank; membership alone
   is insufficient.
4. **Periodic spine-index subset:** the exact spine rank word has minimal period
   `1683`.
5. **Fibonacci selector:** pulling back through `t=F_n` replaces CRT period by
   Pisano period `360`.
6. **Harmonic carrier choice:** recurrence index, spine index/rooted depth,
   and nonlinear label give different convergence laws.

The overlap hostile is particularly useful:

```text
q_20=1683
```

is divisible by all three surviving atoms.  The correct exact rank is four.
This is the rank analogue of the origin sidecar in THM-3454: quotient data
must retain the operation that consumes it, here `min`, not merely the set of
available labels.

## Carrier-dependent asymptotics

The cap-seven set can be viewed through several natural-number carriers.

| carrier | membership density | reciprocal behavior |
|---|---:|---|
| all spine indices `t` | `76/187` | index harmonic sum has coefficient `76/187` |
| actual rooted depths `h=t-1` | `76/187` | after omitting `h=0`, depth harmonic sum has the same coefficient |
| nonlinear labels `q_t` | counting order `sqrt(X)` | reciprocal-label sum converges |
| labelled Fibonacci indices `n` | `43/90` | harmonic sum has coefficient `43/90` |
| Fibonacci spine indices `F_n`, positive rooted depths `F_n-1` (`n>=3`), or labels `q_(F_n)` | exponential sparsity | reciprocal sums converge |

This is the precise form of “every subset of the naturals is a subset of the
harmonic series.”  One must first say which natural-number coordinate is
being used.  The same Boolean selection diverges in the index coordinate and
converges after a nonlinear recurrence embedding.

MISTAKE-209 remains active.  The index spectrum is labelled and legitimately
counts `n=1,2` separately.  The unlabelled value sets identify
`F_1=F_2=1`; start at `n=2` before calling them support subsets.

## Typed connection ledger

| source | target | map | preserved | lost | sidecar | decisive hostile |
|---|---|---|---|---|---|---|
| THM-3453 atoms | `q_t` branch | solve `d|x^2+2` | divisor support | owner masks/time | exact atom rank | `q_20=1683` |
| spine-index residues | Fibonacci indices | `t=F_n` | Boolean membership | uniform residue distribution | Pisano state `(F_n,F_(n+1))` | period `360` |
| index subset | branch labels | `n -> q_(F_n)` | labels if retained | density and divergence | index coordinate | duplicate `n=1,2` |
| literal cover | fixed physical mode | `c=1/(2q)` | mask OR and zero cochain | nonzero current/endpoint/exit | relation-current sidecar | `q_4=83` negative cap-seven control |

The second row explains why looking only at `F_n mod m` one value at a time is
safe for membership but not for proving the period: the faithful recurrence
state is the ordered pair `(F_n,F_(n+1))`.

## What changed for LRC(14)

THM-3453 and THM-3455 close the **support** question on this branch:

```text
which q_t can admit a transverse literal cover by at most seven masks?
```

They do not answer:

```text
can the lawful fixed physical time be transported to the required nonzero
current, owner distinctness in the live row, endpoint typing, and decrement
exit?
```

The first four branch labels give a compact hostile battery:

```text
11 -> rank 6,
27 -> rank 4,
51 -> rank 7,
83 -> no cap-seven literal cover.
```

A fixed-centre bridge succeeds on the first three in their stated ranks and
fails through rank seven at the fourth.  A method that predicts the same
behavior for all four has forgotten atom support; a cap-seven method that succeeds at `83` is not
computing THM-3453's literal predicate.

## New frontiers

### 1. Algebraic-sequence atom sieve

Given a proved divisor antichain `D` and an integral sequence `f(t)`, compute
the local root sets `f(t)=0 mod d`, retain rank priority on overlaps, and only
then transport recurrence or harmonic structure.  THM-3455 is one positive
instance.  A second genuinely different instance is required before this
move belongs in `META-PATTERNS.md`.

### 2. Physical-time lift on the three positive atoms

Pull the explicit THM-3453 witnesses at `9,11,51` through the exact branch
residue classes and test the earliest lawful endpoint-current gate.  The
source must retain the owner residues, not only rank.  The hostile `q_20`
tests whether competing atom witnesses are incorrectly mixed.

### 3. Turnpike preorder atlas

THM-3454 exposed the six distances of four line points.  A separate exact
probe finds `25` realizable labelled weak orders: `10` strict, `10` with one
tie, `4` with two ties, and `1` with four tied pairs.  This is orthogonal to
the atom sieve and is the next clean tournament-with-missing-edges theorem
candidate.

### 4. Full-germ versus one harmonic coefficient

The index density is a useful scalar but does not recover the periodic rank
word.  Its faithful finite sidecar is the complete residue subset modulo the
minimal period; its analytic analogue is the full Dirichlet generating germ.
Test this against THM-3450's margin obstruction before proposing any D5/LRC
current transport.

## Boundary

The exact densities are structural census results, not probabilities of LRC
success.  The fixed common half-twist time and zero cochain do follow after a
labelled witness is restored; no nonzero physical current, decrement, spectral
nonvanishing, or Jacobian map follows.  LRC(14) remains open.
