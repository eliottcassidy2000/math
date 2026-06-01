---
id: HYP-2051
status: OPEN
source: codex-2026-06-01-S550; codex-2026-06-01-S551 addendum
related:
  - HYP-1983
  - HYP-1984
  - HYP-1994
  - HYP-2040
  - HYP-2044
  - HYP-2049
---

# HYP-2051: doubled primes are the one-dyadic lift that lets prime cycles cross parity

**Claim.** The doubled primes `2p` are important because they are the minimal
even objects that still remember a prime odd core.  They sit at dyadic height
`1` in the natural address

```text
N = 2^h * odd_core,
```

with `odd_core=p` prime.  Thus `2p` is not just "an even composite"; it is the
first even sheet over an irreducible odd cycle.

This explains the Goldbach/Lemoine split:

```text
even target:  N = p + q       (prime + prime)
odd target:   N = p + 2q      (prime + doubled prime)
```

Even Goldbach uses two odd prime cycles, so the parity channel is

```text
1 + 1 = 0 mod 2.
```

Lemoine/Levy uses one odd prime cycle and one one-dyadic lift, so the parity
channel is

```text
1 + 0 = 1 mod 2.
```

The doubled prime is the parity bridge that keeps the problem binary.

## Cycle reading

Odd cycles are the primitive parity obstructions: they cannot be 2-colored,
and in the tournament conflict graph they are the objects that enter the
odd-cycle formula for `H`.  Even cycles are compatible with a two-sheet
alternation.  The length `2p` has a special role: the cycle `C_{2p}` is the
bipartite double cover of the prime odd cycle `C_p`.

```text
prime p      = irreducible odd cycle core;
doubled 2p   = bipartite two-sheet lift of that core.
```

So a doubled prime is even enough to live on the parity-balanced sheet, but
still prime enough not to have any nontrivial odd-core factorization.

## Additive reading

Lemoine's form

```text
N = p + 2q
```

is also the diagonal slice of ternary Goldbach:

```text
N = p + q + q.
```

The doubled prime packages two equal primes into one even object.  This is the
arithmetical compression: instead of using three independent odd primes to hit
an odd number, it asks whether one free prime plus a repeated-prime pair is
already dense enough.

Thus `2q` is a diagonal Goldbach gate.  It is the simplest even Goldbach pair
available, namely `(q,q)`, carried as a single term.

## Residue-wheel reading

On the mod-6 wheel, large primes occupy the unit beads

```text
P = {1,5}.
```

Doubling sends them to the even unit-lift beads

```text
2P = {2,4}.
```

Then:

```text
P + P   = {0,2,4} mod 6   (all even residue channels),
P + 2P  = {1,3,5} mod 6   (all odd residue channels).
```

Goldbach and Lemoine are therefore the same convolutional mechanism on the
prime wheel, with `2P` providing the parity-shifted copy of `P`.  For odd
moduli, multiplication by `2` is a unit automorphism, so the doubled-prime set
is not a random sparse set; it is the prime set moved to the other parity
sheet.

## Relation to the repo stack

HYP-1984 says addition is the horizontal `x+2` direction and multiplication is
the vertical `x*2` direction.  The doubled primes are exactly the first vertical
move applied to a prime:

```text
p  --x*2-->  2p.
```

That is why they keep recurring near first-even frontiers such as `n=2p`:
they carry one dyadic seam and one prime branch, with no extra odd factors and
no higher dyadic depth.  In LRC language they are rank-one parity lifts, not
generic composites.

HYP-2049 strengthens this parity-lift reading into the recursion-bridge
language: `2p=p+p=2*p` is the rung where the additive Goldbach diagonal and the
multiplicative doubling bridge coincide.

## n=4 LRC boundary reading

The n=4 character layer makes the gate sharper.  In the mod-4 Fourier
description, the safe-indicator pair correction is controlled by the odd
character `chi4`; even reduced cofactors are pairwise silent.  At the same
time, HYP-2040 and THM-391 record the boundary fact that all-odd triples have
the closed witness `t=1/4`, so any tight n=4 obstruction must involve an even
speed.

Thus a doubled-prime leg can be boundary-active while remaining pairwise
quiet.  It can break the all-odd `t=1/4` witness by spending one dyadic factor,
but because the even reduced cofactor kills the `chi4` pair correction, that
debt need not appear in the easy `R2` ledger.  The natural place to look is
therefore the endpoint, full-support, or higher-order resonance layer.

## Tournament Analysis next object

Do not default to runner vertices for this layer.  The object that preserves
the Goldbach/Lemoine distinction is a bridge cell

```text
b = (N,p,q),
```

with `N=p+2q` for Lemoine bridges and `N=p+q` for Goldbach comparison
bridges.  Useful observables are dyadic height, odd core, representation
abundance in the fixed-`N` fiber, endpoint-debt proxy, `chi4` neutrality, and
endpoint-owner labels when available.

A first switch should orient `b1 -> b2` by lower boundary debt or higher repair
potential, with ties ordered by increasing `(N,q,p)`.  Report score
histograms, directed cycles, SCCs, Hamiltonian-path counts, and edge flips
against the S513 add/multiply denominator tournament.

## Predictions

1. Lemoine representation counts should behave like a cross-correlation
   between the prime set `P` and its doubled image `2P`, not like a new
   independent sequence.
2. Any obstruction to Lemoine should be a wheel-channel thinning of
   `P + 2P`, analogous to the twin-center holes in HYP-1994.
3. Denominators `2p` should be especially clean in LRC/Tournament Analysis:
   they have one dyadic seam and one prime odd core, hence rank-one p-adic
   channel structure.
4. In cycle language, `2p` should be treated as the bipartite lift of the
   prime cycle `p`, not as a featureless even length.
5. The doubled-prime gate should help separate odd/even additive questions:
   `P+P` covers even channels, while `P+2P` covers odd channels.
6. In n=4-style LRC decompositions, doubled-prime legs should be boundary
   active but pairwise quiet, moving debt out of the `chi4` pair ledger and
   into endpoint/full-support bridge-cell fingerprints.

## Status

Open synthesis.  The parity, residue, cycle, and n=4 pairwise-quiet identities
are exact; the unproved content is the explanatory force for Lemoine density
and for LRC first-even denominators.

## Files

`07-reflections/doubled-primes-as-parity-lift-s550.md`.
