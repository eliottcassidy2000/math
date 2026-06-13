---
id: HYP-2157
status: synthesis / proof-program
source: user-2026-06-03; codex-2026-06-03-S604
updates:
  - codex-2026-06-03-S605
related:
  - HYP-2142
  - HYP-2144
  - HYP-2146
  - HYP-2151
  - HYP-2152
  - HYP-2153
  - HYP-2154
  - HYP-2155
  - HYP-2156
  - HYP-2161
  - THM-002
  - THM-401
  - THM-406
---

# HYP-2157: overlap-order calculus unifies Helly, Vitali, correlated residues, and partition functions

## Claim

The common object behind the recent LRC/Helly/OCF threads is an
overlap-order or compatibility-order partition function.

For LRC covering depth, THM-406 gives overlap volumes

```text
S_j = sum_{|I|=j} measure(intersection_{i in I} D_i),
p_0 = sum_j (-1)^j S_j.
```

For OCF, the odd-cycle conflict graph has independence polynomial

```text
I(Omega(T), x) = sum_j alpha_j x^j,
H(T) = I(Omega(T), 2).
```

These are sibling packages.  LRC stores coverage multiplicities and overlap
orders in `{S_j}` or `{p_k}`; OCF stores compatible odd-cycle packets in
`{alpha_j}`.  In both cases, a proof asks how many orders of this package must
be retained before the target predicate is decided.

## Dictionary

```text
Helly number
  = how many overlap/compatibility orders a certificate must retain.

Vitali wall
  = no fixed finite truncation of those orders decides the predicate.

Correlated residue
  = a structured locus where low-order density is blind, so residue/CRT or
    arithmetic packet data replaces first-moment reasoning.

Partition function
  = the sibling world where all compatible orders are already packaged into a
    polynomial or generating function.
```

This refines the slogan "Helly number is proof difficulty."  More precisely,
Helly number is the least order cutoff `h` such that order data up to `h`
contains a certificate.  If no uniform finite `h` works, the problem is at a
Vitali wall and the proof must use a compressed all-orders structure instead.

## LRC Consequence

For LRC, the first moment is blind:

```text
S_1 = E[depth] = 2n delta.
```

THM-406 shows that even `S_2` does not separate additive-chain collapse from
generic rows.  Thus HYP-2153's `p_0=0` family is not a second-moment phenomenon.
It is an all-orders alternating cancellation in `sum (-1)^j S_j`.

The additive-chain and shell labels are candidates for a compressed
all-orders certificate.  They do not replace `{S_j}` by a low-order moment;
they explain why the full alternating sum collapses while the unit boundary
witness floor survives.

## Category / Number-Theory Refinement

The all-orders package becomes canonical through a categorical coimage:

```text
clock circle -> coimage(depth observable) -> depth values.
```

Yoneda gives the recognition principle: if the depth law, factorial moments,
spectral measure, Helly probes, CRT witnesses, and partition-function
evaluations all represent the same answer functor, then the package is not an
arbitrary invariant.  It is the coimage of the observable.

Number theory supplies one of the most important probe families.  At
`delta=1/n`, THM-401/S571 compresses endpoint witnesses through

```text
C = 2n - 1,
V mod C,
antipodal shells {a,-a},
unit clocks a^{-1},
gcd strata for composite C.
```

The `C`-resonance quotient is a compressed all-orders coordinate.  If a unit
shell is missed, a witness exists and `p_0` cannot collapse.  If unit shells are
covered, the remaining cancellation question moves to additive resonances and
nonunit strata.  Thus `2n-1` resonances are not merely a floor trick: they are
the number-theoretic observer maps that expose or block the inclusion-exclusion
cancellation.

The incoming bounded-CRT automaton work through `n=20` (HYP-2142) is the same
kind of move on the two-block side: density and low-order measures are blind,
but correlated residue states can certify emptiness of the live bounded
intersection.

## Collatz / Two-Block Reading

The Collatz rapidity identities and the LRC two-block determinant automata both
live in the "correlated residue where density is blind" row of the dictionary.
The coarse drift, first moment, or pair energy gives the right ambient scale
but misses the residue condition that decides the exceptional behavior.

In the two-block LRC branch, the determinant/CRT state is a compressed
all-orders witness: it remembers enough residue compatibility to decide
bounded emptiness without expanding every overlap order.  This is why Helly
smallness and CRT automaton emptiness belong in the same ledger.

## OCF / Partition-Function Sibling

OCF is the sibling world where the all-orders object is already explicit:

```text
H(T) = I(Omega(T), 2) = sum_j alpha_j 2^j.
```

The coefficients `alpha_j` are not moments of circular arcs, but they play the
same structural role as overlap orders: they count compatible packets of size
`j`.  Truncating OCF at `j <= h` is the exact analogue of keeping only overlap
orders up to `h` in LRC.  The independence polynomial is therefore the packaged
depth distribution for the tournament/OCF world.

This explains why OCF feels easier to state than LRC: the partition function
has already performed the coimage step.  LRC is still trying to discover the
right compressed coordinates for its full overlap sequence.

## Assumption Challenge

The natural vertices are not necessarily runners, arcs, or tournament edges.
For this synthesis the candidate vertex sets are:

```text
overlap orders S_j,
moment truncations,
cover arcs,
CRT residue states,
two-block determinant components,
OCF independent odd-cycle packets,
partition-function coefficients alpha_j,
proof obligations.
```

The chosen quotient should preserve the target predicate: `p_0=0`,
bounded-CRT emptiness, or `H(T)=I(Omega,2)`.  It may destroy phase order,
individual arc geometry, or packet identities.  That loss is acceptable only
when the retained object still decides the predicate or yields a certificate.

## Proof Program

1. For each problem, define the all-orders package `Z(x)=sum a_j x^j`.
2. Identify the target functional: `p_0`, bounded intersection emptiness,
   Hamiltonian path count, or a proof obstruction.
3. Determine whether a finite cutoff `h` decides the target.  If yes, `h` is
   the Helly/order number for that proof route.
4. If no uniform finite cutoff works, find a compressed all-orders coordinate:
   additive-chain shell data, CRT automaton state, determinant residue, or an
   explicit partition function.
5. For LRC `n=14`, test whether the `C=27` shell/CRT structure is a compressed
   all-orders certificate for the `p_0` collapse branch, rather than a low-order
   density estimate.

HYP-2161 sharpens step 5: the `C=2n-1` shell/lift/CRT ledger should be tested as
a Yoneda-style conservative probe family for the coimage.  In that reading,
"`2n-1` resonances are the cancellation" means the nonzero relation-lattice
corrections that cancel the free `p_0` term are first finitely visible through
the THM-401 shell probes.

## See

HYP-2142, HYP-2144, HYP-2146, HYP-2151, HYP-2152, HYP-2153,
HYP-2154, HYP-2155, HYP-2156, HYP-2161, THM-002, THM-406,
`05-knowledge/hypotheses/HYP-2153-lrc-p0-collapse-additive-chains.md`,
`01-canon/theorems/THM-002-ocf.md`,
`01-canon/theorems/THM-406-covering-depth-master-object-factorial-moments-and-spectral-identity.md`.
