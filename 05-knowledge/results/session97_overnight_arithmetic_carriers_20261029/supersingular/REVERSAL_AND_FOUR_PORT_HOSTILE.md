# Reciprocal-choice torsors and a sharp four-port hostile

**Session97, 2026-10-29. Status:** PROVED reversal and multiplier laws;
FINITE-EXACT complete eight-code classification and observation threshold.
This extends [CYCLIC_KERNEL_RECURRENCE.md](CYCLIC_KERNEL_RECURRENCE.md).
The codes are classical objects; this packet's point is their explicit
source-preserving use as arithmetic kernels and joint-information hostiles.

## 1. Reversal acts on the retained choice bits, not on their count

When the cyclic theorem gives `R(n)=2^r`, each reciprocal pair of irreducible
binary cyclotomic factors requires a choice. The bank is naturally a product
of r two-point sets: an **affine Boolean choice space**, not a vector space
with a distinguished zero until one factor from each pair is named.

The coordinate reversal `i→-i mod n`, fixing the extra branch point, sends
each factor to its reciprocal. It therefore complements every choice bit.
It has no fixed kernel for r>0. In particular no positive-genus Lagrangian
for the `(1,n)` branch action is invariant under the full dihedral action
(rotation and reflection). If a self-reciprocal simple component exists,
there is already no rotation-invariant kernel; otherwise reflection swaps
all the choices. This covers every odd n>1.

Forgetting reversal merges pairs of kernels. Recovering one of the two
members needs one bit *after the reversal orbit is known*. Recovering an
arbitrary member of the whole bank needs r bits. These two information
problems must not be conflated.

Unit multipliers `i→ai mod n` act by permutations of the reciprocal pairs,
possibly swapping members. Thus they act by signed permutations on the
choice cube. The action of `(Z/n)^×` has kernel exactly `<2>`: fixing all
choices fixes each factor, including the doubling coset of1, so a must be a
power of2; powers of2 preserve every doubling coset. Hence

`(Z/n)^×/<2> → C2 wr S_r`

is faithful. This is an explicit occurrence of the signed-permutation carrier
from Session96, not an F1 realization of the Jacobian or Sha. The reflection
class -1 is the global bit complement.

For n7 and n23, the bank contains two kernels and reflection swaps them.
For n31, the bank has eight kernels and the multiplier quotient has order6.
Its orbits have sizes2 and6. The complete code certificate below verifies
these orbits, avoiding an inference from cardinality alone.

## 2. Eight length32 codes, identical weight enumerators

The n31 construction produces eight distinct marked extended binary codes
of length32 and dimension16. All have the same weight enumerator

`1+620z^8+13888z^12+36518z^16+13888z^20+620z^24+z^32`.

Let an **octad** mean the support of a weight-eight codeword. Each code has
620 octads. For a named coordinate set S define the weighted port response

`R_t(S)=#{octads B : S subset B}, |S|=t`.

For **every** one of these codes and **every named** S, exact enumeration gives

```
R1(S)=155;
R2(S)=35;
R3(S)=7.
```

Thus all point, pair and triple responses agree, not merely their histograms.
The full weight enumerator also agrees. None of these summaries distinguishes
even the two coordinate-isomorphism classes described next.

## 3. Fourth-order responses recover the missing marked source

For t4 the profile is no longer universal. The histogram records how many
of the `C(32,4)=35960` named four-subsets have each response:

| Code family | response0 | response1 | response2 | response3 | response4 | response7 |
|---|---:|---:|---:|---:|---:|---:|
| two codes | 11160 | 12400 | 7440 | 3720 | 1240 | 0 |
| six codes | 0 | 34720 | 0 | 0 | 0 | 1240 |

The family sizes refer to marked cyclic kernels, not to eight distinct
unmarked code isomorphism classes. Multiplier maps explicitly identify every
code within its two- or six-element orbit. The displayed invariant separates
the orbits under **all** coordinate permutations. Hence these eight marked
codes give **exactly two** S32-equivalence classes.

The generators of the two-element orbit are75035 and110921. Those of the
six-element orbit are73235,83903,92135,102641,118733,129989. Bit i is the
coefficient of x^i in the length31 even-subcode generator; the thirty-second
coordinate is the marked extension point.

A particularly cheap hostile is S={0,1,2,3}. It is contained in no octad for
both generators in the two-element orbit. Every four-subset, including this
one, is contained in an octad for every member of the six-element orbit.
This is an existential four-port obstruction missed by every triple response.

The **full weighted R4 table**, stronger than its histogram or its zero-set,
is distinct for all eight marked codes. Thus, within this exact bank,

* full weighted observations through order3 cannot distinguish any sources;
* full weighted observations through order4 reconstruct the marked source;
* the histogram of order4 answers the coarser two-class decision, not the
  full reconstruction problem.

This is a sharp finite observation threshold, not a universal four-port
carrier theorem. In particular, asking only whether an R4 entry vanishes
cannot distinguish the six zero-free codes. Multiplicity is load-bearing.

## 4. Map and loss ledger

**Source:** one of the eight Frobenius-stable Lagrangians at branch type(1,31),
with the cycle chart and fixed extension point.

**Target:** its extended binary code and its octad-incidence responses.
The code map is injective and transports the source's full symplectic kernel.
Replacing it by the weight enumerator plus all responses through order3
loses every choice among the eight. The operation “count minimum words
through these four named ports” observes the lost data.

The missing carrier can be restored by the full R4 table on this bank, or
more compactly by its three reciprocal-factor choice bits in the declared
chart. The two-class decision needs only the appropriate class bit, but that
bit cannot answer all future marked R4 queries. This is a finite exact
counterpart of the repo's distinction between faithful future-operation
representations and short certificates for a particular decision.

Neither tournament port order nor graceful-tree signs have been mapped
onto these codes here. The common theorem is the observation-loss pattern;
transporting specific operations requires a further explicit map.

## 5. Reproduction

`check_four_port_codes.py` imports the cyclic constructor, enumerates all
codewords in each of the eight length32 codes, obtains every octad, exhausts
all port incidences through order4, and computes all multiplier images.
It performs107 additional checks beyond its imported cyclic checks.
Normal and optimized Python outputs are identical. The output retains the
complete profiles, explicit zero witness and multiplier orbits.

`check_cyclic_frobenius.py` independently verifies the reversal law for every
constructed positive kernel, including all512 kernels at n127. Neither check
uses floating-point arithmetic or assertions that disappear under `-O`.
