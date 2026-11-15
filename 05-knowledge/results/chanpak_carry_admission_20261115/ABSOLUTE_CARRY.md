# Reverse the exact carry: a finite leading-coefficient interface

This is a target-first continuation of the [bounded-checksum interface](../chanpak_crossing_scales_20261115/CHROMATIC_INTERFACE.md) and the [two-core physical gluing](../chanpak_quadratic_anchor_atlas_20261115/README.md). It does not assume the low part is smaller than the gap scale.

For nonempty cores U,W, put
```
u=h(U), h=h(W), x=Q(U)−u, d=x Q(W), t=u h, G=d−x h+2t.
```
Here d is positive even, t positive odd. Mark the transitive spine of the host by bits 1 when the apex beats a vertex. Use a spine of form 1 0^g v, g≥1, with suffix v of width k and value b. Let r be the reversed k-bit value and p(b) the count of an apex over that suffix minus 1. If the suffix is empty, k=b=r=p(0)=0. Leading zero padding of a suffix makes no difference to p(b). With
```
D=(d+t)2^k−t r,    c=G+db+t p(b),
N=2^g D+c,
```
the count formula is exact without a no-carry assumption.

## 1. Normalize the leading coefficient, not the carried quotient

Absorb leading zeros of v into g. In the shortest nonempty suffix the reversed bit integer is odd, so D is odd. For the empty suffix D=d+t is also odd. Denote this canonical odd coefficient by λ. In **every** canonical presentation,
```
1≤c≤λ.
```
In fact the same inequality c≤D holds before normalization. The binary host bound is p(b)≤b. Reverse-complement preserves the apex tournament up to converse, so also p(b)≤2^k−1−r. Therefore
```
c−D ≤ −d(2^k−b−1)+t−x h ≤ 0
```
since x≥u. Equality requires an all-ones suffix (empty counts as all ones) and x=u, which forces U to be a singleton by the path-cut lower bound Q(U)≥(|U|+1)u. Conversely these conditions give c=λ. For canonical empty suffix, this gives the physical diagonal λ=H(W)+Q(W), not an assumed representation of arbitrary λ.

## 2. Invert an exact carried quotient

Given a prescribed odd N and gap g≥1, write M=2^g. Any witness has an odd canonical λ in the explicit backward interval
```
ceil(N/(M+1)) ≤ λ ≤ floor((N−1)/M).
```
Then c=N−Mλ is forced. Its carry into the observed quotient floor(N/M) is exactly floor(c/M), which is not required to vanish. Equivalently if that observed quotient is J, the possible carry j is J−λ, not a freely chosen auxiliary digit. At g≥(log₂ N)/2, the interval has length less than one, so it determines at most one integral candidate. In particular no separate search over low parts or high parts is needed there.

For every fixed λ this is a **finite all-order exact interface**, dual to fixing a checksum in the preceding session. Since
```
λ ≥ d·2^k+t,   d≥(2^|U|−1)2^|W| ≥2^(|U|+|W|−1),
```
the combined core order is at most 1+floor(log₂(λ−1)). From each literal core pair enumerate only widths with d·2^k+t≤λ, and recover
```
r=((d+t)2^k−λ)/t.
```
This must be integral and in [0,2^k−1], with the specified parity for a nonempty canonical suffix. Reverse its bits, compute the exact c, and keep the pair precisely when it matches the forced low part. The algorithm returns the actual host and cores. It works with arbitrary cores up to the forced order bound, not just the small palette. Repeating over 1≤g≤floor(log₂(N−1)) gives a complete decidable inverse for the full binary-gap grammar with unrestricted carries and unrestricted core order. This is not a universal witness theorem, and the core enumeration can be very expensive at small g.

## 3. A fixed-profile reduction removed by the canonical bound, and a literal carried hit

If a suffix starts with zeros, its raw D may be an arbitrarily large multiple of the same λ, while c remains unchanged. Before normalization this can look like an unlimited high-coefficient bank. The exact normalization shows there are finitely many checksums at fixed λ, even **without** a checksum budget or a no-carry restriction. Enlarging the gap introduces no new checksums.

Literal example: U=W=TT1 and suffix v=111111110 give (λ,c)=(1281,1279), using the [first gap-two crossing](../chanpak_crossing_scales_20261115/CHECKSUM_CROSSING.md). At g=1 the output is 3841. The observed quotient is 1920, not 1281; its carry is 639 and the residue is 1. The fixed-g inverse tests λ within [1281,1920] and the physical checksum recovers precisely this carried construction. The example illustrates inversion, not uniqueness of the representation at this small height.

Replay checks the inverse both against the forward signature bank and by counting a small literal example whose carry is nonzero. Neither replay nor this theorem asserts absence outside the stated construction.
