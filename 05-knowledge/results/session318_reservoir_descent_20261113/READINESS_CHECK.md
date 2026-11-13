# Two exact-precision checks after the reservoir change

## The equal-count switch is not a two-slot-neutral move

For the exact switch in [Dyadic Switch](DYADIC_SWITCH.md), use the order of
charts printed there. The full ordered two-slot cover counts satisfy
```
Q_2(first)-Q_2(second) = 12 - 4*3^(k-2).
```
In particular they do not agree for k>=4. This is a concrete loss of a
continuation coordinate in an otherwise exact equal-count switch.

**Proof.** In both binary prefixes, the first 1 is followed by the same
block length l=k−2 of consecutive zeros, followed by a fixed low word
(11 in the first, 1 in the second). The zero block is a transitive
module. All vertices outside the block other than the initial 1 lose
to every vertex of the block. Any path returning to the block after
leaving it must return via that initial 1. The only possible incoming
edge at this first prefix vertex is from apex a; hence there is at
most one return to the block across both slots. Accordingly two
ordered slots give at most three visits to the block. In a collapsed
layout of the two slots containing j visits, partition the l named
block vertices into j ordered nonempty subsets. There are j! S(l,j)
choices, and each subset has one order because the block is transitive.
So Q_2 is a linear combination of 1, 2^l and 3^l for l>=1, by
inclusion-exclusion on the used visits. Direct subset DP gives the
two counts at k=4,5,6:
```
1416 1440
3636 3732
9624 9936
```
The three evaluation points determine the coefficients (distinct bases
give a Vandermonde matrix), yielding the displayed identity for every k.
This concerns one tracked Q coordinate; it does not rule out a different
continuation that does not need Q_2.

## Native residues couple the two read directions

In [Bidirectional Reservoir](BIDIRECTIONAL_RESERVOIR.md), prescribing
the first L **leading** bits of r leaves the coarse reversed end response
free, giving the joint normalized rectangle. Prescribing instead
r modulo 2^L fixes the **last** L chain bits, which are exactly the
leading digits of the end readout. More precisely all cuts with that
residue have E/3^m in a single interval of width at most
2(2/3)^L+2(2/3)^m. This follows by summing the unprescribed geometric
tail and the nonnegative exact end corrections (total at most 2(2/3)^m
after normalization). Thus a deep binary guard already pins the coarse
end response; it is not a compatibility coordinate independent of the
new reservoir. The guard should be checked before composing any proposed
clock refill with this end decoder.
