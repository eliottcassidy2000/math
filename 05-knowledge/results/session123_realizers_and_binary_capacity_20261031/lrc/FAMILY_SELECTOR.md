# A direct polynomial-bit selector for every certified five-tail family

2026-10-31. **PROVED WITH INDEPENDENT PROOF/IMPLEMENTATION AUDIT.** Depends on the completed
Session122 eight-core/five-tail theorem. This is constructive progress, not a
new LRC family or radius.

## Theorem

For each eight-element subset C of{1,...,13}, and every five positive integer
speeds (allowing repetitions), `family_selector.py` returns a rational time
safe for every original speed at radius1/14. Its bit complexity is polynomial
in the total binary encoding of the five arbitrary speeds, with the finite
bounded pool fixed. The implementation expands no unbounded-speed comb.
It also accepts any distinct bounded core of size8,...,13 with exactly13-m
attached speeds, by the same theorem and padding.

## Universal completion forces a positive interval

Fix any prefix of attached speeds, leaving k>=1 arbitrary tails unassigned.
Its complete common safe set is a finite union of rational closed intervals
and isolated rational points. It cannot be empty, since the universal family
theorem supplies a safe time for any remaining tail choice.

If it had no positive-width component, it would be a finite set of rational
points. Choose a positive integer D divisible by every reduced denominator
of those points. Then ||Dt||=0 at every safe point. Choosing each remaining
tail equal to D destroys all possible times, contradicting universal
completion. Therefore the complete prefix cover has positive width whenever
at least one arbitrary tail remains.

This uses the actual universal-completion predicate, not the false claim
that every nonempty rational closed cover has positive width. Point-only
covers are allowed after the last speed has been assigned.

## Algorithm and inherited constraints

Start with any bounded core of size m>=8 and13-m tails. Sort the tails.
Before creating any interval cover, promote every tail v<=13 into the core.
If it is already present, add the smallest missing bounded speed instead;
then remove that redundant tail. Each step strengthens the constraints and
preserves the total of13 core/tail slots. The enlarged core still belongs
to the inherited eight-core family: choose eight of its bounded speeds and
regard its other bounded speeds as attachments. Thus every strengthened
prefix has universal completion. If no tails remain, the core is all of
{1,...,13} and t=1/14 is a witness.

For the remaining sorted tails, all greater than13, maintain the complete
safe cover J of the enlarged core and already assigned prefix. Let k be the
remaining tail count and W the largest component width. The preceding lemma
gives W>0. Use the Session122 constants

    c1=1/7, c2=2/7, c3=3/7, c4=4/7, c5=10/7.

If the smallest remaining speed is at least ceil(c_k/W), select a maximal
safe interval and invoke the polynomial-bit conductor realizer for all
remaining tails. For k=1 use the exact first safe arrival in that interval,
including equality. Otherwise the next speed is below the cutoff: intersect
its complete safe comb with J and continue. After the last finite assignment,
return the first endpoint of the nonempty common cover, even if point-only.
Every returned time is tested against every original speed, not only padded
or quotient constraints.

## Why the finite work is independent of huge input speeds

Suppose every currently assigned speed is at most M. Every component endpoint
has the form integer/(14v) for one assigned speed v (including the ambient
endpoints0 and1). Two distinct endpoints therefore differ by at least

    1/(14 M²).

Indeed their difference has denominator14uv and nonzero integer numerator.
The positive-width lemma gives W>=1/(14M²). Hence any next speed that is
expanded into a finite comb satisfies

    v < ceil(c_k/W) <= 14 c_k M²,

The ceiling bound holds because 14c_k M² is an integer: 14c_k equals
2,4,6,8,20 for k=1,...,5. Starting from M=13, repeated
updates M<-max(M,14c_k M²) bound every expanded speed by an input-independent
integer after at most five steps. This conservative bound is not claimed
practical or optimal; observed small cutoffs in test cases are not promoted
to universal bounds.

Consequently every finite intersection and table-sized operation takes a
fixed finite amount of work independent of the bit lengths of the arbitrary
input tails. Large remaining tails are handled only through the polynomial-
bit conductor realizer. Sorting, comparisons and final original-speed tests
have polynomial bit cost. This proves the stated complexity. No generic
unrestricted n-runner polynomial-time solver is inferred.

The algorithm is a concrete witness-preserving compilation of the earlier
universal theorem. In particular it does not infer the theorem from a finite
sample of successful inputs, nor replace exact common-time composition by
nonempty marginal fibres.
