# Tight intervals, exact correction compilers, and recoverable tournament interfaces

**Session279, 12 November 2026. Chan–Pak remains OPEN.**

The main advance is not another large census. It is a more useful division of
work between **exact resummation, finite correction, and discrete recovery**.
For A000568 the interval can now collapse to one integer through a structurally
proved range. For Chan–Pak, a sufficiently long, specified physical suffix can
make a scalar count determine the root data that composition had seemed to lose.
The two results concern different observables and are not a proof of coverage.

## 1. A000568: remove the common error before approximating what remains

Session278 replaced the difficult gcd interaction of permutation cycle lengths
by a min-kernel envelope, corrected by a signed interaction expansion. We now
integrate any chosen finite backbone of cycle lengths exactly. Then, within
each residual sector, subtract its minimum positive defect from **every** defect.
This exact ground-state extraction removes many edges altogether. Four moments
of the remaining edges give an optimal scalar moment-positivity interval using
explicit rational polynomials, without a floating LP solver.

The certificate needs the actual minimum defect and remaining-edge count as
sidecars. It is not squeezing more information out of the identical old scalar
summary. At order64, sufficient independently supplied residue bits become:

| Exact backbone | Previous signed interval | Ground-state moment interval |
|---|---:|---:|
| {3} | 738 | 362 |
| {3,5} | 493 | 166 |
| {3,5,7} | 310 | 160 |
| {3,5,7,9,11} | 121 | **0** |

Zero means the interval itself determines the integer, not that the answer is
small: A000568(64) has1721 bits. Saved no-reference consumers read only the
moment packets, never a supplied exact answer.

There is an all-size structural explanation for the last row: with backbone
through11, the ground-state degree-four certificate is rationally exact for
**every n<=69**. Its first failure is at70, from the four-cycle core
`(13,15,17,25)`, where five nonzero corrected edges survive. At70 the certificate
requires31 residue bits. This is a proved boundary of this certificate, not
an extrapolation from a lucky order64 experiment.

A complementary impossibility theorem prevents overclaiming: any fixed finite
backbone and fixed moment depth eventually have unbounded absolute ambiguity
in the moment-positivity relaxation. Singleton padding amplifies an explicit
moment-invisible parity perturbation. Additional arithmetic information can
still distinguish it; the theorem is not an impossibility result for all
counting algorithms.

[Proof, packets and exact consumer](cross/PROOF.md).

## 2. Low correction coefficients no longer require unrestricted partitions

The new finite moment-jet recurrence retains weighted occupation moments with

    correction order + occupation degree <= K.

Inserting a new cycle length spends precisely this degree budget. The hierarchy
closes exactly, so every fixed number K of correction coefficients can be
computed in polynomially many arithmetic operations for a fixed finite
backbone. This is a theorem about fixed parameters, not polynomial-time exact
A000568 at unbounded accuracy.

Two independent representations explain the same closure:

- The **moment jet** propagates the finite degree budget in one length sweep.
- The **marked graph compiler** exposes the at most2K cycle instances touched
  by K correction edges, then sums every unmarked cycle through an external-
  field size recurrence. It retains repeated-length multiplicities exactly.

First- and second-order jet implementations are provided. A general marked-
graph prototype also computes C0 through C6 at n64. Since at most four residual
cycles remain outside the backbone through11, six is the maximum edge count.
Their alternating sum is a complete no-reference exact producer, independently
agreeing with the ground-moment consumer. This route computes higher coefficients;
it is not the same workload as the degree-four stratified certificate.

There is also a practical implementation improvement: factorial normalization
and a fixed dyadic scale compile the first jet into integer arithmetic. In one
matched C0/C1 benchmark at n80, it takes about0.56s versus3.49s for the direct
partition coefficient calculation. These are machine-specific low-coefficient
timings, **not a speedup claim over a full optimized exact tournament counter**.

Finally, the primitive-binary-word/Witt divisor identity rewrites the first
correction response into divisor coordinates. Only proper odd gcd indices up
to n/3 are needed between distinct lengths. This gives a concrete connection
between the repository's necklace work and tournament symmetry correction,
not merely a shared generating-function vocabulary.

[Jet, integer compilation and divisor response](enumeration/PROOF.md).
[Marked-graph theorem and complete order64 producer](marked/PROOF.md).

**Integration boundary:** the efficient fixed jet computes unstratified C_j.
The ground-state consumer additionally needs strata `(M,delta,R)`. Its current
producer still enumerates partitions. Efficient generation of those richer
strata is an open next step; the two components are not silently presented as
one already optimized pipeline.

## 3. Chan–Pak: a count can recover its missing interface after a Pisot suffix

The physical rooted state is `(A,W,d)`, where A is the Hamiltonian-path count
and W,d retain endpoint response. Appending a specified triangle block acts by

    M = [[3,3,-1],[3,4,0],[1,0,0]].

This operation is inherited and produces an actual tournament. The new view
uses the exact image of the physical cone `0<=d<=A<=W`. Conditioning that cone
on a prescribed count A=N gives simultaneous rational intervals for W and d.

The matrix has one expanding Pisot eigenvalue and two conjugates of modulus
below one. Its transverse errors therefore contract in **absolute**, not merely
relative, size. An exterior-square bound turns this into an effective theorem:
for **every integer source packet with 1<=A<=W<=96**, after **every k>=25**
triangle blocks the exact count N and depth k uniquely determine both missing
integer coordinates. Four rational bounds plus a four-step contraction prove
all later depths. This is not a finite numerical extrapolation.

An exact inverse then checks divisibility by4 at each step, cone membership,
and finally ownership by an actual physical seed. At the certified depths, the resulting target selector
is necessary and sufficient within the stated bank-and-suffix language.
Ambiguous shorter-depth cases are reported rather than silently rejected. Tests
include210 literal seed words giving84 distinct packets and near-target rejects.
The all-depth theorem covers the entire bounded integer cone, not just those84
measured packets.

The mechanism extends beyond triangles. Every active one-apex block has a cubic
characteristic polynomial f; its two transverse eigenvalues lie inside the unit
circle exactly when `f(1)<0` and `f(-1)<0`. This is an exact symbolic classification
criterion. All63 eligible words through length7 pass. The first failure by word
length is `10101010`, whose transverse eigenvalue lies between1 and2. Its ratios
can converge while the absolute error grows. This is precisely the distinction
an integer-recovery theory must retain.

A fixed finite seed bank with only one repeated suffix remains sparse, producing
at most O(log X) counts below X. It does **not** cover all odd targets. The next
constructive question is whether variable physical digit choices can supply
all targets while preserving the decoding budget. Individual Pisot blocks do
not guarantee contraction under arbitrary switching.

[All-depth theorem, block criterion and physical tests](physical/PROOF.md).

## 4. What the combined viewpoints actually contribute

The useful shared principle is **restore only the information the next operation
can consume, but certify its absolute scale**.

- A degree-K signed correction consumes a finite occupation jet or a graph on
  at most2K marked instances.
- A ground shift removes a common interaction before moment truncation, but
  its minimum and exceptional support must remain visible.
- A physical suffix can contract the missing interface below one lattice unit;
  its inverse still needs radix congruences and an actual seed.

Old LRC factorial-moment duality, Boolean-parity blindness, Witt/necklace divisor
transforms, and Pisot/carry observations each supply a specific lawful piece.
They do not identify the problems with one another. The new connections work
because their maps, lost coordinates, and restoration conditions are explicit.

The remaining central obligations are now sharper: an efficient stratified
moment producer for the tightest A000568 bounds, and a realizable overlapping
absolute-target digit supply compatible with these physical decoders for
Chan–Pak. Neither is replaced by numerical evidence of coverage.

[Independent new-claim review and reproduction receipts](audit/REVIEW.md).
