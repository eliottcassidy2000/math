# A two-port birth law and a Boolean decoder from inactive sources

**PROVED:** conditional-expectation lemma and elementary degree-two guard.
**FINITE-EXACT:** coefficients and constructive witnesses for G25 to G26.
This extends the inherited97-state interface without claiming a smaller
universal attachment carrier.

G_n is the square-sum graph on1,...,n, with distinct vertices adjacent
when their sum is an integer square. For an orientation of its edges,
H_G counts directed Hamiltonian paths using those edges. It can be zero
and need not be odd. Completing the other arcs arbitrarily gives an
actual tournament, but its ordinary Hamiltonian count is a different
function.

The old graph has32 edges and20 directed Hamiltonian paths. Adding26
adds just the two edges (10,26),(23,26), producing34 edges and24 directed
Hamiltonian paths. A bit1 on either new edge points toward26.

The inherited [native birth interface](../eleven_wedge_honeycomb_20261102/NATIVE_26.md)
has97 pre-birth states:38 extending old active sources,58 detecting old
inactive sources that can become active, and one remaining inactive
class. Its two clocks distinguish orders25 and26; the target observer
has71 states. This packet independently reconstructs the literal
orientation sets, rather than estimating their sizes from those states.

## Exact two-port distribution

Choose the32 old orientation bits uniformly. Let the joint distribution
of the two new bits be independent of that old source, but allow the
new bits to be correlated with each other. Write
\[
 r=\Pr(b_{10}\ne b_{23}).
\]
For independent bits of probabilities p,q this is p+q-2pq.
The entire new guarded-count law is
\[
\begin{aligned}
 \Pr(H_G=1)&=(512+4096r)/2^{32},\\
 \Pr(H_G=2)&=464r/2^{32},\\
 \Pr(H_G=3)&=32r/2^{32},\\
 \Pr(H_G=0)&=1-(512+4592r)/2^{32}.
\end{aligned}
\]
Consequently
\[
 \mathbb EH_G=(1+10r)/2^{23},\qquad
 \mathbb EH_G^2=(16+195r)/2^{27}.
\]
Under this specified uniform-old-source law, r is necessary and sufficient
to determine the output distribution, since the mean depends injectively
on r. It is not a sufficient selector for arbitrary old source laws or
later operations. Coupling the new bits to the old source invalidates
the displayed marginal reduction unless that joint coupling is retained.

The two laws with marginals p=q=1/2 but with the new bits always equal or
always different have means1/2^23 and11/2^23. Thus separate port marginals
discard a factor-eleven effect in this actual graph.

### Where the9,184 sources enter

There are4,864 active old orientations. Among the11,232 active new
orientations,2,048 extend an active old source and9,184 extend an inactive
one. Each of the latter has a different old parent and exactly one
successful new-bit assignment. The complete split is:

| New bits | Old inactive sources with H1 | with H2 | with H3 |
| --- | ---: | ---: | ---: |
| 00 | 0 | 0 | 0 |
| 01 | 4240 | 320 | 32 |
| 10 | 4240 | 320 | 32 |
| 11 | 0 | 0 | 0 |

The two middle rows each sum to4592. Hence
\[
 \Pr(\text{old inactive and new active})=287r/2^{28}.
\]
At fair independent bits this recovers9184/2^34.

The guard has an elementary explanation valid for every degree-two
vertex birth. If the new vertex is an endpoint of a directed Hamiltonian
path, deleting it gives an old directed Hamiltonian path. Therefore a
path born from an inactive old source must use the new vertex internally,
requiring one incoming and one outgoing edge. This proves the necessity
of unequal bits. The counts and existence for each unequal assignment
are specific to G25/G26, not a general sufficiency theorem.

First moments are independently counted by the individual path cylinders.
Second moments are counted by ordered compatible pairs of path cylinders,
using their actual shared edge orientations. This checks the coefficients
without the completion-enumeration algorithm that produced the table.

## Constructing a source instead of averaging it away

**Conditional clamping lemma, PROVED.** Let a finite positive probability
law be supported on actual allowed sources, and let F be integer-valued.
If E[(F-h)^2]<1 for an integer h, successive conditioning on coordinates,
choosing a positive-mass branch of no greater conditional energy, ends at
an actual source with F=h. At every step the parent energy is the
population-weighted mean of its conditional energies. No product-law
assumption is needed. Fixed port assignments and any constraint on the
source support are retained.

Apply this to the uniform distribution on **inactive old sources**,
whose population is D=2^32-4864, with either unequal new-bit assignment
held fixed. The energy for target H_G=1 is
\[
 1+\frac{5808-2\cdot4976}{D}
   =1-\frac{4144}{D}
   =\frac{268434893}{268435152}<1.
\]
For a partial old-bit assignment, its population is
2^(number of free bits) minus the number of compatible old active sources.
Its energy correction is the sum of H_G^2-2H_G over compatible
old-inactive/new-active sources. Thus the decoder never enumerates the
huge inactive cube and never discards its mass.

The32 successive conditionings return:

| New-bit integer | Old mask | New mask | Old H_G | New H_G |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 3927609650 | 8222576946 | 0 | 1 |
| 2 | 2449821388 | 11039755980 | 0 | 1 |

For the first source, the unique new path is

    18,7,9,16,20,5,4,21,15,1,8,17,19,6,10,26,23,2,14,22,3,13,12,24,25,11.

Every adjacent sum is an actual square, and every old edge bit is
unchanged. The receipt gives the second path and all conditional masses
and energies for both runs. Assign any orientations to the291 non-guard
arcs to obtain a full26-vertex tournament with exactly this one
square-sum-compatible Hamiltonian path.

This is a genuine finite Boolean inverse under source and port
constraints. It does not realize ordinary H=1 at order26, prove the
unrestricted odd spectrum, or supply a uniform all-size inverse.

