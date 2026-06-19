# LRC14 Seven-Sector AP-Frontier Envelope

The tempting sentence was:

`S7(E)` is maximized by the consecutive block.

That sentence is almost right in the only place where it is dangerous.  The
exact scout found no AP-beater for `k=8,9,10,11` in widened primitive boxes,
and these are the rows where the sector cap has the smallest margins.  But the
sentence is false for `k=12,13`: near-AP tail shapes like
`(0,1,2,3,4,5,6,7,8,9,10,12)` beat the plain AP.

The important thing is where the failure lives.  Those high-`k` rows have huge
slack against `cap_k`; `k=13` is even tautological at this layer because
`cap_13=1`.  So the useful conjecture is not "AP always wins."  It is:

1. AP wins for `k=8..11`.
2. AP-rich tails in `k=12` stay below `6/7`.
3. Everything high-spread or relation-sparse sits far under the frontier.

This is a better proof route than the first guess.  It turns THM-532's
low-height residual into an envelope problem with the hard and easy rows
separated.  It also keeps the error honest: local gap compression is already
known false, so the desired proof is global, probably a rearrangement theorem
for sector-word coimages or a relation-height majorization with a finite
AP-rich cap check.

The lesson for LRC(14): do not ask symmetry to be universally extremal.  Ask it
to own the rows where the proof margin is small.  The rest only needs a coarse
ceiling.

Cross-links: HYP-2604, HYP-2603, THM-532, HYP-2599, HYP-2595, OPEN-Q-108.
