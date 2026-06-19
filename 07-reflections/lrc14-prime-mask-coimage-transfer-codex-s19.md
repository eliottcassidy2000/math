# LRC14 Prime-Mask / Coimage Transfer - Codex S19

The user's recurrence prompt seems right, but not in the naive way.

There is no clean recurrence in runner count hiding here.  The stable object is
a transfer between finite quotients.  THM-539 sees the prime masks of `k-1`:
the killer speed `a(k-1)` annihilates all coarse clocks whose denominators
divide `k-1`.  After MISTAKE-079, this should be read as a mask mechanism,
not as the retracted unbounded-level claim: the `a=3,4` dips survive, while
the natural `a=5` primorial continuation collapses to the floor.  The LRC14
support-six residual sees a different mask: the unit seam
`(Z/14Z)^* -> F_7^*`, followed by the mod-7 coimage tail.

That seam is exact.  The units mod 14 reduce to all nonzero residues mod 7, so
HYP-2617's projective `F_7^*` quotient is literally the coimage of the LRC14
unit action.  This gives the mod-7 atlas a cleaner conceptual status.  It is
not "mod 7 because seven sectors are convenient"; it is mod 7 because the
14-runner unit group forgets exactly the factor that the support-six relation
uses.

The mask computation gave the useful surprise.  For `k=10`, allowing apex
masks `{2,3,5}` does not improve on the empty tracked mask at all: both hit
`73` classes and `72.120496%` signed mass.  Letting the apex be divisible by
`7` jumps to `85` classes and `84.229179%`.  So mod 30 was a true signal in the
spectral family, but it is not the live transfer coordinate in the fixed LRC14
coimage tail.  Here the new coordinate is the prime that becomes invisible in
the mod-7 relation.

The remaining repeated packet is also more arithmetic than it first looked.
For `(1,1,1,1,a,a)`, the largest two classes are `a=2,4`, the nontrivial
quadratic residues mod 7; the three smaller classes are `a=3,5,6`, the
nonresidues.  The `(1,1,1,1,a,b)` packet is finer, but it still collapses to a
short list of multiplicative-character signatures.  That feels like the right
shape of the next theorem: not a generic repeated-residue estimate, but a
repeated-root character-sum estimate.

This does not prove LRC14.  It does improve the proof target.  The final
support-six theorem should now be stated over a small family of repeated-root
Dedekind/cotangent sums, with cases indexed by quadratic characters over
`F_7`, after finite height-2 wall accounting.

Tournament-wise, the useful vertices are again not runners.  They are the
quotients that preserve the proof predicate:

```text
unit seam
> prime mask transfer
> height-2 wall class
> repeated tail packet
> signed Dedekind theorem
> raw support
> raw runner.
```

That quotient throws away timing, but timing is not the scarce information at
this stage.  The scarce information is which finite arithmetic packet survives
long enough to need an analytic bound.
