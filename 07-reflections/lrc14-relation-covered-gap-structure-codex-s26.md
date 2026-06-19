# LRC14 Relation-Covered GAP Structure - Codex S26

**Source:** HYP-2639 / T887.  This is now a refinement of HYP-2637, not a
competing claim number: HYP-2637 gives the weighted relation-fiber/nullity GAP
split; this note says which pieces of that split remain visible to the pinned
LRC observer.

The user nudge was right, but it needed one more bit of typing.

```text
every element in a small relation
```

does smell like additive energy, and additive energy does smell like Freiman.
But the old S578 shifted-AP warning is decisive: AP and shifted AP have the same
additive energy, while AP is floor-tight and shifted AP is very safe.  The
energy did not change.  The location of the collision nodes changed.

So the right object is not raw energy.  It is a labelled relation hypergraph.
The labels are the proof:

```text
pair-sum node C = a+b          summand graph
C divides w?                  multiplicand graph
balanced or observer-coupled  sign/parity type
C in S or hidden above S       visibility shell
```

This also clarifies the even/odd and positive/negative language.  A balanced
relation `a+b=c+d` has coefficient sum zero; it is translation-invariant and
behaves like an even scalar shadow.  A fold `a+b=c` has coefficient sum one; it
is observer-coupled and behaves like an odd marked channel.  Midpoints
`a+c=2b` sit between them: balanced, but carrying the 2-adic midpoint.

The computation gives the useful correction.  The KPS third-pocket examples are
wide and every nonzero element is covered by a small relation motif, but their
`|S+S|` is well above the elementary `3k-4` Freiman window.  So a direct
"relation-covered implies one-dimensional GAP" theorem is probably false in
the form we need.  The replacement target is:

```text
dissociated stranger
or Freiman-small GAP
or relation-covered non-GAP with enough low visible/hidden shell payload.
```

The third case is the interesting one.  It is not tight; the sampled rows have
large margins.  The proof difficulty is not finding slack but proving signed
cancellation over a relation-rich lattice without taking the absolute envelope.
That is exactly where HYP-2634's defect sieve and HYP-2636's block-frequency
transfer should enter.

After reading the incoming S26/S12/S27 work, the picture is cleaner.  HYP-2637's
weighted fibers explain why the KPS examples are not peelable; HYP-2638 claims
the finite Freiman small-excess pocket; KPS S12's sumset-excess scan explains
why the remaining non-AP pockets are numerically loose.  This scout contributes
the typed obstruction between those facts: raw energy cannot distinguish AP
from shifted AP, but visible folds can.

The follow-up Freiman-dimension recursion makes the architecture even more
legible: dimension gives the coarse pocket and the big penalty; visibility/sign
decides which relation channels still have to be controlled in the signed tail.

I think the hidden key is now less "high additive energy" and more:

```text
energy becomes useful only after the summand node is retained.
```

Addition supplies the candidate denominator.  Multiplication tests clearance.
The sign pattern says whether the relation is a harmless balanced shadow or an
observer-coupled fold.  That is the small recursive detail that persists as the
number of runners changes.
