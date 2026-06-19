# LRC Below 14 Modern Reproof Ladder

The below-14 reproof attempt gave a cleaner lesson than I expected.

The AP tight rows are not a nuisance to eliminate.  They are the spine.  For
`N=3..13`, `(1,...,N-1)` has exact `M=1/N` and zero strict safe measure at the
target.  In the one-step AP frontier, AP stays the unique tight row.  But as
soon as the target relaxes from `1/N` to `1/(N+1)`, the strict safe measure
becomes positive.

At the top edge this reproduces a number already central to the LRC14 proof
state: `7/858`.  Among the `91` primitive 12-subsets of `[1,14]`, only
`(1,...,12)` is tight at the LRC(13) target `1/13`, and the minimum relaxed
safe measure at `1/14` is `7/858`, attained by the drop-6 core
`(1,2,3,4,5,7,8,9,10,11,12,13)`.

So the proved range suggests the right proof language for 14 is not "avoid
tightness."  It is "classify tightness, then prove a uniform collar."  The
collar is the object that signed wall transport should explain.

The support-floor ladder says why this stops being routine at 14.  For even
`N=2q`, the half-gap sector quotient has support floor `q-1`:

```text
N=4,6,8,10,12 -> floors 1,2,3,4,5
N=14          -> floor 6
```

The current hard object, THM-538/HYP-2646 support-six signed coset mass, first
exists exactly at `N=14`.  Below that, the same quotient grammar is present but
does not yet expose the support-six conditional convergence problem.

My next proof target from this is an AP-frontier fattening lemma:

```text
AP_N is the unique tight row in the first AP frontier,
and relaxation 1/N -> 1/(N+1) creates a signed wall-transport collar.
```

Then extend from the first frontier to all AP-rich rows by Freiman small-excess.
That would turn the known `N<14` proof into a rehearsal for the exact shape of
the `N=14` proof: classify the AP collar, route structured perturbations
through signed wall transport, and reserve the signed coset factorization for
the first support-six wall.
