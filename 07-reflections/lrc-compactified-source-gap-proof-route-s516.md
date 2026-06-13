# LRC Compactified Source-Gap Proof Route (S516)

THM-382 and THM-383 point to the same correction from opposite sides.
THM-382 says the raw tournament quotient is too coarse until the LRC threshold
is placed into the fiber.  THM-383 says the open circular menu is too thin
until equality walls are compactified.  S516 turns that into a small theorem:
once the observer is marked, LRC safety is just the closed threshold color of
the two adjacent gaps.

That is not the Lonely Runner proof, but it is a cleaner target.  We no longer
need to ask whether a complicated raw A000568 class is safe.  The target fiber
is the single local condition

```text
observer-adjacent gaps = (long, long).
```

The global problem is now forced visitation.  A counterexample would be a
compactified arithmetic walk that never enters `(long,long)`.  That is a much
sharper object than the earlier "bad marked tournament class" language.

The audit script supports the split.  In bounded total `n=3,4,5,6` samples,
the adjacent-gap criterion has zero mismatches and the good fiber is exactly
`(1,1)`.  Initial consecutive systems are wall-only source hits, matching the
old tight-witness story.  The hard n14/n16/n18 ladder rows are not
source-gap-avoiding: they have short open corridors, and the corridor counts
double under gate/double-gate scaling.

The proof route now seems to be:

1. Use THM-384 to replace "lonely" by "visited source-gap fiber."
2. Keep THM-383 boundary compactification.
3. Use THM-369 to eliminate non-sieve-complete systems.
4. For a remaining source-gap-avoiding walk, build the endpoint core.
5. Use THM-380 to demand an owner-compatible pressure cycle.
6. Prove compactified source-gap avoidance makes that labelled pressure graph
   peelable.

That last step is the actual mountain.  But it is a better mountain than the
previous one: it is labelled, local, and has a concrete contradiction target.
The S514 result already warned that coarse pressure SCCs are too lossy.  S516
says exactly where the labels have to live: in the endpoint core over the
compactified source-gap walk.
