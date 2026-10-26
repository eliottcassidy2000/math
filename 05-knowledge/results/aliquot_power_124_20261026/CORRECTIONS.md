# Two inherited scope repairs

2026-10-26. **PROVED correction of a quantifier overreach; FINITE-EXACT
counterexample to a per-class assertion.** The exact removed and replacement
text, source commit and before/after hashes are in
[CORRECTION_PREIMAGES.json](CORRECTION_PREIMAGES.json).

## Distinct-prime Goldbach was an antecedent

The [S57 report, section3](../s57_untouchables_root27_20261026/README.md)
stated without a qualifier that the only odd k whose k-1 lacks a
distinct-prime partition are3,5,7. With k>=3, this is precisely the missing
distinct-prime Goldbach assertion for even numbers>=8. Its surrounding
discussion already kept the relevant global question open.

The repair puts the conjectural antecedent back in that sentence. No
counterexample to Goldbach or to the claim that5 is the only odd untouchable
is asserted. The existing finite census remains finite evidence. The failure
mechanism was turning a sufficient conditional construction
`k-1=p+q, p!=q => s(pq)=k` into an unconditional global coverage statement.
Repeated-prime partitions do not repair it: `s(9)=4`, not7. Nor does the lack
of a distinct-prime partition exclude other preimage shapes. The next
question is whether such shapes can yield a coverage theorem independent of
the unproved antecedent; this packet supplies no such theorem.

## A minimum over a fibre is not a bound on every member

The finite-record bullet in [THM-4748 / two-reversed-arcs-exact-path-counts](../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md)
could be read as saying every class with H in{31,63,127,255} has minimum
feedback-arc number3. The explicit seven-vertex tournament with out-row masks

    (106,60,89,48,97,68,10)

has H=127 but minimum feedback number5. Literal enumeration of all7! orders
checks both statistics independently of the subset path DP; its
odd-cycle profile is `(1,45,9)`, so `H=1+2·45+4·9=127` as a third route.
The integrating root also independently enumerated all5040 orders.
This is a concrete witness; no claim that it is minimal in every possible
parameter is needed.

The corrected sentence takes the minimum over all realizations of each
specified H value. The theorem's fan construction attains3, and its lower
bound excludes smaller values for the stated Mersenne parameters. That
result survives. The failed implication exchanged a fibrewise existence/
minimum statement for a universal assertion about its elements. Retain the
actual tournament as a sidecar to H. A natural next target is the full
feedback-number spectrum inside a fixed-H fibre, not a scalar identification.

Both repairs are recorded near the top of the maintained
[mistakes ledger](../../../01-canon/MISTAKES.md); its header remains at byte0.
