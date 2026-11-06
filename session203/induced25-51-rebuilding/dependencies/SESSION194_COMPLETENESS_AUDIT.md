# Independent audit of the H41/H67 weak-parent certificate

**PASS**, 2026-11-05. No blocking mathematical or computational correction is
needed for the stated strong-child, one-vertex, weak-parent theorem.

## Frozen inputs

The reviewed certificate source hashes are:

- `PROOF.md`:
  `fd0bd87c680b0a90eec9609ec5ddb6d939e2a7bc5debffa34bd924ad3de34e43`
- `verify_completeness.cpp`:
  `270d3f9e75253f196177dadfabe9e8e530b5d0b9fc216e6378a2926703f45648`
- `verify_profiles.py`:
  `79d0fff3a2bdd844ae827a764118aa55b2af48eb3597915721e0fc9db495ec26`
- `inherited_tournaments_through8.txt`:
  `5d5a5f625b4a85e04c0744aeb847ed36f59456f6a39996479aa82b9eb950afe3`

The companion exact weak-parent theorem was independently reviewed at
`STRONGER_RATIO.md` hash
`20a2c73af6bc17350e5a0b7daad526cf9b0ce35f08da435278da9f55231e0985`.
Its current promoted hash is
`ef6407da7ccabcb3a305ef22d045613542cac0a0deaefc7a1283726857bed38f`.
Replacing only its status paragraph by the original reviewed paragraph
reproduces the earlier hash exactly, confirming this is an editorial-only
promotion. Both versions have identical mathematical statements and proofs.

## Why the finite certificate really is complete

The C++ canonicalization is a full graph isomorphism canonicalization, not a
score-sequence proxy. It fixes increasing outdegree-class order and exhausts
all permutations inside each class. Isomorphisms preserve outdegrees, so this
includes every possible isomorphism after class ordering. Equal encoded full
adjacency matrices imply an actual vertex permutation isomorphism. The code
starts each class permutation in sorted vertex order and uses all of them.
Orders are kept in separate sets. At order eight the complete eight-by-eight
matrix fits exactly in 64 unsigned bits; each shift is by eight, not by 64.

Starting from the empty graph, every tournament of the next order is an
extension of some representative of its induced deletion. Exhausting all
attachment masks therefore generates every next-order isomorphism class.
Equality with the supplied next bank, plus the no-duplicate check, proves
completeness inductively. The observed class counts are not used as a
substitute for this upward-extension argument. In total, 62,423 extensions
are checked, with equality at each order one through eight.

The profile script reads a retained completeness receipt. That receipt is
valid evidence here because this audit independently regenerated it on the
exact hash-bound bank and compared the output byte-for-byte. Merely supplying
an unaudited PASS JSON would not establish the completeness claim.

## Endpoint responses and the all-order implication

The subset recurrences correctly count paths in both endpoint directions.
The literal one-break method is also exact: a Hamilton word admits one
two-bin split for each marked endpoint, while a word with one bad adjacency
can be split only at that adjacency. More than one bad adjacency cannot form
two directed path bins.

There are six H41 strong core classes, at orders six and seven, and six H67
classes, all at order seven. This is 83 rooted profiles. Every R and L agrees
with literal counting, and none has surplus one or two. The H5 positive
control does have both surpluses. All four excluded-child target counts also
have verified strong witnesses, so the exclusion is not a support-absence
claim.

The cited monotone Busch minimum bounds every H41 strong core by order seven
and every H67 strong core by order eight. Thus no larger core is omitted.
The separately proved exact weak-parent theorem eliminates arbitrary extra
singleton padding and reduces every candidate to precisely the tested rooted
surplus conditions. Its two exceptional presentations have parent counts
three and nine, so neither applies at 41 or 67. This is a valid finite-to-all-
order deduction; it does not extrapolate the finite bank to arbitrary orders.

The restrictions in the statement are load-bearing: the final child is
strong, the parent is weak, exactly one vertex is added, and old arcs are
unchanged. The certificate does not rule out a strong parent, weak child,
multiple births, changed old arcs, or an unrelated realization of the target.
The prime-atom interpretation for 127 and 203 is correctly scoped to these
particular literal ternary ancestry steps.

## Replays and independent controls

Both `-O3` and `-O0` C++ builds passed and produced identical completeness
receipts, matching the retained producer receipt. Both ordinary and optimized
Python profile runs passed and agreed byte-for-byte with each other and the
retained producer receipt. The four source hashes above were rechecked after
completion and were unchanged.

An additional audit checker, `verify_literal_bins.py`, imports no producer
code. It tests every permutation and every explicit two-bin split for all
twelve core classes, checking adjacency separately inside both bins. Its
10,156 valid ordered covers reproduce all 83 rooted profiles and the excluded
surpluses. Ordinary and optimized outputs also agree.

The retained audit receipts are `completeness.O3.json`,
`completeness.O0.json`, `profiles.normal.json`, `profiles.optimized.json`,
`literal_bins.json`, and `literal_bins.optimized.json`.

The release must retain the companion exact weak-parent theorem and its
dependencies, including the proved single-root-triangle minimum and the
CITED Busch bound. The unrestricted `{7,21}` spectrum remains OPEN.
