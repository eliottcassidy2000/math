# Independent audit of the multi-exit two-birth reset

**PASS — 2026-11-04.** The coordinator independently derived the full path
and old-root endpoint polynomials, checked all injection types, read the
standalone source, and replayed its exact controls. This accepts
[MULTIPORT_RESET.md](MULTIPORT_RESET.md) in its stated all-tournament scope.

The x coefficient is (P−E)(1+A), since the old predecessor must lie outside
B. The joint xy coefficient is P+2PA−EA. At an old root, the corresponding
coefficients are (P−E)Q and PQ. Extracting WP gives both negative cross-port
corrections, with no missing factor of2. The start polynomials of x and y
give gaps [EA] and H−e_B respectively. Empty B and root membership in B are
covered by the same expressions; repeated singleton endpoints recover C2.

The scalar injection has exactly |B|+2 distinct outputs per old path:
one ending x, one ending y, and |B| ending at an old vertex. The latter
family uses the yx block before every B position, not just the first. This
proves the all-order bound and explains why it still applies when F is
allowed to change. H>=5 with at least two exits overshoots both requested
nonzero ternary digits. This is the decisive stopping mechanism.

For the rooted bound, deleting the newborns recovers both ordered old slots.
The three selected lifts are distinct, and a cover meeting B in both slots
has two additional lifts in the other slot. A shortest set-B-to-root path
contains only one B vertex; the complementary Hamilton path contains another
when |B|>=2. If no B vertex reaches the root, the SCC cut has all B after
the root's component, and both H and F factor through the same suffix
Hamilton count. F tripling then forces zero scalar offset. Thus the exact
positive-tripling reduction to singleton B covers weak inputs too.

The small fork is a real boundary, not a counterexample to H>=5: the
multi-exit step starts at H3. Its H13 output has newborn gaps2 and1, and
its literal H43/H41 children retain the stated old-root responses56 and58.
Their reset-menu stopping follows from the inherited all-auxiliary C2
image theorem plus the new multi-exit bound. No other native operation is
excluded by that statement.

The coordinator replay of
[`algebra_multiport_reset.py`](../../../../04-computation/nonzero_preparation_20261103/algebra_multiport_reset.py)
passes **296,258 always-active gates** and matches the frozen
[`receipt`](../../../../04-computation/nonzero_preparation_20261103/algebra_multiport_reset.json)
byte for byte. Producer normal and optimized modes also agree. The source
enumerates actual path objects on both sides through input order4, every
exit subset and retained root, and checks all literal lift images. Its
strong-input scalar probe through5 is complete in its stated universe and
is not used instead of the all-order injection. Source SHA256 is
`b91ca738b48b3fe1aebb3cd11b15b3b52a7ba9b9cd6d5fbf0b6edf1e7d468d08`;
receipt SHA256 is
`e3cdcb0cf11d288cd0361eafe0d8913006b66226286c4401cb90b7044fa4a30a`.

Reproduction commands, the exact finite universes and hostile controls are
in the primary note. Its audited prepromotion SHA256 is
`45e9d914f9a9e061ed4fa952b9d4afa2d5f68012169cdf3034809b7ab24d49ca`.
Status/link-only promotion is permitted. No global spectrum or unbounded
positive reset claim follows.
