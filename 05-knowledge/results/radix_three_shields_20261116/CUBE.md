# Checking the four-independent-arc shortcut

If the first two width-four dials were each implemented as
two independent single-arc flips, their combined affine cube
would have positive increments 2,4,6,12. Even the first
three of those increments would have to coexist as an affine
three-arc cube. Exact subset-DP census through order eight
finds no such cube with the first three prescribed increments.
This checks a stricter requirement than a four-setting
control: a setting may change several arcs, as the successful
two-dial bank does.

The checker scores all labelled tournaments at order six and
seven, and all 6,880 representatives of the archived order-eight
atlas. Each cube has a unique minimum labelled corner among
its settings, with the three respective single-arc successors
at H+2,H+4,H+6. It checks the three two-flip states and one
three-flip state for exact additivity. Relabeling a candidate
to its representative preserves all these tests. Orders below
six are covered by neutral source padding to order six.
It also tests single-arc paths of count increments +2 for
both four-term windows 65…71 and 67…73 around 68: neither
occurs at these orders. This does not refute the literal
three-term 67→69→71 path, or a four-setting construction
made using different arc groups or larger order.

Run search_cube.cpp at arguments 6 and 7. For order eight,
pass the authenticated n8_classes.csv from
[the multiplicity census](../actual_hamilton_multiplicities_20261021/session24-multiplicities/)
to search_cube8.cpp. These are finite negative tests, not
a size independent cube theorem.
