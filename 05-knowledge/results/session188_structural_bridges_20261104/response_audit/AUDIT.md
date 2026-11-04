# Independent audit of named quartic response

**PASS within the stated scopes.** Reviewed Session 188 `response/PROOF.md`
and `response/verify.py`, independently derived the main identities, and ran a
separate literal-query verifier normally and under `python -O`. Both modes give
identical PASS receipts.

## Mathematical findings

* The two-shared-name quartic is zero: its required port incidences would make
  an undirected four-cycle inside a Hamilton path.
* For one shared name i, the only blocks are j,x,i,y,k and k,y,i,x,j. Deleting
  i,x,y gives exactly the two stated nonempty path covers. In particular the
  independence from edges incident to i is correct, including singleton paths.
* The disjoint-pair formula counts the eight orientation/order classes with
  the two boundary edges deliberately untested. The positive quartic sign is
  correct because two outgoing port factors contribute two minus signs.
* For a twin edge u->v, the surviving block j,y,u,v,x,i gives C_ji, while the
  reversed edge gives C_ij. Thus the sign C_ji-C_ij is correct, and root/exterior
  arcs genuinely disappear rather than being averaged out.
* On a transitive m-name exterior, C_first,last=0 for m>=3 and
  C_last,first=2^(m-2). This proves the displayed coefficient for every m>=3
  and every root-incidence pattern. The claimed strong construction contains
  the needed directed cycle in the quotient and remains strong after splitting
  the root in either twin orientation.
* The universal old-order threshold is exact. For n<=4 every quartic is
  orientation-independent: zero for double overlap, two for single overlap,
  four for disjoint pairs at n=4. The printed strong order-five pair has the
  stated adjacency masks and coefficient two.

## Independent computation

`verify.py` imports no code from the author's packet. It obtains each
coefficient by sixteen actual attachment queries, using a memoized Hamilton
path recursion. Separately, it counts the required covers by permutations and
cuts of the old names. This is distinct from the author's full polynomial
expansion and endpoint-cover dynamic program.

Normal and optimized receipts contain:

* 97 named tournaments: every tournament through old order four, 16 retained
  order-five controls, four order-six controls, and the two strong hostiles.
* 5,078 complete named quartic checks: 3,144 single-overlap, 1,284 disjoint,
  and 650 double-overlap zeros.
* 16,728 literal attachment queries.
* 63 independent twin-family checks: all root-incidence patterns for exterior
  sizes two through five, plus strong examples of exterior sizes six, seven,
  and eight. Size two correctly gives coefficient zero.

The author's own exhaustive polynomial census through old order five was
inspected; it is broader than this independent census. Its stated full truth
table check is specific to the retained strong hostile, not a universal full
response reconstruction claim.

## Boundaries

The unbounded claims are supported by the path-block bijections, not by the
finite census. Neither the original packet nor this audit proves an exact
full-response collision, universal named reconstruction, or an exclusion of
all hidden prime-core converses. The theorem does not claim this twin quartic
family captures every possible quartic obstruction.

Reproduce with `python verify.py > normal.json` and
`python -O verify.py > optimized.json`, then compare the two receipts.
