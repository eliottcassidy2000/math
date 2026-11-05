# Independent component-transfer and prime-amplitude addendum audit

2026-11-05. **PASS for the d=9,11,13 unit-source tables, complete core
inventory, and odd prime unit-amplitude classification through 19.**

Frozen source hashes:
- ADDENDUM.md: `6385e6425a5a0e0685a836837bd8c1169df96d7e924d4f990e16ed34dcf3def3`
- verify_addendum.py: `79ee27f074299a50a713ce908d13da54807b60d13b484ce06a07b3d00ead91fe`
- receipt: `2640f42503038bcfac774cc19ea33fd486bee545d64222154c08a0cb2c34836a`

## Component transfer

All seven local quantities have the stated physical definitions. The
four-coordinate update is exact and nonnegative. For the boundary sum,
the four cases are disjoint: both endpoints old, both new, old start/new
end, or old end/new start. Their contributions are respectively hf, Fd,
vU, and ba. In the last case every old vertex must lie in the first path
and every new vertex in the second, explaining why a is the full-prefix
endpoint mass rather than a two-cover response.

The expanded boundary formula is also valid. A reverse boundary between
nonadjacent SCCs is impossible because no path can cover an intervening
component while respecting both designated endpoints. Adjacent SCCs give
exactly the stated ending/starting response product. The empty prefix
(1,0,0,0), singleton updates, and final profile reconstruction all agree.
The beta monotonicity claim follows because every subtracted term is
nonnegative and h>=1; the unit return condition at beta=-1 is correctly
restricted to this fixed-root component-append operation.

## All-order core completeness

The named S195 strong-core theorem was read and checked. Its strict
anchored threshold is 3H(C3)+F2(C3)=15. If a strong target has H<15,
starting from any directed triangle, the canonical mixed-vertex growth
cannot stop at a proper module: the resulting two-vertex cyclic ear would
already force at least 15 paths. Thus every larger target would have a
strong induced order-six member along its singleton chain, with no larger
Hamilton count by induced monotonicity.

The independently complete order-six census has 22,320 labelled strong
tournaments and minimum H=15, not 17. Therefore all strong cores with
H<=13 have order at most five, exactly as the corrected proof states.
There is also an independent published confirmation: Busch's strong bound
gives `13^3=2197 <3125=5^5`, which directly forces order<=5.

The independently generated isomorphism inventory has respectively
1,1,0,2,1,1 strong types for H=3,5,7,9,11,13. This makes the deletion-core
banks complete: d=9 permits either one of two H9 types or two ordered C3
components; d=11 and d=13 each permit their single strong type; singleton
padding is the only other component possibility. No order-six minimum
assumption is silently used at the strict boundary H=15.

## Independent finite tables

`audit.py` obtains core graphs from the auditor's independently complete
H67 census. Local boundary counts come from explicit two-part partitions
and forward-path tables, rather than the producer's root-extension slack
calculation. The global profile uses the sequential four-coordinate
update, rather than the producer's expanded pairwise formula.

Singleton blocks are independently distributed into gaps around the
ordered core bank. All local incidence choices satisfying the outdegree
bound are enumerated, with normalized-prefix and normalized-tail tests.
All examined and matching counts agree exactly:

- d=9, beta=+1: 18,092 examined, 5 matches, bases {20,24,32,44}.
- d=9, beta=-1: 77,271 examined, 15 matches, bases
  {26,27,34,35,41,45,99}.
- d=11, beta=+1: 8,493 examined, 3 matches, bases {29,37,47}.
- d=11, beta=-1: 8,493 examined, no matches.
- d=13, either sign: 8,493 examined, no matches.

The sequential transfer also matches every one of the 5,405 rooted
labelled graphs through order five. Every saved independent witness is
checked by full subset dynamic programming, and 28 actually constructed
K extensions verify the stated counts and full profiles.

The source verifier passes normally and under `-O`, with byte-identical
output equal to the frozen receipt. The power-of-two hostile has the
verified profile (19,1,73,9), hence U=64 and f=46>d=9. No false
power-of-two shortcut is used by the proofs or the classifications.

## Prime-amplitude scope

General two-piece unit stars reduce to the singleton-partner source table
only because the amplitude is an odd prime, by the independently reviewed
physical factorization theorem. The source endpoint budget then gives
`3d<=2p+2beta-1`. At p<=19, all odd possibilities lie at d<=13. An odd
amplitude cannot have a positive dyadic normalization tail, so checking
the undoubled finite bases is complete.

The surviving cases are exactly (3,-1), (5,-1), and (11,+1), with actual
physical witnesses. The extra exclusion (23,-1) correctly uses d<=14 and
oddness to reduce to d<=13. The positive sign at p=23 is explicitly left
undecided because d=15 becomes possible.

Composite amplitudes do not receive the singleton-partner reduction, and
no global Hamilton-count exclusion is inferred from these amplitude
tables. The statement is correctly scoped to indefinitely returning
positive unit binary transitive cloning. H15 classification is separate
and is not included in this report.

No mathematical or implementation blocker was found. Results are in
`AUDIT.json` and `audit.log`. Reproduction uses `python audit.py` and the
sibling independent prime-atoms counting engine and H67 census files.
No author source file was changed.
