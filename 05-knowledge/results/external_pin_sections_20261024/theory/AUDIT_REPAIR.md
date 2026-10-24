# Independent audit of receiver repair and additional unary pins

Status: **VERIFIED proof review and FINITE-EXACT replay.** No correction requested.
Auditor: complementary-trijection/coverage lane, 2026-10-24.

Reviewed `RECEIVER_BANKS.md`, especially Sections 3 and 7, and
`repair_labels.py` with `repair_labels_receipt.json`. This audit also uses the
independently replayed original erasure and actual native-route inputs;
see `../erasure/AUDIT.md` and `../native/AUDIT_ROUTE.md`.

## Exact repair criterion

After subtracting the fixed all-U contribution, adding the two declared
records of one compatible slot changes a grade-r support to

    A_r union (lambda + A_(r-1)).

Consequently its holes are exactly `H_r intersect (lambda + H_(r-1))`.
A surviving hole is equivalent to a literal pair `(q,p)` with
`q in H_r`, `p in H_(r-1)` and `lambda=q-p`. Thus avoiding the union of
these difference sets is necessary AND sufficient, with the boundary and
adjacent-grade indices retained. This is not a criterion for an abstract
label lacking two compatible full native records.

The deficiency bound in Section 3 is also correct: `A intersect (q-B)` has
at least `|A|+|B|-|Q|` entries, each determining one retained witness pair.
Subtracting the actual incompatible-pair count leaves the stated strict
sufficient inequality. Boundary-dependent witness compatibility cannot be
discarded or inferred merely from the sizes of the supports.

## Native certificate and completeness

The inherited erasure receipt contains all 1,140 triples of the twenty
cheap core slots. Exactly six have a hole at grade 8 or 9. The other
triples contribute no forbidden repair label at requested grades 8, 9, 10:
for grade 10, completeness at grade 9 alone already makes the difference
set empty, irrespective of the grade-10 support. The code uses this valid
adjacent-grade argument rather than silently assuming monotone coverage.

The six exact hole lists give 32 forbidden labels in
`Z/128 x (F2)^4`. The checker tests all 2,048 labels against the direct
translated-hole intersection as well as the difference-set description.
It authenticates every old core label against the original electrical
potential vectors, and both path records and the changed clocks of all
eight outside candidates. All eight labels avoid the forbidden set; five
sources have actual switching cost 3 and three have cost 4. The receipt
retains a boundary, grade and literal surviving-hole pair for every
forbidden label. No hypothetical source is supplied by the 2,016 allowed
abstract labels.

## Flexible-pin routing scope and costs

Fix any three cheap core records and any `a=0,...,7` of the eight outside
core records, each to its actual U or V record. Their occupied resources,
counts and charges remain fixed constants. A remaining outside source p
can be reassigned into the receiver without creating or duplicating a
port. The receiver then has seventeen cost-2 slots and p of cost `c=3`
or `4`; the remaining outside set has `e=7-a` slots, each costing at most 4.

All eight possible values of e are included in the exact grade-only cost
replay. For free-core population `8 <= t <= 10+e`, the receiver is complete
at every needed grade 8, 9 and 10, at every actual outside/helper frame.
Its normalization grade is `r0=max(8,t-e)`. The bounds are:

| Event | Maximum actual slot-clock occurrences |
|---|---:|
| Enter the complete band | `46+c` (safe) |
| Change receiver grade by one | `36+c` |
| Exchange two outside choices | `42+c` |
| Change the helper pair along a Johnson edge | `49+c` |

The receiver-refresh bound is `34+c`; a helper edge costs at most 15.
Same-grade receiver changes use at most 18 slots, so adding the two changed
helper/outside slots gives support at most 20. Band entry has support at
most 16, and consecutive-grade changes at most 18. Thus the interior
20-slot and `49+c`-clock claim is valid, including `e=0` where outside
events disappear.

Outside this population interval, either the free V count or the free U
count is at most 7. Two feasible endpoints differ on at most fourteen
free core slots and four free helper slots. Their direct simultaneous
event costs at most `14*4+29=85`, with support at most 18. Taking the
larger interior support bound gives the stated global **20/85** bound for
every nonempty fibre. This does not assert boundary charge coverage.
The population-complement argument bounds distance; it is not an
unguarded physical move and does not negate helper labels.

The native application has private old core owner frames and authenticated
disjoint complete resource stores. Fixing unary actual records therefore
introduces no hidden cross-factor, source capacity, target inverse-port or
converse collision. The externally fixed shared owner remains 4 and both
original population margins remain fixed at every arrival. The additional
pins here are unary record choices, not arbitrary dependent shared-owner
relations. No conclusion is made if all eight repair candidates are frozen
or if new relational pins restrict the receiver. The charge criterion
retains the authenticated paired dyadic aggregate currents, not every
ordinary per-port character or an unspecified global Hamiltonian observer.

## Reproduction

From the repository root, normal and optimized Python each passed:

```text
python 05-knowledge/results/external_pin_sections_20261024/theory/check.py
python -O 05-knowledge/results/external_pin_sections_20261024/theory/check.py
PASS 1,040,572 checks; 42554 full endpoint routes

python 05-knowledge/results/external_pin_sections_20261024/theory/repair_labels.py
python -O 05-knowledge/results/external_pin_sections_20261024/theory/repair_labels.py
PASS 2,294; 32 forbidden labels; all8 native sources repair all1140 erased triples
```

These finite controls accompany the conditional all-size receiver proof;
they do not extend the authenticated native supply beyond its stated
82-source carrier.
