# Independent hostile audit: the source-time P169 high branch is nonredundant

**Status: FINITE-EXACT INDEPENDENT AUDIT ACCEPTS THE SCOPED SOURCE-TIME
INVERSE-BRANCH SIDECAR.  LRC(14) remains OPEN.**  The submitted branch-sheet
implementation was not imported.  The audit reconstructs all thirteen branch
profiles directly from inverse images of the source `E` set, then applies `Q`,
the `P_(13^5)` transfer, and the independently audited common-owner endpoint
coupling.

## Exact coordinate and temporal certificate

Write the left source packet branch as

```text
n=13b+r,  0<=n<169,
b=floor(n/13).
```

This `b` is the high digit of the source-time `P_169` inverse branch.  On the
`91*169=15379` exact temporal cells, direct midpoint arithmetic gives

```text
source guard sheet = floor(n/13)=b.
```

Every digit has exactly `1183=91*13` cells.  This is an equality of two
source-time descriptions on the same cell decomposition, not an identification
by cardinality.

The branch profiles are rebuilt by taking every exact inverse image

```text
y -> (y+n)/169
```

of `E`, grouping only after the individual branch is fixed, intersecting with
`Q`, and then applying `P_(13^5)`.  Summing the thirteen nonnegative profiles
recovers every audited collision-root profile pointwise.  The branch-profile
digest agrees with the candidate.

After `Q`, branches `0` and `12` vanish and branches `1..11` survive.  Every
surviving branch meets each owner-visible Boolean state in exactly one
connected component.  Thus support incidence is eleven identical `1111`
rows bracketed by zero rows; support alone cannot explain the response rank.

## Common-base and marginal gates

On each exact endpoint segment, the audit first decomposes the selected left
source value by branch and only then multiplies by the selected right source
value and endpoint harmonic jump.  Three literal guard controls pass, and the
same-root sector remains zero before integration.

Summing `b` recovers every character row and every inverse-table entry of the
audited Boolean square.  This is the necessary Fubini/normalization gate: the
sidecar refines an existing one-base integrand rather than fitting an unrelated
branch table to its final marginal.

The ranks are

```text
Boolean square against relation:     4,
(branch,state) against relation:      6,
branch marginal against relation:     3.
```

The strict `4 -> 6` jump establishes nonredundancy despite the flat geometric
support matrix.  At fixed relation `(1,0,6)`, the `13 x 4` matrix has rank `3`
and exactly `44=11*4` nonzero entries.

After centering all three axes, the pure interaction has rank `6`; all
`12*3*12=432` possible branch/state/relation modes are nonzero.  As in the
earlier audits, Fourier fullness is interpreted only after marginal recovery
and rank increase.  All eight submitted profile/gamma/table/spectrum digests
agree.  Normal and optimized runs are byte-identical; the independent semantic
digest is `6cac7c90f3ebe3f33a27c1979b74006c545260bb6d41a68f4f76442456114fe6`.

## Three distinct thirteen-valued coordinates

This audit sharpens the coordinate ledger:

```text
b_source = floor(n/13): high digit in source-time P_169 before Q,
s        = u-q:         difference of two collision roots,
r_owner  = a mod 13:    last current-leg inverse digit in P_(13^5).
```

They live at different stages and have different maps.  The present sidecar
marginalizes both `s` and the pointed source tail; the incoming `r_owner`
candidate separately marginalizes `b_source`.  Two rank-six results on these
different quotients do not prove the coordinates equivalent or jointly
independent.

The natural next hostile is a joint table retaining `b_source` with one of the
already audited source coordinates, or retaining `b_source` and `r_owner` on
the full Boolean stalk.  Its first gate must recover both one-coordinate
marginals exactly.  A formal product made after either marginal would not be a
lawful joint ancestry object.

## Connection contract

| field | exact answer |
|---|---|
| source | high digit of the left source `P_169 E` inverse branch |
| target | `F_13(b_source) x V_4(state) x F_13(relation)` |
| map | retain `floor(n/13)` before `Q`, `P_(13^5)`, endpoint multiplication, and integration |
| preserved | source-time sheet, collision-root profiles, common base, Boolean state, literal guards, endpoint phase |
| destroyed | low source digit, source root difference, pointed tail, current-leg owner digit, arrival atom, exact address |
| mandatory marginal | branch sum recovers the audited Boolean-square gamma and table |
| decisive hostile | flat branch/state support versus strict response-rank jump `4 -> 6` |
| survivor | rank-six pure interaction and `432/432` triple modes |

This is not the current-leg `r_owner`, a deep-root label, an arrival atom,
THM-2334 grouped exact address, `U_clock` chronology, physical current,
uniform-row theorem, row exclusion, or LRC(14).

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_source_branch_sheet_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_source_branch_sheet_independent_audit_20260816.py
```
