# Proposed THM-2791 path-sheet addendum

Insert after Section 3, “An exact partial physical germ.”

## 3a. The partial germ preserves the literal rail ancestry sheet

The common rail weight in (15) is not only numerically equal at the two
cylinders.  Before Perron marginalization, write

```text
A_x = {(a,b) in Z/(13^5) x Z/169:
       (x+a)/13^5 lies in Q and
       (((x+a)/13^5)+b)/169 lies in E},

B_x = {e' in Z/(13^5):
       (((x-1/13) mod 1)+e')/13^5 lies in E}.
```

Thus the route-two rail numerator at `x` is the Boolean Cartesian sheet
`A_x x B_x`.  For every `x` in the first selected open cylinder, put
`x'=x+delta mod 1`.  Exact enumeration of the raw `Q`, `E`, and rotated-`E`
walls gives one common contributor chamber, in base-grid coordinates,

```text
[140890500190440,144190879112280).
```

Both the source cylinder and its translated image lie strictly inside this
chamber.  Consequently

```text
A_x=A_x',                    B_x=B_x'
```

by the identity on the literal labels.  Their exact cardinalities are

```text
|A_x|=966606,                |B_x|=28534,
|A_x x B_x|=27581135604.
```

The collision labels are fixed as

```text
(u,v,s,t)=(5,6,1,12),
```

so no ancestry label is lost through `(u,v,a,b,e')`.  For example the same
Boolean copy

```text
(a,b,e')=(59162,26,56658)
```

is present on both cylinders.  Hence (14)--(17) define a
**same-ancestry partial translation germ on the THM-2471 rail sheet**, not
merely an equality of aggregate rail weights.

This refinement stops before endpoint allocation.  It constructs neither a
THM-2779 endpoint origin nor an allocated THM-2625 endpoint atom, and it
does not identify either one with the preserved rail sheet.  Thus the
endpoint-origin, endpoint-current, and row-exclusion invoices of Section 8
remain open.

## Audit sidecar

The assertion-free companion

```text
.scratch/thm2791_hostile_audit/path_sheet_probe.py
```

pins every dependency used in this refinement, directly enumerates both
label factors from interval membership, and has byte-identical ordinary and
optimized transcripts.  Its label-set digest is

```text
15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd.
```
