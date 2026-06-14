# Tournament trace speedups and overlap corrections

The efficient-computation question has a clean moral: trace is not merely a
faster counter. Trace tells us when the structure first needs memory.

For tournaments, `tr(A^k)/k` counts directed `k`-cycles perfectly for
`k=3,4,5`. The reason is tiny but decisive. A non-simple closed walk would
split into two closed directed walks, and tournaments have no directed cycles
of length `1` or `2`. The first possible split has length `3+3=6`.

So `c5=tr(A^5)/5` is not a lucky identity. It is the last point before
overlap enters.

At `k=6`, the first correction is exactly the two-triangle overlap layer:

```text
tr(A^6) = 6*c6 + 3*c3 + 6*p33_meet.
```

The term `3*c3` is a triangle traversed twice. The term `p33_meet` counts
distinct directed triangle pairs with nonempty intersection. Each such pair
gives a length-six figure-eight, and trace sees its rotations.

The tempting shortcut `sum_v tau_v^2` fails because it only counts walks that
return to their chosen start at time `3`. Rotating the same figure-eight
usually starts inside one lobe, so the midpoint is not the start. This is the
whole lesson in miniature: a scalar-looking trace correction secretly asks
for placement data.

The `n=6` exhaustive information audit then gave a nice surprise. The score
sequence alone leaves many Hamiltonian-path buckets mixed. Adding `c5` helps.
Adding `c6` almost closes the ambiguity. And the full low cycle vector

```text
(c3,c4,c5,c6)
```

determines `H` for all labelled `n=6` tournaments in the audit.

The rebase brought in the sharper S5 result, THM-499:

```text
H = 1 + 2(c3+c5) + 4D
```

for `n<=6`, where `D` is the number of vertex-disjoint directed-triangle
pairs. That makes the pattern cleaner. My `c6` correction uses
`p33_meet`, the number of intersecting triangle pairs, and

```text
p33_meet = C(c3,2) - D.
```

So `c6` was not magically approximating `H`; it was carrying the complement
of the first OCF conflict ingredient. The trace boundary and the conflict-graph
boundary are the same wall seen from opposite sides.

I would not overstate that as a theorem beyond `n=6`. It is better as a
pattern: the first corrected trace vector captures exactly the first
cycle-placement obstruction, and that obstruction is also the missing
ingredient for `H` at six vertices.

The next good computation is therefore clear. Do not brute force `c7` or `H`
raw if a corrected trace vector can get there first. Enumerate closed-walk
support types for `tr(A^7)`, separate scalar from placement-sensitive
corrections, and then test the information tournament at `n=7`:

```text
Do (c3,c4,c5,c6,c7_corrected) buckets still determine or sharply compress H,
and which OCF alpha term is each trace correction secretly measuring?
```

That would turn the speedup catalogue into a structural hierarchy: every time
trace fails, the required correction names the next hidden geometry.
