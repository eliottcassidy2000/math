# LRC14 AP84 Color-Grid Bridge

HYP-3470 reconnects the older HYP-2593/HYP-2594 phase-color machinery with the
current AP84 bridge, as the exact canonical AP84 `q=14V` CRT-placement sidecar
under HYP-3459's color-packet legality audit, complementary to HYP-3460's
phase-branch pullback, downstream of HYP-3461's colored-extension gate
carrier, and as the placement sibling of HYP-3458's AP84 coloring-recursion
sidecar.

The main negative result is useful: for

```text
S_m={1,2,...,11,13,84m}
```

the phase-color reservoir is too symmetric to see the HYP-3456 period-`35`
escape clock.  In HYP-2593 form, `P={1,...,11,13}`, `E={0}`, `V=84m`; color
`0` is dead and every live color has the same four intervals:

```text
[15/182,13/154], [29/70,41/98], [57/98,41/70], [141/154,167/182].
```

So mass is the wrong object.  The signal is shifted-grid discrepancy:

```text
A(m+385)-A(m)=5112
```

for the total closed CRT grid count, while individual live colors only become
uniform after period `5005`.  A smaller boundary bonus sits inside this as a
`7`-clock: the closed-open bonus is `0` exactly when `7 | m`, and `2`
otherwise.

This gives a clean division of labor:

```text
HYP-3456: component/corridor escape count, period 35.
HYP-3459: legal AP84 color-packet quotient and discrepancy warning.
HYP-3460: phase-branch color pullback from CRT witnesses to branch/gate carriers.
HYP-3461: colored-extension gate carrier that names AP84 as an endpoint-floor packet.
HYP-3458: coloring-recursion state for the period-35 AP-tail escape count.
HYP-3470: q=14V phase-color CRT placement count, period 385 total and 5005 by color.
```

The proof route should keep these sidecars separate.  If the AP-tail splice
uses branch-union components, use HYP-3454/HYP-3456/HYP-3457.  If it needs
actual `q=14V` CRT witnesses, attach the HYP-3470 color-grid sidecar and do
not scalarize the thirteen live colors to one mass.

Assumption challenge: I considered runners, phase colors, color residues,
fixed live intervals, interval endpoints, boundary hits, AP-tail `m` residues,
HYP-3456 high gaps, Fourier modes, and proof obligations.  The selected
quotient preserves exact CRT placement and color discrepancy, but forgets
branch-cover geometry and component adjacency.  The challenged assumption was
that the phase-color reservoir's mass should reveal the AP-tail clock.  It
does not; the clock is in the grid.

Next hook: prove the four fixed intervals and the floor/ceiling color count
symbolically, then carry the period-`385`/`5005` sidecar only in the AP-tail
subproofs that require actual colored CRT placement.
