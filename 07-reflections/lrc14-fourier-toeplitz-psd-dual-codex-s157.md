# LRC14 Fourier-Toeplitz PSD Dual -- codex S157

The user asked to push the Fourier idea: if the danger arcs cover, then
`C_S(t)-1>=0`, and nonnegative functions have PSD Toeplitz moment matrices.
The useful move was to keep this as a phase-sensitive dual, not just another
count-moment inequality.

S158 had just claimed HYP-2973 for the danger-count distribution, so I
renumbered this route to HYP-2974.  The distinction is clean:

```text
HYP-2973: keep the distribution of C_S(t).
HYP-2974: keep the Fourier locations of C_S(t)-1.
```

The Fourier coefficients are sparse:

```text
c_k = sum_{v|k} sin(pi*(k/v)/7)/(pi*(k/v)).
```

This is the curried-function hook.  The observable is not merely speed `v` or
mode `k`, but the relation `(v divides k)` together with the quotient `k/v`.
That makes low harmonic degree a divisor-fiber probe.  It also explains why
high replacement speeds can be invisible to early modes even while the row is
safe.

The script does not compute true first negative eigenvalues.  Instead it uses
the Fejer vector centered at the midpoint of the largest exact safe component.
That is enough: a negative quadratic form is a real Toeplitz PSD violation.

Main evidence:

```text
HYP-2963 bank rows       21913
zero-safe rows               2  (AP and GW)
positive rows            21911
Fejer PSD-vector hits     21911
misses at degree 512          0
max first degree            280
```

The hardest named row is `P10+GW` at degree `280`; `single swap 6->63` needs
degree `266`.  Many rows certify much earlier: the full-bank degree histogram
is concentrated in degrees `20..49`.

Assumption challenge: the tournament vertices should not be runners or arcs.
For this route they are proof carriers:

```text
exact nonnegative function,
Toeplitz PSD cone,
Fejer center certificate,
divisor-curried Fourier fibers,
HYP-2973 count distribution,
HYP-2972 twist ladder,
endpoint/lift packets,
raw safe interval,
raw runner set.
```

This quotient preserves the cover predicate and a concrete dual certificate.
It destroys endpoint-owner labels except for the safe midpoint used to place
the Fejer vector.  That is acceptable only because HYP-2965/HYP-2968 keep the
endpoint/lift data elsewhere.

The honest gap is formalization.  The negative values are floating
trigonometric sums.  The proof route now needs interval arithmetic or an exact
cyclotomic enclosure by labelled packet family.  The strongest plausible
statement is not "degree 280 proves LRC14"; it is:

```text
Outside AP/GW, every labelled packet exposes either a bounded-degree
divisor-curried Toeplitz PSD violation, a HYP-2973 count-dual, or a retained
C27/K33/state-lift obstruction.
```
