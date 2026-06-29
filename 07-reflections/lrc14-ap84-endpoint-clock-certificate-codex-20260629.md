# LRC14 AP84 Endpoint-Clock Certificate

HYP-3454 is a useful tightening of HYP-3452 rather than a new route.  The AP
tail was already known to enter a rank-one `E:84m/E:84m` phase at `m=5`; the
new contribution is to isolate the exact inequalities that make this phase a
proof object:

```text
I_m=[(14ceil(48m/7)+1)/(588m),(14ceil(48m/7)+13)/(588m)]
```

lies inside the fixed low corridor `[8/49,6/35]`, has length `1/(49m)`, and
is bounded by the moving `E:84m` wall on both sides.  The checked tail
`m=5..70` has no endpoint failures and endpoint rank `1`.

This should be read beside incoming HYP-3453, which gives the bank-level
gate-escape transversal.  HYP-3453 says dead-cover obstructions should expose
low-rank survivor gates; HYP-3454 says that in the AP-tail corridor the
eventual endpoint gate has a closed-form clock.

The remaining work is now better separated.  The four transients `m=1..4`
are finite mixed `E/B1` cases.  The infinite endpoint tail is a closed-form
wall inequality.  HYP-3456 now removes the still-experimental infinite
ingredient by deriving the period-`35` escape count from HYP-3431 fixed low
corridors and the moving high-grid gaps.

For the LRC14 proof route, this should feed HYP-3439 exactly where its AP-tail
rank-`5` descent was vague.  Instead of saying "use HYP-3452 for the tail",
the bridge can require:

```text
finite m=1..4 cases
closed-form endpoint interval for m>=5
HYP-3456 mod-35 boundary count lemma
```

The challenged assumption was that the AP-tail should be organized by raw
dead fraction or by `m` as the tournament vertex.  The useful vertices are
proof obligations: endpoint interval, low-corridor containment, moving
`E:84m` gap, finite transients, and residue clock.  This preserves the LRC
escape predicate while refusing to treat scalar dead fraction as a certificate.
