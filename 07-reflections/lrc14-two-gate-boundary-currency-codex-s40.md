# LRC14 Two-Gate Boundary Currency

The useful reframe from S40 is that the near AP-tail shell story and the far
Delta story are the same proof currency seen in two gauges.

HYP-2661 says the local mouth is protected by shell-1 tower conservation:
damaging `{1,2,4,8}` should already force enough measure.  HYP-2662/HYP-2664
say the far endpoint imbalance is not controlled by a clean trace or QR/NQR
mean alone, but it looks payable by the mass where the base misses exactly one
inner sector, `p1(E')`.  HYP-2665 sharpened this while I was working:
`p1/3` is false for actual `Delta_w^+`, but `3p1/8` remains viable in the
sample.  Those statements should be applied in that order, not averaged into
one scalar from the start.

The S40 scout makes the order concrete.  In the exact bank
`E'={0} union 7-subsets of [1,18]`, `31788` primitive rows, the inequality

```text
p0(E') + (1/7+c)p1(E') <= cap_9
```

has zero violations through `c=3/8`.  The global top row remains AP8,
`(0,1,2,3,4,5,6,7)`, and it misses tower bit `8`.  The shell-1-full stratum is
less dangerous: its maximum is `(0,1,2,3,4,5,6,8)` with `c=3/8` slack
`194613/2242240`, about twice the global slack.

So the shape is:

```text
packet first, scalar last
local packet damage -> shell-1 carry/mouth gate
far endpoint imbalance -> p1 boundary tax
```

This also explains why the new `[16,18]` rows are interesting but not scary.
Some all-damaged `3`-multiple packets jump into the top thirty by `p1` value,
which smells like the old ramified `2n-1=27` story.  But they are not
shell-1-full rows, so the correct routing is not "large p1, panic"; it is
"damaged packet, send to the tower or 3-ramified ledger first."

The next proof should try to state a coarea-style inequality:

```text
Delta_w^+(E') <= p1(E')/3
```

only as the old tempting target; HYP-2665 now refutes it.  The corrected
version is

```text
Delta_w^+(E') <= 3*p1(E')/8
```

or a packet-refined inequality for the rows above `p1/3`, only after the
HYP-2661 gate is applied.  The important bit is that the tax is charged to an
addressed boundary mass, not to a raw absolute endpoint count.

Tournament Analysis: vertices are proof obligations
`shell1_gate, p1_boundary_tax, missing_packet, cap_slack, raw_value`.  The
observable is the cap inequality through `c=3/8`; the switch is the shell-1 packet
split; the Hamiltonian path is
`shell1_gate > p1_boundary_tax > missing_packet > cap_slack > raw_value`.
The challenged assumption is that far Delta and AP-tail shell gates are
separate mechanisms.  They now look like one boundary-currency system.
