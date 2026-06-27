# LRC14 Twist-Ladder Dual Certificate

HYP-2972 tries a different proof angle from the boundary-gap and lift-packet
work: use rational twists as primal certificates, and failed twist ladders as
dual blocker hypergraphs.

For a row `S`, a reduced phase `a/q` is a certificate when

```text
||v a/q|| >= 1/14 for every v in S.
```

If a ladder of candidate twists has no certificate, then each twist is blocked
by at least one speed.  That finite incidence relation is the dual object:

```text
vertices = twists a/q
hyperedge(v) = twists blocked by speed v
```

Default S155 audit on the HYP-2963 bank:

```text
audited rows                 21913
q<=27 certified              21908
q<=27 missed                 5
q<=42 certified              21913
full ladder max first q      41
```

The five `q<=27` misses are exactly the divisor-loaded lcm-tail packet family:

```text
{1,2,3,4,5,6,7,8,9,10,11,13,84m}, m=1..5.
```

All five are rescued by the same witness:

```text
t = 17/41.
```

This is not a claim that `q<=42` proves LRC14 globally.  HYP-2901 already
warns that no fixed finite denominator ladder can do that, because committed
lcm speeds kill all denominators below their wall.  The live theorem is
dynamic: after Moon-core / packet reductions, either a bounded rational ladder
gives a witness, or the blocker hypergraph exposes a committed wall, private
resource, K33/state-lift packet, or next rung.

This route is intentionally not endpoint geometry.  It asks for a global
discrete witness first, and treats failure as finite set cover.
