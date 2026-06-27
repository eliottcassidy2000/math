# LRC14 Dual Cover Certificate

This is a different proof route from the Moon-core / labelled-packet stack.
Instead of asking which residual family a row belongs to, it asks whether the
open danger arcs contain a positive-winding endpoint circulation.

For `delta=1/14`, define

```text
I(s,m) = ((m-delta)/s, (m+delta)/s).
```

For two danger arcs `a=(s,m)` and `b=(r,n)`, the exact integer endpoint credit is

```text
K(a,b) = 14*r*s*(R(a)-L(b))
       = 14*(r*m - s*n) + r + s.
```

Then:

```text
K>0  means overlap,
K=0  means endpoint pinch,
K<0  means strict safe gap.
```

The equality case forces `r+s == 0 mod 14`, which explains why AP/GW boundary
cycles are pair-sum divisible by `14`.

Build `G_open(S)` with vertices the danger arcs and an edge `a -> b` when the
next left endpoint of `b` lies strictly inside `a`, with a wrap label
`epsilon in {0,1}`.  Then:

```text
S is a strict LRC14 counterexample
iff G_open(S) has a directed cycle of total winding 1.
```

The dual certificate is a potential `Phi` satisfying every open edge constraint

```text
Phi(b) + epsilon <= Phi(a).
```

Summing around a cycle forbids positive winding.  This is a Farkas/difference
constraints certificate for LRC14 safety.

AP and GW should sit on the boundary:

```text
G_closed has a unit-winding zero-credit cycle,
G_open has no unit-winding cycle.
```

A quick exact graph check matched that:

```text
AP, GW:                 G_open no, G_closed yes
12->36, 12->84, 12->26: G_open no, G_closed no
```

So the alternate proof target is:

```text
Every primitive 13-speed row either has a dual endpoint potential,
or is AP/GW zero-credit boundary,
or its positive-credit SCC carries the K33/H=7 state-lift obstruction.
```

This does not replace HYP-2964 through HYP-2968.  It turns their positive
boundary gaps, NORK pinch templates, and few-apex lift packets into shadows of
one graph object: the endpoint-credit winding graph.

Tournament analysis uses endpoint transition events as vertices, not runners.
The pair observable is `(winding, sign(K), K33 visibility, pair-sum divisibility
by 14, normalized overlap)`.  A strict counterexample must realize a tournament
class containing a positive-winding SCC with no dual potential and no state-lift
exit.  That is a much sharper falsifier than another scalar low-mass row.

Next exact task: emit `G_open/G_closed` for AP, GW, `12->36`, `12->84`, and the
latest pinch rows; run the winding-cycle / potential check; verify that AP and
GW are the only closed-unit/open-empty atoms in the local bank.
