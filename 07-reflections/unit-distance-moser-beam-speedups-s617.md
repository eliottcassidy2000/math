# Unit Distance Moser-Carrier Counting Speedups, S617

The practical answer to "can we speed up their count?" is yes: do not count
unit distances child-local if the search already lives in a carrier with a
known unit shell.

In the S614 Moser coordinate carrier, the unit shell has `18` directed vectors.
For a parent cluster `S`, build a frontier table

```text
gain(q) = |{u in U : q+u in S}|.
```

Then every child has

```text
E(S union {q}) = E(S) + gain(q).
```

This replaces a pairwise distance recount for every candidate child by one
state-local frontier pass. In the S617 width-`1200` run, the displayed
edge-count-check ratio reaches `211.4x` at `n=22`. The whole script still
spends time canonicalizing states and sorting the beam, but the count itself is
no longer the hot repeated operation.

The beam also gives a useful dense-core extension ledger. It recovers the
known Moser-carrier/lower-bound lane with `57` edges at `n=21` and `60` at
`n=22`. But the retained `57`-edge `21`-cores only accept gain-`3` frontier
points, and retained `56`-edge cores only reach gain `4`. So this carrier beam
explains its own failure to find `61`: it is exploring many children, but not a
gain pattern strong enough to cross the one-bit frontier.

That is exactly the HYP-2176 shape. The graph-only quotient is too coarse; the
Moser carrier is better; and the extension ledger is better still because it
keeps the deletion-core side channel. The next proof-grade object should
combine four pieces:

- Moser or known embedded `21`-cores;
- frontier gain tables for candidate degree-`4/5` extensions;
- automorphism and bitset speedups so the same state is not counted many ways;
- totally-unfaithful obstruction labels before any expensive embedding route.

The late incoming opus S599w-x/HYP-2187 packet is a strong cross-check: it
applies the same state-local frontier-gain idea to LRC, replacing full child
recomputation with survival bitmasks on `Z/(2n-1)`. So the abstraction is not
"Moser beams are special." It is: choose the carrier that preserves the target
predicate, then maintain the child observable by a frontier ledger.

I also challenged the tournament vertex set, per the repo rule. Points are not
the right vertices for this question. The script compares counting tricks and
proof routes instead: `21`-core ledger, unfaithful obstruction library,
frontier-gain count, Moser shell, bitset popcount, graph-only enumeration,
triangular-lattice beam, naive recount. This quotient preserves speed,
geometry retention, side-channel value, `n=22` specificity, proof leverage, and
risk. It destroys individual geometric embedding data; that is why the
obstruction library and core solver outrank the raw incremental counter.

One exploratory side run with target `26` and width `400` continued the same
Moser-carrier lane as `64, 68, 72, 76` edges for `n=23,24,25,26`. I would not
read that as an optimum claim. It is useful because it says the speedup is not
just tuned to `22`: the frontier-gain carrier can keep walking through larger
dense clusters while we develop the proof filters.

S617 therefore pushes the unit-distance program from "search denser graphs" to
"audit extension certificates." The arithmetic-carrier lesson from HYP-2184,
the incoming HYP-2186 equidecomposability/Dehn-volume sharpening of the
pi-over-3 carrier, and the coimage lesson from HYP-2176 line up cleanly: count
on the retained carrier, then ask what side channel the quotient forgot.
