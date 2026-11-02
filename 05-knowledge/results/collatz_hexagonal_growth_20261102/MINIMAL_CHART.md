# The smallest chart for a declared branch-and-clock language

**PROVED all horizons for the specified observer; FINITE-EXACT for the
weighted-polynomial comparison. No full arithmetic-state compression claim.**

Take any odd unit n>1, not only descendants of5. At a remaining horizon
h>0 the observer reports the two SF5 exponents A(n), ordered increasingly,
and accepts either branch selector0 or1. A selection of exponent a emits
the ordinary clock a+1, replaces n by j_a(n), and decrements h. At h0
the observer reports STOP and accepts no further branch selector.
Literal integer labels, arbitrary arithmetic queries and native meeting
labels are not outputs in this language. They require the source register
retained in JOINT_PORTS.md.

## Minimality theorem

At fixed h>=1, the equivalence classes of sources under every permitted
future observation are exactly their unit residues modulo3^(h+1).
Consequently the least autonomous deterministic chart has

    c_h=2*3^h states,       c_0=1.                    (C1)

If the remaining horizon0..H is part of the chart and STOP is observed,
the total is

    1+sum_(h=1)^H 2*3^h = 3^(H+1)-2.                (C2)

These counts ignore no hidden integer queries: permitting such queries
changes the language and generally requires the actual source as a sidecar.
They are not memory lower bounds for an algorithm that can freely reread
an external unbounded integer register.

**Sufficiency and closure.** In a state(r,h), take r modulo3^(h+1).
The two exponent choices depend only on r mod9. The selected successor
is

    ( (2^a r-1)/3 modulo3^h, h-1 ).                  (C3)

The division is integral on the reported branch. Different representatives
give the same reduced residue, because their difference is divisible by
3^(h+1). The integer source remains an odd positive unit under the actual
SF5 branch; the residues alone need not be odd representatives. This is
a closed update precisely because the retained precision decreases.

**Necessity.** The initial menu identifies the residue modulo9:

| n mod9 | A(n) |
|---:|---|
|1|2,4|
|2|3,7|
|4|2,6|
|5|3,5|
|7|4,6|
|8|5,7|

The six menus are distinct. Suppose h>1 and two roots have identical
outputs along the path that always chooses the smaller exponent. Their
first menus agree, so the same a is selected. Inductively their child
residues agree modulo3^h. But

    j_a(n')-j_a(n)=2^a(n'-n)/3,

and2^a is a3-adic unit. Thus the roots agree modulo3^(h+1). The complete
menu trace along this ONE selected path already distinguishes all unit
classes, proving necessity. Every unit residue has an odd representative
greater than1, so none is an empty formal state. Different horizons are
distinguished by the time at which STOP appears. QED.

For H1,2,3 the complete chart sizes are7,25,79. These numbers arise from
the explicit geometric sum(C2); no identification with unrelated groups,
primes, or tournament orders follows from the values alone.

## Counts versus the full clock polynomial

The source-dependent recurrence J5 factors through(C1). chart_checks.py
computes it by memoizing(r,h) and dropping one3-adic digit at each step.
All coefficients are additionally compared to literal native path unions
for EVERY root class at h1..4 (240 root classes). The larger census is:

|h|Full branch/clock classes|Distinct polynomials Z_h|Distinct totals Z_h(1)|
|---:|---:|---:|---:|
|1|6|6|4|
|2|18|18|8|
|3|54|54|12|
|4|162|162|22|
|5|486|486|38|
|6|1458|1458|58|
|7|4374|4374|90|

**FINITE-EXACT:** no two unit classes collide in the full polynomial
through h7. **OPEN:** injectivity of Z_h on these classes at every h.
The minimality proof above uses menus along a named branch, not unlabelled
polynomial coefficients, and therefore does NOT prove this open statement.

Already at h1, root13 (residue4) and root7 (residue7) both have total9,
but their menus are(2,6) and(4,6), with different ordinary clocks. Hence
total expanded size cannot be the exact branch observer. Repeated use
of the total statistic does not recover omitted source or port data.

This answers a finite-horizon minimal-carrier question completely while
retaining its boundary. For one fixed source tree, a direct recurrence may
be cheaper than constructing the full chart: there are2^h odd leaves but
2*3^h possible source classes. No algorithmic speedup is inferred from
memoization without measuring repeated states for the actual workload.
