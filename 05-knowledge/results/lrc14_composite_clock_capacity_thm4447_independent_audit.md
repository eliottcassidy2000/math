# Independent audit of THM-4447

**Verdict: PASS, after provenance and endpoint-wording repairs.** The exact
pointwise count, sharp cap, divisor absorption, arbitrary-tail strengthening,
and symbolic small-clock reduction are theorem-safe. LRC(14) remains open.

## Corrections made before promotion

The clock-\(2,3,4\) residual table is not a new classification:
[the prior effective-clock result](lrc14_effective_clock_empty_core_sep06.md)
already proves it using THM-761/765, hereditary deletion gcds, and pair lcm
constraints. THM-4447 labels it as recovered and independently reproved.

A second wording repair separates low-count chambers from resonance. When
\(m/7\) is nonintegral, the count stays one orbit point below capacity
through a full interval of translates. Resonance is a transition wall, not
the cause of every interior deficit. The exact nonresonant control
\((q,t,y)=(9,1,0)\) has one bad label against capacity two.

## Proof audit

Writing \(q=gm\), \(t=ga\), multiplication by \(a\) permutes the order-\(m\)
grid and each point has \(g\) labelled preimages. Strict tail danger is
bijective to the integers in
\[
 (-mty/q-m/14,\,-mty/q+m/14),
\]
which proves the floor/ceiling formula and the cap
\(g\lceil q/(7g)\rceil\). The referee separately checked the two chamber
levels, literal equality walls, cap attainment, and the equality-partition
iff.

For \(d\mid q\), the physical-row rewrite is literal. Original distinctness
prevents a collision in the absorbed pack, and primitivity prevents all three
tails from being absorbed. Thus it contains at most twelve speeds and the
cited lower-dimensional LRC input is legitimate.

The arbitrary-tail infinite proof was rederived by prime-factor cases. A
prime \(p\ge5\) closes directly. For \(2^a3^b\), the divisors
\(2,3,6,4,8,9\) close every case outside the displayed clock-\(2,3,4\)
signatures. This proof uses no finite cutoff.

## Independent coverage

The clean-room script imports neither the primary audit nor repository
geometry. It checks:

- 36,940 wall/chamber samples through reduced order 140;
- 50,547 literal label counts;
- 1,608 typed physical divisor rewrites;
- all 175,081 primitive arbitrary-tail gcd signatures through clock 512;
- equality partitions, divisor-two covers, resonant walls, and a
  nonresonant low-count control.

Normal and optimized transcripts agree at 254,791 assertion gates.

    python -B 04-computation/lrc14_composite_clock_capacity_thm4447_independent.py
    python -B -O 04-computation/lrc14_composite_clock_capacity_thm4447_independent.py
