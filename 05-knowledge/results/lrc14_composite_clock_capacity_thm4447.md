# THM-4447 result note: exact composite-clock capacity

**Status:** **PROVED ELEMENTARY + PROVED RELATIVE TO CITED LRCUpTo13 +
VERIFIED-EXACT + INDEPENDENTLY AUDITED.** LRC(14) remains **OPEN**.

## Pointwise theorem

Above a pack-safe phase \(y\), label the \(q\) common lifts by
\(x_j=(y+j)/q\). For one tail \(t\), put \(g=\gcd(q,t)\), \(m=q/g\),
\(\delta=ty/q\), and \(A=-m\delta-m/14\). Its exact strict-danger count is
\[
 |K_t(y)|=
 g\left(\left\lceil A+{m\over7}\right\rceil-\lfloor A\rfloor-1\right),
\]
and its sharp uniform capacity is
\[
                         B_q(t)=g\left\lceil{q\over7g}\right\rceil.
\]
The count is either \(B_q(t)\) or \(B_q(t)-g\). For nonintegral \(m/7\),
the deficit occupies a whole phase chamber; equality walls mark transitions
but are not the only deficit phases. If summed capacities equal \(q\), full
coverage occurs exactly when all tails are cap-tight and their bad-label
sets partition the clock.

## Divisor optimization

For every divisor \(d\mid q\), absorb the \(d\)-divisible tails:
\[
 qC\cup T=d\left((q/d)C\cup\{t/d:d\mid t\}\right)
               \cup\{t:d\nmid t\}.
\]
Whenever the absorbed pack has a safe phase, the remaining row closes if
\[
 \sum_{d\nmid t}\gcd(t,d)
       \left\lceil{d\over7\gcd(t,d)}\right\rceil<d.
\]
The operation preserves the physical row and optimizes over its possible
pack/tail decompositions.

## Exact ten-pack consequence

Let \(P\) be ten speeds, \(T\) three further arbitrary speeds, and
\(P\cup T\) a distinct primitive thirteen-speed row. At the exact body clock
\(c=\gcd(P)>1\), cited LRC through thirteen total runners makes every
absorbed pack of at most twelve speeds safe. The divisor theorem closes all
\(c\ge5\). It leaves only:

- \(c=2\): zero or one even tail;
- \(c=3\): all three tails prime to three;
- \(c=4\): exactly one tail has \(2\)-adic valuation one and the other two
  are odd.

This table was already proved in
[the effective-clock result](lrc14_effective_clock_empty_core_sep06.md) by a
different hereditary-gcd argument. It is recorded here as a **recovered,
independently reproved corollary**. The added content is the pointwise
floor/ceiling law, its chamber/equality description, and the general divisor
optimization.

## Verification

The primary path checks 1,140,196 gates, including all 43,991 primitive
gcd-signature triples through clock 210. The clean-room path imports no
primary or geometry code and checks 254,791 gates, including all 175,081
primitive signatures through clock 512. Both paths pass normally and under
optimization.

    python -B 04-computation/lrc14_composite_clock_capacity_thm4447.py
    python -B -O 04-computation/lrc14_composite_clock_capacity_thm4447.py
    python -B 04-computation/lrc14_composite_clock_capacity_thm4447_independent.py
    python -B -O 04-computation/lrc14_composite_clock_capacity_thm4447_independent.py

The finite sweeps audit, but do not supply, the infinite divisor proof.
