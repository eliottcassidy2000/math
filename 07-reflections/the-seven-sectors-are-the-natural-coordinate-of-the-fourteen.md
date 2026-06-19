# The seven sectors are the natural coordinate of the fourteen

**Source:** mac-mini-2026-06-18-S6. Dispatch: push the LRC(14) proof via relation-height
split / sector-word discrepancy / finite low-height + high-height tail. Built on codex's
HYP-2603 (seven-sector net cap). Canon: THM-532.

## One lattice, three windows

The LRC(14) residual has now been written as a positivity question for the **same object**
— the singular-series tail over the cluster's offset relation lattice — tested against three
different windows:

- the **5/7-window** `Σ_gaps(gap−2/7)_+` (HYP-2600a/2601), certified when the absolute tail
  `B(E) < (5/7)^k ≈ 0.013`;
- the **6/7-window** `μ_{1/7}=meas{maxgap>1/7}` (HYP-2602), reduced to "the AP minimises";
- the **seven sectors** `meas(S7)=meas{all [j/7,(j+1)/7) hit}` (HYP-2603/THM-532).

All three are the orbit `{frac(e_i x)}` weighed against a test function, and all three split
the same way: a main term plus a relation-lattice correction whose size is graded by the
**relation height**. Dissociated clusters (high height) sit near the main term; arithmetic
progressions (low height, dense short relations) carry the whole correction. The split the
user named — *finite low-height patterns plus a high-height tail* — is not a trick for one
window; it is the shape of the object.

## Why the seven-sector window is the good coordinate

The threshold is `1/14 = 1/(2·7)`. The `7` is not decoration: it is the modulus. The
sector-Fourier coefficients `â_T(n)` carry exactly the **THM-503 vanishing** — `â_T(7m)=0` —
the same `s(7j)=0` that has organised this problem since the singular series was written
down. The seven sectors are the orbit's coordinates *aligned to that modulus*, and in those
coordinates two things become clean that were murky in the gap measure:

1. **The main term is tiny.** `M7(8) = 20160/823543 ≈ 0.0245`, far below `cap_8 = 0.38`. A
   dissociated 8-cluster fills all seven sectors at once only `2.4%` of the time — its orbit
   is a one-dimensional curve, not a cloud, and a curve threads seven sectors rarely. So the
   high-height certificate has a `~0.357` margin, thirty times the `(5/7)^k` budget the
   5/7-window allows.

2. **The correction is the short-relation weight.** `corr(E) ≈ C·W(E)`, `W(E) = Σ_{triples}
   1/height`, and `W` is maximised — like `meas(S7)` itself — by the consecutive AP. The "AP
   is extremal" statement that resisted proof as a gap inequality becomes, in sector
   coordinates, "the AP is the densest in low-height triples," a pure relation count.

The reframe didn't make the problem smaller; it made the difficulty *visible and graded*. The
hard mass is concentrated in one place — the arithmetic progressions — and labelled by a
single positive number, `W(E)`.

## What it would take to finish

Honest: the crude product bound does not close. `C*·W(consec)` overshoots `corr(consec)`
because the worst correction-per-triple and the most triples live on different shapes; their
product double-counts. So the high-height certificate certifies everything *except* a narrow
band of AP-rich shapes near the single consecutive cluster, and that band still needs a
finite check. The two genuine tasks are now small and concrete: (i) write the absolute bound
`corr ≤ C·W` with the explicit sector-Fourier `C` (a support-3 sum plus a geometric tail —
exactly HYP-2601's calculation in the sector basis); (ii) show `{W > W0}` is finite modulo
scaling and enumerate it. Neither is the old wall. The wall was "infinitely many clusters, no
compactness." What is left is "one arithmetic progression and its immediate neighbours,"
counted in a coordinate the modulus chose for us.

## The recurring lesson, once more

Three sessions running, the same move has paid: the obstruction was real but attached to the
wrong object — the width not the measure (S1), the speed scale not the invariant (S5), and
now the gap not the sector. Each time, choosing the coordinate the problem's own arithmetic
prefers — here the `7` that was in the threshold all along — turned a global extremal
mystery into a graded, mostly-certified, finitely-residual statement. The mathematics keeps
offering its own coordinates; the work is to accept them.
