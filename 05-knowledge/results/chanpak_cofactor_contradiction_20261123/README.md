# A square-root contradiction test for the first missing tournament count

**PROVED scoped / FINITE-EXACT; Chan–Pak remains open.**
This extends the [counterexample-machine packet](../chanpak_tiling_machines_20261123/README.md),
the [graded window endgame](../chanpak_window_endgame_20261123/README.md),
and the [literal factor-seven fork](../chanpak_ternary_buckets_20261122/PHYSICAL-FORKS.md).
The new contradiction tests tolerate extraneous prime factors.
There is no derivation of the needed witnesses from RH or GRH.

## A first hole forces an exact obstruction

The counterexample machine tries allowed odd targets in increasing
order and halts on its first certified absence. Since each query has
a finite bound, this is an unconditional deterministic search, not a
guesser stopping at a timeout. Assume it halts at \(N\). All smaller
allowed targets have witnesses, but unequal known holes 7 and 21
cannot be treated as induction seeds.

Ordinal composition multiplies counts. If \(7\nmid N\) and \(N\) is
composite, two smaller factors witness it. Two or more 7s can be packed
into 49s and, if needed, 343; these are attained. For exactly one 7,
a composite cofactor can be split by giving 7 to a proper divisor
other than 3, or treated with the attained 63 when the cofactor is a
power of 3 of exponent at least two. Thus the only cases to reject are
\(N=p\) with odd prime \(p\ne7\), and \(N=7p\) with odd prime \(p\ne3,7\).
These are the inherited two lanes. Crucially, the [new cofactor tests](COFACTOR-TEST.md)
do not require the total witness height to be an exact prime power
or exact single-seven ray.

For example, fix a suspect prime \(p\), put
\(x=H(T)=p^r m\) with \((m,p)=1\), and let \(K_p\) count
the strongly connected components (SCCs) of \(T\) whose own Hamilton
count is divisible by \(p\). If \(p\) were absent, then
\[
 2K_p\le r+\Omega(m).
\]
Each \(p\)-bearing SCC must either spend a second \(p\) or spend a prime
from \(m\). Counting those costs proves the inequality. Hence among
witnesses with \(m\le x^\varepsilon\) one has the literal ceiling
\[
 p^{K_p}\le x^{\,1/2+C_p\varepsilon},
 \qquad C_p=\frac{\log p}{2\log3}.
\]
A witnessed power excess \(p^{K_p}>x^{1/2+\delta}\) with
\(C_p\varepsilon<\delta\) contradicts the missing prime. The
constant \(C_p\) is fixed before choosing \(\varepsilon\). No
asymptotic estimate for \(K_p\) is assumed; this gives a precise
goal for constructing hosts.

In the \(7p\) lane let \(Z\) count SCCs carrying both \(p\) and
**exactly one** 7, let \(r=\nu_p(H(T))\), and strip all 7s and
all \(p\)'s to form the coprime remainder \(m\). Under absence
of \(7p\) the exact inequality is \(Z\le r-K_p+\Omega(m)\).
It permits a growing signal by accumulating many separated
7-carriers. At any fixed p one witnessed violation forces a
component of count \(7p\), with no numerical division of
Hamilton counts.

## What actually transfers from the analytic endgame

In [Slijepčević's outline, Props. 2.2–2.5](https://hexagonmath.org/pdf/2611.00812v1)
of the verified Dirichlet GRH result, fix the allegedly bad
zero and its character first, choose the loss parameter so that
\(1/2+(e_D+e_L+7)\varepsilon\) stays strictly below its real
part, and arrange the one-sided bound and the forced excursions
on the **same** arbitrarily late windows. The \(\sqrt{x}\) scale
is part of how that proof gets the analytic bound, not a universal
recipe for errors in a graph search. The [independent checks](https://math.mit.edu/~drew/grh/)
concern ordinary RH and Dirichlet GRH, not this problem.

Here the role of the fixed bad object is the first suspect prime
or single-seven cofactor. The rank obstruction is pointwise in
each actual host once that object is fixed. A compatible witness
over the forbidden threshold suffices. What is missing is the
supply of such witnesses for **every** suspect \(p\), at the
same heights at which the cofactor error is small. Knowing that
a mass lives in a nearby residue class or that a fork is in the
same large numeric window does not locate its SCC factors.

## What the supplied fork gives, and where the width is lost

The old \(h=7(40\cdot2^s+1)\) rail is literal and two-sided;
the aligned parent-to-center route is five edits inside
\([h,3h]\), not height-monotone. [Arbitrary fillers](FILLERS.md)
extend its scope: put a nonempty tournament \(M\) in the uniform
growing block. With \(D=H(M)\), \(Q=Q_2(M)\) the two-slot
count including empty slots, the same exterior edits give
\[
 h=7(40Q+D),\qquad 3h-2D,\quad 3h,\quad 3h+2D .
\]
This is a real way to produce additional \(7p\) when \(40Q+D\)
is prime. It is not the old unit fork unless \(M\) is transitive.
The anchored census through order eight gives 284 distinct
prime cofactors at these orders: **283 off the transitive rail**
(seven at order six, 276 at order eight), plus 641 on it.
For example \(D=661,\ Q=8848\) gives \(p=354581\) and the
actual parent \(h=2482067\). These are new outputs for this chart,
not claims that their counts were previously unknown.

[Rigidity](RIGIDITY.md) separately shows that the original two
marked arcs keep exact gaps \(-2,+2\) under arbitrary tagged
ordinal grafts only on the original dyadic chart, up to irrelevant
padding and the \(s=0\) endpoint. Allowing both binary flank
addresses and any of the 28 fixed core or core–apex arcs gives
the same restriction. This leaves genuinely different cores,
nonordinal contacts and edits touching a moving flank open.

Even if all subsequent arithmetic ternary digits were granted
for free, the [ray's inverse dynamics](DYNAMICS.md) covers
less than \(3/280\) of each full ternary shell. It has dense
leading phases yet leaves whole gaps of length
\(\Omega(x/\log x)\), much larger than \(x^{1/2+\delta}\) for
\(\delta<1/2\). Every fixed finite-state superset of this
optimistic guard language instead accepts asymptotically nearly
everything. The same note gives a first-hole countermodel
satisfying ordinal multiplication, all smaller allowed seeds,
and this optimistic ray closure. These are statements about
the chosen rail, not about the real spectrum.

Replays: **python3 replay.py**, **python3 sevenfork-all-arc.py**,
and **c++ -O3 -std=c++20 fillers-census.cpp -o /tmp/fillers-census**
followed by **/tmp/fillers-census 8**. Exact prime lists and
receipts are included. See also [direct-ray separation](COFACTOR-TEST.md#why-near-power-hits-on-the-original-affine-ray-do-not-suffice).
