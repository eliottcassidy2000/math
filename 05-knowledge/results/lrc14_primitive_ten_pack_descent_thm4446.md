# THM-4446 result note: primitive ten-pack descent

**Status:** **PROVED RELATIVE TO CITED LRCUpTo13, THM-3818, THM-4052,
AND THM-4442 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** LRC(14) remains
**OPEN**.

## New closure

For a primitive row \(3C\cup T\), where \(C\) has ten speeds, \(T\) has
three distinct ternary-unit speeds, and \(\gcd(C)>1\), choose a prime
\(p\mid\gcd(C)\). A cited phase for \(C/p\) has \(3p\) common lifts. Each
\(p\)-divisible tail kills at most \(p\) labels, each other tail kills at
most \(\lceil3p/7\rceil\), and primitivity permits at most two tails of the
first kind. Thus
\[
 3p-\left(2p+\left\lceil{3p\over7}\right\rceil\right)>0.
\]
At least one lift is safe. The worst certified margin is one label, at
\(p=2\), and the independent referee found an exact primitive row attaining
it.

The endpoint convention is material: safety is weak, so the danger arc is
open. At order \(7k\), an open arc of length \(1/7\) contains at most \(k\)
grid points.

## Dilation rays

For every \(C_0\in\binom{[13]}{10}\), integer \(h\ge1\), and ternary-unit
tail triple \(T\),
\[
                         G_{3hC_0\cup T}\ne\varnothing.
\]
After dividing the full row by its gcd \(g\), one has \(g\mid h\). Scale
\(h/g=1\) is THM-4442; scale \(h/g>1\) has a nonprimitive ten-pack inside a
primitive full row and is closed by the new descent. This upgrades the 286
bounded THM-4442 bodies to 286 complete integer rays.

## Honest remaining chart

In THM-3818/4004's live \(d=3\), \(11+2\) branch, let \(P\subset C\) be
the two divided pair-shore entries. Any survivor must obey
\[
 \gcd(C)=1,\qquad
 \min_{x\in P,z\notin P}
 \max\left({x\over\gcd(x,z)},{z\over\gcd(x,z)}\right)>91^6,\qquad
 \max T<11\max C.
\]
So a survivor cannot be compressed into the bounded chart by removing a
common body scale. Its obstruction is projective shape plus component
address.

The exact unresolved predicate is
\[
 G_{3C\cup T}=\varnothing\quad\Longleftrightarrow\quad G_C\subseteq F_T,
\]
where \(F_T\) is the set of pack phases whose three ternary sheets are all
spoiled. The missing coordinate is the cyclic tail-sheet event word on a
component of \(G_C\), including endpoint owners. A scalar tail mass or a
numerically oriented tournament loses this address.

## Verification

The primary verifier checks translated order grids through 91, the
\(m=7,21\) endpoint cases, eight primes, all 286 ten-subsets of \([13]\),
and every possible count \(0,1,2\) of prime-divisible tails: 6,292 literal
rows.

The clean-room verifier checks grid orders through 210, primes through 1000,
all admissible small-prime tail event cells, and 292,500 normalization rows:
2,189,781 explicit checks. Normal and optimized runs agree.

    python -B 04-computation/lrc14_primitive_ten_pack_descent_thm4446.py
    python -B -O 04-computation/lrc14_primitive_ten_pack_descent_thm4446.py
    python -B 04-computation/lrc14_primitive_ten_pack_descent_thm4446_independent.py
    python -B -O 04-computation/lrc14_primitive_ten_pack_descent_thm4446_independent.py

The computations audit the proof's sharp finite-orbit mechanism. They are
not used to infer the theorem's all-height quantifier.
