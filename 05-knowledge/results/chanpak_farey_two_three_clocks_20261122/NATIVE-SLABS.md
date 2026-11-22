# Unrestricted exact Farey slabs on the remaining two spines

Let \(N=E_{61}-1=2^{60}(2^{61}-1)-1\). These searches leave both
addresses arbitrary within the stated intervals, not just periodic,
few-run, or synchronized addresses. They are **not** exhaustive at
either whole spine.

Use the canonical Farey charts from the decimal-kernel packet.
Reduction modulo five allows \(x=10X+i\) with \(i=0,1\).
Set \(y=2Y+j\) for \(j=0,1\) and scan the smaller of \(X,Y\),
breaking equality in favor of \(X\). This gives eight disjoint
parity/orientation branches. Their exact affine-plus-read equations
and row fallback are described in [the integer spine-eight
inverse](../chanpak_farey_complementary_exact_20261122/SPINE8-SWEEP.md). The present implementation only covers
a specified middle interval, using rational tangent slabs in
native C++.

## Exactness independent of both floating seeds

After normalization and possibly swapping sides, write the chart
as
\[
 R=QXY+TX+CY+h\big(f(sX+o)+f(uY+p)\big).
\]
The same lower and upper bounds on \(Y\) follow from
\(0\le f(v)\le v\); they are rational functions of \(X\).
For a block beginning at \(X_0\), use any reduced rational
slope \(a/b\), a chosen integer anchor \(Y_0\), and line coordinate
\(z=b(Y-Y_0)+a(X-X_0)\). The program's continued-fraction seed
is a convenience only: every rational slope is valid.

The outward upper bound on \(z\) evaluates both integer
endpoints of the block, by convexity. For the lower bound write
the rational boundary as \((R'-T'X)/(QX+C')\). The sign of its
consecutive difference after adding the line slope is the sign
of
\[
 aD(D+Q)-b(QR'+T'C'),\qquad D=QX+C'.
\]
It is nondecreasing in \(X\). A floating seed proposes the turning
point, then exact integer comparisons move it until its left
difference is negative and its right difference nonnegative
(or until it is an endpoint). Both endpoint line evaluations
and this exact lower evaluation round outward by signed integer
division. Thus neither seed can remove a valid point.

Each line has \(bY=bY_0+z-ad\), \(d=X-X_0\). Modulo the two
prime factors 17,29 (spine four) or 29,257 (spine eight) of the
corresponding \(h\), substitution in \(QXY+TX+CY=R\) gives a
quadratic congruence in \(d\). This remains necessary when the
modulus divides \(b\); no inverse of \(b\) is taken in making
that equation. Find all modular roots by a square-residue table,
including linear, contradictory and identically vanishing cases.
Combine their residue classes with \(ad\equiv z\bmod b\) by
generalized CRT. Every remaining candidate is checked against
the literal original height; modular equality is never used as
a proxy for integer equality.

For arithmetic range, the implementation accepts \(X\ge10^{12}\),
block width at most \(2^{16}\), and target at most \(N\).
The reachable \(QX+C\) is below \(2^{77}\) (since \(Q<2^{32}\),
\(R<2^{121}\) and \(X\le\sqrt{R/Q}\)); the lower denominator
is below \(2^{78}\). The tangent magnitude in this scope is
less than \(2^{16}\) because \(Q\ge45{,}875{,}200\).
Thus \(a\le 2^{32}\), and \(a d(QX+C')\) fits in signed 128 bits.
The remaining line numerator is evaluated as the **anchored**
residual \(b(R'-T'X-Y_0(QX+C'))\), also within signed 128 bits:
the upper-curve residual at the anchor is a remainder below
its denominator, and its displacement changes by
\((T'+QY_0)d\) with \(QY_0\le R/X_0\). The finite-difference
comparison uses four exact 64-bit limbs. The modular quadratic
constant is formed after reducing each term, not by computing a
potentially overflowing intermediate \(bR\). The physical
address words fit in unsigned 64 bits in the stated scope.

## Completed intervals

The coordinate in the table is \(\min(X,Y)\) **after the above
normalization**, not a raw address or a sum.

| Spine | Inclusive interval | Slabs across all eight branches | Literal height tests | Hits |
|---:|---|---:|---:|---:|
| 4 | 150,000,000,000,000 to 150,020,000,000,000 | 3,906,256 | 1,379,231 | 0 |
| 8 | 14,950,000,000,000 to 15,050,000,000,000 | 52,083,336 | 255,009 | 0 |

Native replay:

~~~sh
g++ -O3 -std=c++17 -Wall -Wextra -Werror farey-native-slab.cpp -o farey-native-slab
./farey-native-slab 4 150000000000000 150020000000000 40960
./farey-native-slab 8 14950000000000 15050000000000 15360
python farey-native-check.py
~~~

The Python replay exercises 32 planted long pairs at the two
spines, several parity offsets, and both orientations. The
compiled scan is exact when several distinct subintervals are
distributed across workers, because slab anchors need not be
shared between them. This does not close the whole exponent-61
search.

## A dyadic phase envelope with a controlled additive defect

For a binary word \(n\) of width \(k\), define
\[
 \rho(n)=\sum_{i=0}^{k-1}(1-b_i)2^{-i-1}
        =\frac{2^k-1-\operatorname{rev}_k(n)}{2^k},
 \qquad \epsilon(n)=f(n)-n\rho(n).
\]
Use zero width for the empty prefix. Appending a bit acts on
the reversed phase and additive correction by
\[
 \begin{array}{ll}
 0:& (\rho,\epsilon)\mapsto((\rho+1)/2,\epsilon+[n\text{ odd}]),\\
 1:& (\rho,\epsilon)\mapsto(\rho/2,\epsilon-\rho/2).
 \end{array}
\]
This is an exact substitution in the defining read recurrence.
In particular the additive defect changes by at most 1 upwards
or \(1/2\) downwards per appended bit, not by the magnitude of
the address. If the high prefix \(p\) is fixed and there are
\(g\) following bits, then
\[
 \epsilon(p)-g/2\le\epsilon(n)\le\epsilon(p)+g.
\]
If in addition the low suffix \(w\) is fixed to padded width
\(\ell\), then
\(\rho(n)\in[\rho_\ell(w),\rho_\ell(w)+2^{-\ell})\).
Combining gives an exact affine envelope for \(f(n)\) with
known slope to precision \(2^{-\ell}\) and an error additive
in the free length. This is a dyadic phase/height separation
useful for conditioning the slab search on actual words; it
does not count the resulting branch tree as exhausted.

## Contiguous trillion-row extension at spine eight

**FINITE-EXACT.** The additional inclusive interval
\(15{,}050{,}000{,}000{,}001\) through \(16{,}050{,}000{,}000{,}000\)
also has no hit. All eight parity/direction branches were completed
in all eight disjoint 125-billion-row subintervals, using 16,000-row
slabs. The extension used **500,000,000 slabs**, visited
2,888,867,151 integral line indices, and made **2,316,502 literal
height tests**. It joins the first spine-eight interval without a gap:
the combined interval is \(14{,}950{,}000{,}000{,}000\) through
\(16{,}050{,}000{,}000{,}000\), with \(1{,}100{,}000{,}000{,}001\)
distinct values of the smaller normalized coordinate, independently
of the eight parity/direction assignments.

~~~sh
for j in $(seq 0 7); do
    L=$((15050000000001+j*125000000000))
    U=$((L+125000000000-1))
    ./farey-native-slab 8 "$L" "$U" 16000
done
~~~

The exact full branch output is recorded in the extension transcript.
This is a contiguous new exclusion at the remaining spine, not
the whole spine-eight sweep; even less does it exclude other cores.
