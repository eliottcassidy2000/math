# An integer-only exhaustive inverse for spine length eight

This independently certifies a complete search algorithm and budget,
not completion of the whole search. In the canonical Farey tower use
the lower chart \((h,q,c,t,e,r,k)=(37265,144834560,84839910,57602544,
24360,14906,33838015)\) from the receiver specialization. The upper
orientation adds two. At \(N=E_{61}-1\) reduction modulo five permits
only \(x\equiv0,1\bmod10\). Pin both parities by writing
\(x=10X+i,\ y=2Y+j\), \(i,j\in\{0,1\}\), and partition into \(X\le Y\)
and \(Y<X\). In the second branch exchange normalized variables.
After subtracting constants including parities, the chart is
\[
 R=QXY+TX+CY+h(f(sX+o)+f(uY+p)).
\]
Here \(Q=20q\), \(s,u\) are 10 and 2 in the appropriate order,
and \(T,C,R\) are obtained by expansion. Every solution has
\(0\le X\le Y\) after the chosen orientation and
\(X\le\lfloor\sqrt{R/Q}\rfloor\).
The binary read satisfies \(f(0)=0,\ f(2w+1)=f(w)\) and
\(f(2w)=f(w)+w+[w\text{ odd}]\).

## Integer extrema instead of rounded tangent enclosures

Because \(0\le f(v)\le v\), the opposite coordinate belongs to
\[
 Y_-(X)=\frac{R-h(o+p)-(T+hs)X}{QX+C+hu}
 \ \le Y\le\
 Y_+(X)=\frac{R-TX}{QX+C}.
\]
Choose a block \(X_0,\ldots,X_0+W-1\) and any reduced slope \(a/b\).
The code uses exact continued-fraction approximation of
\((QR+TC)/(QX_0+C)^2\) with \(b\le W\); a poor slope only costs work.
With \(d=X-X_0\) and anchor \(Y_0\), every survivor has integral
index \(z=b(Y-Y_0)+ad\). Minimize the corresponding expression
with \(Y_-\) and maximize it with \(Y_+\) over the integer block.
For either rational boundary \((R'-T'X)/(QX+C')\), the consecutive
difference is
\[
 a-\frac{b(QR'+T'C')}{(QX+C')(Q(X+1)+C')}.
\]
Its sign changes at most once and is tested by exact cross
multiplication. The endpoints and a bisection for that change
give outward ceiling/floor values without floating arithmetic.

## Solve the modular hyperbola along a line

Put \(D=QX_0+C,\ J=bY_0+z\). On a line \(bY=J-ad\).
At either prime \(p=29,257\) dividing \(h\), the absolute equation
therefore implies
\[
 -Qa\,d^2+(bT+QJ-aD)d+bTX_0+DJ-bR\equiv0\pmod p .
\]
Use the exact root set at each prime (including linear, contradictory,
and identically-zero cases), and combine with \(ad\equiv z\bmod b\)
by generalized CRT. Multiplication by \(b\) is used only as a
necessary condition, so the certificate is valid even if \(p\mid b\).
For each offset check the designated orientation and the complete
original integer height on the reconstructed addresses.

Near the axes the code fixes one coordinate and splits off a
low binary suffix of the other using
\[
 f(2^m v+w)
 =f(v)+f(w)+(2^m-1-\operatorname{rev}_m(w))v
        +[w<2^{m-1}][v\text{ odd}].
\]
For a fixed suffix, the linear coefficient in the quotient is
known; the omitted nonlinear part lies between zero and \(h(v+1)\).
The quotient interval is intersected with the necessary slope
congruence modulo \(h\) and the required physical-address residue.
This covers the axes without needing a separate theorem.

## Reproducible finite budget and measured window

Set \(K=\max(1,\lfloor\sqrt[3]{\lfloor R/Q\rfloor}\rfloor)\).
Below \(2K\) use exact axis rows. Above it, a block starting at \(X\)
has width \(\lfloor X/K\rfloor\), truncated at the final endpoint.
All eight E61 branches have \(\max X=30{,}294{,}489{,}443{,}840\)
and \(K=971{,}797{,}457\). The complete schedule is exactly
**1,943,594,914 axis rows and 9,644,631,668 rational slabs per
branch** (15,548,759,312 rows and 77,157,053,344 slabs overall).
These are outer work units, not a running-time estimate; a row
may require many suffixes. The budget is computed by grouping
successive starts of unchanged width.

The inclusive interval
\(15{,}000{,}000{,}000{,}000\le X\le15{,}001{,}000{,}000{,}000\)
for the smaller normalized coordinate gave zero hits across all
eight branches. This excludes only that interval.

~~~sh
python spine8-exact.py --check
python spine8-exact.py --budget
python spine8-exact.py --low 15000000000000 --high 15001000000000
~~~

Controls compare every output with direct enumeration on small
targets and include a planted large tangent point. Inclusive
intervals of the smaller normalized coordinate can run independently.
Removing interval bounds is the complete algorithm; the full sweep
was not run, and the unrestricted exponent-61 question is open.
