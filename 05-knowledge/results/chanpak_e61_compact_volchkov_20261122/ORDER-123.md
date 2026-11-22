# Absolute order-123 receiver screens at exponent 61

This follows the three normalized rooted receiver charts left in
[SELECTOR.md](../chanpak_e61_fibonacci_target_20261122/SELECTOR.md). Let \(E=2^{60}(2^{61}-1)\) and \(N=E-1\).
Each chart tests the exact integer equation
\[
 N=K_{y\bmod 2}+cy+h f(y),
 \quad f(0)=0,\quad
 f(2u+1)=f(u),\quad f(2u)=f(u)+u+[u\ {\rm odd}].
\]
The test at order 123 fixes the binary length of \(y\). The three
charts, in the order of the preceding selector, are:
\[
 (c,h,K_0,K_1,k)=(81,25,43,68,115),\quad
 (411,75,367,392,113),\quad(1021,125,915,940,112).
\]

**Complete exclusion of the triangle receiver.**
For the middle chart, \(c,h\) are both divisible by 3, but
\(N\equiv0\pmod 3\) and neither \(K_0\) nor \(K_1\) is divisible
by 3. This excludes every address, not only addresses of the
required length. The two other charts have \(\gcd(c,h)=1\).
The binary recurrence shows \(f(y)\) is even: in its even branch
the increment \(u+[u\text{ odd}]\) is even. Their exact necessary
address classes are therefore respectively \(y\equiv22,47\pmod{50}\)
and \(y\equiv60,35\pmod{250}\), listed even then odd.

**Complete sublanguage exclusions for the two other charts.**
No address ending in 40 or more equal bits gives the exact target.
Neither does any address with at most nine runs of equal bits.
These are finite exact sublanguage outcomes; the two charts are
not fully excluded.

## Terminal carry with the output parity constraint

For a given terminal run length \(\ell\le60\), put \(A=2^\ell\),
\(P=E/A\), and write \(y=Az+b(A-1)\), \(b\in\{0,1\}\).
For final ones the recurrence gives
\[
 P=c(z+1)+(K_1+1-c+h f(z))/A.
\]
For final zeros it gives
\[
 P=(c+h)z+(K_0+1+h(f(z)-z+[z\text{ odd}]))/A.
\]
An exact solution first bounds \(z\) by the fixed width and
\(cy\le N-K_b\). If \(Z\) denotes the resulting upper bound,
the residual numerator is between \(K_1+1-c\) and
\(K_1+1-c+hZ\) for ones, and between \(K_0+1-hZ\) and \(K_0+1\)
for zeros, using \(0\le f(z)\le z-[z\text{ odd}]\).
It must be a multiple of \(A\). These bounds narrow \(z\) again.
Finally impose \(c(Az+b(A-1))+K_b\equiv N\pmod{2h}\) before
testing any residual, then check the full integer equation.
Every address with a longer homogeneous terminal run belongs
to the \(\ell=40\) test, so one run length suffices for the result.

| Receiver slope \(c\) | Final bit \(b\) | Exact tests after interval and congruence | Hits |
|---:|---:|---:|---:|
| 81 | 0 | 256,117,314 | 0 |
| 81 | 1 | 335,165,867 | 0 |
| 1021 | 0 | 0 | 0 |
| 1021 | 1 | 2,109,494 | 0 |

## Exact elimination of the final two run lengths

Appending a run of \(s\) ones to \(x\) leaves \(f(x)\) unchanged;
appending \(s\) zeros adds \((2^s-1)x+[x\text{ odd}]\).
For a fixed prefix \(x>0\) of digital read \(F=f(x)\),
suppose the remaining \(L\) positions are two alternating runs,
of lengths \(a,L-a\) with \(1\le a<L\). Put \(B=2^L\).
If the runs are zeros then ones, the target equation is
\[
 N-K_1-c(Bx-1)-h(F-x+[x\text{ odd}])
   = hx\,2^a+c\,2^{L-a}.
\]
If they are ones then zeros, it is
\[
 K_0+cB(x+1)+h(F+B(x+1)+2)-N
   = h(x+1)\,2^a+(c+h)\,2^{L-a}.
\]
Thus each case is a positive-integer equation
\(\Delta=U2^a+V2^{L-a}\). If the two summands have distinct
2-adic valuations, their minimum equals \(v_2(\Delta)\);
each possible minimum specifies \(a\). If their valuations
are equal, their equality specifies \(a\) instead.
At most three proposed splits need checking by equality.
This optimization avoids scanning the final run boundary,
and it does not assume any valuation condition is sufficient.

The C++ run search is also replayable in an independent direct
form: append the runs one at a time with these recurrences,
allowing every positive run length, and require exactly \(k\)
bits. Its only numerical pruning uses \(f(x)\) as a lower
bound for every extension and \(f(x)+(2^L-1)(x+1)+1\) as an
upper bound for any \(L\)-bit extension, by the binary splitting
identity and \(f(r)\le r\). The complete direct 7-run scan
reports zero hits on both remaining charts. The supplied
optimized scan uses the same outer enumeration but the two-run
valuation inverse at the last branch and excludes the eighth and ninth
run layers too. The ninth layer tests 5,110,728,658 and 2,988,293,772
split proposals in the singleton and four-vertex charts, respectively,
with zero exact hits. The middle chart is excluded for all odd
exponents since the same modulo-three argument applies to each.

These tests are for the fixed six-vertex core and the order-123
receivers under discussion. They do not decide arbitrary
addresses on the two unresolved charts, smaller order
receivers, or the unrestricted exponent-61 problem.

Replay the searches by compiling the companion C++ sources in
C++17 mode, then run the carry search with argument 40, the direct
run search with argument 7, and the optimized run search with
argument 9. The tests use signed 128-bit integers; all heights and
intermediate terms in the stated scans are below \(2^{127}\).

## A high-prefix inverse in reverse-low coordinates

Fix a positive high prefix \(x\) in \(y=Bx+r\), with \(B=2^t\),
\(0\le r<B\), and put \(R=\operatorname{rev}_t(r)\) (including
padding). By the binary splitting identity the exact equation is
\[
 N=cBx+h(B-1)x+h f(x)-hxR+
    \bigl(K_{r\bmod2}+cr+h f(r)+h[x\text{ odd}][r<B/2]\bigr).
\]
The parenthesized expression lies in
\([K_0,\,K_1+(c+h)(B-1)+h]\). Therefore the compatible \(R\)
lie in one explicitly known integer interval of length at most
\(1+(K_1-K_0+(c+h)(B-1)+h)/(hx)\), intersected with
\([0,B-1]\). For each \(R\) there is only one \(r\), obtained by
reversing again. First impose the congruence modulo \(2h\);
then replay the whole height. This can be much cheaper than
running over every low word \(r\), if the high prefix is long
and prescribed. It does not select the high prefix itself.

For the common prefix of a leading one followed by 44 zeros,
the singleton chart has an empty reverse interval. The
four-vertex chart has a reverse interval of length \(76{,}906{,}758\)
with \(615{,}253\) congruence survivors and no exact match.
Necessarily this also excludes any longer initial zero run.
Compile the separate prefix C++ replay and run with argument 44.

For the four-vertex chart the five final ones prescribed by the order bound give \(y=32z+31\), \(f(y)=f(z)\). Division and modular reduction then force \(z=47\bmod125\) and \(f(z)=142\bmod2042\). For the two singleton parities the conditions are \(y=22\bmod50,\ f(y)=116\bmod162\) and \(y=47\bmod50,\ f(y)=34\bmod162\). These use that \(f\) is even and do not exhaust the remaining integer windows.

A concurrent exact search on [reflected and reverse-complemented halves](../chanpak_e61_fibonacci_target_20261122/REFLECTED-SUFFIX.md) also excludes signed half-differences with binary nonadjacent weight at most five, and gives a distance-multiset (turnpike) inverse. That is a different sublanguage, independent of the run and prefix bounds here. Neither test proves exhaustion of arbitrary suffixes.
