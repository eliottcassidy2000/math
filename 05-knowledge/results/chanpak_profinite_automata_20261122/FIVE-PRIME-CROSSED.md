# Five independent target primes on one strong single-arc rail

Use the [adjacent-opposite crossed-contact construction](../chanpak_hadwiger_planar_folds_20261122/TWO-APEX.md), now
with far-left common tags \(W=10110110\,1^k\,0011\), far-right
reverse complement, and central contacts \(00,10,01,00\).
The two weighted addresses differ by one. Both endpoints are
strong and reverse just the apex arc, with order \(2k+30\).
For \(X=2^k,\ Z=3^k\) their exact midpoint is
\[
 M_k=3412542-982640804X+2328668912Z+71433464766X^2
       -1453657649592XZ+8410552316658Z^2.
\]
The [companion replay](check-five-prime-crossed.py) counts both orientations by the 12-state
transfer at widths 1 through 9. The two growing ports can
each have at most three runs, so the polynomial lies in the
six-dimensional basis \(1,X,Z,X^2,XZ,Z^2\); six widths
certify it by the Vandermonde argument.

**Theorem.** For \(s\ge6,t\ge2,u\ge2,r\ge1,j\ge1\), every
residue \(R\) with
\[
                    R\equiv0\pmod{64},\qquad
                    R\equiv1\pmod{9\cdot25\cdot7\cdot11}
\]
is realized as \(M_k\bmod 2^s3^t5^u7^r11^j\) at infinitely
many widths \(k=22+60v\). In particular this realizes the
five simultaneous local phases of \(E_{61}\) at arbitrary
precision in each of the five primes.

Write \(P(Z)=M(0,Z)\), \(Q(X)=M(X,0)\). For sufficiently
large \(k\), the omitted terms vanish respectively at 2 and
3. At 2, \(Z=3^{22}(3^{60})^v\) runs over the full \(9\bmod16\)
coset at every dyadic precision, since the step has valuation
\(v_2(3^{60}-1)=4\). The polynomial
\[
 P(9+16w)/64=10644932798157+37848067592189w+33642209266632w^2
\]
permutes every dyadic modulus by its odd divided differences.
At 3, \(X=2^{22}(2^{60})^v\) runs over the \(7\bmod9\)
coset, since \(v_3(2^{60}-1)=2\). Here \(Q(7)\equiv1\bmod9\)
and \(Q'(7)\equiv1\bmod3\); it bijects this coset with
\(1\bmod9\) at every precision by the same divided-difference
test. These two inverses select \(v\) modulo \(2^{s-4}\)
and \(3^{t-2}\), respectively.

For the other primes use the base expansion
\(M_{22+60v}=\sum_b c_b b^{22}(b^{60})^v\).
For \(p=5,7,11\), put \(\delta_p=2,1,1\), respectively.
For every base \(b=1,2,3,4,6,9\), one has \(b^{60}\equiv
1\bmod p^{\delta_p}\). Binomial expansion gives
\[
 M_{22+60(v+d p^i)}-M_{22+60v}
       \equiv d p^{\delta_p+i}D_p
                            \pmod{p^{\delta_p+i+1}},
 \qquad D_p\equiv (M_{82}-M_{22})/p^{\delta_p}\pmod p .
\]
The checkpoints are:

| \(p\) | modulus | \(M_{22}\) | \(M_{82}\) | \(D_p\) |
|---:|---:|---:|---:|---:|
| 5 | 125 | 1 | 51 | 2 |
| 7 | 49 | 8 | 36 | 4 |
| 11 | 121 | 56 | 100 | 4 |

Each starting residue is 1 modulo \(p^{\delta_p}\), and
each \(D_p\) is a unit, giving unique successive digits of
\(v\) modulo \(5^{u-2},7^{r-1},11^{j-1}\).
Apply CRT to the five pairwise coprime exponent moduli,
and add arbitrary common periods to make \(k\) large for
the two terms initially discarded. The Euclid center satisfies
the requisite initial phases since the orders of 2 modulo
9, 25, 7 and 11 divide 60, and \(E_{61}=2^{60}(2^{61}-1)\).

**Exact exclusion for this rail.** Modulo 19 the target
is \(15\), but no residue \(k\bmod18\) produces that
midpoint. Thus the five-prime local success is complete
at all indicated precisions while this physical rail has
no absolute solution at any width. The unrestricted
tournament problem is untouched.
