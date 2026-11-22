# An exact inverse that groups the low words by their reversed digits

**PROVED.** After fixing either flank and the parity of the other, an exponent-61 face has the form
\[
                 s y+h f(y)=M,\qquad s,h>0,\quad h\text{ odd}.
\]
Here \(M\) is the target after subtracting the known terms. The following inverse covers unbounded nonnegative addresses. It groups many low words in one phase query and checks the surviving addresses in the original integer equation. The opposite flank and the spine enter only through \(s,h,M\).

### A bounded error in a reversed binary product

Use \(B=2^b\), \(r=y\bmod B\), \(v=\operatorname{rev}_b(r)\), and \(z=\lfloor y/B\rfloor\). The reversed word includes leading zeros. Write bits \(r_i\) from the low end. The recurrence for \(f\) gives weight \(2^{i-j-1}+[i=j+1]\) to each pair \(i>j\) with \(r_i=1,r_j=0\). Expand the ordinary product \(r(B-1-v)/B\). Its terms at \(i>j\) have precisely the first weight. Its terms at \(i=j\) vanish. Its remaining terms sum to at most \(b/2\), since the geometric tail above any given low bit sums to less than one half. The adjacent corrections in \(f\) number at most \(b/2\). Hence
\[
             \left| f(r)-r(B-1-v)/B\right|\le b .
\]
In particular, writing \(\phi_2(y)=\sum_i y_i2^{-1-i}\) for binary radical inversion gives \(f(y)=y(1-\phi_2(y))+O(\log(2+y))\). A low bit fixes a leading digit of the normalised read at arbitrarily large height. [Inherited binary split](../chanpak_e61_fibonacci_target_20261122/SUFFIX-INVERSE.md)

Put \(D(v)=sB+h(B-1-v)\), \(P(v)=M/D(v)\), \(D_{\min}=sB\). The exact split and \(0\le f(z)+[r<B/2][z\text{ odd}]\le z\) give
\[
 -hb/D_{\min}\le P(v)-(z+r/B)
 \le h(M/D_{\min}+b+1)/D_{\min}.
\]
Indeed the difference times \(D(v)\) is \(h\) times the sum of the bounded product error, \(f(z)\), and the single boundary bit. The original nonnegative split gives \(z\le M/D_{\min}\).

### One lookup for a whole interval of reversed words

Let \(g=\gcd(s,h)\). The face is empty unless \(g\mid M\). Otherwise \(y\equiv a\pmod m\) for \(m=h/g\), with the unique residue \(a\) determined by the face equation. Put \(W=2^k\), \(k<b\), group \(v=v_0+i\) with \(W\mid v_0,\ 0\le i<W\), and write \(R=\operatorname{rev}_{b-k}(v_0/W)\), \(j=\operatorname{rev}_k(i)\). Thus \(r=R+(B/W)j\). Define
\[
 z_0=(a-R)B^{-1}\bmod m,\quad E_i=m^{-1}j\bmod W .
\]
The case \(m=1\) uses its sole residue. All allowed quotients at this word are
\[
               z=z_0+(mE_i-j)/W+\ell m,\qquad \ell\in\mathbb Z .
\]
Thus after multiplying \(z+r/B\) by \(W/m\) and removing the common offset \((W/m)(z_0+R/B)\), the entire residue choice is simply \(E_i\bmod W\). In particular the reversed product error did not cost a full integer of uncertainty in \(z\).

Use the increasing derivative \((W/m)P'(v)=WMh/(mD(v)^2)\). On consecutive groups enclose this derivative between \(u/W\) and \((u+d)/W\) by integer comparisons at the ends, subdividing as necessary. Bucket each position \(i\) by
\[
                        E_i-\lfloor u i/W\rfloor\pmod W.
\]
The mean value theorem bounds the unrounded slope error across a group by \(d\); rounding adds less than one. The two product-error bounds above, multiplied by \(W/m\), give the rest of the lookup window around \((W/m)(P(v_0)-z_0-R/B)\). The resulting bucket and its integer translate recover the candidate \(z\) by the displayed congruence. The algorithm checks the original split exactly. It never assumes that successive values of the digital height are monotone. High and low derivative bounds are evaluated with integer multiprecision. The exponent-61 implementation uses integer arithmetic throughout.

For the table width \(W\) the implementation chooses \(B\) at least \(\sqrt{MhW/m}/s\). This keeps the error from the unresolved high read to about one bucket. Large low words can then be processed \(W\) at a time; wider tables cost memory and derivative preprocessing. The table size, the number of blocks, and the number of reported candidates are distinct resources.

### Reproduce

Compile on a platform with 64-bit unsigned long, an unsigned 128-bit compiler type, and GMP:

~~~sh
c++ -O3 -std=c++17 farey_axis.cpp -lgmpxx -lgmp -o /tmp/farey_axis
python3 -B check_farey_axis.py /tmp/farey_axis
python3 -B check_farey_axis.py /tmp/farey_axis --spine 12 --bits 12
python3 -B check_farey_axis.py /tmp/farey_axis --spine 12 --bits 13 --fixed 1 --through-fixed 128
python3 -B check_farey_axis.py /tmp/farey_axis --spine 16 --bits 10
python3 -B check_farey_axis.py /tmp/farey_axis --spine 8 --bits 18
~~~

Homebrew GMP can be located with the compiler flags `-I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib`. The control compares planted solutions and their neighbors against literal enumeration, including a nonunit slope. A shard covers a stated block interval in **each** parity. An empty shard certifies only its own interval. [Exact inverse](farey_axis.cpp) · [Face driver and controls](check_farey_axis.py)

**FINITE-EXACT.** Both coordinate axes at spine 16 are empty. At spine 12 the complete left-only inverse tested 62,914,400 candidates and found none. The other axis fails its initial congruence. Fixing either address to any value 1 through 128 and leaving the other unbounded also gives no hit. The latter operation visits 12,009,472 blocks and checks 84,914,463 candidates with the right flank fixed. The initial congruences already remove the positive fixed left flanks through 128.

At spine 8 the two complete axes are also empty. With the right address fixed to zero, the counts of bucket candidates for the even and odd left addresses were 268,419,115 and 268,408,335. With the left address zero they were 536,908,040 and 536,877,244. The respective numbers of blocks per parity were 33,554,432 and 67,108,864.

At spine 8 the entire strips with either address fixed to any integer from 1 through 128 are empty too. With the right address fixed, this used 287,834,112 blocks and 2,203,565,481 candidates. With the left address fixed, it used 160,432,128 blocks and 1,216,289,690 candidates. Running `check_farey_axis.py` as above with `--spine 8 --bits 16 --fixed 1 --through-fixed 128` reproduces both strips.

**FINITE-EXACT.** On spine 8, even *one* flank having only a single run of ones excludes the target, with arbitrary trailing zeros and no restriction on the opposite word. For each side the preceding strips cover all words of bit length at most seven, as well as zero. For longer words up to the absolute cutoff, write the fixed address as \((2^w-1)2^{j-w}\) and run the exact one-variable inverse for the other address. This adds 4,533 fixed left words and 4,439 fixed right words. In each orientation only 83 and 337, respectively, require the reversed-bucket inverse; direct congruence progressions suffice for the others. Run `python3 -B fixed_run_inverse.py /tmp/farey_axis 8 left 8 123 --bits 15` and repeat with `right`. [One-run and unrestricted-flank search](fixed_run_inverse.py)

At spine 4, every fixed one-run flank of bit length at least seven likewise misses, regardless of its opposite word. This covers 5,131 left words (129 bucket searches) and 5,031 right words (496 bucket searches). Use the same commands with spine 4 and minimum length 7. Continuing the left search through lengths 2–6 and the right search through lengths 5–6 found no hit in 20 and 11 further words, respectively (use `--bits 17`). Thus a nonempty left word with a single one-run would have to be the integer 1, while a right word with a single one-run would have to be below 16. The opposite word is unrestricted throughout.

### How sparse can a locally universal inverse be?

**PROVED.** For fixed positive real numbers \(s,h\), and any bounded parity correction \(k_j\), the number of nonnegative addresses of height at most \(X\) on one face is
\[
 \#\{y:sy+h f(y)+k_{y\bmod2}\le X\}
   =\frac{X}{h}\log(1+h/s)+O_{s,h,k}(\sqrt X).
\]
This counts addresses with multiplicity and so also bounds the size of the distinct image from above.

Indeed, for \(y=O(X)\) the reversed product approximation above gives height \(y[s+h(1-\phi_2(y))]+O(\log(2+X))\). Group by the bottom \(b\) digits with \(B=2^b\) between \(\sqrt X\) and \(2\sqrt X\). In each residue class the reversed fraction lies in an interval \([v/B,(v+1)/B]\), with each \(v\) used once. Replacing the height condition by \(y\le X/[s+h(1-v/B)]\) moves the cutoff by at most \(O(X/B+\log X)\). Counting an arithmetic progression of spacing \(B\) in each interval costs a total \(O(B+X/B+\log X)\), including rounding. The average over \(v\) is a Riemann sum for the integral of \(1/[s+h(1-u)]\) from zero to one, with error \(O(1/B)\). Integrating proves the formula.

Across the union of both axes of the eight charts, the upper density *among odd integers* of exact lower counts is less than **0.000004353**. Indeed \(\log(1+h/s)/h<1/s\), and twice the sum of \(1/c_n+1/t_n\) over the eight chart rows is less than that bound. Already the left-only face of any one chart realizes every odd phase modulo every fixed modulus. Thus these local tests can all succeed even for a family whose absolute image misses more than 99.99956 percent of odd heights. This density estimate makes no prediction for a particular height or for the whole two-flank chart. [All-modulus inverse](LOCAL-WORDS.md) · [Two-flank height frequency](ABSOLUTE-WORDS.md)
