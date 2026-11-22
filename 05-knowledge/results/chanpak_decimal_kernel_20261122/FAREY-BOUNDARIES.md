# Local universality inside prescribed word boundaries and absolute witnesses at order 118

**PROVED.** In each of the eight canonical Farey charts with left tag mask 12, prescribe both a prefix and a suffix of the left address and both a prefix and a suffix of the right address. Prefixes of nonempty canonical binary addresses start with one. For every positive integer \(M\) and odd integer \(T\), there is a pair of finite addresses with exactly those prescribed boundaries whose physical arc reversal has readings \(T,T+2\bmod 2M\). Here a boundary specifies digits and leaves the length of the intervening word free. One can fix the complete right address from its two boundaries alone, before knowing \(M\) or \(T\).

Write \(H(x,y)\) for the lower height and \(h,q,c,t,e,r,k\) for its usual seven readings. Each of the eight charts has \(h\) odd and \(\gcd(t,h)=1\). For the right prefix of value \(p\) and suffix of width \(d\) and value \(s\), insert \(w=\lceil\log_2(h+1)\rceil\) middle bits with value \(v\). The address \(y=(2^wp+v)2^d+s\) can be made divisible by \(h\) since \(2^d\) is invertible modulo \(h\) and the range of \(v\) meets every residue. This \(y\) has the requested terminal digits and \(\gcd(t+qy,h)=1\).

For the left suffix of width \(b\ge1\) and value \(i\), put \(B=2^b\), allow leading zeros in word reversal, and write \(x=Bu+i\). An empty suffix specification may be extended by one arbitrary bit. By binary concatenation, when \(u\) is even,
\[
 H(Bu+i,y)=H(i,y)+s'u+h f(u),\qquad
 s'=B(t+qy)+h(B-1-\operatorname{rev}_b(i)).
\]
The evenness requirement eliminates the single boundary carry in the digital reading. Also \(\gcd(s',h)=1\). The baseline is an odd ordinary tournament count. Consequently there are residues \(a,z\bmod M\) solving
\[
 s'a+hz=(T-H(i,y))/2\bmod M.
\]
From the prescribed left prefix the [register inverse](LOCAL-WORDS.md) supplies \(u=2a\bmod 2M\) and \(f(u)/2=z\bmod M\) simultaneously. Appending the prescribed left suffix completes the pair. This proof preserves both terminal parities, including the case of two long runs of terminal ones. Repeated odd primes need no separate compatibility test. At powers of two the address residues forced by the suffixes remain fixed and the \(h f(u)\) register supplies the missing height phase. The proof uses \(O(\log^2(2M))\) filler bits on the left and \(O(\log(2h))\) on the right. The physical height can far exceed \(T\).

**FINITE-EXACT.** For \(E=2^{60}(2^{61}-1)\), the following two unit arc flips in the length-four Farey chart lie on opposite sides of \(E\). If the indicated midpoint offset is \(v\), the two actual ordinary Hamilton path counts are \(E+v-1\) and \(E+v+1\).

| Left address | Right address | Tournament order | Midpoint minus \(E\) |
| ---: | ---: | ---: | ---: |
| 999325040886328 | 1159777636980293 | 118 | \(-5\,378\,156\,638\) |
| 1000669119751650 | 1158219846718188 | 118 | \(+464\,575\,890\) |

For both witnesses the lower chart is \((2465,2293760,1344270,712944,1740,986,423385)\), and reversing the marked middle arc adds two. The order is \(17+\operatorname{bitlength}(x)+\operatorname{bitlength}(y)\). To find them, partition both coordinates by their two least significant bits, apply the [rational tangent search](FAREY-HYPERBOLA.md) to each resulting hyperbola, and compare candidate readings with \(E-1\) in integer arithmetic. These are individual certificates. Whether this chart attains midpoint \(E\) remains open.

Run the command below from this directory. It reconstructs both physical charts, verifies both endpoints of each flip, tests every odd target modulo 4, 8, 16 and 32 with four choices of boundary words in every chart, and tests repeated odd primes both separately and together. [Integer verifier and boundary-word constructor](farey_boundaries.py)

~~~sh
python3 -B farey_boundaries.py
~~~
