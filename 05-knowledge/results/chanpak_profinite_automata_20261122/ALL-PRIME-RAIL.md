# One fixed physical seam rail passes every modulus and misses the exact target

**Proved.** With \(T=E_{61}-1\), where \(E_{61}=2^{60}(2^{61}-1)\), there is a single fixed literal one-arc seam family \((H_n,H_n+2)\), \(n\ge1\), such that for **every modulus** \(M\), arbitrarily large \(n\) satisfy \(H_n\equiv T\pmod M\). No \(n\) satisfies \(H_n=T\). This is one run-length-compressible family, not a different choice of physical source for each prime and not a family obtained by discarding the target. It does not exclude other realizations of the exact target.

## Physical rail and exact non-hit

Use the six-vertex balanced pair with outgoing masks \((4,17,42,3,13,27)\) and \((4,17,34,7,13,27)\), tag mask 22, the changed arc at \((2,3)\), and left address \(x=2\). The [established physical star formula](../chanpak_gap_star_ports_20261122/README.md) gives the lower count \(493+769y+29p(y)\) on even positive right addresses \(y\), and its literal companion is two larger. Both tournaments are strong. Here \(p(2z+1)=p(z)\), \(p(2z)=p(z)+z+[z\text{ odd}]\) for \(z>0\), and \(p(0)=0\).

Start with \(y_0=73\) and append \(10\) exactly \(n\) times, so \(y_{n+1}=4y_n+2\). This gives \(p(y_{n+1})=p(y_n)+2y_n+2\) and \(p(73)=32\). Solving the two recurrences gives
\[
 9H_n=A4^n+174n+C,\qquad A=522665,\quad C=-4643.
\]
This expression is an integer formula after dividing by 9 at every \(n\ge0\); the physical even-address seam uses only \(n\ge1\). Its order is \(16+2n\). The lower counts increase strictly, and direct arithmetic gives
\[
 H_{52}=1177878401721400035792487788156173005
 <T=2658455991569831744654692615953842175
 <H_{53}=4711513606885600143169951152624690571 .
\]
Thus the family never hits the target.

## The prime-by-prime lifts

At 2, require \(n\) large enough to kill \(4^n\bmod 2^s\). Then prescribing \(H_n\equiv T\pmod{2^s}\) is the affine congruence \(174n\equiv9T-C\pmod{2^s}\). It specifies one residue modulo \(2^{s-1}\), since \(174=2\cdot87\) and the right side is even. At precision at least 2 this residue is odd.

At 3, \(H_{15}\equiv T\equiv0\pmod3\). For \(j\ge1\) one has
\[
 \frac{4^{3^j}-1}{3^{j+1}}\equiv7\pmod9.
\]
The base \(j=1\) is \(63/9=7\); cubing a number \(1+3^{j+1}a\) shows the residue remains unchanged. Therefore for \(n\equiv0\pmod3\) and any digit \(d=0,1,2\),
\[
 9(H_{n+3^j d}-H_n)
 \equiv 3^{j+1}(7A4^n+58)d \pmod{3^{j+3}}.
\]
Here \(4^n\equiv1\pmod9\) and \(7A+58\equiv6\pmod9\). Dividing the displayed congruence by 9 gives \(H_{n+3^j d}-H_n\equiv2\cdot3^j d\pmod{3^{j+1}}\). This is a unit derivative on the selected residue and lifts it to arbitrary ternary precision. It explicitly accounts for the denominator 9.

At an odd prime \(p\ne3,29\), put \(h=\operatorname{ord}_p(4)\). For any torsion phase \(t\bmod h\), the map \(k\mapsto H_{t+hk}\) is a bijection modulo \(p^s\) for every \(s\). In fact for \(n\ge0\), \(j\ge0\), and any nonnegative digit \(d\), binomial lifting gives
\[
 (4^h)^{p^j d}-1\equiv0\pmod{p^{j+1}},\qquad
 H_{n+hp^j d}-H_n\equiv (174h/9)p^j d\pmod{p^{j+1}}.
\]
The coefficient is always a unit, independent of \(n\) and the already chosen digits. In particular at \(p=7\), \(h=3\). Use the phase \(15\equiv0\pmod3\) selected at 3, and observe \(H_{15}\equiv T\equiv0\pmod7\). The leading \(n\)-digit is then \(15\equiv1\pmod7\), whichever higher 7-digits are chosen.

At 29 use \(h=14\) and \(n\equiv15\pmod{14}\). It has \(H_n\equiv H_{15}\equiv T\equiv2\pmod{29}\) because the drift is divisible by 29. Put \(a=(4^{14}-1)/29\). The derivative after removing the already supplied factor 29 is
\[
 \frac{A4^n a+6\cdot14}{9}\equiv2\pmod{29},
\]
independent of \(n\) in this torsion phase: \(4^n\equiv4^{15}\pmod{29}\). One direct control is \((H_{29}-H_{15})/29\equiv2\pmod{29}\). Expanding \((4^{14})^{29^j d}\) modulo \(29^{j+2}\) gives \(H_{n+14\cdot29^j d}-H_n\equiv2\cdot29^{j+1}d\pmod{29^{j+2}}\). This is a unit lift of every higher 29-digit.

## Exact cross-prime compatibility

First impose at least binary precision 2 and moduli 3, 7, and 29, even if not requested. Solve separately at each of these four primes to the desired higher precision. The base constraints agree: the 29 torsion condition \(n\equiv15\pmod{14}\) just requires odd \(n\) and \(n\equiv1\pmod7\); these are the locked phases at 2 and 7. The 7 torsion \(n\equiv0\pmod3\) is the locked ternary phase. CRT therefore combines these four lifts.

Process all other requested primes \(p\) in increasing order. Before processing \(p\), fix a torsion phase modulo \(\operatorname{ord}_p(4)\) compatible with every old congruence, e.g. the current CRT representative modulo that order. Its prime divisors are all strictly smaller than \(p\). The current CRT modulus has no \(p\)-factor: prior nonexceptional primes were smaller and the exceptional clauses have only 2, 3, 7 and 29 as their factors. The good-prime lift works in **every** torsion phase, so it sets the desired new \(p^s\) condition without changing any old constraint. Combine by generalized CRT. Finally increase \(n\) by any large multiple of the final CRT period to meet the nonunit binary size condition and the physical requirement \(n\ge1\).

This proves one family approximates the same integer at all primes jointly, not merely at every prime separately. It also provides a coherent profinite sequence by requesting all prime powers up to growing thresholds and choosing the representatives increasingly. No local residue, regardless of the number of primes combined, can witness the absolute miss on this particular family.

## The whole open target class

The same proof applies to any profinite lower height in the class
147 modulo 2436. This is the class of the target T, and encodes the
four starting phases 3 mod 4, 0 mod 3, 0 mod 7, and 2 mod 29.
At every other prime the construction accepts every residue. The
linear drift is essential to the proof: at a good prime its first
difference is a unit times the lifted prime power, while the
torsion-preserving first difference of the exponential term is
divisible by one more power of that prime.

## Compressed exact witnesses

For the fifteen usual Ogg primes \(2,3,5,7,11,13,17,19,23,29,31,41,47,59,71\), the replay finds the **least** width
\[
 n=5683352354064687
\]
in this one fixed rail, with total shift period \(9713789940653138340\). For the product of their squares, the least width is \(47754852698807371455951618379739667\), with shift period \(90380953075798947196807652194614900\). These are rail-restricted minima, not component-order bounds for other physical witnesses.

To certify minimality, enumerate the free lower torsion digits rather than using only one arbitrary CRT branch. The orders of 4 at \(5,7,11,13,17,19,23,29,31,41,47,59,71\) are respectively \(2,3,5,6,4,9,11,14,5,10,23,29,35\). The 29 phase is uniquely \(n\equiv1\bmod14\) (the geometric coefficient is a unit). The ternary root is uniquely \(n\equiv0\bmod3\), as direct evaluation at \(0,1,2\) shows. Every other prime's lift is unique in its selected torsion phase. At first powers the only free choices are \(n\bmod4\in\{1,3\}\), \(n\bmod9\in\{0,3,6\}\), and \(n\bmod29\) arbitrary. These give 174 complete periodic classes and the minimum is taken over all their first positive representatives. At squared precision requirements fix the \(n\bmod9\) and \(n\bmod29\) digits, leaving exactly 2 classes from \(n\bmod4\). Both cases exclude \(n=1,2\) already modulo 3; no binary exceptional small index is lost. The physical orders are enormous, but the word description is prefix \(1001001\) and one repetition instruction \((10)^n\). Verification uses modular exponentiation, never construction of those vertex lists.

For arbitrary \(M\), the constructive CRT period satisfies \(L\le1218 M^2\): the extra initial powers of \(2,3,7,29\) multiply the requested modulus by at most \(2\cdot3\cdot7\cdot29=1218\), and the product of the new torsion orders is at most the product of their primes, which divides \(M\). A binary safety cutoff requires at most one more period. Thus the compressed witness has \(O(\log M)\) bits in its repetition exponent; this is not logarithmic physical order. The converse height observation is immediate but useful: when \(M>T\), every non-hit congruence witness in this rail must have \(H_n\ge T+M\), since a lower positive number congruent to \(T\) would be negative. Refining the modulus therefore forces genuine escape in the absolute path count.

Replay: e61-all-primes-check.py, Python standard library. It checks the binary recurrence, derivative controls, exact non-hit bracket, and all requested simultaneous phases.
