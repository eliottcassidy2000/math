# The last-four-digit table is an exact 5-adic fold

Write \(E_p=2^{p-1}(2^p-1)\) whether or not \(2^p-1\) is prime. The [GIMPS list](https://www.mersenne.org/primes/) supplies the known Mersenne exponents. The following is not a frequency assertion about those exponents; it proves the governing rule for **all** odd exponent residues.

For \(k\ge1\) and \(p\ge k+1\), \(E_p=0\bmod2^k\). On the other CRT factor put \(x=2^p\) and write
\[
 E_p=(x^2-x)/2,\quad
 E_p-E_r=(2^p-2^r)(2^p+2^r-1)/2,\quad
 8E_p+1=(2^{p+1}-1)^2 .
\]
The multiplicative order of 2 modulo \(5^k\) is \(4\cdot5^{k-1}\) by LTE from \(2^4-1=15\). Unlike modulo \(5^k\) alone, the powers modulo the even modulus \(10^k\) have a **preperiod**. The above period holds after stabilization.

* If \(p\equiv1\bmod4\), then \(x\equiv2\bmod5\); the derivative of \((x^2-x)/2\) is a unit modulo 5. As the exponent runs through this sheet modulo \(4\cdot5^{k-1}\), the center runs exactly once through every class \(E_p\equiv1\bmod5\) modulo \(5^k\), by simple Hensel lifting. CRT gives the decimal endings in \(16\bmod20\) for \(k\ge2\). The last digit is 6 and the tens digit is odd.
* If \(p\equiv3\bmod4\), then \(x\equiv3\bmod5\) and \(2x-1\) runs through all multiples of 5 modulo \(5^k\). The centers are exactly the fold of this ball by the square test \(8E_p+1\in(5\mathbb Z/5^k)^2\). Moreover
  \[
      v_5(8E_p+1)=2\left(1+v_5((p+1)/4)\right).
  \]
  At two decimal places, \(E_p\equiv3\bmod25\) and \(E_p=0\bmod4\), hence the suffix is always **28**. The very small exponents \(p=2,3\) check directly.

For \(k\ge2\), the period \(4\cdot5^{k-1}\) is least even on the stabilized odd-exponent sequence: a smaller even period must preserve the \(p\equiv1\bmod4\) sheet, whose last digit distinguishes it from the other, and injectivity of that sheet then forces the full period.

At \(k=4\), the \(p=3\bmod4\) sheet has precisely eleven possible numerical suffixes: the zero square modulo 625 and the ten units squares modulo 25 multiplied by 25, each transported through the CRT rule above. The \(p=1\bmod4\) sheet has 125 possible numerical suffixes (they also vanish modulo 16). **The four-digit parity string itself has least even exponent period 100, not 500.** The parity of digit \(j\) is \(\lfloor E_p/10^j\rfloor\bmod2\), hence depends only on \(E_p\bmod 2\cdot10^j\). Four digits therefore need only \(E_p\bmod2000\), where the odd factor is \(125\). This gives period 100. Test the proper-divisor shifts 20 and 50 (every proper divisor of 100 divides one of those), with exact witnesses in the replay. Its six possible words are \(EOOE,EOEE,EEOE,OOOE,OOEE,OEOE\). The possible words do not claim equal frequencies, nor a distribution theorem among Mersenne primes. To restrict to large prime exponents one may remove the residue classes \(p=0\bmod5\) from the period table (100 suffices for parity and 500 for digits).

At \(p=61\) the center ends **2176**, with parity \(EOOE\), and the endpoint suffixes are 2175 and 2177. A true physical seam has to meet that decimal guard, but matching it does not solve the absolute suffix equation. The useful split is between a regular sheet admitting exact inverses and a singular but explicitly parameterized square image. This is a numerical analogue of retaining a contact's side information rather than treating the observed parity pattern as its full state.

For arbitrary trailing length \(d\ge2\), the same argument gives **least** parity period \(4\cdot5^{d-2}\) on eventual odd exponents. Indeed the parity word depends only on \(E_p\bmod2^d5^{d-1}\), giving this upper bound. Its factor 4 is necessary from the tens digit. For \(d\ge3\) take \(p=1\bmod4\) and shift repeatedly by \(4\cdot5^{d-3}\). Across five consecutive shifts, \(x=2^p\) runs through all lifts of its fixed class modulo \(5^{d-2}\) to modulo \(5^{d-1}\). The unit derivative on this sheet makes the five \(E_p\) classes all the lifts too. CRT with the shared residue at \(2^d\) yields five equally spaced residues modulo \(2\cdot10^{d-1}\). They occupy both halves, so the parity of the highest retained digit changes. This kills the period divided by 5. The argument for \(d=2\) needs only the tens digit. This strengthens the four-digit verification above without a distribution assumption about primes.

The concurrent [decimal phase packet](../chanpak_decimal_kernel_20261122/DECIMAL.md) independently obtains these periods and goes further: it gives a Hensel/discrete-log inverse with at most two exponent progressions for any numerical suffix, as well as the complete table for the known 52 Mersenne exponents.
