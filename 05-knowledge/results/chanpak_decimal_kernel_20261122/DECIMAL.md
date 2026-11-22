# The decimal fibers of a perfect-shaped center

Let \(E_p=2^{p-1}(2^p-1)\) for any integer \(p\ge2\). Primality of \(2^p-1\) makes this an even perfect number. It is unnecessary for the following congruences. Write a decimal suffix with leading zero padding. For \(p\ge k+1\), its last \(k\) digits are zero modulo \(2^k\), so they are determined by \(E_p\bmod 5^k\) via CRT. The parity of digit \(j\), starting at \(j=0\) on the right, is determined by the residue modulo \(2\cdot10^j\). Hence the parities of all \(k\) digits require only \(E_p\bmod 5^{k-1}\) and the known zero modulo \(2^k\).

**PROVED.** For odd \(p\ge k+1\) and \(k\ge2\) the least positive **even** exponent shift preserving the exact \(k\)-digit suffix for all such \(p\) is \(4\cdot5^{k-1}\). For the vector of digit parities the least such shift is \(4\cdot5^{k-2}\).

Indeed 2 has order \(4\cdot5^{j-1}\) modulo \(5^j\), from \(v_5(16^d-1)=1+v_5(d)\). These orders give the upper bounds. For sharpness put \(x=2^{p-1}\) and \(F(x)=x(2x-1)\). On \(p\equiv1\pmod4\) the variable \(x\) traverses every unit congruent to 1 modulo 5; for two such units,
\[
 F(x)-F(y)=(x-y)(2x+2y-1),\qquad v_5(F(x)-F(y))=v_5(x-y).
\]
Thus this branch traces exactly the whole residue class 1 modulo 5 at every higher precision, preserving all pairwise valuations. A putative smaller period divides the displayed period after taking a gcd on the odd-exponent progression. Halving the displayed period interchanges the two classes modulo 4. It changes even the tens parity, since \(E_p\equiv16\pmod{20}\) on the first branch and \(E_p\equiv8\pmod{20}\) on the second. Dividing the exact period by five loses a unit of 5-adic precision on the first branch. For the parity period, when \(k\ge3\), put \(d=4\cdot5^{k-3}\). Five exponents spaced by \(d\) on the first branch visit all five different lifts modulo \(5^{k-1}\) of the same lower residue. By CRT their values modulo \(2\cdot10^{k-1}\) are five equally spaced points on that residue circle. Both halves of the circle are met, so their highest digit parities cannot all agree. This eliminates every proper divisor in both period claims.

## Why 28 and 1328 recur

Put \(U=2^{p+1}-1\). Completing the square gives
\[
                      8E_p+1=U^2 .
\]
For \(p\equiv3\pmod4\) one has \(5\mid U\) and thus \(E_p\equiv3\pmod{25}\). Together with \(4\mid E_p\) this forces the last two digits **28**. On \(p\equiv1\pmod4\) the congruence \(E_p\equiv16\pmod{20}\) forces a last digit 6 and an odd tens digit. This proves the observation for every odd exponent \(p\ge3\), with \(p=2\) giving 6 separately. More precisely, residues \(p=1,5,9,13,17\bmod20\) give last two digits 76, 96, 16, 36, 56 respectively. All residues \(3\bmod4\), including 11 modulo 20, give 28.

For \(p\equiv3\pmod4\), \(v_5(U)=1+v_5(p+1)\). Therefore if \(5^{\lceil k/2\rceil-1}\mid p+1\), the critical branch satisfies \(8E_p+1\equiv0\pmod{5^k}\); the same valuation identity proves necessity. At \(k=4\) the CRT lift of \(-1/8\pmod{625}\) which is zero modulo 16 is **1328**. Hence, for odd \(p\ge5\), the last four digits are **1328 if and only if \(p\equiv19\pmod{20}\)**. On this critical branch even the *full* four digits have period 100 and their parity vector has period 20. For period 100, shifting by 100 changes \(U\) by a multiple of 125, and \(U+U'\) is a multiple of five. Shifting by 20 can change 2528 to 4528 (use \(p\equiv3,23\pmod{100}\) above the transient). A shift by four can change the parity of 8128 to that of 1328 after iterating it. These eliminate all smaller branch-preserving divisors.

The complete critical image at precision \(k\) consists of the squares of all multiples of five for \(8E_p+1\bmod5^k\). Indeed the residues \(16^{(p+1)/4}-1\) traverse every multiple of five modulo \(5^k\). This yields exactly
\[
       1+\sum_{1\le a<k/2}2\cdot 5^{k-2a-1}
\]
values, by classifying the square according to its even valuation \(2a\). At four digits the branch gives 11 strings: 0128, 0528, 1328, 2128, 2528, 4128, 4528, 6128, 6528, 8128, 8528. The regular branch supplies another 125 strings. In particular, on odd \(p\ge5\), **2176 occurs exactly when \(p\equiv61\pmod{500}\)**, since it is the regular value for \(p=61\).

## Inverting an arbitrary exact decimal suffix at any precision

**PROVED.** Given \(k\ge1\) and \(0\le r<10^k\), all odd \(p\ge k+1\) having suffix \(r\) form the union of **at most two arithmetic progressions**. The following procedure returns them without enumerating a full period.

1. Require \(2^k\mid r\). Set \(D=8r+1\). If \(5^k\mid D\), use the single progression \(p\equiv-1\pmod{4\cdot5^{\lceil k/2\rceil-1}}\), by the square identity and the valuation of \(U\).
2. Otherwise \(v_5(D)<k\) must be even, say \(2a\). Compute the two square roots of the unit \(D/5^{2a}\) modulo \(5^{k-2a}\), if they exist. Unit roots modulo five lift uniquely one base-five digit at a time because their derivatives are units.
3. For each unit root \(w\), put \(W=5^a w\pmod{5^{k-a}}\). Retain only \(W\equiv0,3\pmod5\), the two possibilities for odd \(p\). Then solve \(2^{p+1}\equiv1+W\pmod{5^{k-a}}\). Since 2 generates this unit group, each retained root yields one progression of modulus \(4\cdot5^{k-a-1}\). Its discrete logarithm lifts one base-five digit at a time with five trials per digit.

These tests are necessary and sufficient. Dividing the original square relation by \(5^{2a}\) shows that \(U\equiv W\pmod{5^{k-a}}\) is exactly the residue needed, including all singular root lifts in one progression. The first test then restores the binary component of the decimal suffix. For \(k=4,r=2176\) the sole progression is \(61\bmod500\). For \(r=2528\) the two progressions are \(3,35\bmod100\); for \(r=1328\) they coalesce into \(19\bmod20\). Run `python3 -B check.py --suffix 0128` or supply any other zero-padded suffix to apply this inverse.

## Four parity digits by exponent modulo 100

In each row add the column heading to the row heading to get \(p\bmod100\). The domain is prime exponents \(p\ge7\) (the same table holds for every odd \(p\ge5\) coprime to five). E and O are read from thousands to units. A classification on exponent residue classes gives no frequency assertion about the Mersenne primes.

| \(p\bmod20\) | +0 | +20 | +40 | +60 | +80 |
|---:|---|---|---|---|---|
| 1 | OOOE | EOOE | EOOE | EOOE | OOOE |
| 3 | EOEE | EOEE | EOEE | EOEE | EOEE |
| 7 | EOEE | EOEE | EOEE | EOEE | EOEE |
| 9 | EEOE | EEOE | EEOE | OEOE | OEOE |
| 11 | EOEE | EOEE | EOEE | EOEE | EOEE |
| 13 | EOOE | OOOE | OOOE | OOOE | EOOE |
| 17 | OEOE | EEOE | EEOE | OEOE | OEOE |
| 19 | OOEE | OOEE | OOEE | OOEE | OOEE |

The three small primes give 0006, 0028 and 0496. For the 52 exponents on the GIMPS list retrieved 2026-11-22, the pattern counts are EOEE 15, EOOE 9, OEOE 9, EEOE 7, OOOE 7, OOEE 3 and EEEE 2. Among these, 33 centers end in 6 including the initial 6; the other 32 have the asserted odd tens digit. The rank of the latest listed Mersenne prime is provisional. [GIMPS list of known Mersenne primes](https://www.mersenne.org/primes/)

## Applying the projection to the three receiver equations

The full target \(N=E_{61}-1\) has final digits 2175, and the other endpoint has 2177. In the inherited receiver formula \(N=K_j+c y+h f(y)\), the three maximal-order charts have data shown below. The last column discards the unknown nonlinear term exactly by reducing modulo \(h\). [Inherited target selector](../chanpak_e61_fibonacci_target_20261122/SELECTOR.md)

| Rooted receiver rows | \(c,h,K_0,K_1\) | Necessary address condition |
|---|---|---|
| (0) | 81,25,43,68 | \(y\equiv22\pmod{25}\) |
| (2,4,1) | 411,75,367,392 | Neither parity is possible modulo 3 |
| (2,12,9,1) | 1021,125,915,940 | \(y\equiv60\pmod{125}\) if even, \(35\pmod{125}\) if odd |

The triangle receiver is excluded at **every** odd exponent at least three and any suffix width with this core: \(E_p\equiv1\pmod3\) whereas \(K_0,K_1\) are 1 and 2 modulo three. At order 123 the four-vertex chart already forces five terminal ones, so \(y\equiv31\pmod{32}\) and the table gives \(y\equiv1535\pmod{4000}\). Equivalently \(y=32(z+1)-1\) with \(z\equiv47\pmod{125}\). For the singleton at every \(p\equiv3\pmod4\) the last-two-digit law forces \(y\equiv14\pmod{25}\). These are necessary filters on physical constructions. A modular hit retains its exact suffix equation as a separate obligation.
