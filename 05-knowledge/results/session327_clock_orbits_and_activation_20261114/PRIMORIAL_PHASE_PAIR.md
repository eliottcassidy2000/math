# The 30,030 primorial and the exact 99/100 paired sieve

**Proved:** exact finite-wheel periods/densities, and a composite lift in every surviving cylinder. The two inputs are the literal strong tournament rails from Session 326, so this is an arithmetic sieve on physically attained counts. It is not a primality criterion, an infinitude theorem, or an absolute odd-spectrum proof.

## 1. Which 30,000 boundary?

The [least-resistance Paley note](../paley_least_resistance_runs_one_modular_form_20261014.md) separates a proved series/parallel reduction from experimentally found modular formulas. Its [adjacent six-vertex terminal](../paley_k6p6_adjoint_weight3_level32_20261014.md) uses 30,000 as a finite scan cutoff for a modular-form identification, not a primorial or a proven change of behavior. The genuine primorial here is \(13\#=30030\). The [old ABC-gap note](../abc_difference_families_20261011c.md) already has a fixed-wheel composite progression \(9450+30030t\) for additive neighbors. The construction below transports that **failure mechanism**, not the additive progression: its native progression is in the exponent, because powers of two propagate by multiplicative orders.

A small explicit check distinguishes the absolute gap questions: \(29989,30011\) are consecutive primes on either side of 30000 (gap 22), whereas \(30029,30047\) are consecutive around 30030 (gap 18). Symmetric primes \(30030\pm17\) are not consecutive because 30029 intervenes. The two aligned rail targets \(X_m=99\cdot2^m+23\) and \(Y_m=100\cdot2^m+23\) both first exceed 30000 at \(m=9\), and then both fail primality despite surviving the primorial wheel:
\[
X_9=50711=17\cdot2983,\quad Y_9=51223=181\cdot283.
\]
The corresponding tournament order is 18. This is a falsifiable first-boundary prediction, not the claim that prime gaps control these rails. The relative gap tends to \(1/100\), while the exact difference is \(2^m\); doubling both coefficients only retimes the common exponent.

## 2. Their shared \(13\#\)-roughness has period 60

For each odd prime \(p\) with \(p\nmid A\), a hit on \(A\cdot2^m+23\) requires \(2^m\equiv-23A^{-1}\pmod p\). There is at most one such class modulo \(\mathrm{ord}_p(2)\). A zero target on the right is impossible for powers of two; if \(p\mid A\) check the constant 23 separately. Both rails are odd for \(m\ge1\).

| rail | active forbidden exponent classes for \(p\le13\) |
|---|---|
| \(X\) | \(3\bmod4\) from 5; \(1\bmod12\) from 13 |
| \(Y\) | \(0\bmod2\) from 3; \(5\bmod10\) from 11; \(8\bmod12\) from 13 (already even) |

Combined survival is exactly
\[
m\bmod60\in\{9,17,21,29,33,41,53,57\}.
\]
Indeed \(m\) is odd and \(1\bmod4\); remove \(1\bmod12\), leaving \(5,9\bmod12\); finally remove \(5\bmod10\). Period 60 is least: shifts \(60/2,60/3,60/5\) send respective survivors \(9,41,33\) to nonsurvivors.

Each one of the eight surviving classes contains an explicitly forced simultaneous-composite progression. Impose additionally \(m\equiv1\bmod8\) and \(m\equiv2\bmod23\). The former gives \(17\mid X_m\), compatible because every survivor is \(1\bmod4\); the latter gives \(47\mid Y_m\) and its modulus is coprime to 120. Their common stride is 2760. In the order of survivors displayed above, representative lifts modulo 2760 are
\[
2049,\ 2417,\ 2601,\ 209,\ 393,\ 761,\ 1313,\ 1497.
\]
The factors are proper for these positive exponents, and remain proper indefinitely. This is a repeatable failure certificate for **every** shared primorial survivor, with the same two auxiliary primes for all eight.

## 3. Exact first-25-prime comparison

The first 25 primes are precisely those below 100. The table uses least periods and counts residues in one such period. “Joint” means **the same integer \(m\)** passes both sieves, not independent choices of an exponent.

| wheel | least period | survivors | density |
|---|---:|---:|---:|
| \(X\) | 61,209,720 | 32,353,278 | \(5392213/10201620\) |
| \(Y\) | 169,909,740 | 39,935,280 | \(8644/36777\) |
| joint | 9,854,764,920 | 184,192,008 | \(187/10005\) |

The 99-only row agrees with the earlier independent first-25-prime sieve. The other rows are new comparisons. The joint density is about 1.869%, versus about 12.423% for the product of the individual densities. The two survival tests are sharply correlated against each other on the shared exponent; a formal product of densities is false.

### Exact computation without scanning billions of exponents

Maintain the survivor indicator as a signed sum of congruence-class indicators. Begin with the full class modulo 1. Each forbidden congruence multiplies by its complementary indicator; intersect the old classes with it using generalized CRT, subtract the intersection contributions, and coalesce identical classes. For a known period \(M\), a class of modulus \(d\mid M\) contributes exactly \(M/d\) times its coefficient. This supplies all counts without probability or scan extrapolation.

For clarity, the joint density also admits a small direct CRT derivation. First apply the active forbidden moduli dividing 360, leaving
\[
29,69,77,141,149,173,213,293,317,357\pmod{360}.
\]
Among these ten base classes, the counts for \((m\equiv2\bmod3,\ m\equiv2\bmod5)\) are respectively \(3,1,4,2\) in the order \((0,0),(0,1),(1,0),(1,1)\). Outside the base the extra 11-coordinate loses a fraction \(1/11\) only on the first tagged group (the 67-hit on \(Y\)); the extra 7-coordinate loses \(1/7\) only on the second tagged group (the 71-hit on \(Y\)). The 29-hits cannot occur on the base classes. The independent extra coordinates contribute \(22/23\) (47), \(28/29\) (59), \(39/41\) (83, two distinct hits), and \(11/13\) (53, two distinct hits). Precisely:
\[
\frac{3+1(6/7)+4(10/11)+2(10/11)(6/7)}{360}
\frac{22}{23}\frac{28}{29}\frac{39}{41}\frac{11}{13}
=\frac{187}{10005}.
\]
Here all extra prime coordinates are coprime to 360 and pairwise independent by CRT; the explicit tags retain exactly the low-prime dependencies.

**Least-period certificate.** A period is the lcm of remaining moduli after deleting individually dominated forbidden classes. The 99-side dominated primes are 29,37; the 100-side ones are 13,29,41,59. This produces the first two displayed candidate periods; their lcm is a candidate for the joint sieve. If any candidate \(M\) were not least, then \(M/p\) would also be a period for some prime divisor \(p\mid M\), since least periods divide every period. In each following row, \(p:k\) means \(k\) survives but \(k+M/p\) fails. These are finite witnesses using only the active prime congruences.

| wheel | witnesses \(p:k\) |
|---|---|
| \(X\) | 2:5, 3:21, 5:36, 11:62, 13:29, 29:149, 41:69 |
| \(Y\) | 2:3, 3:1, 5:7, 7:27, 11:131, 13:33, 23:61, 41:267 |
| joint | 2:69, 3:69, 5:141, 7:4397, 11:533, 13:149, 23:1469, 29:3309, 41:509 |

## 4. A universal failure mechanism without factoring

Fix **any** finite prime wheel, and any exponent \(m\ge1\) for which both rails survive it. Let \(a=X_m\), \(b=Y_m\). They are odd and \(\gcd(a,b)=1\), since their difference is \(2^m\). The integer 2 is a unit modulo both. Choose a positive stride \(D\) divisible by both multiplicative orders and by the finite sieve period. Then
\[
a\mid X_{m+jD},\qquad b\mid Y_{m+jD}\qquad(j\ge1)
\]
by exact substitution, while these new exponents stay in the original surviving cylinder. The divisors \(a,b\) are proper by strict growth. No factorization of \(a,b\) or primitive-divisor theorem is needed; for existence one can take the wheel period times \((ab)!\), which is a multiple of both orders. Thus **every surviving cylinder of every finite wheel** contains infinitely many pairs with both strong tournament counts composite. This does not rule out infinitely many prime pairs there as well. It rules out proving even pair primality from survival in any fixed finite wheel.

This provides a distinct version of the native continuation warning: two physically valid channels and finitely many exact modular guards do not make an absolute prime selector. The full wheel is still useful as a falsifiable filter, with its actual shared phase and joint density rather than an independence estimate.

Replay: \`python3 replay_primorial_wheel.py\`. Normal and optimized receipts agree.
