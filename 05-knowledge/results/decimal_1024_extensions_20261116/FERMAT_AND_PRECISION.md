# From 16 to 1024 to the fifth Fermat number

## An affine rule that produces a Fermat divisor

Suppose 2^m=b^r+c with r odd and bc+1=2^ℓ. Then
`D=b·2^m+1=b^(r+1)+2^ℓ`.
Modulo D, b·2^m=−1, so raising to the even exponent r+1 and substituting proves
`D | 2^((r+1)m+ℓ)+1`.
If the exponent on the right is a power of two and D>1 is odd, every prime factor p of D has ord_p(2) exactly twice that exponent.

For normalized 16=10+6 we have 8=5+3; for normalized 1024=1000+24 we have 128=125+3. Both satisfy 5·3+1=16, so the same rule gives:

| Normalized seed | New divisor | Divides |
|---|---:|---:|
| 8=5+3 | 41 | 2^10+1=25·41 |
| 128=125+3 | 641 | 2^32+1=641·6700417 |

The first identity produces a divisor of **1024+1**, not a derivation of the second identity itself. The second produces the famous divisor of the fifth Fermat number. This packages the classical 641 congruence into one rule that explains why 41 and 641 arise from these two matching decimal arrows. The middle synchronized n=2 gives 32=25+7 and fails the pure-power condition 5·7+1=36. So the first two stages are not evidence for an infinite homogeneous chain.

## Two unbounded continuations that do not preserve equality

The equation `2^k−5^n=3` has the two previously classified exact matches (k,n)=(3,1),(7,3). Nevertheless:

* **Concurrent local fitting.** Given arbitrary r≥2 and s≥1 there are arbitrarily large k,n with 2^k−5^n≡3 mod 2^r·5^s. For k≥r and n≥s, independently solve 5^n≡−3 mod 2^r and 2^k≡3 mod 5^s. The first has one exponent class modulo 2^(r−2) because 5 generates the 1 mod 4 group (v₂(5^(2^j)−1)=j+2); the second has one class modulo 4·5^(s−1) because 2 generates mod 5 and v₅(2^4−1)=1. Each class lifts one digit at a time, and can be moved to be arbitrarily large. For the second congruence the least k at s=1,2,3,4,5,6 is 3,7,7,107,607,8107. The cubic match k=7 is an exact local plateau because 2^7−3=125. It breaks exactly at one further digit.
* **Real relative closeness.** For each n choose k=ceil(n log₂5), the smallest binary exponent above 5^n. The relative excess is 2^(k−n log₂5)−1. Irrational rotation gives arbitrarily small positive excesses, which is unsurprising once the problem is projected onto a circle. Exact integer comparisons give the first upper records (n,k)=(1,3),(2,5),(3,7),(31,72),(59,137). At n=31 the relative excess 0.01412… beats 3/125, but its signed integer defect is 65753609792252635571 rather than 3. Each type of approximation discards precisely what the other needs.

The reusable distinction is between more matching digits, smaller normalized error, and exact equality of two actual integers. This applies to extrapolating a residue controller into a prescribed tournament count; the note does not construct such a tournament. For wider context on controlling gaps between powers rather than assuming relative closeness gives an integer match, see Yann Bugeaud, “B′” (https://arxiv.org/abs/2209.00275).

The common second carry 5·3+1=16 is itself the ordinary Collatz arrow 5→16. Swapping 3 and 5 in the cubic seed gives 2^5=3³+5, the other sign/mode of the cubic near-power, and D=3·32+1=97 divides 2^24+1 by the same rule. This is the base-6 identity 256=216+5·8; 3→8 is not an ordinary Collatz arrow.

There is a classification that unifies exactly these two radices. Require the same correction at n=1 and n=3, under K(n)=Au^n−v^n with consecutive v=u+1 and u≥2. The match forces A=(v³−v)/(u³−u)=(u+2)/(u−1). Hence (u,v,A)=(2,3,4) or (4,5,2). When the unnormalized radix is twice v, these are exactly `16=6+5·2, 256=6³+5·8` and `16=10+3·2, 1024=10³+3·8`. The 1-and-3 collision comes from the ratio of two tetrahedral numbers, binom(v+1,3)/binom(u+1,3), rather than an arbitrary comparison of near powers.

## One variable can match arbitrarily many decimal digits

The actual arrow defect K(n)=2·4^n−5^n has a stronger property than independent local fitting. For r≥3 restrict to odd n≥r and write n=1+2t. Modulo 2^r, K≡−5·25^t (the positive binary term has vanished); modulo 5^r, K≡8·16^t. The first expression hits the entire class 3 mod 8 with t of period 2^(r−3), the second the entire class 3 mod 5 with period 5^(r−1), by the valuation of 25−1 and 16−1. These periods are coprime. Thus K hits **every residue 3 mod 40** modulo 10^r, once per odd-source exponent period 2^(r−2)5^(r−1). Exact K=3 never recurs beyond n=3. The first next solutions to K≡3 after enforcing n≥r are n=303, 2803, 22803, 272803 for r=4,5,6,7 respectively. In particular the arrow 303→910 ends in the requested phase 0003 in this signed normalization without a remotely small absolute defect.

## Turn the gap into a ternary clock

Let q=128/125=1+3/125. The integer sequence U_n=(128^n−125^n)/3 satisfies U_0=0, U_1=1 and U_(n+2)=253U_(n+1)−16000U_n. Since 128 and 125 are coprime, the Euclidean algorithm applied to their power differences proves gcd(U_m,U_n)=U_gcd(m,n) (each smaller power of 125 is a unit modulo the relevant difference). Also v₃(U_n)=v₃(n) by binomial lifting. Normalize C_n=U_n/125^(n−1). Then
`C_(m+n)=C_m+C_n+(3/125)C_m C_n,quad v₃(C_m−C_n)=v₃(m−n)`.
So the seed's one normalized gap of 3 generates a coherent bijection modulo every power of 3, and an isometry and bijection of the 3-adic integers by passage to the limit. The clock does not encode parity modulo any 3^r; an actual Collatz transport must retain that coordinate separately. In the variable y=1+(3/125)C_n=q^n, the transported map acts as sqrt(y) for even n and qy³ for odd n.

A small MRB-inspired analytic analogue uses the *same* q in the alternating series Σ_(n≥1)(−1)^n(q^(1/n)−1). Writing L=log(q) and expanding, the series equals −Σ_(j≥1)L^j η(j)/j!, where η(j)=Σ_(n≥1)(−1)^(n−1)n^(−j). Treat the j=1 term by alternating convergence; the j≥2 double sum is absolutely convergent, since ζ(j)≤ζ(2). With 0<η(j)≤1 the error of the linear approximation −L log 2 has magnitude at most q−1−L. This is a controlled analogue, not an identification with the original MRB constant, which has the moving logarithm log(n).

A second MRB bridge comes from making the fitting radix real: B_n=8(2−3/4^n)^(1/n) is the unique real radix for the synchronized +3 equation (B₁=B₃=10). Replacing B_n by 8·2^(1/n) makes an error O(4^(−n)/n) by the mean value theorem. Thus its alternating radix-repair series has an exponentially convergent correction to the fixed-log slice of
`F(a,b)=Σ_(n≥1)(−1)^(n−1)(exp((a+b log n)/n)−1)`.
The fixed ratio above is F(log q,0) with the opposite sign; the moving-root MRB slice is F(0,1), also up to sign convention. Expanding gives Σ_(j≥1)[(a−b∂_s)^j η(s)]_(s=j)/j!: the j=1 terms converge separately by alternation, and the remaining double sum is bounded absolutely by a constant times Σ(1+log n)²/n².

The base-six sibling has the same local/global split. For odd n sufficiently large, K₆(n)=4·2^n−3^n ranges over exactly 5 mod 24 modulo 6^r for r≥3: the two odd-source clocks are −3·9^t mod 2^r and 8·4^t mod 3^r, of coprime orders. The exact coefficient 5 occurs only at n=1,3 (the next term is −17 and the sequence decreases past n=2). The minus-Collatz branch, by contrast, has stable odd-source K₋(n)=2^(2n−1)−5^n≡27 mod 40, so its lone exact −3 match cannot be reproduced even mod 40 by large odd n.

There is also a positive counting interpretation of the ternary clock's U_n. Use an alphabet of 128 symbols, with 125 ordinary and three marked. Cyclic rotation of the three marked symbols is free on words containing at least one marked symbol. Counting orbits by the first marked position gives U_n=Σ_(j=0)^(n−1)125^j128^(n−1−j), with no cancellation. In particular the recurrence with a negative term describes a family with a positive, first-marked-position decomposition.

The totient side is separate from the additive identity but exact: φ(1024)=512, and its totient fiber has eleven inputs. An odd prime dividing a preimage has exponent one and p−1 a power of two; primality forces a Fermat exponent and the candidates ≤513 are 3,5,17,257. Their p−1 exponents are 1,2,4,8. For each s=0,…,9 take the unique subset with exponent sum s, multiply its primes by 2^(10−s), and additionally take the odd preimage 3·257=771. The same two inputs at s=9 differ by a factor of two without changing the totient. This is a particular explicit fiber, not a claim about Carmichael’s conjecture in general.
