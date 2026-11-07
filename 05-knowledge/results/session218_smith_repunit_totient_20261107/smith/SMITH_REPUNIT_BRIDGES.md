# Smith numbers, repunit primes and moving-base maps

Status: PROVED (elementary statements below); FINITE-EXACT (receipt only). No fixed-base infinitude or literature priority claim.
Source snapshot: `session216-stage`, remote commit `88ff4b0afd38315a4ad96db6164b80b55488e649`; source checkout was read only. Parent handles repository startup/publishing; this lane writes only new scratch artifacts.

## Inheritance and scope

Closest proved mechanisms: THM-4473, `01-canon/theorems/THM-4473-collatz-digit-chains-rotation-repunits.md` (multidigit rotation fixed points); THM-4749, `01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md` (factorization controls tournament degree); THM-4746 Camion pseudoprimes, whose exact filename should be retained by the integrating parent. Corrected near misses: MISTAKE-651 (omitted nonsquarefree cases), MISTAKE-653 (repunit survey), and the THM-4473 one-digit boundary correction. Smith normal forms elsewhere in the repo are NOT Smith numbers. No theorem here identifies them.

Board: base-dependent digit sum; prime factors with multiplicity; repunit length; carries/block grouping; Camion exponent scaling. Anchor: exact Smith/repunit bridge. Niche: finite all-base sieve. Wildcard: long-prime seed by separated digit boundaries. No intrinsic antisymmetric relation has been identified, so no cosmetic tournament is attached.

Write `s_b(n)` for base-b digit sum, `A_b(n)=sum_p v_p(n)s_b(p)`, and `Delta_b(n)=s_b(n)-A_b(n)`. A base-b Smith number is COMPOSITE and has Delta=0. Define `F(n)=sum_p v_p(n)p`. Set `R_b(k)=(b^k-1)/(b-1)`. Bases b>=2 and lengths k>=2 throughout. Multiplicities are essential: 4 has factor sum 2+2, not 2.

## 1. Finite all-base sieve and the unique unbounded exception (PROVED)

If composite n!=4 is base-b Smith, then `b-1 | n-F(n)` and `b<=n`. Thus all Smith bases for a fixed n!=4 can be found by testing `b=d+1` for positive divisors d of n-F(n). For n=4 the Smith bases are exactly b>=5.

Proof. Digit sums are congruent to their input modulo b-1, so Smith equality implies the divisibility. For integers x,y>=2, x+y<=xy, with equality only x=y=2. Induction on a prime factorization shows F(n)<=n, with equality for a composite only n=4. When b>n all inputs n and p are single digits, so equality is impossible unless n=4. For n=4 bases 2,3,4 have digit sums 1,2,1 versus prime-factor sums 2,4,4 respectively; b>=5 gives 4=2+2. This proves necessity and the exhaustive sieve, not sufficiency of the congruence.

Hostile: n=22 gives n-F(n)=9 and candidate bases 2,4,10, but only 10 is Smith (check exact receipt). A modular filter is not a Smith characterization.

## 2. Moving-base constructions and their exact boundary (PROVED)

For every prime p>=5, both `2p` and `p^2` are Smith in base b=p-1. Indeed p has digits11, 2p has digits22, and p² has digits121, so both composites have digit sum4, matching the multiplicity-counted factor sum4. These digits require b>=3. The prime p=3 gives b=2, where neither 6 nor9 is Smith; p=2 gives the invalid base1.

Consequently 14,26,34,38 are Smith in bases 6,12,16,18 respectively. This is an explicit moving-base embedding of twice-prime integers, not a characterization of the Camion exponent sequence, and not a fixed-base Smith infinitude result. It preserves the factor p as the length-two repunit R_(p-1)(2); it loses any assertion about digit sums in a common base. E.g.14 is not decimal Smith (digit sum5, factor sum9).

There are infinitely many repunit primes when the base may vary: every prime p>=3 is R_(p-1)(2). Euclid's proof supplies infinitely many such p: a finite list cannot contain a prime divisor of one plus its product. This says nothing about infinitely many repunit primes at one fixed base.

**Exact one-digit multiplier classification.** If R_b(k) is prime and 2<=c<b, then cR_b(k) is Smith iff c=k=2. Proof: cR has k digits all c; its prime-factor sum is A_b(c)+k. Since c<b, A_b(c)=F(c)<=c (also when c is prime). Equality requires (c-1)k=F(c)<=c. For k>=2 this forces c=2, then k=2; the converse follows. This includes all bases and explains why longer primes need a multiplier with multiple digits.

**Carry-free powers.** If the coefficients of `(1+x+...+x^(k-1))^a` are all <b and R_b(k) is prime, a>=2, then R_b(k)^a is Smith iff a=k=2. Evaluation at x=1 gives digit sum k^a, while factors contribute ak. The equation k^(a-1)=a has only a=k=2: for a>=3, 2^(a-1)>a. Without the no-carry hypothesis this statement is not claimed.

## 3. Long repunit prime seeds, with a binary obstruction (PROVED)

Let `c=(b-1)t+1`, with t>=1, and choose m>=1 with `t<=R_b(m)`. For every k>=m,

`s_b(c R_b(k)) = k + s_b(t) + s_b(R_b(m)-t) - m`.

Proof: cR_b(k)=t b^k+R_b(k)-t. The low term lies in [0,b^k), and
`R_b(k)-t=b^m R_b(k-m)+(R_b(m)-t)` with the second term in [0,b^m). Both sums have separated blocks, proving the formula; R_b(0)=0.

Hence if `A_b(c)=s_b(t)+s_b(R_b(m)-t)-m`, then every prime R_b(k) with k>=max(m,2) yields a Smith number cR_b(k). Prime factors shared by c and R cause no problem because A_b is completely additive.

Concrete decimal universal seed: c=1540, t=171, m=4. Here 1540=2²*5*7*11, A_10(c)=18, s_10(171)=9, s_10(1111-171)=s_10(940)=13, and 9+13-4=18. Therefore **1540 R_10(k) is Smith whenever k>=4 and R_10(k) is prime**. This is conditional on repunit primality, not a claim that all these products are Smith. Hostile k=6: R=111111=3*7*11*13*37, so A_10(R)=3+7+2+4+10=26, not6. Therefore 1540R has digit sum24 and factor sum44 and fails. At k=4 the composite repunit1111=11*101 happens itself to have factor weight4, so the seed also succeeds there: primality is sufficient, not necessary.

Binary hostile: for any 1<=c<=2^k,
`s_2(c(2^k-1))=k`.
Indeed c(2^k-1)=(c-1)2^k+(2^k-c); the low block is the k-bit complement of c-1. If R_2(k) is prime and c>=2, its prime-factor sum is k+A_2(c)>k. Thus no multiplier 2<=c<=2^k can seed a binary Smith number from that prime. Large multipliers genuinely can: k=3,c=34 gives238=2*7*17, binary digit sum6=factor weight1+3+2. This is not mere base relabelling.

## 4. Carry and base-grouping transport (PROVED)

For positive integers x,y, when multiplying in base b, raw digit convolution has total coefficient sum s_b(x)s_b(y). If the carry out of column j is q_j, the carrying equation gives
`s_b(xy)=s_b(x)s_b(y)-(b-1)sum_j q_j`.
Thus the total carry mass C is a nonnegative integer and
`Delta_b(xy)=Delta_b(x)+Delta_b(y)+s_b(x)s_b(y)-s_b(x)-s_b(y)-(b-1)C`.
In particular Smith is not generally multiplicatively closed. It is never preserved by multiplying a Smith number by b^t for t>=1: digit sum stays constant whereas factor weight grows by t A_b(b)>0.

For an integer h>=1 and B=b^h, write n=sum_i d_i b^i and define
`G_h(n)=sum_(i=hj+r) d_i (b^r-1)/(b-1)`.
Then `s_B(n)=s_b(n)+(b-1)G_h(n)` and
`Delta_B(n)=Delta_b(n)+(b-1)(G_h(n)-sum_p v_p(n)G_h(p))`.
Grouping preserves the integer and prime factorization but destroys Smith equality unless the displayed boundary defect vanishes. This is the necessary sidecar for base changes.

## 5. Natural scaling of Camion repunits is a hostile, not a Smith lift (PROVED)

For every b>=2, b²-1 is NOT a base-b Smith number. For b=2 it is prime. For b>=3 its digits are (b-1)(b-1), sum2b-2. Write n=(b-1)(b+1).
For a composite integer m>=4, F(m)<=m/2+2: split m=uv with u,v>=2, use F(m)<=u+v and (u-2)(v-2)>=0. If b>=8 and b+1 is composite, A_b(b+1)=F(b+1)<b-1, while A_b(b-1)<=b-1. If b+1 is prime, its base-b digit sum is2<b-1 for b>=4. Thus A_b(n)<2b-2 for b>=8. Direct remaining b=3,4,5,6,7 give A_b(n)=6,5,9,7,11 versus digit sums4,6,8,10,12.

Consequently `2^(2r)-1=(2^r)^2-1` is never Smith in its natural grouped base B=2^r. This applies to every even exponent, including the Camion exponents 14,26,34,38, independently of the tournament criterion. In the simultaneous-prime family B-1 and (B+1)/3 prime, r odd>=3, factors are 3, B-1, (B+1)/3, all below B, so A_B=(4B+7)/3 whereas s_B=2B-2. Equality would force B=13/2.

## 6. Repunit factorization and prime order (PROVED)

`R_b(uv)=R_b(u)R_(b^u)(v)`. Thus primality of R_b(k) forces k prime. Also gcd(R_b(k),b-1)=gcd(k,b-1), since each of its k terms is1 modulo b-1. If k=l is prime and R_b(l) prime, l cannot divide b-1 (otherwise l divides R, while R>l). More generally for prime q|R_b(l), either q=l and b=1 mod q, or ord_q(b)=l and l|(q-1). Proof: q cannot divide b; b^l=1 mod q so its order divides l; order1 implies b=1 and then q|l. Conversely the order-l condition implies q|R because b-1 is invertible modulo q. The order divides q-1 because multiplication by b permutes the nonzero residues and the finite subgroup it generates partitions them into equal cosets.

The composite-length factorization preserves arithmetic value and explicit factors; it does not preserve prime status, Smith equality, or tournament degree without their individual sidecars. It aligns with THM-4749's factorization mechanism, rather than declaring a new tournament equivalence.

## Reproduction and audit

Run `python3 session218-smith/check_smith_repunit.py` from the workspace parent (or run the script by its published relative path). The dependency-free script prints the small JSON receipt. This includes positive controls, boundary failures, explicit bounded universes and complete prime factorization by trial division. An independent reviewer separately checked the proof statements and key controls; see the parent session audit. No deterministic primality assertion for enormous decimal repunits is needed by the conditional seed theorem.
