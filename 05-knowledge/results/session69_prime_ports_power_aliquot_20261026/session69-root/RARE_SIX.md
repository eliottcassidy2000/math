# Why the self-return with a prime-square companion is unique

**PROVED, elementary consequence of the classical Euclid–Euler structure.** No novelty claim for that classical structure. This strengthens the statement “six is first” to uniqueness for a precisely defined conjunction, without assuming that odd perfect numbers do not exist or that there are infinitely many Mersenne primes.

Let s(n)=sigma(n)-n for positive integers. Then

    s(N)=N and s(p²)=N with p prime
       if and only if N=6 and p=5.

So six is the unique perfect target whose inverse aliquot fibre contains a prime square; its companion source is25. This is not uniqueness of all multiple fibres, all self-return notions, or all perfect numbers.

## Proof

For prime p, s(p²)=p+1, hence N=p+1. If p=2 then N=3 is not perfect. Otherwise N is even. Write N=2^a m, with a>=1 and m odd. Put q=2^(a+1)-1. Perfection gives

    q sigma(m)=2^(a+1)m.

Since gcd(q,2^(a+1))=1, q divides m. Write m=qt; then

    sigma(m)=(q+1)t=m+t.

If t>1, the three distinct divisors1,t,m already give sigma(m)>=1+t+m, a contradiction. Therefore t=1, m=q, and sigma(q)=q+1, forcing q to be prime. The exponent r=a+1 must then be prime too: a composite exponent would factor2^r-1 nontrivially. This is the classical even-perfect-number argument, included to make the implication explicit.

For r=2 we obtain N=6, and p=N-1=5. If r is odd, then modulo3,

    2^(r-1)=1,   2^r-1=1,   N=1.

Thus p=N-1 is divisible by3. Primality forces p=3, giving N=4, which is not perfect. Hence r=2 is the only possibility. Conversely s(6)=6 and s(25)=6 by direct divisor sums. QED.

## Microcosm/macrocosm boundary

The twin-prime mechanism s(2p)=s((p+2)²)=p+3 explains the first two finite multiple fibres at6 and8. Its p=3 instance is uniquely self-returning within that family. The theorem above is stronger: it excludes ANY other perfect target with a prime-square companion, regardless of whether its other source has the form2p.

This does not make the earlier power/doubling event equally rare. That event has density1/2 under its declared odd-root convention. Nor is the self-return an exceptional primitive ABC triple: normalizing N+s(N)=sigma(N) always yields1+1=2. The integer N and its removed gcd are essential data.

## Finite controls

`check_rare_six.py` independently checks all sources up to10000 by a divisor sieve, confirms the perfect numbers6,28,496,8128 in that window, and tests the prime-square-companion predicate. It additionally checks all prime Mersenne exponents through19 plus31 against the modular proof. Those tests corroborate the argument; the all-order uniqueness comes from the proof above, not the cutoff.
